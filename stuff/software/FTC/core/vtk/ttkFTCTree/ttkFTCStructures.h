#ifndef TTKFTCSTRUCTURES_H
#define TTKFTCSTRUCTURES_H

#include <FTCTree.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkCharArray.h>

struct LocalFTC {
   ftc::FTCTree tree;
   ftc::idNode  offset;
};

struct WrapperData {
   template<typename vtkArrayType>
   inline vtkSmartPointer<vtkArrayType> initArray(const char* fieldName, size_t nbElmnt)
   {
      vtkSmartPointer<vtkArrayType> arr = vtkSmartPointer<vtkArrayType>::New();
      arr->SetName(fieldName);
      arr->SetNumberOfComponents(1);
      arr->SetNumberOfTuples(nbElmnt);

#ifndef TTK_ENABLE_KAMIKAZE
      if (!arr) {
         cerr << "[ttkFTCTree] Error, unable to allocate " << fieldName
              << " the program will likely crash" << endl;
      }
#endif
      return arr;
   }

   inline static ftc::NodeType getNodeType(ftc::FTCTree_MT& tree, const ftc::idNode nodeId,
                                           ftc::Params params)
   {
      const ftc::Node* node = tree.getNode(nodeId);
      int upDegree{};
      int downDegree{};
      if (params.treeType == ftc::TreeType::Join or params.treeType == ftc::TreeType::Contour) {
         upDegree   = node->getNumberOfUpSuperArcs();
         downDegree = node->getNumberOfDownSuperArcs();
      } else {
         downDegree = node->getNumberOfUpSuperArcs();
         upDegree   = node->getNumberOfDownSuperArcs();
      }
      int degree = upDegree + downDegree;

      // saddle point
      if (degree > 1) {
         if (upDegree == 2 and downDegree == 1)
            return ftc::NodeType::Saddle2;
         else if (upDegree == 1 and downDegree == 2)
            return ftc::NodeType::Saddle1;
         else if (upDegree == 1 and downDegree == 1)
            return ftc::NodeType::Regular;
         else
            return ftc::NodeType::Degenerate;
      }
      // local extremum
      else {
         if (upDegree)
            return ftc::NodeType::Local_minimum;
         else
            return ftc::NodeType::Local_maximum;
      }
   }
};

struct ArcData : public WrapperData {
   vector<vtkIdType>               point_ids;
   vtkSmartPointer<vtkCharArray>   point_regularMask;
   vtkSmartPointer<vtkFloatArray>  point_scalar;
   vtkSmartPointer<vtkIntArray>    cell_ids;
   vtkSmartPointer<vtkIntArray>    cell_sizeArcs;
   vtkSmartPointer<vtkDoubleArray> cell_spanArcs;

   inline int init(vector<LocalFTC>& ftcTree, ftc::Params params)
   {
      ftc::idSuperArc nbArcs       = 0;
      ftc::idSuperArc nbNodes      = 0;
      ftc::idSuperArc samplePoints = 0;
      ftc::idVertex   nbVerts      = 0;

      for (auto& t : ftcTree) {
         ftc::FTCTree_MT* tree = t.tree.getTree(params.treeType);
         nbArcs += tree->getNumberOfSuperArcs();
         nbNodes += tree->getNumberOfNodes();
         samplePoints += params.samplingLvl >= 0
                             ? tree->getNumberOfNodes() + (nbArcs * params.samplingLvl)
                             : tree->getNumberOfVertices();
         nbVerts += tree->getNumberOfVertices();
      }

      point_ids.resize(nbVerts, ftc::nullVertex);
      cell_ids          = initArray<vtkIntArray>("SegmentationId", samplePoints);
      point_regularMask = initArray<vtkCharArray>("RegularMask", samplePoints);
      point_scalar      = initArray<vtkFloatArray>("Scalar", samplePoints);

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs = initArray<vtkIntArray>("RegionSize", samplePoints);
         }
         cell_spanArcs = initArray<vtkDoubleArray>("RegionSpan", samplePoints);
      }

      return 0;
   }

   inline bool hasPoint(const ftc::idVertex vertId)
   {
      return point_ids[vertId] != ftc::nullVertex;
   }

   inline void addPoint(const ftc::idVertex globalId, const vtkIdType id, const float scalar,
                        const bool reg)
   {
      point_ids[globalId] = id;
      setPoint(id, scalar, reg);
   }

   inline void setPoint(const vtkIdType id, const float scalar, const bool reg)
   {
      point_scalar->SetTuple1(id, scalar);
      point_regularMask->SetTuple1(id, reg);
   }

   inline void fillArrayCell(const vtkIdType pos, const ftc::idSuperArc arcId, LocalFTC& ftcTree,
                             Triangulation* triangulation, ftc::Params params)
   {
      const ftc::idNode idOffset = ftcTree.offset;
      ftc::FTCTree_MT*  tree     = ftcTree.tree.getTree(params.treeType);
      ftc::SuperArc*    arc      = tree->getSuperArc(arcId);

      if (params.normalize) {
         cell_ids->SetTuple1(pos, idOffset + arc->getNormalizedId());
      } else {
         cell_ids->SetTuple1(pos, idOffset + arcId);
      }

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs->SetTuple1(pos, tree->getArcSize(arcId));
         }

         float               downPoints[3];
         const ftc::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ftc::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ftc::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ftc::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         cell_spanArcs->SetTuple1(pos, Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkUnstructuredGrid* skeletonArcs, ftc::Params params)
   {
      // Some arcs might have been less sampled than the desired value, if they have not enought
      // regular vertices. Here we ensur that we will no keep noise in these arrays.
      const size_t nbPoints = skeletonArcs->GetNumberOfPoints();
      const size_t nbCells  = skeletonArcs->GetNumberOfCells();

      cell_ids->SetNumberOfTuples(nbCells);
      skeletonArcs->GetCellData()->AddArray(cell_ids);

      if (params.advStats) {
         if (params.segm) {
            cell_sizeArcs->SetNumberOfTuples(nbCells);
            skeletonArcs->GetCellData()->AddArray(cell_sizeArcs);
         }
         cell_spanArcs->SetNumberOfTuples(nbCells);
         skeletonArcs->GetCellData()->AddArray(cell_spanArcs);
      }

      point_scalar->SetNumberOfTuples(nbPoints);
      skeletonArcs->GetPointData()->AddArray(point_scalar);
      point_regularMask->SetNumberOfTuples(nbPoints);
      skeletonArcs->GetPointData()->AddArray(point_regularMask);

      point_ids.clear();
   }
};

struct NodeData : public WrapperData{
   vtkSmartPointer<vtkIntArray> ids;
   vtkSmartPointer<vtkIntArray> vertIds;
   vtkSmartPointer<vtkIntArray> type;
   vtkSmartPointer<vtkIntArray> regionSize;
   vtkSmartPointer<vtkIntArray> regionSpan;

   inline int init(vector<LocalFTC>& ftcTree, ftc::Params params)
   {
      ftc::idNode numberOfNodes = 0;
      for (auto& t : ftcTree) {
         ftc::FTCTree_MT* tree = t.tree.getTree(params.treeType);
         numberOfNodes += tree->getNumberOfNodes();
      }

      ids     = initArray<vtkIntArray>("NodeId", numberOfNodes);
      vertIds = initArray<vtkIntArray>("VertexId", numberOfNodes);
      type    = initArray<vtkIntArray>("NodeType", numberOfNodes);

      if (params.advStats) {
         if (params.segm) {
            regionSize = initArray<vtkIntArray>("RegionSize", numberOfNodes);
         }
         regionSpan = initArray<vtkIntArray>("RegionSpan", numberOfNodes);
      }

      return 0;
   }

   inline void fillArrayPoint(vtkIdType arrIdx, const ftc::idNode nodeId, LocalFTC& ftcTree,
                              vtkDataArray* idMapper, Triangulation* triangulation,
                              ftc::Params params)
   {
      const ftc::idNode   idOffset   = ftcTree.offset;
      ftc::FTCTree_MT*    tree       = ftcTree.tree.getTree(params.treeType);
      const ftc::Node*    node       = tree->getNode(nodeId);
      const ftc::idVertex l_vertexId = node->getVertexId();
      const ftc::idVertex g_vertexId = idMapper->GetTuple1(l_vertexId);

      ids->SetTuple1(arrIdx, idOffset + nodeId);
      vertIds->SetTuple1(arrIdx, g_vertexId);
      type->SetTuple1(arrIdx, static_cast<int>(getNodeType(*tree, nodeId, params)));

      if (params.advStats) {
         ftc::idSuperArc saId = getAdjSa(node);
         if (params.segm) {
            regionSize->SetTuple1(arrIdx, tree->getArcSize(saId));
         }

         ftc::SuperArc* arc = tree->getSuperArc(saId);

         float               downPoints[3];
         const ftc::idVertex downNodeId   = tree->getLowerNodeId(arc);
         const ftc::idVertex downVertexId = tree->getNode(downNodeId)->getVertexId();
         triangulation->getVertexPoint(downVertexId, downPoints[0], downPoints[1], downPoints[2]);

         float               upPoints[3];
         const ftc::idVertex upNodeId   = tree->getUpperNodeId(arc);
         const ftc::idVertex upVertexId = tree->getNode(upNodeId)->getVertexId();
         triangulation->getVertexPoint(upVertexId, upPoints[0], upPoints[1], upPoints[2]);

         regionSpan->SetTuple1(arrIdx, Geometry::distance(downPoints, upPoints));
      }
   }

   inline void addArray(vtkPointData* pointData, ftc::Params params)
   {
      pointData->AddArray(ids);
      pointData->AddArray(vertIds);
      pointData->AddArray(type);
      if (params.advStats) {
         if (params.segm) {
            pointData->AddArray(regionSize);
         }
         pointData->AddArray(regionSpan);
      }
   }

  private:

   ftc::idSuperArc getAdjSa(const ftc::Node* node)
   {
      if (node->getNumberOfDownSuperArcs() == 1) {
         return node->getDownSuperArcId(0);
      }

      if (node->getNumberOfUpSuperArcs() == 1) {
         return node->getUpSuperArcId(0);
      }

      // Degenerate case, arbitrary choice
      if (node->getNumberOfDownSuperArcs()) {
         return node->getDownSuperArcId(0);
      }

      if (node->getNumberOfDownSuperArcs()) {
         return node->getDownSuperArcId(0);
      }

      // Empty node
      cerr << "[ttkFTCTree]: node without arcs:" << node->getVertexId() << endl;
      return ftc::nullSuperArc;
   }
};

struct VertData: public WrapperData {
   vtkSmartPointer<vtkIntArray>    ids;
   vtkSmartPointer<vtkIntArray>    sizeRegion;
   vtkSmartPointer<vtkDoubleArray> spanRegion;
   vtkSmartPointer<vtkCharArray>   typeRegion;

   inline int init(vector<LocalFTC>& ftcTrees, ftc::Params params)
   {
      if (!params.segm)
         return 0;

      ftc::idVertex numberOfVertices = 0;

      for (auto& t : ftcTrees) {
         ftc::FTCTree_MT* tree = t.tree.getTree(params.treeType);
         numberOfVertices += tree->getNumberOfVertices();
      }

      ids        = initArray<vtkIntArray>("SegmentationId", numberOfVertices);
      typeRegion = initArray<vtkCharArray>("RegionType", numberOfVertices);

      if (params.advStats) {
         sizeRegion = initArray<vtkIntArray>("RegionSize", numberOfVertices);
         spanRegion = initArray<vtkDoubleArray>("RegionSpan", numberOfVertices);
      }

      return 0;
   }

   void fillArrayPoint(const ftc::idSuperArc arcId, LocalFTC& l_tree, Triangulation* triangulation,
                       vtkDataArray* idMapper, ftc::Params params)
   {
      if (!params.segm)
         return;

      ftc::FTCTree_MT* tree = l_tree.tree.getTree(params.treeType);
      const ftc::idNode idOffset = l_tree.offset;
      ftc::SuperArc* arc = tree->getSuperArc(arcId);

      const ftc::idNode   upNodeId     = arc->getUpNodeId();
      const ftc::Node*    upNode       = tree->getNode(upNodeId);
      const ftc::idVertex l_upVertexId = upNode->getVertexId();
      const ftc::idVertex g_upVertexId = idMapper->GetTuple1(l_upVertexId);
      const ftc::NodeType upNodeType   = getNodeType(*tree, upNodeId, params);
      float               coordUp[3];
      triangulation->getVertexPoint(l_upVertexId, coordUp[0], coordUp[1], coordUp[2]);

      const ftc::idNode   downNodeId     = arc->getDownNodeId();
      const ftc::Node*    downNode       = tree->getNode(downNodeId);
      const int           l_downVertexId = downNode->getVertexId();
      const int           g_downVertexId = idMapper->GetTuple1(l_downVertexId);
      const ftc::NodeType downNodeType   = getNodeType(*tree, downNodeId, params);
      float               coordDown[3];
      triangulation->getVertexPoint(l_downVertexId, coordDown[0], coordDown[1], coordDown[2]);

      const int    regionSize = tree->getSuperArc(arcId)->getNumberOfRegularNodes();
      const double regionSpan = Geometry::distance(coordUp, coordDown);

      ftc::idSuperArc nid = arc->getNormalizedId();

      ftc::ArcType regionType;
      // RegionType
      if (upNodeType == ftc::NodeType::Local_minimum && downNodeType == ftc::NodeType::Local_maximum)
         regionType = ftc::ArcType::Min_arc;
      else if (upNodeType == ftc::NodeType::Local_minimum || downNodeType == ftc::NodeType::Local_minimum)
         regionType = ftc::ArcType::Min_arc;
      else if (upNodeType == ftc::NodeType::Local_maximum || downNodeType == ftc::NodeType::Local_maximum)
         regionType = ftc::ArcType::Max_arc;
      else if (upNodeType == ftc::NodeType::Saddle1 && downNodeType == ftc::NodeType::Saddle1)
         regionType = ftc::ArcType::Saddle1_arc;
      else if (upNodeType == ftc::NodeType::Saddle2 && downNodeType == ftc::NodeType::Saddle2)
         regionType = ftc::ArcType::Saddle2_arc;
      else
         regionType = ftc::ArcType::Saddle1_saddle2_arc;

      // fill extrema and regular verts of this arc

      // critical points
      if(params.normalize){
         ids->SetTuple1(g_upVertexId   , idOffset + nid);
         ids->SetTuple1(g_downVertexId , idOffset + nid);
      } else{
         ids->SetTuple1(g_upVertexId   , idOffset + arcId);
         ids->SetTuple1(g_downVertexId , idOffset + arcId);
      }

      if (params.advStats) {
         sizeRegion->SetTuple1(g_upVertexId   , regionSize);
         sizeRegion->SetTuple1(g_downVertexId , regionSize);
         spanRegion->SetTuple1(g_upVertexId   , regionSpan);
         spanRegion->SetTuple1(g_downVertexId , regionSpan);
      }
      typeRegion->SetTuple1(g_upVertexId   , static_cast<char>(regionType));
      typeRegion->SetTuple1(g_downVertexId , static_cast<char>(regionType));

      // regular nodes
      for (const ftc::idVertex l_vertexId : *arc) {
         const ftc::idVertex g_vertexId = idMapper->GetTuple1(l_vertexId);
         if (params.normalize) {
            ids->SetTuple1(g_vertexId, idOffset + nid);
         } else {
            ids->SetTuple1(g_vertexId, idOffset + arcId);
         }
         if (params.advStats) {
            sizeRegion->SetTuple1(g_vertexId, regionSize);
            spanRegion->SetTuple1(g_vertexId, regionSpan);
         }
         typeRegion->SetTuple1(g_vertexId, static_cast<char>(regionType));
      }

   }

   void addArray(vtkPointData* pointData, ftc::Params params)
   {
      if (!params.segm)
         return;

      pointData->AddArray(ids);

      if (params.advStats) {
         pointData->AddArray(sizeRegion);
         pointData->AddArray(spanRegion);
      }
      pointData->AddArray(typeRegion);
   }
};

#endif /* end of include guard: TTKFTCSTRUCTURES_H */
