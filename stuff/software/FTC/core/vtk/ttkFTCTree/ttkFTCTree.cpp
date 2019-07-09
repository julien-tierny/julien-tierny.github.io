#include <ttkFTCTree.h>

// only used on the cpp
#include <vtkConnectivityFilter.h>
#include <vtkDataObject.h>
#include <vtkThreshold.h>

using namespace ftc;

vtkStandardNewMacro(ttkFTCTree);

int ttkFTCTree::FillInputPortInformation(int port, vtkInformation* info)
{
   if (port == 0)
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
   return 1;
}

int ttkFTCTree::FillOutputPortInformation(int port, vtkInformation* info)
{
   switch (port) {
      case 0:
      case 1:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
         break;

      case 2:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
         break;
   }

   return 1;
}

int ttkFTCTree::addCompleteSkeletonArc(const ftc::idSuperArc arcId, const int cc, vtkPoints* points,
                                       vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   FTCTree_MT*   tree     = ftcTree_[cc].tree.getTree(GetTreeType());
   vtkDataArray* idMapper = connected_components_[cc]->GetPointData()->GetArray("VertexIdentifier");
   SuperArc*     arc      = tree->getSuperArc(arcId);
   float         point[3];
   vtkIdType     pointIds[2];

   const idVertex downNodeId     = tree->getLowerNodeId(arc);
   const idVertex l_downVertexId = tree->getNode(downNodeId)->getVertexId();
   const idVertex g_downVertexId = idMapper->GetTuple1(l_downVertexId);
   triangulation_[cc]->getVertexPoint(l_downVertexId, point[0], point[1], point[2]);
   const double scalarMin = inputScalars_[cc]->GetTuple1(l_downVertexId);

   // Get or create first point of the arc
   vtkIdType nextPointId;
   if (!arcData.hasPoint(g_downVertexId)) {
      nextPointId = points->InsertNextPoint(point);
      arcData.addPoint(g_downVertexId, nextPointId, scalarMin, false);
   } else {
      nextPointId = arcData.point_ids[g_downVertexId];
   }

   pointIds[0] = nextPointId;

   for (const idVertex vertexId : *arc) {
      triangulation_[cc]->getVertexPoint(vertexId, point[0], point[1], point[2]);
      pointIds[1] = points->InsertNextPoint(point);
      const double scalar = inputScalars_[cc]->GetTuple1(vertexId);
      arcData.setPoint(pointIds[1], scalar, true);

      const vtkIdType nextCell = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
      arcData.fillArrayCell(nextCell, arcId, ftcTree_[cc], triangulation_[cc], params_);

      pointIds[0] = pointIds[1];
   }

   const idVertex upNodeId     = tree->getUpperNodeId(arc);
   const idVertex l_upVertexId = tree->getNode(upNodeId)->getVertexId();
   const idVertex g_upVertexId = idMapper->GetTuple1(l_upVertexId);
   triangulation_[cc]->getVertexPoint(l_upVertexId, point[0], point[1], point[2]);
   const double scalarMax = inputScalars_[cc]->GetTuple1(l_upVertexId);

   // Get or create last point of the arc
   if (!arcData.hasPoint(g_upVertexId)) {
      nextPointId = points->InsertNextPoint(point);
      arcData.addPoint(g_upVertexId, nextPointId, scalarMax, false);
   } else {
      nextPointId = arcData.point_ids[g_upVertexId];
   }

   pointIds[1] = nextPointId;

   const vtkIdType nextCell = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
   arcData.fillArrayCell(nextCell, arcId, ftcTree_[cc], triangulation_[cc], params_);

   return 0;
}

int ttkFTCTree::addDirectSkeletonArc(const idSuperArc arcId, const int cc, vtkPoints* points,
                                     vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   FTCTree_MT* tree = ftcTree_[cc].tree.getTree(GetTreeType());
   vtkDataArray* idMapper = connected_components_[cc]->GetPointData()->GetArray("VertexIdentifier");
   SuperArc* arc = tree->getSuperArc(arcId);
   float     point[3];
   vtkIdType pointIds[2];

   const idVertex downNodeId     = tree->getLowerNodeId(arc);
   const idVertex l_downVertexId = tree->getNode(downNodeId)->getVertexId();
   const idVertex g_downVertexId = idMapper->GetTuple1(l_downVertexId);
   triangulation_[cc]->getVertexPoint(l_downVertexId, point[0], point[1], point[2]);
   const double scalarMin = inputScalars_[cc]->GetTuple1(l_downVertexId);
   // Get or create first point of the arc
   if (!arcData.hasPoint(g_downVertexId)) {
      const vtkIdType nextPointId = points->InsertNextPoint(point);
      pointIds[0]                 = nextPointId;
      arcData.addPoint(g_downVertexId, nextPointId, scalarMin, false);
   } else {
      pointIds[0] = arcData.point_ids[g_downVertexId];
   }

   const idVertex upNodeId     = tree->getUpperNodeId(arc);
   const idVertex l_upVertexId = tree->getNode(upNodeId)->getVertexId();
   const idVertex g_upVertexId = idMapper->GetTuple1(l_upVertexId);
   triangulation_[cc]->getVertexPoint(l_upVertexId, point[0], point[1], point[2]);
   const double scalarMax = inputScalars_[cc]->GetTuple1(l_upVertexId);
   // Get or create last point of the arc
   if (!arcData.hasPoint(g_upVertexId)) {
      const vtkIdType nextPointId = points->InsertNextPoint(point);
      pointIds[1]                 = nextPointId;
      arcData.addPoint(g_upVertexId, nextPointId, scalarMax, false);
   } else {
      pointIds[1] = arcData.point_ids[g_upVertexId];
   }

   const vtkIdType nextCell = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
   arcData.fillArrayCell(nextCell, arcId, ftcTree_[cc], triangulation_[cc], params_);

   return 0;
}

int ttkFTCTree::addSampledSkeletonArc(const idSuperArc arcId, const int cc, vtkPoints* points,
                                      vtkUnstructuredGrid* skeletonArcs, ArcData& arcData)
{
   FTCTree_MT* tree = ftcTree_[cc].tree.getTree(GetTreeType());
   vtkDataArray* idMapper = connected_components_[cc]->GetPointData()->GetArray("VertexIdentifier");
   SuperArc* arc = tree->getSuperArc(arcId);
   float     point[3];
   vtkIdType pointIds[2];

   const idVertex downNodeId     = tree->getLowerNodeId(arc);
   const idVertex l_downVertexId = tree->getNode(downNodeId)->getVertexId();
   const idVertex g_downVertexId = idMapper->GetTuple1(l_downVertexId);
   triangulation_[cc]->getVertexPoint(l_downVertexId, point[0], point[1], point[2]);
   const double scalarMin = inputScalars_[cc]->GetTuple1(l_downVertexId);

   // Get or create first point of the arc
   vtkIdType nextPointId;
   if (!arcData.hasPoint(g_downVertexId)) {
      nextPointId = points->InsertNextPoint(point);
      arcData.addPoint(g_downVertexId, nextPointId, scalarMin, false);
   } else {
      nextPointId = arcData.point_ids[g_downVertexId];
   }

   pointIds[0] = nextPointId;

   const idVertex upNodeId     = tree->getUpperNodeId(arc);
   const idVertex l_upVertexId = tree->getNode(upNodeId)->getVertexId();
   const idVertex g_upVertexId = idMapper->GetTuple1(l_upVertexId);
   triangulation_[cc]->getVertexPoint(l_upVertexId, point[0], point[1], point[2]);
   const double scalarMax = inputScalars_[cc]->GetTuple1(l_upVertexId);

   const double delta       = (scalarMax - scalarMin) / (params_.samplingLvl + 1);
   double       scalarLimit = scalarMin + delta;
   double       scalarAvg   = 0;

   // Get or create last point of the arc
   if (!arcData.hasPoint(g_upVertexId)) {
      nextPointId = points->InsertNextPoint(point);
      arcData.addPoint(g_upVertexId, nextPointId, scalarMax, false);
   } else {
      nextPointId = arcData.point_ids[g_upVertexId];
   }

   int       c = 0;
   float     sum[3]{0, 0, 0};
   for (const idVertex vertexId : *arc) {
      triangulation_[cc]->getVertexPoint(vertexId, point[0], point[1], point[2]);
      const double scalarVertex = inputScalars_[cc]->GetTuple1(vertexId);

      if (scalarVertex < scalarLimit) {
         sum[0] += point[0];
         sum[1] += point[1];
         sum[2] += point[2];
         scalarAvg += scalarVertex;
         ++c;
      } else {
         if (c) {
            sum[0] /= c;
            sum[1] /= c;
            sum[2] /= c;
            scalarAvg /= c;

            pointIds[1] = points->InsertNextPoint(sum);
            arcData.setPoint(pointIds[1], scalarAvg, true);
            const vtkIdType nextCell = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
            arcData.fillArrayCell(nextCell, arcId, ftcTree_[cc], triangulation_[cc], params_);

            pointIds[0] = pointIds[1];
         }

         scalarLimit += delta;
         sum[0]    = 0;
         sum[1]    = 0;
         sum[2]    = 0;
         scalarAvg = 0;
         c         = 0;
      }
   }

   // The up id
   pointIds[1] = nextPointId;

   const vtkIdType nextCell = skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
   arcData.fillArrayCell(nextCell, arcId, ftcTree_[cc], triangulation_[cc], params_);

   return 0;
}

int ttkFTCTree::doIt(vector<vtkDataSet*>& inputs, vector<vtkDataSet*>& outputs)
{
   Memory m;

   vtkUnstructuredGrid* outputSkeletonNodes = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
   vtkUnstructuredGrid* outputSkeletonArcs  = vtkUnstructuredGrid::SafeDownCast(outputs[1]);
   vtkDataSet*          outputSegmentation  = outputs[2];

   if(inputs[0]->IsA("vtkUnstructuredGrid")){
      // This data set may have several connected components,
      // we need to apply the FTC Tree for each one of these components
      // We then reconstruct the global tree using an offest mecanism
      vtkSmartPointer<vtkUnstructuredGrid> input = vtkSmartPointer<vtkUnstructuredGrid>::New();
      input->ShallowCopy(inputs[0]);
      identify(input);

      vtkSmartPointer<vtkConnectivityFilter> connectivity =
          vtkSmartPointer<vtkConnectivityFilter>::New();
      connectivity->SetInputData(input);
      connectivity->SetExtractionModeToAllRegions();
      connectivity->ColorRegionsOn();
      connectivity->Update();

      nbCC_ = connectivity->GetOutput()->GetCellData()->GetArray("RegionId")->GetRange()[1] + 1;
      connected_components_.resize(nbCC_);

      if (nbCC_ > 1) {
         // Warning, in case of several connected components, the ids seen by
         // the base code will not be consistent with those of the original
         // mesh
         for (int cc = 0; cc < nbCC_; cc++) {
            vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
            threshold->SetInputConnection(connectivity->GetOutputPort());
            threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
            threshold->ThresholdBetween(cc, cc);
            threshold->Update();
            connected_components_[cc] = ttkUnstructuredGrid::New();
            connected_components_[cc]->ShallowCopy(threshold->GetOutput());
         }
      } else {
         connected_components_[0] = ttkUnstructuredGrid::New();
         connected_components_[0]->ShallowCopy(input);
      }
   } else {
      nbCC_ = 1;
      connected_components_.resize(nbCC_);
      connected_components_[0] = ttkImageData::New();
      connected_components_[0]->ShallowCopy(inputs[0]);
      identify(connected_components_[0]);
   }

   // now proceed for each triangulation obtained.

   if (setupTriangulation()) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTCTree] Error : wrong triangulation." << endl;
      return -1;
#endif
   }

   // Fill the vector of scalar/offset, cut the array in pieces if needed
    getScalars();
    getOffsets();

   if(debugLevel_) {
       cout << "Launch on field : " << ScalarField << endl;
   }

   ftc::idNode acc_nbNodes = 0;

   // Build tree
   for (int cc = 0; cc < nbCC_; cc++) {
      ftcTree_[cc].tree.setScalars(inputScalars_[cc]->GetVoidPointer(0));
      ftcTree_[cc].tree.setVertexSoSoffsets(offsets_[cc].data());
      ftcTree_[cc].tree.setTreeType(GetTreeType());
      ftcTree_[cc].tree.setSegmentation(GetWithSegmentation());
      ftcTree_[cc].tree.setNormalizeIds(GetWithNormalize());

      switch (inputScalars_[cc]->GetDataType()) {
        vtkTemplateMacro({ if(ftcTree_[cc].tree.build<VTK_TT>()) return -1; });
      }

      ftcTree_[cc].offset = acc_nbNodes;
      acc_nbNodes += ftcTree_[cc].tree.getTree(GetTreeType())->getNumberOfNodes();
   }

   UpdateProgress(0.50);

   // Construct output
   if (getSkeletonNodes(outputSkeletonNodes)) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTCTree] Error : wrong properties on skeleton nodes." << endl;
      return -7;
#endif
   }

   if (getSkeletonArcs(outputSkeletonArcs)) {
#ifndef TTK_ENABLE_KAMIKAZE
      cerr << "[ttkFTCTree] Error : wrong properties on skeleton arcs." << endl;
      return -8;
#endif
   }

   if (GetWithSegmentation()) {
      outputSegmentation->ShallowCopy(inputs[0]);
      if (getSegmentation(outputSegmentation)) {
#ifndef TTK_ENABLE_KAMIKAZE
         cerr << "[ttkFTCTree] Error : wrong properties on segmentation." << endl;
         return -9;
#endif
      }
   }

   UpdateProgress(1);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   printCSVStats();
#endif

   {
      stringstream msg;
      msg << "[ttkFTCTree] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 0;
}

int ttkFTCTree::getOffsets()
{
   offsets_.resize(nbCC_);
   for (int cc = 0; cc < nbCC_; cc++) {
      vtkDataArray* inputOffsets;
      if (OffsetFieldId != -1) {
         inputOffsets = connected_components_[cc]->GetPointData()->GetArray(OffsetFieldId);
         if (inputOffsets) {
            InputOffsetScalarFieldName = inputOffsets->GetName();
            UseInputOffsetScalarField  = true;
         }
      }

      const idVertex numberOfVertices = connected_components_[cc]->GetNumberOfPoints();

      if (UseInputOffsetScalarField and InputOffsetScalarFieldName.length()) {
         inputOffsets =
             connected_components_[cc]->GetPointData()->GetArray(InputOffsetScalarFieldName.data());
         offsets_[cc].resize(numberOfVertices);
         for (int i = 0; i < numberOfVertices; i++) {
            offsets_[cc][i] = inputOffsets->GetTuple1(i);
         }
      } else {
         if (hasUpdatedMesh_ and offsets_[cc].size()) {
            // don't keep an out-dated offset array
            offsets_[cc].clear();
         }

         if (offsets_[cc].empty()) {
            offsets_[cc].resize(numberOfVertices);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) \
    schedule(static, max(1, numberOfVertices / threadNumber_))
#endif
            for (int i = 0; i < numberOfVertices; i++) {
               offsets_[cc][i] = i;
            }
         }
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if (offsets_[cc].empty()) {
         cerr << "[ttkFTCTree] Error : wrong input offset scalar field for " << cc << endl;
         return -1;
      }
#endif
   }

   return 0;
}

int ttkFTCTree::getScalars()
{
   inputScalars_.resize(nbCC_);
   for (int cc = 0; cc < nbCC_; cc++) {
      vtkPointData* pointData = connected_components_[cc]->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
      if (!pointData) {
         cerr << "[ttkFTCTree] Error : input cc " << cc << " has no point data." << endl;
         return -1;
      }
#endif

      if (ScalarField.length()) {
         inputScalars_[cc] = pointData->GetArray(ScalarField.data());
      } else {
         inputScalars_[cc] = pointData->GetArray(ScalarFieldId);
         if (inputScalars_[cc])
            ScalarField = inputScalars_[cc]->GetName();
      }

#ifndef TTK_ENABLE_KAMIKAZE
      if (!inputScalars_[cc]) {
         cerr << "[ttkFTCTree] Error : input scalar " << cc << " field pointer is null." << endl;
         return -3;
      }
#endif
   }

   return 0;
}

int ttkFTCTree::getSegmentation(vtkDataSet* outputSegmentation)
{
   VertData vertData;
   vertData.init(ftcTree_, params_);

   for (int cc = 0; cc < nbCC_; cc++) {
      FTCTree_MT*   tree = ftcTree_[cc].tree.getTree(GetTreeType());
      vtkDataArray* idMapper =
          connected_components_[cc]->GetPointData()->GetArray("VertexIdentifier");
      const idSuperArc numberOfSuperArcs = tree->getNumberOfSuperArcs();
      // #pragma omp for
      for (idSuperArc arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
         vertData.fillArrayPoint(arcId, ftcTree_[cc], triangulation_[cc], idMapper, params_);
      }
   }

   vtkPointData* pointData = outputSegmentation->GetPointData();
   vertData.addArray(pointData, params_);

   return 0;
}

int ttkFTCTree::getSkeletonArcs(vtkUnstructuredGrid* outputSkeletonArcs)
{
   vtkSmartPointer<vtkUnstructuredGrid> skeletonArcs = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints>           points       = vtkSmartPointer<vtkPoints>::New();

   ArcData arcData;
   arcData.init(ftcTree_, params_);

   const int       samplingLevel = params_.samplingLvl;
   for (int cc = 0; cc < nbCC_; cc++) {
      FTCTree_MT* tree = ftcTree_[cc].tree.getTree(GetTreeType());

      const idVertex numberOfSuperArcs = tree->getNumberOfSuperArcs();
#ifndef TTK_ENABLE_KAMIKAZE
      if (!numberOfSuperArcs) {
         cerr << "[ttkFTCTree] Error : tree has no super arcs." << endl;
         return -2;
      }
#endif

      for (idVertex arcId = 0; arcId < numberOfSuperArcs; ++arcId) {
         const int numberOfRegularNodes = tree->getArcSize(arcId);
         if (numberOfRegularNodes > 0 and samplingLevel > 0) {
            addSampledSkeletonArc(arcId, cc, points, skeletonArcs, arcData);
         } else if (samplingLevel == -1) {
            addCompleteSkeletonArc(arcId, cc, points, skeletonArcs, arcData);
         } else {
            addDirectSkeletonArc(arcId, cc, points, skeletonArcs, arcData);
         }
      }
   }

   skeletonArcs->SetPoints(points);
   arcData.addArray(skeletonArcs, params_);
   outputSkeletonArcs->ShallowCopy(skeletonArcs);

  // const idVertex p_size = points->GetNumberOfPoints();
  // const idVertex s_size = tree->getNumberOfVertices();
  // cout << "arcs points " << p_size << endl;
  // cout << "scal points " << s_size << endl;
  // cout << "nb arcs     " << tree->getNumberOfSuperArcs()<< endl;
  // if(p_size != s_size){
  //    exit(3);
  // }

  return 0;
}

int ttkFTCTree::getSkeletonNodes(vtkUnstructuredGrid* outputSkeletonNodes)
{
   vtkSmartPointer<vtkUnstructuredGrid> skeletonNodes = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints>           points        = vtkSmartPointer<vtkPoints>::New();

   NodeData nodeData;
   nodeData.init(ftcTree_, params_);

   for (int cc = 0; cc < nbCC_; cc++) {
      FTCTree_MT* tree = ftcTree_[cc].tree.getTree(GetTreeType());
      vtkDataArray* idMapper =
          connected_components_[cc]->GetPointData()->GetArray("VertexIdentifier");

      const idNode numberOfNodes = tree->getNumberOfNodes();
#ifndef TTK_ENABLE_KAMIKAZE
      if (!numberOfNodes) {
         cerr << "[ttkFTCTree] Error : tree has no nodes." << endl;
         return -2;
      }
#endif

      for (idNode nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
         const Node* node = tree->getNode(nodeId);
#ifndef TTK_ENABLE_KAMIKAZE
         if (!node) {
            cerr << "[ttkFTCTree] Error : node " << nodeId << " is null." << endl;
            return -7;
         }
#endif
         const idVertex local_vertId  = node->getVertexId();
         float          point[3];
         triangulation_[cc]->getVertexPoint(local_vertId, point[0], point[1], point[2]);
         const vtkIdType nextPoint = points->InsertNextPoint(point);
         nodeData.fillArrayPoint(nextPoint, nodeId, ftcTree_[cc], idMapper, triangulation_[cc],
                                 params_);
      }
   }

   skeletonNodes->SetPoints(points);
   vtkPointData* pointData = skeletonNodes->GetPointData();
   nodeData.addArray(pointData, params_);
   outputSkeletonNodes->ShallowCopy(skeletonNodes);

   return 0;
}

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
void ttkFTCTree::printCSVStats()
{
   switch (GetTreeType()) {
      case ftc::TreeType::Join:
         cout << "JT" << endl;
         printCSVTree(ftcTree_[0].tree.getTree(ftc::TreeType::Join));
         break;
      case ftc::TreeType::Split:
         cout << "ST" << endl;
         printCSVTree(ftcTree_[0].tree.getTree(ftc::TreeType::Split));
         break;
      default:
         printCSVExtr(ftcTree_[0].tree.getTree(ftc::TreeType::Contour));
         cout << "JT" << endl;
         printCSVTree(ftcTree_[0].tree.getTree(ftc::TreeType::Join), false);
         cout << "ST" << endl;
         printCSVTree(ftcTree_[0].tree.getTree(ftc::TreeType::Split), false);
         break;
   }
}

void ttkFTCTree::printCSVExtr(const ftc::FTCTree_MT* const tree) const
{
   vector<tuple<float, idVertex>> stats;

   stats.clear();
   idVertex nbExtrStats = tree->getNbExtrStats();
   stats.resize(nbExtrStats);

   float offsetTime = tree->getOffsetTime();
   float startTime  = tree->getExtrStats(0).begin;

   for (idVertex i = 0; i < nbExtrStats; i++) {
      if (tree->getExtrStats(i).end > 0) {
         stats[i] = make_tuple(tree->getExtrStats(i).end, 0);
      } else {
         stats.pop_back();
      }
      if (tree->getExtrStats(i).begin < startTime) {
         startTime = tree->getExtrStats(i).begin;
      }
   }

   // note: really bad perfs for the sort, but here we don't mind.
   stats.emplace_back(startTime, 0);

   sort(begin(stats), end(stats),
        [](tuple<float, int> a, tuple<float, int> b) { return get<0>(a) < get<0>(b); });

   nbExtrStats = stats.size();

   cout << "extr_time; ";
   for (idVertex i = 0; i < nbExtrStats; i++) {
      cout << get<0>(stats[i]) - offsetTime << "; ";
   }
   cout << endl << "extr_tasks; ";
   for (idVertex i = 0; i < nbExtrStats; i++) {
      cout << (nbExtrStats - i - (1 + (i == 0))) << "; ";
   }
   cout << endl;
}

void ttkFTCTree::printCSVTree(const ftc::FTCTree_MT* const tree, const bool printExtract) const
{
   if (printExtract) {
      printCSVExtr(tree);
   }

   vector<tuple<float, idVertex>> stats;
   float offsetTime = tree->getOffsetTime();
   float startTime;

   stats.clear();

   idSuperArc nbArcGrowthStats = tree->getNbArcGrowthStats();
   stats.resize(nbArcGrowthStats);

   startTime = 1000000000;

   for (idSuperArc i = 0; i < nbArcGrowthStats; i++) {
      if (tree->getArcGrowthStats(i).end > 0) {
         stats[i] =
             make_tuple(tree->getArcGrowthStats(i).end, tree->getArcGrowthStats(i).remain);
      } else {
         stats.pop_back();
      }

      if (tree->getArcGrowthStats(i).begin != -1 && tree->getArcGrowthStats(i).begin < startTime) {
         startTime = tree->getArcGrowthStats(i).begin;
      }
   }


   // note: really bad perfs for the sort, but here we don't mind.
   stats.emplace_back(startTime, stats.size() - 1);

   sort(begin(stats), end(stats), [](tuple<float, int> a, tuple<float, int> b) {
      return get<0>(a) < get<0>(b);
   });

   nbArcGrowthStats = stats.size();

   // print
   cout << "arc_time; ";
   for (idSuperArc i = 0; i < nbArcGrowthStats; i++) {
      cout << get<0>(stats[i]) - offsetTime << "; ";
   }
   cout << endl << "arc_tasks; ";
   for (idSuperArc i = 0; i < nbArcGrowthStats; i++) {
      cout << get<1>(stats[i]) << "; ";
   }
   cout << endl;

   stats.clear();

   idVertex nbArcTrunkStats = tree->getNbArcTrunkStats();
   stats.resize(nbArcTrunkStats);

   startTime = tree->getArcTrunkStats(0).begin;

   for (idVertex i = 0; i < nbArcTrunkStats; i++) {
      if (tree->getArcTrunkStats(i).end > 0) {
         stats[i] = make_tuple(tree->getArcTrunkStats(i).end, 0);
      } else {
         stats.pop_back();
      }
      if (tree->getArcTrunkStats(i).begin < startTime) {
         startTime = tree->getArcTrunkStats(i).begin;
      }
   }

   // note: really bad perfs for the sort, but here we don't mind.
   stats.emplace_back(startTime, 0);

   sort(begin(stats), end(stats), [](tuple<float, int> a, tuple<float, int> b) {
      return get<0>(a) < get<0>(b);
   });

   nbArcTrunkStats = stats.size();

   cout << "trunk_arc_time; ";
   for (idVertex i = 0; i < nbArcTrunkStats; i++) {
      cout << get<0>(stats[i]) - offsetTime << "; ";
   }
   cout << endl << "trunk_arc_tasks; ";
   for (idVertex i = 0; i < nbArcTrunkStats; i++) {
      cout << (nbArcTrunkStats - i - (1 + (i == 0))) << "; ";
   }
   cout << endl;

   stats.clear();

   idVertex nbTrunkStats = tree->getNbTrunkStats();
   stats.resize(nbTrunkStats);

   startTime = tree->getTrunkStats(0).begin;

   for (idVertex i = 0; i < nbTrunkStats; i++) {
      if (tree->getTrunkStats(i).end > 0) {
         stats[i] = make_tuple(tree->getTrunkStats(i).end, 0);
      } else {
         stats.pop_back();
      }
      if (tree->getTrunkStats(i).begin < startTime) {
         startTime = tree->getTrunkStats(i).begin;
      }
   }

   // note: really bad perfs for the sort, but here we don't mind.
   stats.emplace_back(startTime, 0);

   sort(begin(stats), end(stats), [](tuple<float, int> a, tuple<float, int> b) {
      return get<0>(a) < get<0>(b);
   });

   nbTrunkStats = stats.size();

   cout << "trunk_time; ";
   for (idVertex i = 0; i < nbTrunkStats; i++) {
      cout << get<0>(stats[i]) - offsetTime << "; ";
   }
   cout << endl << "trunk_tasks; ";
   for (idVertex i = 0; i < nbTrunkStats; i++) {
      cout << (nbTrunkStats - i - (1 + (i == 0))) << "; ";
   }
   cout << endl;

   stats.clear();

  idVertex nbSegmStats = tree->getNbSegmStats();
  stats.resize(nbSegmStats);

  startTime = tree->getSegmStats(0).begin;

   for (idVertex i = 0; i < nbSegmStats; i++) {
      if (tree->getSegmStats(i).end > 0) {
         stats[i] = make_tuple(tree->getSegmStats(i).end, 0);
      } else {
         stats.pop_back();
      }
      if (tree->getSegmStats(i).begin < startTime) {
         startTime = tree->getSegmStats(i).begin;
      }
   }

   // note: really bad perfs for the sort, but here we don't mind.
   stats.emplace_back(startTime, 0);

   sort(begin(stats), end(stats),
        [](tuple<float, int> a, tuple<float, int> b) { return get<0>(a) < get<0>(b); });

   nbSegmStats = stats.size();

   cout << "segm_time; ";
   for (idVertex i = 0; i < nbSegmStats; i++) {
      cout << get<0>(stats[i]) - offsetTime << "; ";
   }
   cout << endl << "segm_tasks; ";
   for (idVertex i = 0; i < nbSegmStats; i++) {
      cout << (nbSegmStats - i - (1 + (i == 0))) << "; ";
   }
   cout << endl;

}
#endif

int ttkFTCTree::setupTriangulation()
{
   triangulation_.resize(nbCC_);
   ftcTree_.resize(nbCC_);

   for (int cc = 0; cc < nbCC_; cc++) {
      triangulation_[cc] = ttkTriangulation::getTriangulation(connected_components_[cc]);
#ifndef TTK_ENABLE_KAMIKAZE
         if (!triangulation_[cc]) {
            cerr << "[ttkFTCTree] Error : ttkTriangulation::getTriangulation() is null." << endl;
            return -1;
      }
#endif

      triangulation_[cc]->setWrapper(this);
      ftcTree_[cc].tree.setDebugLevel(debugLevel_);
      ftcTree_[cc].tree.setThreadNumber(threadNumber_);
      ftcTree_[cc].tree.setupTriangulation(triangulation_[cc]);

      hasUpdatedMesh_ = ttkTriangulation::hasChangedConnectivity(triangulation_[cc], connected_components_[cc], this);

#ifndef TTK_ENABLE_KAMIKAZE
      if (triangulation_[cc]->isEmpty()) {
         cerr << "[ttkFTCTree] Error : ttkTriangulation on connected component" << cc << " allocation problem." << endl;
         return -1;
      }
#endif
   }
   return 0;
}

// protected

ttkFTCTree::ttkFTCTree()
    : ScalarField{},
      UseInputOffsetScalarField{},
      InputOffsetScalarFieldName{},
      ScalarFieldId{},
      OffsetFieldId{-1},
      params_{},
      triangulation_{},
      inputScalars_{},
      offsets_{},
      hasUpdatedMesh_{}
{
   SetNumberOfInputPorts(1);
   SetNumberOfOutputPorts(3);
}

ttkFTCTree::~ttkFTCTree()
{
}

void ttkFTCTree::identify(vtkDataSet* ds) const
{
   vtkSmartPointer<vtkIntArray> identifiers=vtkSmartPointer<vtkIntArray>::New();
   const vtkIdType nbPoints = ds->GetNumberOfPoints();
   identifiers->SetName("VertexIdentifier");
   identifiers->SetNumberOfComponents(1);
   identifiers->SetNumberOfTuples(nbPoints);

   for (int i = 0; i < nbPoints; i++) {
      identifiers->SetTuple1(i, i);
   }

   ds->GetPointData()->AddArray(identifiers);
}
