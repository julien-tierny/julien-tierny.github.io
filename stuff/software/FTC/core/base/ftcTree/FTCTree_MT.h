/// \ingroup base
/// \class ttk::ftc::FTCTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the
/// sublevel set tree of scalar data and more
/// (data segmentation
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#ifndef FTCTREE_MT_H
#define FTCTREE_MT_H

#include "AtomicUF.h"
#include "AtomicVector.h"
#include "DataTypes.h"
#include "FTCCommon.h"
#include "Node.h"
#include "Scalars.h"
#include "SuperArc.h"

#include <functional>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk
{
namespace ftc
{
   using UF = AtomicUF *;

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   static DebugTimer _launchGlobalTime;
#endif

   // ( 1 per tree )
   struct MT_Data {
      TreeType treeType;

      // components : tree / nodes / extrema
      AtomicVector<SuperArc> *superArcs;
      AtomicVector<Node> *    nodes;
      AtomicVector<idNode> *  roots;
      AtomicVector<idNode> *  leaves;

      // vertex 2 node / superarc
      vector<idCorresp> *             vert2tree;
      vector<idVertex> *              visitOrder;
      vector<list<vector<idVertex>>> *trunkSegments;
      vector<idVertex> *              trunkVerts;

      // Track informations
      vector<UF> *ufs, *propagationUfs;
      AtomicVector<Propagation> *propagations;
      // valences
      vector<valence> *valences;
      // opened nodes
      vector<valence> *openNodeOnVert;

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
      vector<ActiveTask> *extrStats;
      vector<ActiveTask> *arcGrowthStats;
      vector<ActiveTask> *arcTrunkStats;
      vector<ActiveTask> *trunkStats;
      vector<ActiveTask> *segmStats;
      float               startTime;
#endif

      // current nb of tasks
      idNode activeTasks;

      // Segmentation, stay empty for Contour tree as
      // they are created by Merge Tree
      Segments segments_;

#ifndef TTK_DISABLE_FTC_PRIORITY
      bool isPrior;
#endif
   };

   // also one per tree
   struct Comparison {
      VertCompFN vertLower, vertHigher, vertEqLower;
   };


   class FTCTree_MT : virtual public Debug
   {
     public:
      // global
      Params *const  params_;
      Triangulation *mesh_;
      Scalars<float> *const scalars_;

      // local
      MT_Data    mt_data_;
      Comparison comp_;

      using sortedVertIt = vector<idVertex>::iterator;

     public:

      // -----------
      // CONSTRUCT
      // -----------

      // Tree with global data and partition number
      FTCTree_MT(Params *const params, Triangulation *mesh, Scalars<float> *const scalars, TreeType type);

      virtual ~FTCTree_MT();

      // --------------------
      // Init
      // --------------------

      void initComp(void)
      {
         if (isST()) {
            comp_.vertEqLower = [this](const idVertex a, const idVertex b) -> bool const {
               return this->scalars_->isEqHigher(a, b);
            };
            comp_.vertLower = [this](const idVertex a, const idVertex b) -> bool const {
               return this->scalars_->isHigher(a, b);
            };
            comp_.vertHigher = [this](const idVertex a, const idVertex b) -> bool const {
               return this->scalars_->isLower(a, b);
            };
         } else {
            comp_.vertEqLower = [this](const idVertex a, const idVertex b) -> bool const {
               return this->scalars_->isEqLower(a, b);
            };
            comp_.vertLower = [this](const idVertex a, const idVertex b) -> bool const {
               return this->scalars_->isLower(a, b);
            };
            comp_.vertHigher = [this](const idVertex a, const idVertex b) -> bool const {
               return this->scalars_->isHigher(a, b);
            };
         }
      }

      bool compLower(const idVertex a, const idVertex b)
      {
         return comp_.vertLower(a, b);
      }

      /// \brief if sortedVertices_ is null, define and fill it
      /// Also fill the mirror vector
      template <typename ScalarType>
      void sortInput(void);

      // -------------------
      // Process
      // -------------------

      /// \brief clear local data for new computation
      virtual void alloc(const bool allocScalars = true)
      {
         const auto size = scalars_->getSize();

         if (allocScalars) {
            scalars_->alloc();
         }

         createAtomicVector<SuperArc>(mt_data_.superArcs);
         mt_data_.superArcs->reserve(size / 2);

         // Stats alloc

         createAtomicVector<Node>(mt_data_.nodes);
         mt_data_.nodes->reserve(size / 2);

         createAtomicVector<idNode>(mt_data_.roots);
         mt_data_.roots->reserve(10);

         createAtomicVector<idNode>(mt_data_.leaves);
         mt_data_.leaves->reserve(size / 3);

         // Known size

         createVector<idCorresp>(mt_data_.vert2tree);
         mt_data_.vert2tree->resize(size);

         createVector<list<vector<idVertex>>>(mt_data_.trunkSegments);

         createVector<idVertex>(mt_data_.trunkVerts);
         mt_data_.trunkVerts->reserve(size / 10);

         createVector<idVertex>(mt_data_.visitOrder);
         mt_data_.visitOrder->resize(size);

         createVector<UF>(mt_data_.ufs);
         mt_data_.ufs->resize(size);

         createVector<UF>(mt_data_.propagationUfs);
         mt_data_.propagationUfs->resize(size);

         createVector<valence>(mt_data_.valences);
         mt_data_.valences->resize(size);

         createVector<valence>(mt_data_.openNodeOnVert);
         mt_data_.openNodeOnVert->resize(size);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
         createVector<ActiveTask>(mt_data_.extrStats);
         mt_data_.extrStats->resize(size / 2);
         createVector<ActiveTask>(mt_data_.arcGrowthStats);
         mt_data_.arcGrowthStats->resize(size / 2);
         createVector<ActiveTask>(mt_data_.arcTrunkStats);
         mt_data_.arcTrunkStats->reserve(size / 2);
         createVector<ActiveTask>(mt_data_.trunkStats);
         mt_data_.trunkStats->reserve(threadNumber_ * 50);
         createVector<ActiveTask>(mt_data_.segmStats);
         mt_data_.segmStats->reserve(size / 2);
#endif

         // allocVectPropagations(size / 10);

         mt_data_.segments_.clear();
      }

      virtual void init(const bool initScalars = true)
      {
         if (initScalars) {
            scalars_->init();
         }

         initVector<idCorresp>(mt_data_.vert2tree, nullCorresp);
         initVector<idVertex>(mt_data_.visitOrder, nullVertex);
         initVector<UF>(mt_data_.ufs, nullptr);
         initVector<UF>(mt_data_.propagationUfs, nullptr);
         initVector<valence>(mt_data_.valences, -1);
         initVector<valence>(mt_data_.openNodeOnVert, 0);
      }

      void allocVectPropagations(const idVertex nbLeaves)
      {
         if (!mt_data_.propagations) {
            mt_data_.propagations = new AtomicVector<Propagation>(nbLeaves, comp_.vertHigher);
         }
         mt_data_.propagations->clear();
         mt_data_.propagations->reserve(nbLeaves);
      }

      /// \brief Compute the merge
      virtual void build(const bool computeCT);

      virtual int leafSearch();

      void leafGrowth();

      void arcGrowth(const idVertex startVert, const idVertex orig);

      tuple<bool, bool> propage(Propagation &currentState, UF curUF);

      idVertex calculValence(const idVertex idVert);

      void closeAndMergeOnSaddle(const idVertex saddleVert);

      void closePendingTrunkNodes(idVertex saddleVert);

      void closeArcsUF(idNode closeNode, UF uf);

      // This trunk deal with list of vertices. This is because the arc growth  step create a node
      // when all incoming tasks have reached it. So the nodes in
      // the trunk are not created yet
      idVertex trunk(const bool ct);

      virtual idVertex trunkSegmentation(const vector<idVertex> &pendingNodesVerts,
                                         const idVertex begin,
                                         const idVertex stop);

      // fill treedata_.trunkSegments
      idVertex trunkCTSegmentation(const vector<idVertex> &pendingNodesVerts,
                                   const idVertex begin,
                                   const idVertex stop);

      // segmentation

      /// \brief use vert2tree to compute the segmentation of the fresh builded merge tree.
      void buildSegmentation(const idVertex begin = nullVertex, const idVertex stop = nullVertex);

      // Create the segmentation of all arcs by operating the pending operations
      void finalizeSegmentation(void);

      void normalizeIds();

      // -------------
      // ACCESSOR
      // ------------

      // Tree info for wrapper

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME

      void setStartTime(void) {
         mt_data_.startTime = _launchGlobalTime.getElapsedTime();
      }

      float getOffsetTime(void) const {
         return mt_data_.startTime;
      }

      const ActiveTask& getExtrStats(const idSuperArc taskId) const
      {
         return (*mt_data_.extrStats)[taskId];
      }

      idVertex getNbExtrStats(void) const
      {
         return mt_data_.extrStats->size();
      }

      const ActiveTask& getArcGrowthStats(const idSuperArc taskId) const
      {
         return (*mt_data_.arcGrowthStats)[taskId];
      }

      idSuperArc getNbArcGrowthStats(void) const
      {
         return mt_data_.arcGrowthStats->size();
      }

      const ActiveTask& getArcTrunkStats(const idVertex taskId) const
      {
         return (*mt_data_.arcTrunkStats)[taskId];
      }

      idVertex getNbArcTrunkStats(void) const
      {
         return mt_data_.arcTrunkStats->size();
      }

      const ActiveTask& getTrunkStats(const idVertex taskId) const
      {
         return (*mt_data_.trunkStats)[taskId];
      }

      idVertex getNbTrunkStats(void) const
      {
         return mt_data_.trunkStats->size();
      }

      const ActiveTask& getSegmStats(const idVertex taskId) const
      {
         return (*mt_data_.segmStats)[taskId];
      }

      idVertex getNbSegmStats(void) const
      {
         return mt_data_.segmStats->size();
      }
#endif

      inline idVertex getArcSize(const idSuperArc arcId)
      {
          return getSuperArc(arcId)->size();
      }

      inline bool isJT(void) const
      {
         return mt_data_.treeType == TreeType::Join;
      }

      inline bool isST(void) const
      {
         return mt_data_.treeType == TreeType::Split;
      }

#ifndef TTK_DISABLE_FTC_PRIORITY
      void setPrior(void)
      {
         mt_data_.isPrior = true;
      }

      bool isPrior(void) const
      {
         return mt_data_.isPrior;
      }
#endif

      // global
      // called for the tree used by the wrapper (only).
      // On this implementation, the warpper communicate with ContourForest
      // A child class of this one.

      inline void setupTriangulation(Triangulation *m, const bool preproc = true)
      {
         mesh_ = m;
         if (mesh_ && preproc) {
            // propage through vertices (build)
            mesh_->preprocessVertexNeighbors();
         }
      }

      inline void setTreeType(const int local_treeType)
      {
         params_->treeType = static_cast<TreeType>(local_treeType);
      }

      inline void setSegmentation(const bool segm)
      {
          params_->segm = segm;
      }

      inline void setNormalizeIds(const bool normalize)
      {
          params_->normalize = normalize;
      }

      // scalar

      inline void setScalars(const void *local_scalars)
      {
         scalars_->setScalars((float*)local_scalars);
      }

      inline float getValue(const idVertex v)
      {
         return scalars_->getVal(v);
      }

      // offset

      inline void setVertexSoSoffsets(idVertex* sos)
      {
         scalars_->setOffsets(sos);
      }

      // arcs

      inline idSuperArc getNumberOfSuperArcs(void) const
      {
         return mt_data_.superArcs->size();
      }

      inline SuperArc *getSuperArc(idSuperArc i)
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)i >= mt_data_.superArcs->size()) {
            cout << "[Merge Tree] get superArc on bad id :" << i;
            cout << " / " << mt_data_.superArcs->size() << endl;
            return nullptr;
         }
#endif
         return &((*mt_data_.superArcs)[i]);
      }

      inline const SuperArc *getSuperArc(idSuperArc i) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)i >= mt_data_.superArcs->size()) {
            cout << "[Merge Tree] get superArc on bad id :" << i;
            cout << " / " << mt_data_.superArcs->size() << endl;
            return nullptr;
         }
#endif
         return &((*mt_data_.superArcs)[i]);
      }

      // nodes

      inline idNode getNumberOfNodes(void) const
      {
         return mt_data_.nodes->size();
      }

      inline Node *getNode(idNode nodeId)
      {
         return &((*mt_data_.nodes)[nodeId]);
      }

      inline void setValence(const idVertex v, const idVertex val)
      {
          (*mt_data_.valences)[v] = val;
      }

      // leaves / root

      inline idNode getNumberOfLeaves(void) const
      {
         return mt_data_.leaves->size();
      }

      inline AtomicVector<idNode> &getLeaves(void) const
      {
         // break encapsulation...
         return (*mt_data_.leaves);
      }

      inline idNode getLeave(const idNode id) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)id > (mt_data_.leaves->size())) {
            stringstream msg;
            msg << "[MergTree] getLeaves out of bounds : " << id << endl;
            err(msg.str(), fatalMsg);
            return (*mt_data_.leaves)[0];
         }
#endif
         return (*mt_data_.leaves)[id];
      }

      inline const vector<idNode> &getRoots(void) const
      {
         // break encapsulation...
         return (*mt_data_.roots);
      }

      // vertices

      inline idVertex getNumberOfVertices(void) const
      {
         return scalars_->getSize();
      }

      // vert2tree

      inline void setVert2Tree(decltype(mt_data_.vert2tree) const vect2tree)
      {
         mt_data_.vert2tree = vect2tree;
      }

      // --------------------
      // VERT 2 TREE Special functions
      // --------------------

      // test vertex correpondance

      inline bool isCorrespondingArc(const idVertex val) const
      {
         return !isCorrespondingNull(val) && (*mt_data_.vert2tree)[val] >= 0;
      }

      inline bool isCorrespondingNode(const idVertex val) const
      {
         return (*mt_data_.vert2tree)[val] < 0;
      }

      inline bool isCorrespondingNull(const idVertex val) const
      {
         return (*mt_data_.vert2tree)[val] == nullCorresp;
      }

      // Get vertex info

      inline idNode getCorrespondingNodeId(const idVertex val) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!isCorrespondingNode(val)) {
            stringstream debug;
            debug << "[FTCTree_MT] : getCorrespondingNode, ";
            debug << "Vertex :" << val << " is not a node :";
            debug << (*mt_data_.vert2tree)[val] << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return corr2idNode(val);
      }

      inline idSuperArc getCorrespondingSuperArcId(const idVertex val) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!isCorrespondingArc(val)) {
            stringstream debug;
            debug << "[FTCTree_MT] : getCorrespondingSuperArcId, ";
            debug << "Vertex :" << val << " is not on an arc :";
            debug << (*mt_data_.vert2tree)[val] << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return (*mt_data_.vert2tree)[val];
      }

      // Get corresponding elemnt

      inline SuperArc *vertex2SuperArc(const idVertex vert)
      {
         return &((*mt_data_.superArcs)[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const idVertex vert)
      {
         return &((*mt_data_.nodes)[getCorrespondingNodeId(vert)]);
      }

      // Update vertex info


      inline void updateCorrespondingArc(const idVertex vert, const idSuperArc arc)
      {
         (*mt_data_.vert2tree)[vert] = arc;
      }

      inline void updateCorrespondingNode(const idVertex vert, const idNode node)
      {
         (*mt_data_.vert2tree)[vert] = idNode2corr(node);
      }

      inline idCorresp idNode2corr(const idNode id) const
      {
         // transform idNode to special value for the array : -idNode -1
         return -(idCorresp)(id + 1);
      }

      inline idNode corr2idNode(const idCorresp corr) const
      {
         return -(idNode)((*mt_data_.vert2tree)[corr] + 1);
      }

      // --------------------------------
      // Arcs and node manipulations
      // --------------------------------
      // SuperArcs

      idSuperArc openSuperArc(idNode downNodeId);

      idSuperArc makeSuperArc(idNode downNodeId, idNode upNodeId, const bool withNodes = true);

      void closeSuperArc(idSuperArc superArcId, idNode upNodeId);

      // Nodes

      vector<idNode> sortedNodes(const bool parallel = false);

      void sortLeaves(const bool parallel = false);

      idNode makeNode(idVertex vertexId);

      idNode makeNode(const Node *const n);

      idSuperArc insertNode(Node *node, const bool segm = true);

      // get node starting / ending this arc
      // orientation depends on Join/Split tree
      Node *getDownNode(const SuperArc *a);
      Node *getUpNode(const SuperArc *a);
      idNode getDownNodeId(const SuperArc *a);
      idNode getUpNodeId(const SuperArc *a);

      // get node above / below this arc
      // in term of scalar value
      Node *getLowerNode(const SuperArc *a);
      Node * getUpperNode(const SuperArc *a);
      idNode getLowerNodeId(const SuperArc *a);
      idNode getUpperNodeId(const SuperArc *a);

      idNode getParent(const idNode n)
      {
         return getSuperArc(getNode(n)->getUpSuperArcId(0))->getUpNodeId();
      }

      void delNode(idNode node);

      // ---------------------------
      // Operators : clone / move & print
      // ---------------------------

      FTCTree_MT *clone() const;

      void move(FTCTree_MT *mt);

      // Print
      string printArc(idSuperArc a);

      string printNode(idNode n);

      void printTree2(void);

      int printTime(DebugTimer &t, const string &s, idVertex nbScalars = -1,
                    const int debugLevel = 2) const;

      void printLeavesStats(void);

     protected:

      // -----
      // Tools
      // -----

      idNode getNodeFromVertInc(const vector<idVertex> &range,
                                const idVertex          v,
                                const idNode            last = 0) const;

      idNode getNodeFromVertDicho(const vector<idVertex> &range, const idVertex v) const;

      tuple<idVertex, idVertex> getBoundsFromVerts(const vector<idVertex> &nodes) const;

      idSuperArc upArcFromVert(const idVertex v)
      {
         return getNode(getCorrespondingNodeId(v))->getUpSuperArcId(0);
      }

      inline idVertex getChunkSize(const idVertex nbVerts = -1, const idVertex nbtasks = 100,
                                   idVertex minSize = 10000) const
      {
         const idVertex s = (nbVerts == -1) ? scalars_->getSize() : nbVerts;
#ifndef NDEBUG
         // Debug mode
         minSize = 1;
#endif
         return max(minSize, 1 + (s / (nbtasks * 1)));
      }

      // Carefull: these parameters are linked to those of the getChunkSize function.
      inline idVertex getChunkCount(const idVertex nbVerts = -1, const idVertex nbTasks = 100,
                                    const idVertex minSize = 10000) const
      {
         const idVertex s = (nbVerts == -1) ? scalars_->getSize() : nbVerts;
         return 1 + (s / getChunkSize(s, nbTasks, minSize));
      }

      void sortUpArcs(const idNode nid)
      {
         auto comp = [&](const idSuperArc a, const idSuperArc b) -> bool {
            return comp_.vertLower(getUpperNode(getSuperArc(a))->getVertexId(),
                                   getUpperNode(getSuperArc(b))->getVertexId());
         };

         getNode(nid)->sortUpArcs(comp);
      }

      void sortDownArcs(const idNode nid)
      {
         auto comp = [&](const idSuperArc a, const idSuperArc b) -> bool {
            return comp_.vertHigher(getUpperNode(getSuperArc(a))->getVertexId(),
                                    getUpperNode(getSuperArc(b))->getVertexId());
         };

         getNode(nid)->sortDownArcs(comp);
      }

      template <typename type>
      void createVector(vector<type> *&ptr)
      {
         if(!ptr)
            ptr = new vector<type>;
         ptr->clear();
      }

      template <typename type>
      void createAtomicVector(AtomicVector<type> *&ptr)
      {
         if(!ptr)
            ptr = new AtomicVector<type>;
         ptr->clear();
      }

      template <typename type>
      void initVector(vector<type> *&vect, const type val)
      {
         int s = vect->size();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(static, max(1, s / threadNumber_))
#endif
         for (int i = 0; i < s; i++) {
            (*vect)[i] = val;
         }
      }

   };

   ostream &operator<<(ostream &o, Node const &n);
   ostream &operator<<(ostream &o, SuperArc const &a);

}
}

#endif /* end of include guard: MERGETREE_H */
