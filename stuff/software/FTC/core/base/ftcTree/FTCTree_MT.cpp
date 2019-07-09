/// \ingroup base
/// \class ttk::ftc::FTCTree_MT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date Dec 2016.
///
///\brief TTK processing package that efficiently computes the
/// contour tree of scalar data and more
/// (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#include "FTCTree_MT.h"

#include <stack>

#ifdef TTK_DISABLE_FTC_PRIORITY
#define PRIOR(x)
#else
#define PRIOR(x) priority(x)
#endif

#ifdef __INTEL_COMPILER
#define HIGHER
#endif

#ifndef TTK_ENABLE_OPENMP
#define HIGHER
#endif

using namespace ftc;

FTCTree_MT::FTCTree_MT(Params *const         params,
                       Triangulation *       mesh,
                       Scalars<float> *const scalars,
                       TreeType              type)
    : params_(params), mesh_(mesh), scalars_(scalars)
{
   mt_data_.treeType = type;

   mt_data_.superArcs      = nullptr;
   mt_data_.nodes          = nullptr;
   mt_data_.roots          = nullptr;
   mt_data_.leaves         = nullptr;
   mt_data_.trunkVerts     = nullptr;
   mt_data_.vert2tree      = nullptr;
   mt_data_.trunkSegments  = nullptr;
   mt_data_.visitOrder     = nullptr;
   mt_data_.ufs            = nullptr;
   mt_data_.propagationUfs = nullptr;
   mt_data_.propagations   = nullptr;
   mt_data_.valences       = nullptr;
   mt_data_.openNodeOnVert = nullptr;

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   mt_data_.extrStats      = nullptr;
   mt_data_.arcGrowthStats = nullptr;
   mt_data_.arcTrunkStats  = nullptr;
   mt_data_.trunkStats     = nullptr;
   mt_data_.segmStats      = nullptr;
#endif

   mt_data_.activeTasks = 0;

#ifndef TTK_DISABLE_FTC_PRIORITY
   mt_data_.isPrior = false;
#endif
}

FTCTree_MT::~FTCTree_MT()
{
   // remove allocated union finds
   if (mt_data_.ufs) {
      sort(mt_data_.ufs->begin(), mt_data_.ufs->end());
      auto it = unique(mt_data_.ufs->begin(), mt_data_.ufs->end());
      mt_data_.ufs->resize(distance(mt_data_.ufs->begin(), it));
      for(auto* addr : *mt_data_.ufs) if (addr) delete addr;
   }

   // remove containers
   if (mt_data_.superArcs) {
      delete mt_data_.superArcs;
      mt_data_.superArcs = nullptr;
   }
   if (mt_data_.nodes) {
      delete mt_data_.nodes;
      mt_data_.nodes = nullptr;
   }
   if (mt_data_.roots) {
      delete mt_data_.roots;
      mt_data_.roots = nullptr;
   }
   if (mt_data_.leaves) {
      delete mt_data_.leaves;
      mt_data_.leaves = nullptr;
   }
   if (mt_data_.trunkVerts) {
      delete mt_data_.trunkVerts;
      mt_data_.trunkVerts = nullptr;
   }
   if (mt_data_.vert2tree) {
      delete mt_data_.vert2tree;
      mt_data_.vert2tree = nullptr;
   }
   if (mt_data_.trunkSegments) {
      delete mt_data_.trunkSegments;
      mt_data_.trunkSegments = nullptr;
   }
   if (mt_data_.visitOrder) {
      delete mt_data_.visitOrder;
      mt_data_.visitOrder = nullptr;
   }
   if (mt_data_.ufs) {
      delete mt_data_.ufs;
      mt_data_.ufs = nullptr;
   }
   if (mt_data_.propagationUfs) {
      delete mt_data_.propagationUfs;
      mt_data_.propagationUfs = nullptr;
   }
   if (mt_data_.propagations) {
      delete mt_data_.propagations;
      mt_data_.propagations = nullptr;
   }
   if (mt_data_.valences) {
      delete mt_data_.valences;
      mt_data_.valences = nullptr;
   }
   if (mt_data_.openNodeOnVert) {
      delete mt_data_.openNodeOnVert;
      mt_data_.openNodeOnVert = nullptr;
   }

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   if (mt_data_.extrStats) {
      delete mt_data_.extrStats;
      mt_data_.extrStats = nullptr;
   }
   if (mt_data_.arcGrowthStats) {
      delete mt_data_.arcGrowthStats;
      mt_data_.arcGrowthStats = nullptr;
   }
   if (mt_data_.trunkStats) {
      delete mt_data_.trunkStats;
      mt_data_.trunkStats = nullptr;
   }
   if (mt_data_.segmStats) {
      delete mt_data_.segmStats;
      mt_data_.segmStats = nullptr;
   }
#endif
}

void FTCTree_MT::arcGrowth(const idVertex startVert, const idVertex orig)
{
   // current task id / propag

   // local order (ignore non regular verts)
   idVertex localOrder = -1;
   UF startUF = (*mt_data_.ufs)[startVert]->find();
   // get or recover states
   Propagation *localPropagation;
   if (startUF->getNbPropagation()) {
      localPropagation = startUF->getCurrentPropagation();
   } else {
      const std::size_t currentPropId= mt_data_.propagations->getNext();
      localPropagation = &(*mt_data_.propagations)[currentPropId];
      localPropagation->setStartVert(startVert);
      startUF->addPropagation(localPropagation);
   }

   localPropagation->addNewVertex(startVert);

   // ARC OPENING
   const idNode     startNode  = getCorrespondingNodeId(startVert);
   const idSuperArc currentArc = openSuperArc(startNode);
   startUF->addArcToClose(currentArc);
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   (*mt_data_.arcGrowthStats)[currentArc].begin  = _launchGlobalTime.getElapsedTime();
   (*mt_data_.arcGrowthStats)[currentArc].origin = orig;
#endif

   // TASK PROPAGATION
   while (!localPropagation->empty()) {
      // Next vertex

      idVertex currentVert = localPropagation->getNextMinVertex();

      if (!isCorrespondingNull(currentVert) && !isCorrespondingNode(currentVert)) {
         // already seen
         continue;
      } else if (currentVert == startVert && localOrder != -1) {
         // first node can be duplicate after merge
         continue;
      }

      // local order to avoid sort
      (*mt_data_.visitOrder)[currentVert] = localOrder++;

      // Saddle & Last detection + propagation
      bool isSaddle, isLast;
      tie(isSaddle, isLast) = propage(*localPropagation, startUF);

      // regular propagation
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write seq_cst
#endif
      (*mt_data_.ufs)[currentVert] = startUF;

      // Saddle case
      if (isSaddle) {

         // need a node on this vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write seq_cst
#endif
         (*mt_data_.openNodeOnVert)[currentVert] = 1;

         // If last close all and merge
         if (isLast) {
            // finish works here
            closeAndMergeOnSaddle(currentVert);

            // last task detection
            idNode remainingTasks;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read seq_cst
#endif
            remainingTasks = mt_data_.activeTasks;

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.arcGrowthStats)[currentArc].end    = _launchGlobalTime.getElapsedTime();
            (*mt_data_.arcGrowthStats)[currentArc].remain = remainingTasks;
#endif

            if (remainingTasks == 1) {
               // only trunk remaining
               return;
            }

            // made a node on this vertex
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write seq_cst
#endif
            (*mt_data_.openNodeOnVert)[currentVert] = 0;

            // recursively continue
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskyield
#endif
            arcGrowth(currentVert, orig);
         } else {
            // Active tasks / threads
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update seq_cst
#endif
            mt_data_.activeTasks--;

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            // In case of saddle not last, we also do stats
            idNode remainingTasks;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read seq_cst
#endif
            remainingTasks = mt_data_.activeTasks;

            (*mt_data_.arcGrowthStats)[currentArc].end    = _launchGlobalTime.getElapsedTime();
            (*mt_data_.arcGrowthStats)[currentArc].remain = remainingTasks;
#endif
         }

         // stop at saddle
         return;
      }

      if (currentVert != startVert) {
         updateCorrespondingArc(currentVert, currentArc);
      }
      getSuperArc(currentArc)->setLastVisited(currentVert);

      // Intermediate trunk detection
#ifndef TTK_DISABLE_FTC_EARLY_TRUNK
      if (localOrder % TRUNKCHECKCHUNK == 0) {


         idNode remainingTasks;
# ifdef TTK_ENABLE_OPENMP
# pragma omp atomic read seq_cst
# endif
         remainingTasks = mt_data_.activeTasks;
         if (remainingTasks == 1) {
// Force the down vertex of the arc to be opened
# ifdef TTK_ENABLE_OPENMP
# pragma omp atomic write seq_cst
# endif
            (*mt_data_.openNodeOnVert)[startVert] = 1;

# ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.arcGrowthStats)[currentArc].end    = _launchGlobalTime.getElapsedTime();
            (*mt_data_.arcGrowthStats)[currentArc].remain = 1;
# endif
            // go trunk
            return;
         }
      }
#endif

   }  // end wile propagation

   // close root
   const idVertex closeVert      = getSuperArc(currentArc)->getLastVisited();
   bool           existCloseNode = isCorrespondingNode(closeVert);
   idNode closeNode = (existCloseNode) ? getCorrespondingNodeId(closeVert) : makeNode(closeVert);
   closeSuperArc(currentArc, closeNode);
   // getSuperArc(currentArc)->decrNbSeen();
   idNode rootPos             = mt_data_.roots->getNext();
   (*mt_data_.roots)[rootPos] = closeNode;

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   (*mt_data_.arcGrowthStats)[currentArc].end = _launchGlobalTime.getElapsedTime();
#endif
}

void FTCTree_MT::build(const bool computeCT)
{
    string treeString;
   // Comparator init (template)
   initComp();
   switch(mt_data_.treeType){
       case TreeType::Join:
           treeString = "JT";
           break;
       case TreeType::Split:
           treeString = "ST";
           break;
       default:
           treeString = "computeCT";
           break;
   }

   // Build Merge treeString using tasks
   DebugTimer precomputeTime;
   bool fromCT = leafSearch();
   if(!fromCT)
   printTime(precomputeTime, "[FTC] leaf Search " + treeString, scalars_->getSize(), 3);

   DebugTimer buildTime;
   leafGrowth();
   int nbProcessed = 0;
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
   // count process
   for (int i = 0; i < scalars_->getSize(); i++) {
       if((*mt_data_.vert2tree)[i] != nullCorresp)
           ++nbProcessed;
   }
#endif
   printTime(buildTime, "[FTC] leaf Growth "+treeString, nbProcessed, 3);

   DebugTimer bbTime;
   idVertex bbSize = trunk(computeCT);
   printTime(bbTime, "[FTC] trunk "+treeString, bbSize, 3);

   // Segmentation
   if (computeCT && params_->segm) {
      DebugTimer segmTime;
      buildSegmentation();
      printTime(segmTime, "[FTC] segment " + treeString, scalars_->getSize(), 3);
   }
}

void FTCTree_MT::buildSegmentation(const idVertex begin, const idVertex stop)
{
   // const idVertex startVert = begin == nullVertex ? (*scalars_->mirrorVertices)[0] : begin;
   // const idVertex endVert =
   //     begin == nullVertex ? (*scalars_->mirrorVertices)[scalars_->size - 1] : stop;
   const idSuperArc nbArcs = mt_data_.superArcs->size();

   // Make reserve
   // SuperArc i correspond to segment i,
   // one arc correspond to one segment
   vector<idVertex> sizes(nbArcs);

#ifdef TTK_FTC_GRAINSIZE_SEGM_ARCS
   const int grainSizeSegmArcs = TTK_FTC_GRAINSIZE_SEGM_ARCS;
#else
   const int grainSizeSegmArcs = 50;
#endif

   // get the size of each segments_
   const idSuperArc arcChunkSize = getChunkSize(nbArcs, scalars_->getSize(), grainSizeSegmArcs);
   const idSuperArc arcChunkNb   = getChunkCount(nbArcs, scalars_->getSize(), grainSizeSegmArcs);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   mt_data_.segmStats->resize((idVertex)arcChunkNb + arcChunkNb);
#endif

   for(idSuperArc arcChunkId = 0; arcChunkId < arcChunkNb; ++arcChunkId){
       // WHY shared(sizes) is needed ??
       // Here, using priority clause caused segfault on the taskwait
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(arcChunkId) shared(sizes)
#endif
       {
           const idSuperArc lowerBound = arcChunkId*arcChunkSize;
           const idSuperArc upperBound = min(nbArcs, (arcChunkId+1)*arcChunkSize );
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
           (*mt_data_.segmStats)[arcChunkId].begin  = _launchGlobalTime.getElapsedTime();
           (*mt_data_.segmStats)[arcChunkId].origin = lowerBound;
#endif
           for (idSuperArc a = lowerBound; a < upperBound; ++a) {
              sizes[a] = max(0, (*mt_data_.superArcs)[a].getNbVertSeen() - 1);
           }
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
           (*mt_data_.segmStats)[arcChunkId].end = _launchGlobalTime.getElapsedTime();
#endif
       }
   }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

   // change segments size using the created vector
   mt_data_.segments_.resize(sizes);

   DebugTimer segmentsSet;
   // Fill segments using vert2tree
   // current status of the segmentation of this arc
   vector<idVertex> posSegm(nbArcs, 0);

#ifdef TTK_FTC_GRAINSIZE_SEGM_VERTS
   const int grainSizeSegmVerts = TTK_FTC_GRAINSIZE_SEGM_VERTS;
#else
   const int grainSizeSegmVerts = 400000;
#endif

   // Segments are connex region of geometrie forming
   // the segmentation (sorted in ascending order)
   const idVertex nbVert    = scalars_->getSize();
   const idVertex chunkSize = getChunkSize(nbVert, nbVert, grainSizeSegmVerts);
   const idVertex chunkNb   = getChunkCount(nbVert, nbVert, grainSizeSegmVerts);
   for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(posSegm)
#endif
       {
          const idVertex lowerBound = chunkId * chunkSize;
          const idVertex upperBound = min(nbVert, (chunkId+1)*chunkSize);
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
          (*mt_data_.segmStats)[arcChunkNb + chunkId].begin  = _launchGlobalTime.getElapsedTime();
          (*mt_data_.segmStats)[arcChunkNb + chunkId].origin = lowerBound;
#endif
          for (idVertex i = lowerBound; i < upperBound; ++i) {
             const auto vert = scalars_->getSortedVert(i);
             if (isCorrespondingArc(vert)) {
                idSuperArc sa = getCorrespondingSuperArcId(vert);
                idVertex   vertToAdd;
                if((*mt_data_.visitOrder)[vert] != nullVertex){
                   // Opposite order for Split Tree
                   vertToAdd = (*mt_data_.visitOrder)[vert];
                   if(isST()) vertToAdd = getSuperArc(sa)->getNbVertSeen() - vertToAdd -2;
                   mt_data_.segments_[sa][vertToAdd] = vert;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                    posSegm[sa]++;
                } else if (mt_data_.trunkSegments->size() == 0){
                    // MT computation
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
                   vertToAdd = posSegm[sa]++;

                   mt_data_.segments_[sa][vertToAdd] = vert;
                }

             }  // end is arc
          } // end for
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
          (*mt_data_.segmStats)[arcChunkNb + chunkId].end  = _launchGlobalTime.getElapsedTime();
#endif
       } // end task
   }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

   printTime(segmentsSet, "[FTC] segmentation set vertices", -1, 4);

   if (mt_data_.trunkSegments->size() == 0) {
      // sort arc that have been filled by the trunk
      // only for MT
      DebugTimer segmentsSortTime;
      for (idSuperArc a = 0; a < nbArcs; ++a) {
         if (posSegm[a]) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(a)
#endif
            mt_data_.segments_[a].sort(scalars_);
         }
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
      printTime(segmentsSortTime, "[FTC] segmentation sort vertices", -1, 4);
   } else {
       // Contour tree: we create the arc segmentation for arcs in the trunk
       DebugTimer segmentsArcTime;
       for (idSuperArc a = 0; a < nbArcs; ++a) {
          // CT computation, we have already the vert list
          if ((*mt_data_.trunkSegments)[a].size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(a)
#endif
             mt_data_.segments_[a].createFromList(scalars_, (*mt_data_.trunkSegments)[a],
                                                  mt_data_.treeType == TreeType::Split);
          }
       }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
      printTime(segmentsArcTime, "[FTC] segmentation arcs lists", -1, 4);
   }

   // Update SuperArc region
   // ST have a segmentation wich is in the reverse-order of its build
   // ST have a segmentation sorted in ascending order as JT
   for(idSuperArc arcChunkId = 0; arcChunkId < arcChunkNb; ++arcChunkId){
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(arcChunkId)
#endif
      {
         const idSuperArc lowerBound = arcChunkId * arcChunkSize;
         const idSuperArc upperBound = min(nbArcs, (arcChunkId + 1) * arcChunkSize);
         for (idSuperArc a = lowerBound; a < upperBound; ++a) {
            // avoid empty region
            if (mt_data_.segments_[a].size()) {
               (*mt_data_.superArcs)[a].concat(mt_data_.segments_[a].begin(),
                                               mt_data_.segments_[a].end());
            }
         }

      }
   }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
}

idVertex FTCTree_MT::calculValence(const idVertex idVert)
{
   const auto &neighNumb = mesh_->getVertexNeighborNumber(idVert);
   valence     val       = 0;

   for (valence n = 0; n < neighNumb; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(idVert, n, neigh);
      comp_.vertLower(neigh, idVert) && ++val;
   }

   return val;
}

FTCTree_MT *FTCTree_MT::clone() const
{
   FTCTree_MT *newMT = new FTCTree_MT(params_, mesh_, scalars_, mt_data_.treeType);

   newMT->mt_data_.superArcs = mt_data_.superArcs;
   newMT->mt_data_.nodes     = mt_data_.nodes;
   newMT->mt_data_.leaves    = mt_data_.leaves;
   newMT->mt_data_.roots     = mt_data_.roots;
   newMT->mt_data_.vert2tree = mt_data_.vert2tree;

   return newMT;
}

void FTCTree_MT::closeAndMergeOnSaddle(const idVertex saddleVert)
{
   const idNode closeNode = makeNode(saddleVert);

   // Union of the UF coming here (merge propagation and closing arcs)
   const auto &nbNeigh = mesh_->getVertexNeighborNumber(saddleVert);
   for (valence n = 0; n < nbNeigh; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(saddleVert, n, neigh);

      if (comp_.vertLower(neigh, saddleVert)) {
         if ((*mt_data_.ufs)[neigh]->find() != (*mt_data_.ufs)[saddleVert]->find()) {

            (*mt_data_.ufs)[saddleVert] =
                AtomicUF::makeUnion((*mt_data_.ufs)[saddleVert], (*mt_data_.ufs)[neigh]);
         }
      }
   }

   // close arcs on this node
   closeArcsUF(closeNode, (*mt_data_.ufs)[saddleVert]);

   (*mt_data_.ufs)[saddleVert]->find()->mergeStates();
   (*mt_data_.ufs)[saddleVert]->find()->setExtremum(saddleVert);
}

void FTCTree_MT::closeArcsUF(idNode closeNode, UF uf)
{
   for (const idSuperArc sa : uf->find()->getOpenedArcs()) {
      closeSuperArc(sa, closeNode);
   }
   uf->find()->clearOpenedArcs();
}

void FTCTree_MT::closePendingTrunkNodes(idVertex saddleVert)
{
   idNode closeNode = makeNode(saddleVert);

   const auto &nbNeigh = mesh_->getVertexNeighborNumber(saddleVert);
   for (valence n = 0; n < nbNeigh; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(saddleVert, n, neigh);

      if (comp_.vertLower(neigh, saddleVert)) {
         UF neighUF = (*mt_data_.ufs)[neigh];
         if (neighUF && neighUF->find()->getNbOpenedArcs()) {
            for (const auto arc : neighUF->find()->getOpenedArcs()) {
               if (getSuperArc(arc)->getUpNodeId() == nullNodes) {
                  closeSuperArc(arc, closeNode);
               }
            }

            neighUF->find()->clearOpenedArcs();
         }
      }
   }
}

void FTCTree_MT::closeSuperArc(idSuperArc superArcId, idNode upNodeId)
{
#ifndef TTK_ENABLE_KAMIKAZE

   if ((size_t)superArcId >= getNumberOfSuperArcs()) {
      cout << "[Merge Tree] closeSuperArc on a inexisting arc !" << endl;
      return;
   }

   if ((size_t)upNodeId >= getNumberOfNodes()) {
      cout << "[Merge Tree] closeOpenedArc on a inexisting node !" << endl;
      return;
   }

#endif
   (*mt_data_.superArcs)[superArcId].setUpNodeId(upNodeId);
   (*mt_data_.nodes)[upNodeId].addDownSuperArcId(superArcId);
}

void FTCTree_MT::delNode(idNode node)
{
   Node *mainNode = getNode(node);

   if (mainNode->getNumberOfUpSuperArcs() == 0 && mainNode->getNumberOfDownSuperArcs() == 0) {
      return;
   }

   if (mainNode->getNumberOfUpSuperArcs() == 0) {

      // Root: No Superarc
#ifndef TTK_ENABLE_KAMIKAZE
      if (mainNode->getNumberOfDownSuperArcs() != 1) {
         // Root with several childs: impossible /\ .
         cout << endl << "[FTCTree_MT]:delNode won't delete ";
         cout << mainNode->getVertexId() << " (root) with ";
         cout << static_cast<unsigned>(mainNode->getNumberOfDownSuperArcs()) << " down ";
         cout << static_cast<unsigned>(mainNode->getNumberOfUpSuperArcs()) << " up " << endl;
         return;
      }
#endif

      idSuperArc downArc  = mainNode->getDownSuperArcId(0);
      Node *     downNode = getNode((*mt_data_.superArcs)[downArc].getDownNodeId());

      downNode->removeUpSuperArcId(downArc);
      mainNode->clearDownSuperArcs();

   } else if (mainNode->getNumberOfDownSuperArcs() < 2) {
      // Have one up arc

      // We delete the upArc of this node,
      // if there is a down arc, we reattach it to the upNode

      idSuperArc upArc  = mainNode->getUpSuperArcId(0);
      idNode     upId   = (*mt_data_.superArcs)[upArc].getUpNodeId();
      Node *     upNode = getNode(upId);

      upNode->removeDownSuperArcId(upArc);
      mainNode->clearUpSuperArcs();

      if (mainNode->getNumberOfDownSuperArcs()) {
         // Have one down arc

         // Reconnect
         idSuperArc downArc = mainNode->getDownSuperArcId(0);
         (*mt_data_.superArcs)[downArc].setUpNodeId(upId);
         upNode->addDownSuperArcId(downArc);
         mainNode->clearDownSuperArcs();

         // Segmentation
         (*mt_data_.superArcs)[downArc].concat((*mt_data_.superArcs)[upArc]);
      }
   }
#ifndef TTK_ENABLE_KAMIKAZE
   else
      cerr << "delete node with multiple childrens " << endl;
#endif
}

void FTCTree_MT::finalizeSegmentation(void)
{
   for (auto &arc : *mt_data_.superArcs) {
      arc.createSegmentation(scalars_);
   }
}

tuple<idVertex, idVertex> FTCTree_MT::getBoundsFromVerts(const vector<idVertex> &trunkVerts) const
{
   idVertex begin, stop;

   if (isST()) {
      begin = 0;
      stop  = scalars_->getMirror(trunkVerts[0]);
   } else if (isJT()) {
      begin = scalars_->getMirror(trunkVerts[0]);
      stop  = scalars_->getSize();
   } else {
      begin = scalars_->getMirror(trunkVerts[0]);
      stop  = scalars_->getMirror(trunkVerts[trunkVerts.size() - 1]);
   }

   return make_tuple(begin, stop);
}

Node *FTCTree_MT::getDownNode(const SuperArc *a)
{
   return &((*mt_data_.nodes)[a->getDownNodeId()]);
}

idNode FTCTree_MT::getDownNodeId(const SuperArc *a)
{
   return a->getDownNodeId();
}

Node *FTCTree_MT::getLowerNode(const SuperArc *a)
{
   if (isST())
      return getUpNode(a);

   return getDownNode(a);
}

idNode FTCTree_MT::getLowerNodeId(const SuperArc *a)
{
   if (isST())
      return getUpNodeId(a);

   return getDownNodeId(a);
}

Node *FTCTree_MT::getUpNode(const SuperArc *a)
{
   return &((*mt_data_.nodes)[a->getUpNodeId()]);
}

idNode FTCTree_MT::getUpNodeId(const SuperArc *a)
{
   return a->getUpNodeId();
}

Node *FTCTree_MT::getUpperNode(const SuperArc *a)
{
   if (isST())
      return getDownNode(a);

   return getUpNode(a);
}

idNode FTCTree_MT::getUpperNodeId(const SuperArc *a)
{
   if (isST())
      return getDownNodeId(a);

   return getUpNodeId(a);
}

idNode FTCTree_MT::getNodeFromVertInc(const vector<idVertex> &range,
                                     const idVertex          v,
                                     const idNode            last) const
{
    idNode idRes = last;
    // last element of rande is maximum, cannot be lower.
    // no need to check size
    while (comp_.vertLower(range[idRes + 1], v)) {
       ++idRes;
    }
    return idRes;
}

idNode FTCTree_MT::getNodeFromVertDicho(const vector<idVertex> &range, const idVertex v) const
{
   idNode rangeBegin = 0;
   // Last node is virtual, avoid it
   idNode rangeEnd   = range.size()-1;

   while (rangeEnd - rangeBegin > 1) {
     idNode rangeMiddle = (rangeBegin + rangeEnd) / 2;
     if (comp_.vertHigher(range[rangeMiddle], v)) {
        rangeEnd = rangeMiddle;
     } else {
        rangeBegin = rangeMiddle;
     }
   }
   return rangeBegin;
}

idSuperArc FTCTree_MT::insertNode(Node *node, const bool segm)
{
   // Normal insert : existing arc stay below inserted (JT example)
   //  *   - <- upNodeId
   //  | \ |   <- newSA
   //  |   * <- newNodeId
   //  |   |   <- currentSA
   //  - - -
   // already present
   if (isCorrespondingNode(node->getVertexId())) {
      Node *myNode = vertex2Node(node->getVertexId());
      // If it has been hidden / replaced we need to re-make it
      idSuperArc correspondingArcId = myNode->getUpSuperArcId(0);
      updateCorrespondingArc(myNode->getVertexId(), correspondingArcId);
   }

   idNode     upNodeId, newNodeId;
   idSuperArc currentSA, newSA;

   // Create new node
   currentSA = getCorrespondingSuperArcId(node->getVertexId());
   upNodeId  = (*mt_data_.superArcs)[currentSA].getUpNodeId();
   newNodeId = makeNode(node);

   // Connectivity
   // Insert only node inside the partition : created arc don t cross
   newSA = makeSuperArc(newNodeId, upNodeId);

   (*mt_data_.superArcs)[currentSA].setUpNodeId(newNodeId);
   (*mt_data_.nodes)[upNodeId].removeDownSuperArcId(currentSA);
   (*mt_data_.nodes)[newNodeId].addDownSuperArcId(currentSA);

   // cut the vertex list at the node position and
   // give each arc its part.
   if (segm) {
      if (mt_data_.treeType == TreeType::Split) {
         (*mt_data_.superArcs)[newSA].concat(
             get<1>((*mt_data_.superArcs)[currentSA].splitBack(node->getVertexId(), scalars_)));
      } else {
         (*mt_data_.superArcs)[newSA].concat(
             get<1>((*mt_data_.superArcs)[currentSA].splitFront(node->getVertexId(), scalars_)));
      }
   }

   return newSA;
}

void FTCTree_MT::leafGrowth()
{
   const auto &nbLeaves = getNumberOfLeaves();

   // printLeavesStats();

   // memory allocations here
   allocVectPropagations(nbLeaves + 2);

   // elevation: trunk only
   if (nbLeaves == 1) {
      const idVertex v            = (*mt_data_.nodes)[0].getVertexId();
      (*mt_data_.openNodeOnVert)[v] = 1;
      return;
   }

   mt_data_.activeTasks = nbLeaves;

   // gcc uses a stack (then a queue) for the tasks,
   // but our algorithm is really better (in sequential at least)
   // with an arc growth haveing the same order as the corresponding tree
   // That is why we still have the sort step
   auto comp = [this](const idNode a, const idNode b) {
#ifdef HIGHER
      return this->comp_.vertHigher(this->getNode(a)->getVertexId(),
                                    this->getNode(b)->getVertexId());
#else
      return this->comp_.vertLower(this->getNode(a)->getVertexId(),
                                   this->getNode(b)->getVertexId());
#endif
   };

   sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), comp);

   for (idNode n = 0; n < nbLeaves; ++n) {
      const idNode l = (*mt_data_.leaves)[n];
      int          v = getNode(l)->getVertexId();

#ifdef TTK_ENABLE_OPENMP
#pragma omp task untied PRIOR(isPrior())
#endif
      arcGrowth(v, n);
// #pragma omp taskyeil
   }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
}

int FTCTree_MT::leafSearch()
{
   int ret = 0;
   // if not already computed by CT
   if (getNumberOfNodes() == 0) {
      const auto nbScalars = scalars_->getSize();
      const auto chunkSize = getChunkSize();
      const auto chunkNb   = getChunkCount();

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
      mt_data_.extrStats->resize(chunkNb);
#endif

      // Extrema extract and launch tasks
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) PRIOR(isPrior())
#endif
         {
            const idVertex lowerBound = chunkId * chunkSize;
            const idVertex upperBound = min(nbScalars, (chunkId + 1) * chunkSize);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.extrStats)[chunkId].begin = _launchGlobalTime.getElapsedTime();
#endif

            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const auto neighNumb = mesh_->getVertexNeighborNumber(v);
               bool       val       = false;

               for (valence n = 0; n < neighNumb; ++n) {
                  idVertex neigh;
                  mesh_->getVertexNeighbor(v, n, neigh);
                  if (comp_.vertLower(neigh, v)) {
                     val = true;
                     break;
                  }
               }

               if (!val) {
                  const idNode nextLeaf        = mt_data_.leaves->getNext();
                  (*mt_data_.leaves)[nextLeaf] = makeNode(v);
               }
            }

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.extrStats)[chunkId].end = _launchGlobalTime.getElapsedTime();
#endif
         }
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

      for (const idNode n : *mt_data_.leaves) {
         const idVertex v   = getNode(n)->getVertexId();
         (*mt_data_.ufs)[v] = new AtomicUF(v);
      }
      mt_data_.activeTasks += getNumberOfLeaves();
   } else {
      // Here The CT has already done this job
      // ! Not compatible with this new version of the function
      ret = 1;
   }

   if (debugLevel_ >= 4) {
      cout << "- [FTC] found " << mt_data_.leaves->size() << " leaves" << endl;
   }

   return ret;
}

idNode FTCTree_MT::makeNode(idVertex vertexId)
{
#ifndef TTK_ENABLE_KAMIKAZE
   if (vertexId < 0 || vertexId >= scalars_->getSize()) {
      cout << "[Merge Tree] make node, wrong vertex :" << vertexId << " on " << scalars_->getSize()
           << endl;
      return -1;
   }
#endif

   if (isCorrespondingNode(vertexId)) {
      return getCorrespondingNodeId(vertexId);
   }

   idNode newNodeId = mt_data_.nodes->getNext();
   (*mt_data_.nodes)[newNodeId].setVertexId(vertexId);
   updateCorrespondingNode(vertexId, newNodeId);

   return newNodeId;
}

idNode FTCTree_MT::makeNode(const Node *const n)
{
   return makeNode(n->getVertexId());
}

idSuperArc FTCTree_MT::makeSuperArc(idNode downNodeId, idNode upNodeId, const bool withNodes)

{
   idSuperArc newSuperArcId = mt_data_.superArcs->getNext();
   (*mt_data_.superArcs)[newSuperArcId].setDownNodeId(downNodeId);
   (*mt_data_.superArcs)[newSuperArcId].setUpNodeId(upNodeId);

   if (withNodes) {
      (*mt_data_.nodes)[downNodeId].addUpSuperArcId(newSuperArcId);
      (*mt_data_.nodes)[upNodeId].addDownSuperArcId(newSuperArcId);
   }

   return newSuperArcId;
}

void FTCTree_MT::move(FTCTree_MT *mt)
{
   // we already have common data
   mt_data_.superArcs     = mt->mt_data_.superArcs;
   mt->mt_data_.superArcs = nullptr;
   mt_data_.nodes         = mt->mt_data_.nodes;
   mt->mt_data_.nodes     = nullptr;
   mt_data_.leaves        = mt->mt_data_.leaves;
   mt->mt_data_.leaves    = nullptr;
   mt_data_.roots         = mt->mt_data_.roots;
   mt->mt_data_.roots     = nullptr;
   mt_data_.vert2tree     = mt->mt_data_.vert2tree;
   mt->mt_data_.vert2tree = nullptr;
}

void FTCTree_MT::normalizeIds(void)
{
   DebugTimer normTime;
   sortLeaves(true);

   auto getNodeParentArcNb = [&](const idNode curNode, const bool goUp) -> idSuperArc {
      if (goUp) {
         return getNode(curNode)->getNumberOfUpSuperArcs();
      }

      return getNode(curNode)->getNumberOfDownSuperArcs();
   };

   auto getNodeParentArc = [&](const idNode curNode, const bool goUp, idSuperArc i) -> idSuperArc {
      if (goUp) {
         return getNode(curNode)->getUpSuperArcId(i);
      }

      return getNode(curNode)->getDownSuperArcId(i);
   };

   auto getArcParentNode = [&](const idSuperArc curArc, const bool goUp) -> idNode {
      if (goUp) {
         return getSuperArc(curArc)->getUpNodeId();
      }

      return getSuperArc(curArc)->getDownNodeId();
   };

   std::queue<tuple<idNode, bool>> q;
   std::stack<tuple<idNode, bool>> qr;
   for (const idNode n : *mt_data_.leaves) {
      bool goUp = isJT() || isST() || getNode(n)->getNumberOfUpSuperArcs();
      if (goUp)
         q.emplace(make_tuple(n, goUp));
      else
         qr.emplace(make_tuple(n, goUp));
   }

   while (!qr.empty()) {
      q.emplace(qr.top());
      qr.pop();
   }

   // Normalized id
   idSuperArc nIdMin = 0;
   idSuperArc nIdMax = getNumberOfSuperArcs() - 1;

   vector<bool> seenUp(getNumberOfSuperArcs(), false);
   vector<bool> seenDown(getNumberOfSuperArcs(), false);

   while (!q.empty()) {
      bool   goUp;
      idNode curNodeId;
      tie(curNodeId, goUp) = q.front();
      q.pop();

      if (goUp)
         sortUpArcs(curNodeId);
      else
         sortDownArcs(curNodeId);

      // Assign arc above
      const idSuperArc nbArcParent = getNodeParentArcNb(curNodeId, goUp);
      for (idSuperArc pid = 0; pid < nbArcParent; pid++) {
         const idSuperArc currentArcId = getNodeParentArc(curNodeId, goUp, pid);
         if (goUp) {
            if (getSuperArc(currentArcId)->getNormalizedId() == nullSuperArc) {
               getSuperArc(currentArcId)->setNormalizeIds(nIdMin++);
            }
            if (!seenUp[currentArcId]) {
               q.emplace(make_tuple(getArcParentNode(currentArcId, goUp), goUp));
               seenUp[currentArcId] = true;
            }
         } else {
            if (getSuperArc(currentArcId)->getNormalizedId() == nullSuperArc) {
               getSuperArc(currentArcId)->setNormalizeIds(nIdMax--);
            }
            if (!seenDown[currentArcId]) {
               q.emplace(make_tuple(getArcParentNode(currentArcId, goUp), goUp));
               seenDown[currentArcId] = true;
            }
         }
      }
   }

#ifndef TTK_ENABLE_KAMIKAZE
   if (std::abs((long)nIdMax - (long)nIdMin) > 1) {
      cout << "[FTC] error during normalize, tree compromized: " << nIdMin << " " << nIdMax << endl;
   }
#endif

   printTime(normTime, "[FTC] normalize ids", -1, 4);
}

idSuperArc FTCTree_MT::openSuperArc(idNode downNodeId)
{
#ifndef TTK_ENABLE_KAMIKAZE
   if ((size_t)downNodeId >= getNumberOfNodes()) {
      cout << "[Merge Tree] openSuperArc on a inexisting node !" << endl;
      return -2;
   }
#endif

   idSuperArc newSuperArcId = mt_data_.superArcs->getNext();
   (*mt_data_.superArcs)[newSuperArcId].setDownNodeId(downNodeId);
   (*mt_data_.nodes)[downNodeId].addUpSuperArcId(newSuperArcId);

   return newSuperArcId;
}

string FTCTree_MT::printArc(idSuperArc a)
{
   const SuperArc *sa = getSuperArc(a);
   stringstream    res;
   res << a;
   res << " : ";
   res << getNode(sa->getDownNodeId())->getVertexId() << " -- ";
   if (sa->getUpNodeId() != nullNodes) {
      res << getNode(sa->getUpNodeId())->getVertexId();
   } else {
       res << "X";
   }

   res.seekg(0, ios::end);
   while (res.tellg() < 25) {
      res << " ";
      res.seekg(0, ios::end);
   }
   res.seekg(0, ios::beg);

   res << "segm #" << sa->regionSize() << " / " << scalars_->getSize();  // << " -> ";

   res.seekg(0, ios::end);

   while (res.tellg() < 45) {
      res << " ";
      res.seekg(0, ios::end);
   }
   res.seekg(0, ios::beg);

   res << sa->printReg();
   res << " (" << sa->getNbVertSeen() <<  ")";
   return res.str();
}

string FTCTree_MT::printNode(idNode n)
{
   const Node * node = getNode(n);
   stringstream res;
   res << n;
   res << " : (";
   res << node->getVertexId() << ") \\ ";

   for (idSuperArc i = 0; i < node->getNumberOfDownSuperArcs(); ++i) {
      res << "+";
      res << node->getDownSuperArcId(i) << " ";
   }

   res << " / ";

   for (idSuperArc i = 0; i < node->getNumberOfUpSuperArcs(); ++i) {
      res << "+";
      res << node->getUpSuperArcId(i) << " ";
   }

   return res.str();
}


int FTCTree_MT::printTime(DebugTimer &t, const string &s, idVertex nbScalars, const int debugLevel) const
{
   if (nbScalars == -1) {
      nbScalars = scalars_->getSize();
   }

   if (debugLevel_ >= debugLevel) {
      stringstream st;
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
      int          speed = nbScalars / t.getElapsedTime();
#endif
      for (int i = 3; i < debugLevel; i++)
         st << "-";
      st << s << " in ";
      st.seekg(0, ios::end);
      while (st.tellg() < 25) {
         st << " ";
         st.seekg(0, ios::end);
      }
      st.seekg(0, ios::beg);
      st << t.getElapsedTime();

#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
      st.seekg(0, ios::end);
      while (st.tellg() < 35) {
         st << " ";
         st.seekg(0, ios::end);
      }
      st.seekg(0, ios::beg);
      st << " at " << speed << " vert/s";
#endif
      cout << st.str() << endl;
   }
   return 1;
}

void FTCTree_MT::printTree2()
{
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
   {
      cout << "Nodes----------" << endl;
      for (idNode nid = 0; nid < getNumberOfNodes(); nid++) {
         cout << printNode(nid) << endl;
      }

      cout << "Arcs-----------" << endl;
      for (idSuperArc said = 0; said < getNumberOfSuperArcs(); ++said) {
         cout << printArc(said) << endl;
      }

      cout << "Leaves" << endl;
      for (const auto &l : *mt_data_.leaves)
         cout << " " << (*mt_data_.nodes)[l].getVertexId();
      cout << endl;

      cout << "Roots" << endl;
      for (const auto &r : *mt_data_.roots)
         cout << " " << (*mt_data_.nodes)[r].getVertexId();
      cout << endl;
   }
}

void FTCTree_MT::printLeavesStats(void)
{
   const auto nbLeaves = getNumberOfLeaves();

   // AVG

   float avgScal = 0;
   float avgPos[3] = {0,0,0};

   for (const idNode lid : *mt_data_.leaves) {
      const idVertex lvert = getNode(lid)->getVertexId();
      const float scal = scalars_->getVal(lvert);
      float pos[3];
      mesh_->getVertexPoint(lvert, pos[0], pos[1], pos[2]);

      avgScal += scal;
      avgPos[0] += pos[0];
      avgPos[1] += pos[1];
      avgPos[2] += pos[2];
   }

   avgScal /= nbLeaves;
   avgPos[0] /= nbLeaves;
   avgPos[1] /= nbLeaves;
   avgPos[2] /= nbLeaves;

   // STDEV

   float stddevScal = 0;
   float stddevPos[3] = {0,0,0};

   for (const idNode lid : *mt_data_.leaves) {
      const idVertex lvert = getNode(lid)->getVertexId();
      const float scal = scalars_->getVal(lvert);
      float pos[3];
      mesh_->getVertexPoint(lvert, pos[0], pos[1], pos[2]);

      stddevScal += pow(scal - avgScal, 2);
      stddevPos[0] += pow(pos[0] - avgPos[0], 2);
      stddevPos[1] += pow(pos[1] - avgPos[1], 2);
      stddevPos[2] += pow(pos[2] - avgPos[2], 2);
   }

   stddevScal /= nbLeaves;
   stddevPos[0] /= nbLeaves;
   stddevPos[1] /= nbLeaves;
   stddevPos[2] /= nbLeaves;

   stddevScal = sqrt(stddevScal);
   stddevPos[0] = sqrt(stddevPos[0]);
   stddevPos[1] = sqrt(stddevPos[1]);
   stddevPos[2] = sqrt(stddevPos[2]);

   // Print

   cout << "JT:" << isJT() << " nb:" << nbLeaves;
   cout << " avg:" << avgScal << " - " << avgPos[0] << " " << avgPos[1] << " " << avgPos[2];
   cout << " stddev:" << stddevScal << " - " << stddevPos[0] << " " << stddevPos[1] << " " << stddevPos[2] << endl;
}

tuple<bool, bool> FTCTree_MT::propage(Propagation &localPropagation, UF curUF)
{
   bool        becameSaddle = false, isLast = false;
   const auto &nbNeigh = mesh_->getVertexNeighborNumber(localPropagation.getCurrentMinVertex());
   valence decr = 0;

   // once for all
   auto* curUFF = curUF->find();

   // propagation / is saddle
   for (valence n = 0; n < nbNeigh; ++n) {
      idVertex neigh;
      mesh_->getVertexNeighbor(localPropagation.getCurrentMinVertex(), n, neigh);

      if (comp_.vertLower(neigh, localPropagation.getCurrentMinVertex())) {
         UF neighUF = (*mt_data_.ufs)[neigh];

         // is saddle
         if (!neighUF || neighUF->find() != curUFF) {
            becameSaddle = true;
         } else if (neighUF) {
             ++decr;
         }

      } else {
         if (!(*mt_data_.propagationUfs)[neigh] ||
             (*mt_data_.propagationUfs)[neigh]->find() != curUFF) {
            localPropagation.addNewVertex(neigh);
            (*mt_data_.propagationUfs)[neigh] = curUFF;
         }
      }
   }

   if (becameSaddle) {
      // is last
      valence oldVal;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
      {
         oldVal = (*mt_data_.valences)[localPropagation.getCurrentMinVertex()];
         (*mt_data_.valences)[localPropagation.getCurrentMinVertex()] -= decr;
      }
      if (oldVal == -1) {
         idVertex val = calculValence(localPropagation.getCurrentMinVertex());
         valence  newVal;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
         {
            newVal = (*mt_data_.valences)[localPropagation.getCurrentMinVertex()];
            (*mt_data_.valences)[localPropagation.getCurrentMinVertex()] += val + 1;
         }
         // Check if all other threads reaching this points have finished
         // during valence computation: use oldVal.
         oldVal = newVal + (val + 1) + decr;
      }
      if (oldVal == decr) {
         isLast = true;
      }
   }

   return make_tuple(becameSaddle, isLast);
}

void FTCTree_MT::sortLeaves(const bool para)
{
   auto indirect_sort = [&](const idNode a, const idNode b) {
      return comp_.vertLower(getNode(a)->getVertexId(), getNode(b)->getVertexId());
   };

   if (para) {
#ifdef __clang__
      std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#else
#ifndef _MSC_VER
#ifdef TTK_ENABLE_OPENMP
      __gnu_parallel::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#else
      std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#endif
#else
      std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
#endif
#endif
   } else {
      std::sort(mt_data_.leaves->begin(), mt_data_.leaves->end(), indirect_sort);
   }
}

vector<idNode> FTCTree_MT::sortedNodes(const bool para)
{
   vector<idNode> sortedNodes(mt_data_.nodes->size());
   std::iota(sortedNodes.begin(), sortedNodes.end(), 0);

   auto indirect_sort = [&](const idNode a, const idNode b) {
      return comp_.vertLower(getNode(a)->getVertexId(), getNode(b)->getVertexId());
   };

   if (para) {
#ifdef __clang__
      std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#else
#ifndef _MSC_VER
#ifdef TTK_ENABLE_OPENMP
      __gnu_parallel::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#else
      std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#endif
#else
      std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
#endif
#endif
   } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single
#endif
      {
         std::sort(sortedNodes.begin(), sortedNodes.end(), indirect_sort);
      }
   }

   return sortedNodes;
}

idVertex FTCTree_MT::trunk(const bool ct)
{
   DebugTimer bbTimer;
   const auto &nbScalars = scalars_->getSize();

   for (idVertex v = 0; v < nbScalars; ++v) {
      if ((*mt_data_.openNodeOnVert)[v]) {
         mt_data_.trunkVerts->emplace_back(v);
      }
   }

   const idNode nbNodes = mt_data_.trunkVerts->size();
   if (!nbNodes) {
      return 0;
   }

   // TODO: Make parallel
   sort(mt_data_.trunkVerts->begin(), mt_data_.trunkVerts->end(), comp_.vertLower);

   // If we stopped during an arc growth and not at saddle, we need to
   // close this arc first
   if (isCorrespondingNode((*mt_data_.trunkVerts)[0])) {
      const idNode downArcNode = getCorrespondingNodeId((*mt_data_.trunkVerts)[0]);
      if (getNode(downArcNode)->getNumberOfUpSuperArcs()) {
         const idSuperArc pendingArc = getNode(downArcNode)->getUpSuperArcId(0);
         if (getSuperArc(pendingArc)->getUpNodeId() == nullNodes) {
            const idNode upArcNode = makeNode((*mt_data_.trunkVerts)[1]);
            closeSuperArc(pendingArc, upArcNode);
         }
      }
   }

   // Close arcs pending on nodes of the trunk
   for (const idVertex v : *mt_data_.trunkVerts) {
      closePendingTrunkNodes(v);
   }


#ifdef TTK_FTC_GRAINSIZE_TRUNK_ARCS
   const int grainSizeTrunkArcs = TTK_FTC_GRAINSIZE_TRUNK_ARCS;
#else
   const int grainSizeTrunkArcs = 650;
#endif

   const idNode arcChunkSize = getChunkSize(nbNodes, nbScalars, grainSizeTrunkArcs);
   const idNode arcChunkNb   = getChunkCount(nbNodes, nbScalars, grainSizeTrunkArcs);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   mt_data_.arcTrunkStats->resize(arcChunkNb);
#endif

   // Create arcs of the monotone path called trunk
   for (idNode arcChunkId = 0; arcChunkId < arcChunkNb; arcChunkId++) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(arcChunkId) PRIOR(isPrior())
#endif
      {
         const idNode lowerBound = arcChunkId * arcChunkSize;
         const idNode upperBound = min(nbNodes-1, (arcChunkId + 1) * arcChunkSize);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
         (*mt_data_.arcTrunkStats)[arcChunkId].begin = _launchGlobalTime.getElapsedTime();
#endif

         for (idNode n = lowerBound + 1; n <= upperBound; ++n) {
            const idVertex downArcNode = getCorrespondingNodeId((*mt_data_.trunkVerts)[n - 1]);
            const idVertex upArcNode   = getCorrespondingNodeId((*mt_data_.trunkVerts)[n]);
            if (getNode(downArcNode)->getNumberOfUpSuperArcs()) {
               const idSuperArc pendingArc = getNode(downArcNode)->getUpSuperArcId(0);
               if (getSuperArc(pendingArc)->getUpNodeId() == nullNodes) {
                  closeSuperArc(pendingArc, upArcNode);
               }
            } else {
               const idSuperArc newArc = makeSuperArc(downArcNode, upArcNode);
               getSuperArc(newArc)->setLastVisited((*mt_data_.trunkVerts)[n]);
            }
         }

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
         (*mt_data_.arcTrunkStats)[arcChunkId].end = _launchGlobalTime.getElapsedTime();
#endif
      }
   }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

   const idSuperArc lastArc = openSuperArc(getCorrespondingNodeId((*mt_data_.trunkVerts)[nbNodes - 1]));

   // Root (close last arc)
   // if several CC still the backbone is only in one.
   // But the root may not be the max node of the whole dataset: TODO
   const idNode rootNode = makeNode(scalars_->getSortedVert(isST() ? 0 : nbScalars - 1));
   closeSuperArc(lastArc, rootNode);
   getSuperArc(lastArc)->setLastVisited(getNode(rootNode)->getVertexId());

   printTime(bbTimer, "[FTC] trunk seq.", -1, 4);
   bbTimer.reStart();

   // Segmentation
   idVertex begin, stop, processed;
   tie(begin, stop) = getBoundsFromVerts(*mt_data_.trunkVerts);

   // Add fake last vertex to the mt_data_.trunkVerts to avoid size checking
   if (isST()) {
      mt_data_.trunkVerts->emplace_back(scalars_->getSortedVert(0));
   } else {
      mt_data_.trunkVerts->emplace_back(scalars_->getSortedVert(nbScalars - 1));
   }

   if(ct){
       processed = trunkCTSegmentation(*mt_data_.trunkVerts, begin, stop);
   } else {
       processed = trunkSegmentation(*mt_data_.trunkVerts, begin, stop);
   }
   printTime(bbTimer, "[FTC] trunk para.", -1, 4);

   return processed;
}

idVertex FTCTree_MT::trunkCTSegmentation(const vector<idVertex> &trunkVerts,
                                         const idVertex begin,
                                         const idVertex stop)
{
#ifdef TTK_FTC_GRAINSIZE_TRUNK_VERTS
   const int grainSizeTrunkVerts = TTK_FTC_GRAINSIZE_TRUNK_VERTS;
#else
   const int grainSizeTrunkVerts = 50000;
#endif
   const auto sizeBackBone  = abs(stop - begin);
   const auto chunkSize     = getChunkSize(sizeBackBone, scalars_->getSize(), grainSizeTrunkVerts);
   const auto chunkNb       = getChunkCount(sizeBackBone, scalars_->getSize(), grainSizeTrunkVerts);
   // si pas efficace vecteur de la taille de node ici a la place de acc
   mt_data_.trunkSegments->resize(getNumberOfSuperArcs());

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   mt_data_.trunkStats->resize(chunkNb);
#endif

   // CAUTION: Duplicate code to avoid comp_.vertLower indirection
   if(isST()){
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(trunkVerts) PRIOR(isPrior())
#endif
         {
            idNode   currentNode = 0;
            idVertex nextNodeThresold = 0;
            idVertex nextNodeValue    = 0;

            const idVertex lowerBound = begin + chunkId * chunkSize;
            const idVertex upperBound = min(stop, (begin + (chunkId + 1) * chunkSize));
            idSuperArc thisArc = 0;
            if (lowerBound != upperBound) {
               const idVertex pos =  upperBound - 1;
               currentNode        = getNodeFromVertDicho(trunkVerts, scalars_->getSortedVert(pos));
               nextNodeThresold   = trunkVerts[currentNode + 1];
               nextNodeValue      = scalars_->getMirror(nextNodeThresold);
               thisArc            = upArcFromVert(trunkVerts[currentNode]);
            }

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.trunkStats)[chunkId].begin  = _launchGlobalTime.getElapsedTime();
            (*mt_data_.trunkStats)[chunkId].origin = lowerBound;
#endif

            vector<idVertex> regularList;
            if (params_->segm) {
               regularList.reserve(25);
            }

            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const idVertex s = scalars_->getSortedVert(lowerBound + upperBound - 1 - v);
               if (isCorrespondingNull(s)) {
                  const idNode previousNode = currentNode;
                  if (nextNodeValue > scalars_->getMirror(s)) {
                     currentNode      = getNodeFromVertInc(trunkVerts, s, previousNode + 1);
                     nextNodeThresold = trunkVerts[currentNode + 1];
                     nextNodeValue    = scalars_->getMirror(nextNodeThresold);
                     thisArc          = upArcFromVert(trunkVerts[currentNode]);
                  }
                  updateCorrespondingArc(s, thisArc);

                  if (params_->segm) {
                     if (previousNode == currentNode) {
                        regularList.emplace_back(s);
                     } else {
                        // accumulated to have only one atomic update when needed
                        const idSuperArc oldArc = upArcFromVert(trunkVerts[previousNode]);
                        if (regularList.size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
                           {
                              // Swap to avoid a linear copy
                              (*mt_data_.trunkSegments)[oldArc].emplace_back(vector<idVertex>());
                              (*mt_data_.trunkSegments)[oldArc].back().swap(regularList);
                              // This clear should be useless after the swap
                              // regularList.clear();
                           }
                        }
                        // hand.vtu, sequential: 28554
                        regularList.emplace_back(s);
                     }
                  }
               }
            }
            // force increment last arc
            const idNode     baseNode = getCorrespondingNodeId(trunkVerts[currentNode]);
            const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
            if (regularList.size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
               {
                  (*mt_data_.trunkSegments)[upArc].emplace_back(vector<idVertex>());
                  (*mt_data_.trunkSegments)[upArc].back().swap(regularList);
               }
            }
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.trunkStats)[chunkId].end = _launchGlobalTime.getElapsedTime();
#endif
         }
      }
   } else {
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(trunkVerts) PRIOR(isPrior())
#endif
         {
            idNode   currentNode = 0;
            idVertex nextNodeThresold = 0;
            idVertex nextNodeValue    = 0;

            const idVertex lowerBound = begin + chunkId * chunkSize;
            const idVertex upperBound = min(stop, (begin + (chunkId + 1) * chunkSize));
            idSuperArc     thisArc    = 0;
            if (lowerBound != upperBound) {
               const idVertex pos = lowerBound;
               currentNode        = getNodeFromVertDicho(trunkVerts, scalars_->getSortedVert(pos));
               nextNodeThresold   = trunkVerts[currentNode + 1];
               nextNodeValue      = scalars_->getMirror(nextNodeThresold);
               thisArc            = upArcFromVert(trunkVerts[currentNode]);
            }

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.trunkStats)[chunkId].begin  = _launchGlobalTime.getElapsedTime();
            (*mt_data_.trunkStats)[chunkId].origin = lowerBound;
#endif

            vector<idVertex> regularList;
            if (params_->segm) {
               regularList.reserve(25);
            }

            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const idVertex s = scalars_->getSortedVert(v);
               if (isCorrespondingNull(s)) {
                  const idNode previousNode = currentNode;
                  if (nextNodeValue < scalars_->getMirror(s)) {
                     currentNode      = getNodeFromVertInc(trunkVerts, s, previousNode + 1);
                     nextNodeThresold = trunkVerts[currentNode + 1];
                     nextNodeValue    = scalars_->getMirror(nextNodeThresold);
                     thisArc          = upArcFromVert(trunkVerts[currentNode]);
                  }
                  updateCorrespondingArc(s, thisArc);

                  if (params_->segm) {
                     if (previousNode == currentNode) {
                        regularList.emplace_back(s);
                     } else {
                        // accumulated to have only one atomic update when needed
                        const idSuperArc oldArc = upArcFromVert(trunkVerts[previousNode]);
                        if (regularList.size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
                           {
                              // Swap to avoid a linear copy
                              (*mt_data_.trunkSegments)[oldArc].emplace_back(vector<idVertex>());
                              (*mt_data_.trunkSegments)[oldArc].back().swap(regularList);
                              // This clear should be useless after the swap
                              // regularList.clear();
                           }
                        }
                        // hand.vtu, sequential: 28554
                        regularList.emplace_back(s);
                     }
                  }
               }
            }
            // force increment last arc
            const idNode     baseNode = getCorrespondingNodeId(trunkVerts[currentNode]);
            const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
            if (regularList.size()) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
               {
                  (*mt_data_.trunkSegments)[upArc].emplace_back(vector<idVertex>());
                  (*mt_data_.trunkSegments)[upArc].back().swap(regularList);
               }
            }
#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.trunkStats)[chunkId].end = _launchGlobalTime.getElapsedTime();
#endif
         }
      }
   }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif

   // count added
   idVertex tot = 0;
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
   for (const auto& l : *mt_data_.trunkSegments) {
       idVertex arcSize = 0;
       for (const auto& v: l){
          arcSize += v.size();
       }
       tot += arcSize;
   }

   cout << "Trunk vertices: " << tot << endl;
#endif
   return tot;
}

idVertex FTCTree_MT::trunkSegmentation(const vector<idVertex> &trunkVerts,
                                       const idVertex begin,
                                       const idVertex stop)
{
#ifdef TTK_FTC_GRAINSIZE_TRUNK_VERTS
   const int grainSizeTrunkVerts = TTK_FTC_GRAINSIZE_TRUNK_VERTS;
#else
   const int grainSizeTrunkVerts = 50000;
#endif
   const auto sizeBackBone  = abs(stop - begin);
   const auto chunkSize     = getChunkSize(sizeBackBone, scalars_->getSize(), grainSizeTrunkVerts);
   const auto chunkNb       = getChunkCount(sizeBackBone, scalars_->getSize(), grainSizeTrunkVerts);

   // si pas efficace vecteur de la taille de node ici a la place de acc
   idVertex tot             = 0;
   // CAUTION: Duplicate code to avoid comp_.vertLower indirection
   if(isST()){
       for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(trunkVerts, tot) PRIOR(isPrior())
#endif
         {
            idVertex acc              = 0;
            idNode   currentNode      = 0;
            idVertex nextNodeThresold = 0;
            idVertex nextNodeValue    = 0;
            // Warning, this seems not to work without the tasks,
            // may bug if no omp
            const idVertex lowerBound = begin + chunkId * chunkSize;
            const idVertex upperBound = min(stop, (begin + (chunkId + 1) * chunkSize));

            currentNode = getNodeFromVertDicho(trunkVerts, scalars_->getSortedVert(upperBound - 1));
            nextNodeThresold   = trunkVerts[currentNode + 1];
            nextNodeValue      = scalars_->getMirror(nextNodeThresold);
            idSuperArc thisArc = upArcFromVert(trunkVerts[currentNode]);

            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const idVertex s = scalars_->getSortedVert(lowerBound + upperBound - 1 - v);
               if (isCorrespondingNull(s)) {
                  const idNode previousNode = currentNode;
                  // if we overcome the current arc, update
                  if(nextNodeValue > scalars_->getMirror(s)){
                     currentNode      = getNodeFromVertInc(trunkVerts, s, previousNode + 1);
                     nextNodeThresold = trunkVerts[currentNode + 1];
                     nextNodeValue    = scalars_->getMirror(nextNodeThresold);
                     thisArc          = upArcFromVert(trunkVerts[currentNode]);
                  }

                  updateCorrespondingArc(s, thisArc);

                  if (params_->segm) {
                     if (previousNode == currentNode) {
                        ++acc;
                     } else {
                        // accumulated to have only one atomic update when needed
                        const idSuperArc oldArc = upArcFromVert(trunkVerts[previousNode]);
                        getSuperArc(oldArc)->atomicIncVisited(acc);
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                        tot += acc;
#endif
                        acc = 1;
                     }
                  }
               }
            }
            // force increment last arc
            if (acc) {
               const idNode     baseNode = getCorrespondingNodeId(trunkVerts[currentNode]);
               const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
               getSuperArc(upArc)->atomicIncVisited(acc);
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
               tot += acc;
#endif
            }
         }  // end task
      }
  } else {
      for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId) shared(trunkVerts, tot) PRIOR(isPrior())
#endif
         {
            idVertex acc              = 0;
            idNode   currentNode      = 0;
            idVertex nextNodeThresold = 0;
            idVertex nextNodeValue    = 0;
            // Warning, this seems not to work without the tasks,
            // may bug if no omp
            const idVertex lowerBound = begin + chunkId * chunkSize;
            const idVertex upperBound = min(stop, (begin + (chunkId + 1) * chunkSize));

            currentNode        = getNodeFromVertDicho(trunkVerts, scalars_->getSortedVert(lowerBound));
            nextNodeThresold   = trunkVerts[currentNode + 1];
            nextNodeValue      = scalars_->getMirror(nextNodeThresold);
            idSuperArc thisArc = upArcFromVert(trunkVerts[currentNode]);

            for (idVertex v = lowerBound; v < upperBound; ++v) {
               const idVertex s = scalars_->getSortedVert(v);
               if (isCorrespondingNull(s)) {
                  const idNode previousNode = currentNode;
                  // if we overcome the current arc, update
                  if(nextNodeValue < scalars_->getMirror(s)){
                     currentNode      = getNodeFromVertInc(trunkVerts, s, previousNode + 1);
                     nextNodeThresold = trunkVerts[currentNode + 1];
                     nextNodeValue    = scalars_->getMirror(nextNodeThresold);
                     thisArc          = upArcFromVert(trunkVerts[currentNode]);
                  }

                  updateCorrespondingArc(s, thisArc);

                  if (params_->segm) {
                     if (previousNode == currentNode) {
                        ++acc;
                     } else {
                        // accumulated to have only one atomic update when needed
                        const idSuperArc oldArc = upArcFromVert(trunkVerts[previousNode]);
                        getSuperArc(oldArc)->atomicIncVisited(acc);
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                        tot += acc;
#endif
                        acc = 1;
                     }
                  }
               }
            }
            // force increment last arc
            if (acc) {
               const idNode     baseNode = getCorrespondingNodeId(trunkVerts[currentNode]);
               const idSuperArc upArc    = getNode(baseNode)->getUpSuperArcId(0);
               getSuperArc(upArc)->atomicIncVisited(acc);
#ifdef TTK_ENABLE_FTC_TREE_PROCESS_SPEED
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
               tot += acc;
#endif
            }
         }  // end task
      }
   }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
   return tot;
}

ostream &ttk::ftc::operator<<(ostream &o, SuperArc const &a)
{
   o << a.getDownNodeId() << " <>> " << a.getUpNodeId();
   return o;
}

ostream &ttk::ftc::operator<<(ostream &o, Node const &n)
{
   o << n.getNumberOfDownSuperArcs() << " .-. " << n.getNumberOfUpSuperArcs();
   return o;
}
