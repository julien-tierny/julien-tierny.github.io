/// \ingroup base
///
/// \class ttk::ftc::FTCTree_CT
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

#include "FTCTree_CT.h"

#include <iterator>
#include <string>

using namespace ftc;

FTCTree_CT::FTCTree_CT(Params *const params, Triangulation *mesh, Scalars<float> *const scalars)
    : FTCTree_MT(params, mesh, scalars, TreeType::Contour),
      jt_(new FTCTree_MT(params, mesh, scalars, TreeType::Join)),
      st_(new FTCTree_MT(params, mesh, scalars, TreeType::Split))
{
   ct_data_.cloned = false;
}

FTCTree_CT::~FTCTree_CT()
{
   if (jt_) {
      delete jt_;
      jt_ = nullptr;
   }
   if (st_) {
      delete st_;
      st_ = nullptr;
   }
}

void FTCTree_CT::build(TreeType tt)
{

   const bool bothMT = tt == TreeType::Contour || tt == TreeType::Join_Split;

   initComp();

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   setStartTime();
   jt_->setStartTime();
   st_->setStartTime();
#endif

   if(bothMT){
      DebugTimer precomputeTime;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
         {
            leafSearch();
         }
      }

#ifndef TTK_DISABLE_FTC_PRIORITY
      // Set priority
      if(st_->getNumberOfLeaves() < jt_->getNumberOfLeaves())
         st_->setPrior();
      else
         jt_->setPrior();
#endif

      printTime(precomputeTime, "[FTC] leafSearch", -1, 3);
   }

   DebugTimer mergeTreesTime;

   // JT & ST

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
   {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
       {
          if (tt == TreeType::Join || bothMT) {

#ifndef TTK_DISABLE_FTC_MIX_MT
# ifdef TTK_ENABLE_OPENMP
# pragma omp task untied if(threadNumber_ > 1)
# endif
#endif
             jt_->build(tt == TreeType::Contour);
          }
          if (tt == TreeType::Split || bothMT) {
#ifndef TTK_DISABLE_FTC_MIX_MT
# ifdef TTK_ENABLE_OPENMP
# pragma omp task untied if(threadNumber_ > 1)
# endif
#endif
             st_->build(tt == TreeType::Contour);
          }
       }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
  }

   printTime(mergeTreesTime, "[FTC] merge trees ", -1, 3);

   // Combine

   if (tt == TreeType::Contour) {

      DebugTimer combineFullTime;
      insertNodesInMT();
      insertNodesInCT();

      DebugTimer combineTime;
#ifndef TTK_DISABLE_FTC_PARA_COMBINE
      combineParallel();
#else
      combineSequential();
#endif
      printTime(combineTime, "[FTC] combine trees", -1, 4);
      printTime(combineFullTime, "[FTC] combine full", -1, 3);
   }

   // Debug

   if (debugLevel_ > 3) {
      cout << "- [FTC] final number of nodes :";
      switch (tt) {
         case TreeType::Join:
            cout << jt_->getNumberOfNodes();
            break;
         case TreeType::Split:
            cout << st_->getNumberOfNodes();
            break;
         case TreeType::Join_Split:
            cout << jt_->getNumberOfNodes() + st_->getNumberOfNodes();
            break;
         default:
            cout << getNumberOfNodes();
      }
      cout << endl;
   }
}

// Try old combine, be carefull of the old createCTArcSegmentation

#ifndef TTK_DISABLE_FTC_PARA_COMBINE
void FTCTree_CT::combineParallel()
{
   DebugTimer stepTime;
   AtomicVector<pair<bool, idNode>> *growingNodes, *remainingNodes;
   const bool DEBUG = false;

   // Reserve
   mt_data_.superArcs->reset();
   mt_data_.superArcs->reserve(jt_->getNumberOfSuperArcs()+2);

   ct_data_.seenDownVal.clear();
   ct_data_.seenDownVal.resize(getNumberOfNodes(), 0);
   ct_data_.seenUpVal.clear();
   ct_data_.seenUpVal.resize(getNumberOfNodes(), 0);

   growingNodes = new AtomicVector<pair<bool,idNode>>;
   growingNodes->reserve(jt_->getNumberOfLeaves()+st_->getNumberOfLeaves());

   remainingNodes = new AtomicVector<pair<bool,idNode>>;
   remainingNodes->reserve(jt_->getNumberOfLeaves()+st_->getNumberOfLeaves());

   // Add JT & ST Leaves

   // Add leves to remaining nodes
   const auto& nbSTLeaves = st_->getNumberOfLeaves();
   if (nbSTLeaves > 1) {
      for (idNode n = 0; n < nbSTLeaves; ++n) {
         const idNode   nId    = st_->getLeave(n);
         const idVertex id     = remainingNodes->getNext();
         (*remainingNodes)[id] = {false, nId};
      }
   } else {
      ct_data_.cloned = true;
      move(jt_);
      // avoid dual delete, ugly
      delete growingNodes;
      delete remainingNodes;
      return;
   }

   // count how many leaves can be added, if more than one : ok!
   const auto& nbJTLeaves = jt_->getNumberOfLeaves();
   for (idNode n = 0; n < nbJTLeaves; ++n) {
      const auto &   nId    = jt_->getLeave(n);
      const idVertex id     = remainingNodes->getNext();
      (*remainingNodes)[id] = {true, nId};
   }

   if (DEBUG) {
      cout << "leaves for combination : " << remainingNodes->size()
           << " in : " << stepTime.getElapsedTime() << endl;
   }

   // Warning, have a reserve here, can't make it at the begnining, need build output
   mt_data_.leaves->reserve(jt_->getLeaves().size() + st_->getLeaves().size());
   mt_data_.superArcs->reserve(jt_->getNumberOfSuperArcs()+st_->getNumberOfSuperArcs());
   mt_data_.nodes->reserve(jt_->getNumberOfNodes());

   if (remainingNodes->empty()) {
      cout << "[FTCTree_CT::combine ] Nothing to combine" << endl;
   }

      // Here we avoid the use of hyperthreads, which is nefast for this step
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
   {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
      {
         do {
            if (!remainingNodes->empty()) {
               auto *tmpNodes = growingNodes;
               growingNodes   = remainingNodes;
               remainingNodes = tmpNodes;
            }

            idNode curNbLeaves = growingNodes->size();

            // Debug print
            if (DEBUG) {
               cout << endl;
               for (idNode n = 0; n < curNbLeaves; ++n) {
                  idNode currentNodeId;
                  bool   isJT;
                  tie(isJT, currentNodeId) = (*growingNodes)[n];
                  FTCTree_MT *xt           = (isJT) ? jt_ : st_;
                  cout << isJT << ":" << xt->getNode(currentNodeId)->getVertexId() << " ";
               }
               cout << endl;
            }

#ifdef TTK_DISABLE_FTC_TRUNK_COMBINE
#define TRUNK_COMBINE false
#else
#define TRUNK_COMBINE true
#endif

            if (TRUNK_COMBINE and threadNumber_ > 1 and curNbLeaves == 2) {
               remainingNodes->clear();
               printTime(stepTime, "[FTC] parallel combine", -1, 3);
               DebugTimer trunkTime;
               trunkCombine();
               printTime(trunkTime, "[FTC] trunk combine", -1, 3);
            } else {
               idNode curLeaf = 0;

               while(curLeaf < curNbLeaves){
                  idVertex accSegmentation = 0;
#ifdef TTK_FTC_GRAINSIZE_COMBINE_MAIN
                  const idVertex grainSizeCombineMain = TTK_FTC_GRAINSIZE_COMBINE_MAIN;
#else
                  const idVertex grainSizeCombineMain = 10000;
#endif

                  idNode lowerBound = curLeaf;

                  // get leaves to have at least grainSizeCombineMain vertices to process

                  do {
                     // we have at leat one curLeaf
                     idNode currentNodeId;
                     bool   isJT;
                     tie(isJT, currentNodeId)   = (*growingNodes)[curLeaf];
                     FTCTree_MT *     xt        = (isJT) ? jt_ : st_;
                     const idSuperArc leafArcId = xt->getNode(currentNodeId)->getUpSuperArcId(0);
                     accSegmentation += xt->getSuperArc(leafArcId)->getNbVertSeen();
                     ++curLeaf;
                  } while (curLeaf < curNbLeaves && accSegmentation < grainSizeCombineMain);

                  idNode upperBound = curLeaf;

                  // process these leaves in a task
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(lowerBound, upperBound) if(accSegmentation > 1000)
#endif
                  {

                     for (idNode n = lowerBound; n < upperBound; ++n) {
                        idNode currentNodeId;
                        bool   isJT;

                        // INFO QUEUE

                        tie(isJT, currentNodeId) = (*growingNodes)[n];

                        FTCTree_MT *xt = (isJT) ? jt_ : st_;
                        FTCTree_MT *yt = (isJT) ? st_ : jt_;

                        // INFO JT / ST

                        // i <- Get(Q)
                        const Node *currentNode = xt->getNode(currentNodeId);

                        if (DEBUG) {
                           if (xt == jt_)
                              cout << endl << "JT ";
                           else
                              cout << endl << "ST ";
                           cout << "node : " << currentNode->getVertexId() << endl;
                        }

                        // "choose a non-root leaf that is not a split in ST" so we ignore such
                        // nodes
                        if (currentNode->getNumberOfUpSuperArcs() == 0) {
                           if (DEBUG) {
                              cout << " ignore already processed" << endl;
                           }
                           continue;
                        }

                        idNode correspondingNodeId =
                            yt->getCorrespondingNodeId(currentNode->getVertexId());

                        if (yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {
                           if (DEBUG) {
                              cout << "put remain:" << isJT << "::" << xt->printNode(currentNodeId)
                                   << endl;
                              cout << " which is in yt : " << yt->printNode(correspondingNodeId)
                                   << endl;
                           }
                           const idVertex id     = remainingNodes->getNext();
                           (*remainingNodes)[id] = {isJT, currentNodeId};
                           continue;
                        }

                        // NODES IN CT

                        const idVertex curVert = currentNode->getVertexId();
                        const idNode   node1   = getCorrespondingNodeId(curVert);

                        // j <- GetAdj(XT, i)
                        idSuperArc  curUpArc   = currentNode->getUpSuperArcId(0);
                        idNode      parentId   = xt->getSuperArc(curUpArc)->getUpNodeId();
                        const Node *parentNode = xt->getNode(parentId);

                        if (DEBUG) {
                           cout << " parent node :" << parentNode->getVertexId() << endl;
                        }

                        idVertex     parVert = parentNode->getVertexId();
                        const idNode node2   = getCorrespondingNodeId(parVert);

                        // CREATE ARC

                        idSuperArc processArc = currentNode->getUpSuperArcId(0);

                        // create the arc in in the good direction
                        // and add it to crossing if needed
                        idSuperArc createdArc;
                        if (scalars_->isLower(
                                currentNode->getVertexId(),
                                parentNode->getVertexId())) {  // take care of the order
                           createdArc = makeSuperArc(node1, node2, false);
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                           ++ct_data_.seenUpVal[node1];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                           ++ct_data_.seenDownVal[node2];
                        } else {
                           createdArc = makeSuperArc(node2, node1, false);
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                           ++ct_data_.seenUpVal[node2];
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
                           ++ct_data_.seenDownVal[node1];
                        }

                        // replace the down node visit
                        getSuperArc(createdArc)->atomicIncVisited();

                        // Segmentation
                        if (params_->segm) {
                           createCTArcSegmentation(createdArc, isJT, processArc);
                        }

                        if (DEBUG) {
                           cout << "create arc : " << printArc(createdArc) << endl;
                        }

                        // DEL NODES

                        // DelNode(XT, i)
                        {
                           if (DEBUG) {
                              cout << curVert << " delete xt (" << (xt == jt_) << ") ";
                              cout << "node :" << xt->printNode(currentNodeId) << endl;
                           }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
                           {
                              xt->delNode(currentNodeId);
                           }
                        }

                        // DelNode(YT, i)
                        {
                           if (DEBUG) {
                              cout << curVert << " delete yt (" << isJT << ") node :";
                              cout << yt->printNode(correspondingNodeId) << endl;
                           }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
                           {
                              yt->delNode(correspondingNodeId);
                           }
                        }

                        // PROCESS QUEUE

                        if (parentNode->getNumberOfDownSuperArcs() == 0 &&
                            parentNode->getNumberOfUpSuperArcs()) {

                           idVertex oldVert2tree;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
                           {
                              oldVert2tree = (*xt->mt_data_.vert2tree)[parVert];
                              (*xt->mt_data_.vert2tree)[parVert] = nullVertex;
                           }

                           // add only once
                           if (oldVert2tree != nullVertex) {
                              const idVertex id     = remainingNodes->getNext();
                              (*remainingNodes)[id] = {isJT, parentId};

                              if (DEBUG) {
                                 cout << curVert << " will see : " << parentNode->getVertexId()
                                      << endl;
                              }
                           }
                        }
                     }
                  } // end task
               } // end while launching tasks
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
            }

            growingNodes->reset();
         } while (!remainingNodes->empty());
      }
   }

   if (DEBUG) {
       printTree2();
   }

   delete growingNodes;
   delete remainingNodes;
   printTime(stepTime, "[FTC] Combine end", -1, 3);
}

#else

int FTCTree_CT::combineSequential()
{
   DebugTimer stepTime;
   queue<pair<bool, idNode>> growingNodes, remainingNodes;

   const bool DEBUG = false;

   // Reserve
   mt_data_.nodes->reserve(jt_->getNumberOfNodes());
   mt_data_.superArcs->reserve(jt_->getNumberOfSuperArcs()+2);
   mt_data_.leaves->reserve(jt_->getNumberOfLeaves()+st_->getNumberOfLeaves());

   // Add JT & ST Leaves to growingNodes

   // Add leves to growing nodes
   const auto& nbSTLeaves = st_->getNumberOfLeaves();
   if (nbSTLeaves > 1) {
      for (idNode n = 0; n < nbSTLeaves; ++n) {
         const auto &nId = st_->getLeave(n);
         growingNodes.emplace(false, nId);
      }
   } else {
      move(jt_);
      return 0;
   }

   // count how many leaves can be added, if more than one : ok!
   const auto& nbJTLeaves = jt_->getNumberOfLeaves();
   if (nbJTLeaves > 1) {
      for (idNode n = 0; n < nbJTLeaves; ++n) {
         const auto &nId = jt_->getLeave(n);
         growingNodes.emplace(true, nId);
      }
   } // else can't clone, not same up and down

   if (DEBUG) {
      cout << "growingNodes : " << growingNodes.size() << " in : " << stepTime.getElapsedTime()
           << endl;
   }

   // Warning, have a reserve here, can't make it at the begnining, need build output
   mt_data_.leaves->reserve(jt_->getLeaves().size() + st_->getLeaves().size());
   mt_data_.superArcs->reserve(jt_->getNumberOfSuperArcs());
   mt_data_.nodes->reserve(jt_->getNumberOfNodes());

   if (growingNodes.empty()) {
      cout << "[FTCTree_CT::combine ] Nothing to combine" << endl;
   }

#ifdef TTK_ENABLE_FTC_TREE_DUAL_QUEUE_COMBINE
   do {
      while (!remainingNodes.empty()) {
         bool       isJT;
         idNode     currentNodeId;
         FTCTree_MT *xt;

         tie(isJT, currentNodeId) = remainingNodes.front();
         remainingNodes.pop();
         if (isJT) {
            // node come frome jt
            xt = jt_;
         } else {
            // node come from st
            xt = st_;
         }
         if (xt->getNode(currentNodeId)->getNumberOfUpSuperArcs() == 1) {
            growingNodes.emplace(isJT, currentNodeId);
            if (DEBUG) {
               cout << "repush in growing:" << isJT << "::" << xt->printNode(currentNodeId) << endl;
            }
         }
      }
#endif

      while (!growingNodes.empty()) {
         idNode currentNodeId;
         bool   isJT;

         // INFO QUEUE

         tie(isJT, currentNodeId) = growingNodes.front();
         growingNodes.pop();

         FTCTree_MT *xt = (isJT) ? jt_ : st_;
         FTCTree_MT *yt = (isJT) ? st_ : jt_;

         // INFO JT / ST

         // i <- Get(Q)
         const Node *currentNode = xt->getNode(currentNodeId);

         if (DEBUG) {
            if (xt == jt_)
               cout << endl << "JT ";
            else
               cout << endl << "ST ";
            cout << "node : " << currentNode->getVertexId() << endl;
         }

         // "choose a non-root leaf that is not a split in ST" so we ignore such nodes
         if (currentNode->getNumberOfUpSuperArcs() == 0) {
            if (DEBUG) {
               cout << " ignore already processed" << endl;
            }
            continue;
         }

         idNode correspondingNodeId = yt->getCorrespondingNodeId(currentNode->getVertexId());

         if (yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {
            if (DEBUG) {
               cout << "put remain:" << isJT << "::" << xt->printNode(currentNodeId) << endl;
               cout << " which is in yt : " << yt->printNode(correspondingNodeId) << endl;
            }
#ifdef TTK_ENABLE_FTC_TREE_DUAL_QUEUE_COMBINE
            remainingNodes.emplace(isJT, currentNodeId);
#else
            growingNodes.emplace(isJT, currentNodeId);
#endif
            continue;
         }

         // NODES IN CT

         idNode   node1, node2;
         idVertex curVert = currentNode->getVertexId();
         // NODE1
         if (isCorrespondingNode(curVert)) {
            // already a node in the tree
            node1 = getCorrespondingNodeId(curVert);
         } else {
            // create a new node
            node1 = makeNode(currentNode);

            // check if leaf
            if (!currentNode->getNumberOfDownSuperArcs())
               mt_data_.leaves->emplace_back(node1);
            else if (!currentNode->getNumberOfUpSuperArcs())
               mt_data_.leaves->emplace_back(node1);
         }

         // j <- GetAdj(XT, i)
         idSuperArc curUpArc = currentNode->getUpSuperArcId(0);
         idNode      parentId   = xt->getSuperArc(curUpArc)->getUpNodeId();
         const Node *parentNode = xt->getNode(parentId);

         if (DEBUG) {
            cout << " parent node :" << parentNode->getVertexId() << endl;
         }

         idVertex parVert = parentNode->getVertexId();
         // NODE2
         if (isCorrespondingNode(parVert)) {
            // already a node in the tree
            node2 = getCorrespondingNodeId(parVert);
         } else {
            // create a new node
            node2 = makeNode(parentNode);
            if (!parentNode->getNumberOfUpSuperArcs())
               mt_data_.leaves->emplace_back(node2);
         }

         // CREATE ARC

         idSuperArc processArc = currentNode->getUpSuperArcId(0);

         // create the arc in in the good direction
         // and add it to crossing if needed
         idSuperArc createdArc;
         if (scalars_->isLower(currentNode->getVertexId(),
                               parentNode->getVertexId())) {  // take care of the order
            createdArc = makeSuperArc(node1, node2);
         } else {
            createdArc = makeSuperArc(node2, node1);
         }

         // Segmentation
         if (params_->segm) {
            createCTArcSegmentation(createdArc, isJT, processArc);
         }

         if (DEBUG) {
            cout << "create arc : " << printArc(createdArc) << endl;
         }

         // DEL NODES

         // DelNode(XT, i)
         {
            if (DEBUG) {
               cout << " delete xt (" << (xt == jt_) << ") ";
               cout << "node :" << xt->printNode(currentNodeId) << endl;
            }

            xt->delNode(currentNodeId);
         }

         // DelNode(YT, i)
         {
            if (DEBUG) {
               cout << " delete yt (" << isJT << ") node :";
               cout << yt->printNode(correspondingNodeId) << endl;
            }

            yt->delNode(correspondingNodeId);
         }


         // PROCESS QUEUE

         if (parentNode->getNumberOfDownSuperArcs() == 0 && parentNode->getNumberOfUpSuperArcs()) {
            growingNodes.emplace(isJT, parentId);

            if (DEBUG) {
               cout << "will see : " << parentNode->getVertexId() << endl;
            }
         }

      }
#ifdef TTK_ENABLE_FTC_TREE_DUAL_QUEUE_COMBINE
   } while (!remainingNodes.empty());
#endif

   if (DEBUG) {
       printTree2();
   }

   return 0;
}

#endif

#ifndef TTK_DISABLE_FTC_PARA_COMBINE

void FTCTree_CT::createCTArcSegmentation(idSuperArc ctArc, const bool isJT, idSuperArc xtArc)
{
   const FTCTree_MT *xt = (isJT) ? jt_ : st_;

   const list<Region> &xtRegions = xt->getSuperArc(xtArc)->getRegions();
   for (const Region &reg : xtRegions) {
      segm_it        cur        = reg.segmentBegin;
      segm_it        end        = reg.segmentEnd;
      const idVertex regionSize = end - cur;
#ifdef TTK_FTC_GRAINSIZE_COMBINE_MAIN
      const idVertex grainSizeCombineMain = TTK_FTC_GRAINSIZE_COMBINE_MAIN;
#else
      const idVertex grainSizeCombineMain = 10000;
#endif

      // Parallel process only if large enought
      // Sequential otherwise
      if (regionSize > grainSizeCombineMain) {
         const idVertex segmChunkSize =
             getChunkSize(regionSize, scalars_->getSize(), grainSizeCombineMain);
         const idVertex segmChunkNb =
             getChunkCount(regionSize, scalars_->getSize(), grainSizeCombineMain);

         for (idVertex segmChunkId = 0; segmChunkId < segmChunkNb; ++segmChunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(segmChunkId)
#endif
            {
               auto       curVert    = cur + segmChunkId * segmChunkSize;
               const auto upperBound = min(end, cur + (segmChunkId + 1) * segmChunkSize);
               for (; curVert != upperBound; ++curVert) {
                  if (isCorrespondingNull(*curVert)) {
                     updateCorrespondingArc(*curVert, ctArc);
                     getSuperArc(ctArc)->atomicIncVisited();
                  }
               }
            }
         }
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
      } else {
         for (; cur != end; ++cur) {
            if (isCorrespondingNull(*cur)) {
               updateCorrespondingArc(*cur, ctArc);
               getSuperArc(ctArc)->incVisited();
            }
         }
      }
   }
}

#else

void FTCTree_CT::createCTArcSegmentation(idSuperArc ctArc, const bool isJT, idSuperArc xtArc)
{
   const FTCTree_MT *xt = (isJT) ? jt_ : st_;

   /*Here we prefere to create lots of small region, each arc having its own segmentation with
    * no overlap instead of having a same vertice in several arc and using vert2tree to decide
    * because we do not want to maintain vert2tree information during the whole computation*/
   const list<Region> &xtRegions = xt->getSuperArc(xtArc)->getRegions();
   for (const Region &reg : xtRegions) {
      segm_it cur    = reg.segmentBegin;
      segm_it end    = reg.segmentEnd;
      segm_it tmpBeg = reg.segmentBegin;
      // each element inside this region
      for (; cur != end; ++cur) {
         if (isCorrespondingNull(*cur)) {
            updateCorrespondingArc(*cur, ctArc);
         } else {
            // already set, we finish a region
            if (cur != tmpBeg) {
               getSuperArc(ctArc)->concat(tmpBeg, cur);
            }
            // if several contiguous vertices are discarded
            // cur will be equals to tmpBeg and we will not create empty regions
            tmpBeg = cur + 1;
         }
      }
      // close last region
      if (cur != tmpBeg) {
         getSuperArc(ctArc)->concat(tmpBeg, cur);
      }
   }
}

#endif

void FTCTree_CT::finalizeSegmentation(void)
{
   DebugTimer  finSegmTime;
   const auto &nbArc = getNumberOfSuperArcs();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
   for (idSuperArc i = 0; i < nbArc; i++) {
      getSuperArc(i)->createSegmentation(scalars_);
   }

   printTime(finSegmTime, "[FTC] post-process segm", -1, 4);
}

void FTCTree_CT::insertNodesInCT(void)
{
   const idNode nbNodes = jt_->getNumberOfNodes();
   mt_data_.nodes->reserve(nbNodes);
   mt_data_.leaves->reserve(jt_->getNumberOfLeaves()+st_->getNumberOfLeaves());
   ct_data_.initDownVal.resize(nbNodes, 0);
   ct_data_.initUpVal.resize(nbNodes, 0);

   for (idNode nid = 0; nid < nbNodes; nid++) {
       Node* jnode = jt_->getNode(nid);
       const idVertex v = jnode->getVertexId();
       const idNode nodeid = makeNode(v);

       // valence
       const valence jtVal = jnode->getNumberOfDownSuperArcs();
       const valence stVal = st_->getNode(st_->getCorrespondingNodeId(v))->getNumberOfDownSuperArcs();

       getNode(nodeid)->reserveDownArcs(jtVal);
       getNode(nodeid)->reserveUpArcs(stVal);

       //leaves
       if (!jnode->getNumberOfDownSuperArcs() || !jnode->getNumberOfUpSuperArcs()) {
          mt_data_.leaves->push_back(nodeid);
          ct_data_.initDownVal[nodeid] = nullValence;
          ct_data_.initUpVal[nodeid]   = nullValence;
       } else {
          ct_data_.initDownVal[nodeid] = jtVal;
          ct_data_.initUpVal[nodeid]   = stVal;
       }
   }
}

void FTCTree_CT::insertNodesInMT(void)
{
   vector<idNode> sortedJTNodes = jt_->sortedNodes(true);
   vector<idNode> sortedSTNodes = st_->sortedNodes(true);

#ifdef TTK_ENABLE_OPENMP
#pragma omp task
#endif
   for (const idNode& t : sortedSTNodes) {

      idVertex vertId = st_->getNode(t)->getVertexId();
      if (jt_->isCorrespondingNode(vertId)) {
          continue;
      }
      jt_->insertNode(st_->getNode(t));
   }

#ifdef TTK_ENABLE_OPENMP
#pragma omp task
#endif
   for (const idNode& t : sortedJTNodes) {

      idVertex vertId = jt_->getNode(t)->getVertexId();
      if (st_->isCorrespondingNode(vertId)) {
          continue;
      }
      st_->insertNode(jt_->getNode(t));
   }

#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
#endif
}

int FTCTree_CT::leafSearch()
{
   const auto  nbScalars = scalars_->getSize();
   const auto  chunkSize = getChunkSize();
   const auto  chunkNb   = getChunkCount();

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
      mt_data_.extrStats->resize(chunkNb);
#endif

   // Extrema extract and launch tasks
   for (idVertex chunkId = 0; chunkId < chunkNb; ++chunkId) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(chunkId)
#endif
      {
         const idVertex lowerBound = chunkId * chunkSize;
         const idVertex upperBound = min(nbScalars, (chunkId + 1) * chunkSize);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
            (*mt_data_.extrStats)[chunkId].begin = _launchGlobalTime.getElapsedTime();
#endif

         for (idVertex v = lowerBound; v < upperBound; ++v) {
            const auto &neighNumb   = mesh_->getVertexNeighborNumber(v);
            bool        upval       = false;
            bool        downval     = false;

            for (valence n = 0; n < neighNumb; ++n) {
               idVertex neigh;
               mesh_->getVertexNeighbor(v, n, neigh);
               bool compareNV = scalars_->isLower(neigh, v);
               if (!downval && compareNV) {
                  downval = true;
               }
               if (!upval && !compareNV) {
                  upval = true;
               }
               if(upval && downval)
                  break;
            }

            if (!downval) {
               const idNode nextLeaf             = jt_->mt_data_.leaves->getNext();
               (*jt_->mt_data_.leaves)[nextLeaf] = jt_->makeNode(v);
            }

            if (!upval) {
               const idNode nextLeaf             = st_->mt_data_.leaves->getNext();
               (*st_->mt_data_.leaves)[nextLeaf] = st_->makeNode(v);
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

   for (const idNode n : *jt_->mt_data_.leaves) {
      const idVertex v        = jt_->getNode(n)->getVertexId();
      (*jt_->mt_data_.ufs)[v] = new AtomicUF(v);
   }
   jt_->mt_data_.activeTasks += jt_->getNumberOfLeaves();

   for (const idNode n : *st_->mt_data_.leaves) {
      const idVertex v        = st_->getNode(n)->getVertexId();
      (*st_->mt_data_.ufs)[v] = new AtomicUF(v);
   }
   st_->mt_data_.activeTasks += st_->getNumberOfLeaves();

   return 0;
}

void FTCTree_CT::postProcessCombine()
{
   const idSuperArc nbSa = getNumberOfSuperArcs();
   //TODO parallel test here
   for (idSuperArc said = 0; said < nbSa; ++said) {
      const SuperArc *sa       = getSuperArc(said);
      const idNode    downNode = sa->getDownNodeId();
      const idNode    upNode   = sa->getUpNodeId();
      // we do not want to insert trunk arc in duplicate
      // Trick, these arc can only be in the first position
      // of the node's arc list
      if (!getNode(downNode)->getNumberOfUpSuperArcs() ||
          getNode(downNode)->getUpSuperArcId(0) != said) {
         getNode(downNode)->addUpSuperArcId(said);
      }
      if (!getNode(upNode)->getNumberOfDownSuperArcs() ||
          getNode(upNode)->getDownSuperArcId(0) != said) {
         getNode(upNode)->addDownSuperArcId(said);
      }
   }
}


void FTCTree_CT::trunkCombine(void)
{
   DebugTimer stepTime;
   const idNode nbNodes = getNumberOfNodes();
   vector<idVertex> trunkVerts;
   trunkVerts.reserve(nbNodes / 2);

   for (idNode i = 0; i < nbNodes; ++i) {
      const valence initDownVal = ct_data_.initDownVal[i];
      const valence curDownVal  = ct_data_.seenDownVal[i];
      const valence initUpVal   = ct_data_.initUpVal[i];
      const valence curUpVal    = ct_data_.seenUpVal[i];

      if ((initDownVal != nullValence && curDownVal != initDownVal) ||
          (initUpVal != nullValence && curUpVal != initUpVal)) {
         trunkVerts.emplace_back(getNode(i)->getVertexId());
      }
   }

   if(!trunkVerts.size()) return;

   sort(trunkVerts.begin(), trunkVerts.end(), comp_.vertLower);

   const idNode nbNodesTrunk = trunkVerts.size();

   for (idNode i = 1; i < nbNodesTrunk; ++i) {
      const idSuperArc aid = makeSuperArc(getCorrespondingNodeId(trunkVerts[i - 1]),
                                          getCorrespondingNodeId(trunkVerts[i]));
      getSuperArc(aid)->incVisited();
   }

   printTime(stepTime, "[FTC] Combine : trunk skeleton", -1, 4);

   if (params_->segm) {
      idVertex beginTrunk;
      idVertex stopTrunk;
      tie(beginTrunk, stopTrunk) = getBoundsFromVerts(trunkVerts);

      // Add fake last vertex to the trunkVerts to avoid size checking
      trunkVerts.emplace_back(scalars_->getSortedVert(scalars_->getSize() - 1));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
#endif
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp single nowait
#endif
         {
            trunkSegmentation(trunkVerts, beginTrunk, stopTrunk);
            // already have a taskwait
         }
      }
   }

   printTime(stepTime, "[FTC] Combine : trunk end", -1, 4);
}
