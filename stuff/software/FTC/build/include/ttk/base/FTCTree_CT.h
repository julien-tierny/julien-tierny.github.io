/// \ingroup base
/// \class ttk::ftc::FTCTree_CT
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef FTCTREE_CT_H
#define FTCTREE_CT_H

#include "DataTypes.h"
#include "FTCTree_MT.h"

#include <queue>
#include <set>

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk
{
namespace ftc
{
   struct CT_Data {
      // Initial valence of each node
      vector<valence> initDownVal, initUpVal;
      // Number of time the first step of combineParallel touch this node
      vector<valence> seenDownVal, seenUpVal;
      bool            cloned;
   };

   class FTCTree_CT : public FTCTree_MT
   {
     protected:
      FTCTree_MT *jt_, *st_;
      CT_Data     ct_data_;

     public:
      // -----------------
      // Constructors
      // -----------------

      FTCTree_CT(Params* const params, Triangulation* mesh, Scalars<float>* const scalars);
      virtual ~FTCTree_CT();

      // -----------------
      // ACCESSOR
      // -----------------

      inline FTCTree_MT* getJoinTree(void) const
      {
         return jt_;
      }

      inline FTCTree_MT* getSplitTree(void) const
      {
         return st_;
      }

      inline FTCTree_MT* getTree(const TreeType tt)
      {
         switch (tt) {
            case TreeType::Split:
               return getSplitTree();
               break;
            case TreeType::Join:
               return getJoinTree();
               break;
            case TreeType::Contour:
               return this;
               break;
            default:
               return this;
               break;
         }
         return this;
      }

      inline void setupTriangulation(Triangulation* m, const bool preproc = true)
      {
         FTCTree_MT::setupTriangulation(m, preproc);
         jt_->setupTriangulation(m, false);
         st_->setupTriangulation(m, false);
      }

      inline int setDebugLevel(const int d)
      {
         Debug::setDebugLevel(d);
         jt_->setDebugLevel(d);
         st_->setDebugLevel(d);
         return 0;
      }

      inline int setThreadNumber(const int n)
      {
         Debug::setThreadNumber(n);
         jt_->setThreadNumber(n);
         st_->setThreadNumber(n);
         return 0;
      }

      // -----------------
      // PROCESS
      // -----------------

      virtual void alloc()
      {
         switch (params_->treeType) {
            case TreeType::Contour:
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections if (threadNumber_ > 1)
#endif
            {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  FTCTree_MT::alloc();
               }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  jt_->alloc(false);
               }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  st_->alloc(false);
               }
            }
            break;

            case TreeType::Join_Split:
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections if (threadNumber_ > 1)
#endif
            {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  jt_->alloc();
               }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  st_->alloc(false);
               }
            }
            break;

            case TreeType::Join:
            jt_->alloc();
            break;

            case TreeType::Split:
            st_->alloc();
            break;

            default:
            cerr << "[FTCTree_CT]: error unsupported tree type for allocation" << endl;
            break;
         }
      }

      virtual void init()
      {
         switch (params_->treeType){

            case TreeType::Contour:
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections if (threadNumber_ > 1)
#endif
            {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  FTCTree_MT::init();
               }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  jt_->init(false);
               }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  st_->init(false);
               }
            }
            break;

            case TreeType::Join_Split:
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel sections if (threadNumber_ > 1)
#endif
            {
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  jt_->init();
               }
#ifdef TTK_ENABLE_OPENMP
#pragma omp section
#endif
               {
                  st_->init(false);
               }
            }
            break;

            case TreeType::Join:
            jt_->init();
            break;

            case TreeType::Split:
            st_->init();
            break;

            default:
            cerr << "[FTCTree_CT]: error unsupported tree type for initalization" << endl;
            break;
         }
      }

      virtual void build(TreeType tt);

      virtual int leafSearch();

      // add missing node from JT in ST and vice versa
      void insertNodesInMT();

      // create all the nodes of the CT using MT
      void insertNodesInCT();

#ifndef TTK_DISABLE_FTC_PARA_COMBINE
      void combineParallel();
#else
      int combineSequential();
#endif

      void createCTArcSegmentation(idSuperArc ctArc, const bool isJT, idSuperArc xtArc);

      void finalizeSegmentation(void);

      void postProcessCombine();

      void trunkCombine(void);
   };

}
}

#endif  // CONTOURTREE_H
