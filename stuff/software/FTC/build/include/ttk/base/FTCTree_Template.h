/// \ingroup base
/// \class ttk::ftc::FTCTree
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
///
/// \sa vtkFTCTree.cpp %for a usage example.

#ifndef FTCTREE_TPL_H
#define FTCTREE_TPL_H

#include "FTCTree.h"

#include <type_traits>


namespace ttk
{
namespace ftc
{
   // -------
   // PROCESS
   // -------

   template <typename scalarType>
   int ftc::FTCTree::build(void)
   {

      if (!std::is_same<scalarType, float>::value) {
         std::cerr << "FTC ERROR: Scalar field must be of type float" << std::endl;
         return 1;
      }

      // -----
      // INPUT
      // -----

      params_->printSelf(debugLevel_, threadNumber_);

#ifdef TTK_ENABLE_OPENMP
      omp_set_num_threads(threadNumber_);
      omp_set_nested(1);
#endif

      // ----
      // INIT
      // ----

      setDebugLevel(debugLevel_);
      scalars_->setSize(mesh_->getNumberOfVertices());

      // Alloc / reserve
      DebugTimer initTime;
      alloc();
      printTime(initTime, "[FTC] alloc", -1, 3);

#ifdef TTK_ENABLE_FTC_COUNT_INIT
      DebugTimer startTime;
#endif

      // init values
      DebugTimer setTimer;
      scalars_->removeNaN();
      init();
      printTime(setTimer, "[FTC] init", -1, 3);

#ifndef TTK_ENABLE_FTC_COUNT_INIT
      DebugTimer startTime;
#endif

#ifdef withStatsTime
      getJoinTree()->setStartTime();
      getSplitTree()->setStartTime();
#endif

      // for fast comparison
      // and regions / segmentation
      DebugTimer sortTime;
      scalars_->sort();
      printTime(sortTime, "[FTC] sort step", -1, 3);

      // -----
      // BUILD
      // -----

      DebugTimer buildTime;
      FTCTree_CT::build(params_->treeType);
      printTime(buildTime, "[FTC] build tree", -1, 3);

      printTime(startTime, "[FTC] Total ", -1, 1);

#ifdef TTK_ENABLE_FTC_EARLY_EXIT
      cout << "early exit" << endl;
      exit(0);
#endif

      // Build the list of regular vertices of the arc
      if (params_->segm) {
         switch (params_->treeType) {
            case TreeType::Join:
               getJoinTree()->buildSegmentation();
               getJoinTree()->finalizeSegmentation();
               break;
            case TreeType::Split:
               getSplitTree()->buildSegmentation();
               getSplitTree()->finalizeSegmentation();
               break;
            case TreeType::Join_Split:
               getJoinTree()->buildSegmentation();
               getSplitTree()->buildSegmentation();
               getJoinTree()->finalizeSegmentation();
               getSplitTree()->finalizeSegmentation();
               break;
            case TreeType::Contour:
               if (!ct_data_.cloned) {
#ifndef TTK_DISABLE_FTC_PARA_COMBINE
                  postProcessCombine();
                  buildSegmentation();
#endif
               }
               finalizeSegmentation();
               break;
            default:
               break;
         }
      }

      // Normalization
      if (params_->normalize) {
         switch (params_->treeType) {
            case TreeType::Join:
               getJoinTree()->normalizeIds();
               break;
            case TreeType::Split:
               getSplitTree()->normalizeIds();
               break;
            case TreeType::Join_Split:
               getJoinTree()->normalizeIds();
               getSplitTree()->normalizeIds();
               break;
            case TreeType::Contour:
               normalizeIds();
               break;
            default:
               break;
         }
      }

      if (debugLevel_ > 4) {
         switch (params_->treeType) {
            case TreeType::Join:
               jt_->printTree2();
               break;
            case TreeType::Split:
               st_->printTree2();
               break;
            case TreeType::Join_Split:
               jt_->printTree2();
               st_->printTree2();
               break;
            default:
               printTree2();
         }
      }
      return 0;
   }
}
}

#endif /* end of include guard: FTCTREE_TPL_H */
