/// \ingroup base
//
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-02-10
///
///\brief TTK Shared structures and classes for FTC Tree.

#pragma once

#include "AtomicVector.h"
#include "DataTypes.h"

#include <iostream>

namespace ttk
{
namespace ftc
{
   // Compute parameters (global)
   struct Params {
      TreeType treeType;
      bool     segm        = true;
      bool     normalize   = true;
      bool     advStats    = true;
      int      samplingLvl = 0;

      void printSelf(int debugLvl, idThread threadNumber)
      {
         if (debugLvl > 1) {
            if (debugLvl > 2) {
               std::cout << "------------" << std::endl;
            }
            std::cout << "[FTC] number of threads : " << (int)threadNumber << std::endl;
            if (debugLvl > 2) {
               std::cout << "* debug lvl  : " << debugLvl << std::endl;
               std::cout << "* tree type  : ";
               if (treeType == TreeType::Contour) {
                  std::cout << "Contour";
               } else if (treeType == TreeType::Join) {
                  std::cout << "Join";
               } else if (treeType == TreeType::Split) {
                  std::cout << "Split";
               } else if (treeType == TreeType::Join_Split) {
                  std::cout << "Join + Split";
               }
               std::cout << std::endl;
               std::cout << "------------" << std::endl;
            }
         }
      }
   };

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   struct ActiveTask {
      float    begin  = -1;
      float    end    = -1;
      idVertex origin = nullVertex;
      idVertex remain = -1;
   };
#endif

}
}

