/// \ingroup base
///
/// \class ttk::ftc::FTCTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date December 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa vtkFTCTree.cpp %for a usage example.

#ifndef FTCTREE_H
#define FTCTREE_H

#include "FTCTree_CT.h"
#include "DataTypes.h"

// base code includes
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

namespace ttk
{
namespace ftc
{

   class FTCTree : public FTCTree_CT
   {
     public:

      // -----------------
      // CONSTRUCTORS
      // -----------------

      FTCTree();
      virtual ~FTCTree();

      // -------
      // PROCESS
      // -------

      // Initialize structures then build tree
      // Need triangulation, scalars and all params set before call
      template <typename scalarType>
      int build(void);

   };


}
}

#include "FTCTree_Template.h"

#endif  // TASKEDTREE_H
