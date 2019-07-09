/// \ingroup base
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK FTC package managing data types.
/// These are only alias on existing primitive type,
/// but allow semantic in the code.
/// This can also help changing these type to support larger data sets

#pragma once

#include <functional>
#include <limits>
#include <set>
#include <tuple>

namespace ttk
{
namespace ftc
{
   // Types
   // --------

   /// \brief SuperArc index in vect_superArcs_
   using idSuperArc = long unsigned int;
   /// \brief Node index in vect_nodes_
   using idNode = unsigned int;
   /// \brief Vertex index in scalars_
   using idVertex = int;
   /// \brief Edge index in vect_edgeList_
   using idEdge = int;
   /// \brief Cell index in vect_cellList_
   using idCell = int;

   /// \brief type used to recover Node/Arc in vert2tree SIGNED ONLY
   // Warning, in long long int the max super arc is -1, might not be able to deal with
   // too large data
   using idCorresp = long long int;

   /// \brief for the segmentation, we have an array of segment containing area of the mesh
   using idSegment = idSuperArc;

   /// \brief type use to store threads related numbers
   using numThread = unsigned char;

   /// \brief type stored by UnionFind
   using ufDataType = long int;

   /// \brief manage number of threads
   using idThread = unsigned char;

   /// \brief for task identifiers
   using idTask = idNode;

   /// \brief for vertex up/down valence
   using valence = signed char;

   // For tasks:
   // Set using scalar value comparison
   using VertCompFN     = std::function<bool(const idVertex, const idVertex)>;
   using SetPropagation = std::set<idVertex, VertCompFN>;
   using SetCompFN      = std::function<bool(const SetPropagation &, const SetPropagation &)>;

   // Special values for types
   // --------------------------

   // QUESTION impact on performance using max (0 would be faster alloacted)
   static const idSuperArc     nullSuperArc  = std::numeric_limits<idSuperArc>::max();
   static const idNode         nullNodes     = std::numeric_limits<idNode>::max();
   static const idVertex       nullVertex    = std::numeric_limits<idVertex>::max();
   static const idEdge         nullEdge      = std::numeric_limits<idEdge>::max();
   static const idCell         nullCell      = std::numeric_limits<idCell>::max();
   static const idCorresp      nullCorresp   = std::numeric_limits<idCorresp>::max();
   static const idSegment      nullSegment   = std::numeric_limits<idSegment>::max();
   static const ufDataType     nullUfData    = std::numeric_limits<ufDataType>::max();
   static constexpr ufDataType specialUfData = std::numeric_limits<ufDataType>::max() - 1;
   static const idThread       nullThread    = std::numeric_limits<idThread>::max();
   static const valence        nullValence   = std::numeric_limits<valence>::max();

   // Enum data
   // ----------

   enum TreeType : char { Join = 0, Split = 1, Contour = 2, Join_Split = 3 };

   enum SimplifMethod : char { Persist = 0, Span = 1, NbVert = 2, NbArc = 3 };

   enum ComponentState : char { Visible, Hidden, Pruned, Merged };

   enum class TreeComponent { Arc = -1, Local_minimum, Saddle1, Saddle2, Local_maximum };

   enum class ArcType:char { Min_arc = 0, Max_arc, Saddle1_arc, Saddle2_arc, Saddle1_saddle2_arc };

   enum class NodeType { Local_minimum = 0, Saddle1, Saddle2, Degenerate, Local_maximum, Regular };
}
}

// Tests

// When detecting the trunk early, we check for it regularly.
// all the TRUNKCHECKCHUNK vertices.
#ifndef NDEBUG
#define TRUNKCHECKCHUNK  1
#else
#define TRUNKCHECKCHUNK  10000
#endif

// Stats time are not compatible with CT, only MT
// #define TTK_ENABLE_FTC_TREE_STATS_TIME 1
