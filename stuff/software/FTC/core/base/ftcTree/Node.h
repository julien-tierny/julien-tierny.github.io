/// \ingroup base
/// \class ttk::ftc::Node
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK container representing a node of the FTCTree_MT

#pragma once

#include "DataTypes.h"

#include <Debug.h>

#include <functional>
#include <vector>

namespace ttk
{
namespace ftc
{
   class Node
   {
      friend class FTCTree_MT;

     private:
      // mesh vertex where this node is
      idVertex vertexId_;
      // link with superArc above and below
      std::vector<idSuperArc> vect_downSuperArcList_, vect_upSuperArcList_;
      // Id to manage size by hand
      idSuperArc nbDownArcs, nbUpArcs;

     public:

      // -----------------
      // CONSTRUCTOR
      // -----------------

      // This node will need to receive a vertex id before being printed
      Node()
          : vertexId_(nullVertex),
            vect_downSuperArcList_(),
            vect_upSuperArcList_(),
            nbDownArcs(0),
            nbUpArcs(0)
      {
      }

      Node(idVertex id)
          : vertexId_(id),
            vect_downSuperArcList_(),
            vect_upSuperArcList_(),
            nbDownArcs(0),
            nbUpArcs(0)
      {
      }

      // -----------------
      // ACCESSOR
      // ------------------


      // Vertex id

      inline idVertex getVertexId() const
      {
         return vertexId_;
      }

      inline void setVertexId(idVertex vertexId)
      {
         vertexId_ = vertexId;
      }

      // vector arcs

      inline idSuperArc getNumberOfDownSuperArcs() const
      {
         return nbDownArcs;
      }

      inline idSuperArc getNumberOfUpSuperArcs() const
      {
         return nbUpArcs;
      }

      inline idSuperArc getNumberOfSuperArcs() const
      {
         return nbDownArcs + nbUpArcs;
      }

      inline void reserveDownArcs(const idSuperArc nbDown)
      {
          vect_downSuperArcList_.reserve(nbDown);
      }

      inline void reserveUpArcs(const idSuperArc nbUp)
      {
          vect_upSuperArcList_.reserve(nbUp);
      }

      inline idSuperArc getDownSuperArcId(idSuperArc neighborId) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if ((size_t)neighborId >= nbDownArcs) {
            cerr << "[Merge Tree:Node] get down on bad neighbor !";
            cerr << endl;
            return 0;
         }
#endif
         return vect_downSuperArcList_[neighborId];
      };

      inline idSuperArc getUpSuperArcId(idSuperArc neighborId) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (neighborId >= nbUpArcs) {
            cerr << "[FTCTree_MT:Node] No SuperArc to access " << static_cast<unsigned>(neighborId);
            cerr << endl;
         }
#endif
         return vect_upSuperArcList_[neighborId];
      }

      inline void addDownSuperArcId(idSuperArc downSuperArcId)
      {
         vect_downSuperArcList_.emplace_back(downSuperArcId);
         ++nbDownArcs;
      }

      inline void addUpSuperArcId(idSuperArc upSuperArcId)
      {
         vect_upSuperArcList_.emplace_back(upSuperArcId);
         ++nbUpArcs;
      }

      // Find and remove the arc
      inline void removeDownSuperArcId(idSuperArc idSa)
      {
         for (idSuperArc i = 0; i < vect_downSuperArcList_.size(); ++i) {
            if (vect_downSuperArcList_[i] == idSa) {
               const idSuperArc last     = --nbDownArcs;
               vect_downSuperArcList_[i] = vect_downSuperArcList_[last];
               vect_downSuperArcList_.pop_back();

               return;
            }
         }
      }

      // Find and remove the arc
      inline void removeUpSuperArcId(idSuperArc idSa)
      {
         for (idSuperArc i = 0; i < vect_upSuperArcList_.size(); ++i) {
            if (vect_upSuperArcList_[i] == idSa) {
               const idSuperArc last   = --nbUpArcs;
               vect_upSuperArcList_[i] = vect_upSuperArcList_[last];
               vect_upSuperArcList_.pop_back();

               return;
            }
         }
      }

      inline idSuperArc clearDownSuperArcs(void)
      {
         idSuperArc s = vect_downSuperArcList_.size();
         vect_downSuperArcList_.clear();
         nbDownArcs = 0;
         return s;
      }

      inline idSuperArc clearUpSuperArcs(void)
      {
         idSuperArc s = vect_upSuperArcList_.size();
         vect_upSuperArcList_.clear();
         nbUpArcs = 0;
         return s;
      }

      void sortUpArcs(std::function<bool(const idSuperArc, const idSuperArc)> comp)
      {
         sort(vect_upSuperArcList_.begin(), vect_upSuperArcList_.end(), comp);
      }

      void sortDownArcs(std::function<bool(const idSuperArc, const idSuperArc)> comp)
      {
         sort(vect_downSuperArcList_.begin(), vect_downSuperArcList_.end(), comp);
      }

   };

}
}

