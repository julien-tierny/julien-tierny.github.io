/// \ingroup base
//
/// \class ttk::ftc::Propagation
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-02-11
///
///\brief TTK class for the local propagation in the merge tree

#pragma once

#include "DataTypes.h"

#include <boost/heap/fibonacci_heap.hpp>

namespace ttk
{
namespace ftc
{
   class Propagation
   {
     private:
      // Current vertex_ stay in cache
      idVertex vertex_;
      // priority queue to recover next min
      boost::heap::fibonacci_heap<idVertex, boost::heap::compare<VertCompFN>> propagation_;

     public:
      Propagation(idVertex startVert, VertCompFN vertComp)
          : vertex_(startVert), propagation_(vertComp)
      {
      }

      Propagation(VertCompFN vertComp)
          : vertex_(nullVertex), propagation_(vertComp)
      {
         // need to use setStartVert before using the propagation
      }

      void setStartVert(const idVertex v)
      {
         vertex_ = v;
      }

      idVertex getCurrentMinVertex(void) const
      {
         return vertex_;
      }

      idVertex getNextMinVertex(void)
      {
         vertex_ = propagation_.top();
         propagation_.pop();
         return vertex_;
      }

      void addNewVertex(const idVertex v)
      {
         propagation_.emplace(v);
      }

      void merge(Propagation& other)
      {
         propagation_.merge(other.propagation_);
         vertex_ = propagation_.top();
      }

      bool empty()
      {
         return propagation_.empty();
      }

      // DEBUG ONLY
      bool find(idVertex v)
      {
         return std::find(propagation_.begin(), propagation_.end(), v) != propagation_.end();
      }
   };
}
}
