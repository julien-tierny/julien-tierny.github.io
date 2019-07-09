/// \ingroup base
//
/// \class ttk::ftc::SharedData
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2018-02-11
///
///\brief TTK container to share data between threads

#pragma once

#include "DataTypes.h"
#include "AtomicVector.h"
#include "Propagation.h"

namespace ttk
{
namespace ftc
{
   class SharedData {
     private:
      idVertex                   extremum_;
      AtomicVector<Propagation*> propagations_;
      AtomicVector<idSuperArc>   openedArcs_;

     public:
      explicit SharedData(idVertex e) : extremum_(e), propagations_(50), openedArcs_(50)
      {
      }

      idVertex getExtremum(void) const
      {
         return extremum_;
      }

      void setExtremum(const idVertex v)
      {
         extremum_ = v;
      }

      idNode getNbPropagation(void) const
      {
         return propagations_.size();
      }

      Propagation* getPropagation(const idNode p)
      {
         return propagations_[p];
      }

      bool onePropagation(void) const
      {
         return propagations_.size() == 1;
      }

      void addPropagation(Propagation* curProp)
      {
         const idThread& thisTask = propagations_.getNext();
         propagations_[thisTask]  = curProp;
      }

      void clearPropagation(void)
      {
         propagations_.reset();
      }

      void keepOnlyFirstProp(void)
      {
         propagations_.reset(1);
      }

      void addArc(const idSuperArc arc)
      {
         idSuperArc thisArc  = openedArcs_.getNext();
         openedArcs_[thisArc] = arc;
      }

      AtomicVector<idSuperArc>& getOpenedArcs(void)
      {
         return openedArcs_;
      }

      idNode getNbOpenedArcs(void) const
      {
         return openedArcs_.size();
      }

      void clearOpenedArcs(void)
      {
         openedArcs_.reset();
      }

      void merge(SharedData& other)
      {
         for (auto* prop : other.propagations_) {
            addPropagation(prop);
         }

         for (auto& arc : other.openedArcs_) {
            addArc(arc);
         }
      }

      void reserve(const size_t& s)
      {
         propagations_.reserve(s);
         openedArcs_.reserve(s);
      }
   };
}
}
