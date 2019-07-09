/// \ingroup base
///
/// \class ttk::ftc::AtomicUF
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK FTC Tree Union Find management.
/// This class suport concurrent call on the Find function,
/// but not on the Union function.

#pragma once

#include "DataTypes.h"
#include "SharedData.h"

#include <memory>
#include <vector>

namespace ttk
{
namespace ftc
{

   class AtomicUF
   {
     private:
      unsigned   rank_;
      AtomicUF * parent_;
      SharedData data_;

     public:
      inline explicit AtomicUF(idVertex extrema = nullVertex) : rank_(0), data_(extrema)
      {
         parent_ = this;
      }

      // heavy recursif
      AtomicUF *find()
      {
         if (parent_ == this)
            return this;
         else {
            decltype(parent_) tmp = parent_->find();
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
            parent_ = tmp;

            return parent_;
         }
      }

      // Shared data get/set

      inline idVertex getExtremum(void) const
      {
         return data_.getExtremum();
      }

      inline Propagation* getCurrentPropagation(void)
      {
#ifndef TTK_ENABLE_KAMIKAZE
          if(!data_.onePropagation()) {
             std::cout << "AtomicUF :: getCurrentPropagation : nb state 1 !=  "
                       << data_.getNbPropagation() << std::endl;
          }
#endif
          return data_.getPropagation(0);
      }

      inline size_t getNbPropagation(void) const
      {
         return data_.getNbPropagation();
      }

      inline void clearOpenedArcs(void)
      {
         data_.clearOpenedArcs();
      }

      inline void setExtremum(idVertex v)
      {
         data_.setExtremum(v);
      }

      inline void reserveData(const size_t& s)
      {
          data_.reserve(s);
      }

      inline void addPropagation(Propagation *s)
      {
         data_.addPropagation(s);
      }

      inline void clearPropagation(void)
      {
         data_.clearPropagation();
      }

      inline void addArcToClose(idSuperArc a)
      {
         data_.addArc(a);
      }

      idNode getNbOpenedArcs(void) const
      {
         return data_.getNbOpenedArcs();
      }

      AtomicVector<idSuperArc>& getOpenedArcs(void)
      {
         return data_.getOpenedArcs();
      }

      inline Propagation *mergeStates(void)
      {
         Propagation * p           = data_.getPropagation(0);
         const valence nbPropLocal = data_.getNbPropagation();

         for (valence i = 1; i < nbPropLocal; ++i) {
            p->merge(*data_.getPropagation(i));
         }

         data_.keepOnlyFirstProp();

         return p;
      }

      // UF get / set

      inline int getRank() const
      {
         return rank_;
      }

      inline void setRank(const int &rank)
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
         rank_ = rank;
      }

      inline void setParent(AtomicUF *parent)
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
         parent_ = parent;
      }

      static inline AtomicUF *makeUnion(AtomicUF *uf0, AtomicUF *uf1)
      {
         uf0 = uf0->find();
         uf1 = uf1->find();

         if (uf0 == uf1) {
            return uf0;
         } else if (uf0->getRank() > uf1->getRank()) {
            uf1->setParent(uf0);
            uf0->data_.merge(uf1->data_);
            return uf0;
         } else if (uf0->getRank() < uf1->getRank()) {
            uf0->setParent(uf1);
            uf1->data_.merge(uf0->data_);
            return uf1;
         } else {
            uf1->setParent(uf0);
            uf0->setRank(uf0->getRank() + 1);
            uf0->data_.merge(uf1->data_);
            return uf0;
         }

         return NULL;
      }

      static inline AtomicUF *makeUnion(std::vector<AtomicUF *> &sets)
      {
         AtomicUF *n = NULL;

         if (!sets.size())
            return NULL;

         if (sets.size() == 1)
            return sets[0];

         for (int i = 0; i < (int)sets.size() - 1; i++)
            n = makeUnion(sets[i], sets[i + 1]);

         return n;
      }

      inline bool operator<(const AtomicUF &other) const
      {
         return rank_ < other.rank_;
      }

      inline bool operator>(const AtomicUF &other) const
      {
         return rank_ > other.rank_;
      }
   };

}
}
