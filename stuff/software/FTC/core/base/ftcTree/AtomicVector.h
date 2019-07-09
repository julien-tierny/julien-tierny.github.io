/// \ingroup base
///
/// \class ttk::ftc::AtomicVector
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date 2017-02-09
///
///\brief TTK processing package that manage a paralle vecrion of vector

#pragma once

#ifdef TTK_ENABLE_OPENMP
#include <omp.h>
#endif // TTK_ENABLE_OPENMP
#include <iterator>
#include <vector>

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#include <typeinfo>
#endif

namespace ttk
{
   template <typename type>
   class AtomicVector : public std::vector<type>
   {
     private:
      std::size_t nextId;
      // For new value initialization
      const type defaultValue;

     public:
      AtomicVector(const std::size_t initSize = 1, const type &dv = type{})
          : std::vector<type>{}, nextId{0}, defaultValue{dv}
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!initSize) {
            std::cout << "Caution, Atomic vector need a non-0 init size !" << std::endl;
            std::vector<type>::resize(initSize, defaultValue);
         } else
#endif
         {
            std::vector<type>::resize(initSize, defaultValue);
         }
      }

      // copy constructor
      AtomicVector(const AtomicVector &other) : std::vector<type>(other), nextId(other.nextId)
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!std::vector<type>::size()) {
            reserve(1);
         }
#endif
      }

      AtomicVector(AtomicVector &&other) = default;

      virtual ~AtomicVector() = default;

      // ---
      // STL
      // ---

      void reserve(const std::size_t &newSize)
      {
         if (newSize > std::vector<type>::size()) {
#ifndef TTK_ENABLE_KAMIKAZE
#ifdef TTK_ENABLE_OPENMP
            if (omp_in_parallel()) {
               // WARNING: In parallel we do not want to make reserve as it can lead to
               // data race, we should not enter here
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(AtomicUFReserve)
#endif
               {
                  std::vector<type>::resize(newSize, defaultValue);
               }

            } else
#endif
#endif
            {
               std::vector<type>::resize(newSize, defaultValue);
            }
         }
      }

      void clear(void)
      {
         reset();

         // Remove old content
         std::size_t oldSize = std::vector<type>::size();
         std::vector<type>::clear();
         reserve(oldSize);
      }


      std::size_t size(void) const
      {
         return nextId;
      }

      bool empty(void) const
      {
         return nextId == 0;
      }

      void push_back(const type &elmt)
      {
         const auto &curPos = getNext();
         (*this)[curPos]    = elmt;
      }

      void emplace_back(const type &elmt)
      {
         // not really an emplace back as we still have a copy
         const auto &curPos = getNext();
         (*this)[curPos]    = elmt;
      }

      void swap(AtomicVector<type>& other)
      {
         const auto tmpId = other.nextId;
         std::vector<type>::swap(other);
         other.nextId = nextId;
         nextId = tmpId;
      }

      // -----
      // Other
      // -----

      void reset(const std::size_t &nId = 0)
      {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
         nextId = nId;
      }

      std::size_t getNext(void)
      {
         std::size_t resId;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif
         resId = nextId++;

         if (resId == std::vector<type>::size()) {
// #ifndef TTK_ENABLE_KAMIKAZE
//             type tmp;
//             std::cout << "RE-RESERVE VECTOR : " << nextId << " " << typeid(tmp).name() << std::endl;
// #endif
            reserve(std::vector<type>::size() * 2);
         }

         return resId;
      }



      // --------
      // OPERATOR
      // --------

      AtomicVector<type> &operator=(const AtomicVector<type> &other)
      {
         std::vector<type>::operator=(other);
         nextId                     = other.nextId;
      }

      // ---------
      // ITERATORS
      // ---------
      // allow foreach on the vector

      typedef typename std::vector<type>::iterator       iterator;
      typedef typename std::vector<type>::const_iterator const_iterator;

      iterator end()
      {
         return this->begin() + nextId;
      }

      const_iterator cend() const
      {
         return this->cbegin() + nextId;
      }

      typedef typename std::vector<type>::reverse_iterator       riterator;
      typedef typename std::vector<type>::const_reverse_iterator const_riterator;

      riterator rbegin()
      {
         return this->rend() - (nextId - 1);
      }

      const_riterator crbegin() const
      {
         return this->crend() - (nextId - 1);
      }
   };
}
