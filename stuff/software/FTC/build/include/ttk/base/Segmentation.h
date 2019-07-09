/// \ingroup base
//
/// \class ttk::ftc::Segment
/// \class ttk::ftc::Segments
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date September 2016.
///
///\brief TTK classes for the segmentation of the contour tree
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).

#pragma once

#include "DataTypes.h"
#include "Scalars.h"
#include "FTCCommon.h"

#ifndef TTK_ENABLE_KAMIKAZE
#include <iostream>
#endif
#include <sstream>

#include <list>
#include <tuple>
#include <vector>

using namespace std;

namespace ttk
{
namespace ftc
{

   using segm_it           = std::vector<idVertex>::iterator;
   using segm_rev_it       = std::vector<idVertex>::reverse_iterator;
   using segm_const_it     = std::vector<idVertex>::const_iterator;
   using segm_const_rev_it = std::vector<idVertex>::const_reverse_iterator;

   // Segmentation data
   struct Region {
      // inverted in case of split tree
      segm_it segmentBegin;
      segm_it segmentEnd;
   };


   // one segment: like a vector<idVertex>
   // have a fixed size
   class Segment
   {
     private:
      vector<idVertex> vertices_;

     public:
      Segment(idVertex size);

      void sort(const Scalars<float>* s);
      void createFromList(const Scalars<float>* s, list<vector<idVertex>>& regularsList,
                          const bool reverse);
      void increase(const idVertex newSize);

      segm_const_it begin(void) const;
      segm_const_it end(void) const;
      segm_it       begin(void);
      segm_it       end(void);
      idVertex      size(void) const;

      idVertex operator[](const size_t& idx) const;
      idVertex& operator[](const size_t& idx);
   };

   // All the segments of the mesh, like a vector<Segment>
   class Segments
   {
     private:
      vector<Segment> segments_;

     public:
      Segments();

      Segments(const Segment&) = delete;

      // add one vertex
      tuple<segm_it, segm_it> addLateSimpleSegment(idVertex v);

      // sort all segment (in parallel?)
      void sortAll(const Scalars<float>* s);

      // callable once
      void resize(const vector<idVertex>& sizes);

      // vector like
      idSegment size(void) const;
      void clear(void);
      Segment& operator[](const size_t& idx);
      const Segment& operator[](const size_t& idx) const;

      // print
      inline string print(void) const
      {
         stringstream res;
         res << "{" << endl;
         for (const auto& s : segments_) {
            // if (s.size()) {
            //    res << s[0] << " " << s[s.size() - 1] << " : " << s.size() << endl;
            // }

            for (const auto& it : s) {
               res << it << " ";
            }
         }
         res << "}" << endl;
         return res.str();
      }

   };

   // The segmentation of one arc is a list of segment
   class ArcRegion
   {
     private:
      // list of segment composing this segmentation and for each segment
      // the begin and the end inside it (as a segment may be subdivided)
      list<Region> segmentsIn_;
      // when and how to compact ?
      vector<idVertex> segmentation_;

#ifndef TTK_ENABLE_KAMIKAZE
      // true if the segmentation have been sent to the segmentation_ vector
      bool segmented_;
#endif

     public:
      // create a new arc region with the associated Segment in Segments
      ArcRegion();

      ArcRegion(const segm_it& begin, const segm_it& end);

      // vertex used for the split (not in segmentation anymore), remaining region
      tuple<idVertex, ArcRegion> splitFront(idVertex v, const Scalars<float>* s);

      // vertex used for the split (not in segmentation anymore), remaining region
      tuple<idVertex, ArcRegion> splitBack(idVertex v, const Scalars<float>* s);

      idVertex findBelow(idVertex v, const Scalars<float>* s,
                         const vector<idCorresp>& vert2treeOther = vector<idCorresp>()) const;

      void concat(const segm_it& begin, const segm_it& end);

      void concat(const ArcRegion& r);

      // if contiguous with current region, merge them
      // else return false
      bool merge(const ArcRegion& r);

      void clear(void)
      {
         segmentsIn_.clear();
         segmentation_.clear();
      }

      // Put all segments in one vector in the arc
      // Suppose that all segment are already sorted
      // see Segments::sortAll
      // For Contour Tree segmentaion, you have to precise the current arc since
      // a segment can contain vertices for several arcs
      void createSegmentation(const Scalars<float>* s);

      inline idVertex count(void) const
      {
         size_t res = 0;
         for (const auto& reg : segmentsIn_) {
            res += abs(distance(reg.segmentBegin, reg.segmentEnd));
         }
         return res;
      }

      inline string print(void) const
      {
         stringstream res;
         res << "{";
         for (const auto& reg : segmentsIn_) {
            // res << " " << *reg.segmentBegin;
            // res << "-" << *(reg.segmentEnd - 1);

            auto it = reg.segmentBegin;
            for(; it != reg.segmentEnd; ++it){
               res << *it << ", ";
            }
         }
         res << " }";
         return res.str();
      }

      // Direct access to the list of region
      const decltype(segmentsIn_) & getRegions(void) const
      {
         return segmentsIn_;
      }

      decltype(segmentsIn_) & getRegions(void)
      {
         return segmentsIn_;
      }

      // vector like manip

      inline idVertex size(void) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!segmented_)
            std::cerr << "Needs to create segmentation before size" << std::endl;
#endif
         return segmentation_.size();
      }

      idVertex operator[](idVertex v) const
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!segmented_)
            cerr << "Needs to create segmentation before getting segmentation" << endl;
#endif
         return segmentation_[v];
      }

      idVertex& operator[](idVertex v)
      {
#ifndef TTK_ENABLE_KAMIKAZE
         if (!segmented_)
            cerr << "Needs to create segmentation before getting segmentation" << endl;
#endif
         return segmentation_[v];
      }

      decltype(segmentation_)::iterator begin(void)
      {
         return segmentation_.begin();
      }

      decltype(segmentation_)::iterator end(void)
      {
         return segmentation_.end();
      }
   };

}
}
