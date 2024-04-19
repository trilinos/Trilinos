// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STKSEARCH_KDTREE_BoundingBoxHierarchy_h_
#define STKSEARCH_KDTREE_BoundingBoxHierarchy_h_

#include <vector>                       // for vector
         // for AxisAlignedBB
#include <stk_search/kdtree/KDTree_BoundingBox.hpp>
#ifdef _OPENMP
#include "omp.h" 
#endif

namespace stk {
  namespace search {

enum { MAX_TREE_LEVELS = 100 };

/**
 *  Class for general hierarchy of bounding boxes.  This class is used for fast overlap searches.  Note, the 
 *  calculations for creating the hierarchy and searching the hierarchy have been HEAVILY optimized for both
 *  reduction of operation count as well as for cache efficency.
 *
 *  NKC NOTE, possible further optimizations include
 *    2)  Implement a more efficent search when search for multiple objects.  Currently each object is searched for
 *        independantly.  Might be able to get some better performance by taking into account search is searching
 *        some set of objects against another set of objects.
 */
template <typename BoxType>
class ObjectBoundingBoxHierarchy_T {
 public:
  /**
   * Create a hierarchy from a list of object bounding boxes.  The top level of the hierarchy will contain all boxes.  In the next level
   *  each leave of the tree will contain approximiatly half of the boxes.  The hierarchy continues down until the tree node contains only 
   *  a single box
   */
  KOKKOS_FORCEINLINE_FUNCTION void set_right_child_offset(const int right_child_offset_) {right_child_offset = right_child_offset_;}
  KOKKOS_FORCEINLINE_FUNCTION const BoxType& GetBox() const {return m_box;}
  KOKKOS_FORCEINLINE_FUNCTION BoxType& GetBox() {return m_box;}

  //
  //  Right child offset stores one of two things.
  //    If the offset is <= 0 the current object is a terminal node of the tree.  The value is the negative of the 
  //  the object number associated with the object represented by the terminal node.  If the tree is created from
  //  inputBoxes that are not std::search::OblectBoundingBoxes, then the object number is the offset into the
  //  vector of inputBoxes.
  //    If the value of right_child_offset is positive it is the offset from the current object to the objects
  //  right child.  Note that the offset to the left child is always one.  Thus for a given object the left child
  //  can be found at this[1] and the right child at this[right_child_offset]
  //
 public:
  int right_child_offset = 0;
  BoxType m_box;
};


template <typename RangeBoxType>
class ProximitySearchTree_T {
 public:
  // Create an empty tree
  ProximitySearchTree_T() = default;

  // Create and initialize a search tree from a set of input boxes.  Assumed ident based on original box order

  ProximitySearchTree_T(const std::vector<RangeBoxType>& inputBoxes);
  
  ProximitySearchTree_T(std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> >& inputBoxes);

  // Recreate new search tree from scratch
  void InitializeSearch(std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> >& inputBoxes);

  // Update exiting tree with (subtly) differnt coordinates
  //
  //  This is faster than 'InitializeSearch' but must been done on the existing box set
  //  and could lead to a less ideal tree than 'InitializeSearch' if motions were large.
  //
  int UpdateSearch(std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> >& boxes);

  //
  // Get the a box that bounds the whole search tree
  //
  inline const RangeBoxType BoundingBox() const {
    return ( empty() ?  RangeBoxType() : m_tree[0].GetBox() );
  }

  //
  // Find all overlaps with a given search object
  //
  template<typename DomainBox>
  inline void SearchForOverlap(const DomainBox& searchObject,
                                                    std::vector<int>& returnList) const;  

  template<typename DomainBox>
  inline void SearchForOverlap(const DomainBox& searchObject,
                                                    std::vector<int>& returnIndexList,
                                                    std::vector<RangeBoxType>& returnBoxList) const;

  // Find if there is any overlap of the with the given box and the search tree
  template<typename DomainBox>
  inline bool AnyOverlap(const DomainBox& searchObject) const;

  // Check if a tree is empty, thus any search would return no boxes
  bool empty() const {return m_tree.empty();}

  // Empty out everything in the tree
  void clear() {m_tree.clear();}


  std::vector<ObjectBoundingBoxHierarchy_T<RangeBoxType> > &data()             { return m_tree; }
  const std::vector<ObjectBoundingBoxHierarchy_T<RangeBoxType> > &data() const { return m_tree; }

 protected:
  std::vector<ObjectBoundingBoxHierarchy_T<RangeBoxType> > m_tree;
};

  }
}

#include<stk_search/kdtree/KDTree_impl.hpp>


#endif // STKSEARCH_KDTREE_BoundingBoxHierarchy_h_

