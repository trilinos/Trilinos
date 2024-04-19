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


#ifndef STKSEARCH_KDTREE_IMPL_H
#define STKSEARCH_KDTREE_IMPL_H

// #######################  Start Clang Header Tool Managed Headers ########################
#include <algorithm>                                   // for min
#include <cassert>                                     // for assert
#include <cmath>                                       // for fabs
#include <stk_math/StkMath.hpp>                        // for stk::math::max, stk::math::min
#include <stk_search/kdtree/KDTree_Threaded_Sort.hpp>  // for ThreadedSort
#include <utility>                                     // for pair, make_pair

#include "KDTree_BoundingBox.hpp"

// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
  namespace search {

///< Set output box to contain both input boxes
template <typename BoxType>
KOKKOS_INLINE_FUNCTION void UnionBoxes(const BoxType& inBox1, const BoxType& inBox2, BoxType& outputBox) {
  outputBox.set_box(stk::math::min(inBox1.get_x_min(), inBox2.get_x_min()),
                    stk::math::min(inBox1.get_y_min(), inBox2.get_y_min()),
                    stk::math::min(inBox1.get_z_min(), inBox2.get_z_min()),
                    stk::math::max(inBox1.get_x_max(), inBox2.get_x_max()),
                    stk::math::max(inBox1.get_y_max(), inBox2.get_y_max()),
                    stk::math::max(inBox1.get_z_max(), inBox2.get_z_max()));
}


inline void adjust_moment_split(const int num_boxes, const float global_center_x, const float global_center_y, const float global_center_z, float& moment_x, float& moment_y, float& moment_z) {
  //
  //  Adjust the moment to be taken about the global_center location
  //  Note:  Here we are actually calculating num_boxes * moment.  Only need to compare relative values
  //  of moment so multiplying by num_boxes is fine.  This is done to avoid any division operations.
  //
  //  Models are often symmetric making the global x, y, or z moment near identical.  This can make the splitting direction,
  //  and thus sometimes downstream results sensitive to infintesimal roundoff.  Apply the 1.0012345 and 1.0000531 scaling
  //  factor to slightly bias splitting to x then y directions.  Note, the cliff roundoff ambigouity is still there, it is just moved
  //  to a place that is somewhat less likely to occur in practice.
  //
  moment_x *= num_boxes;
  moment_x -= global_center_x * global_center_x;
  moment_x *= 1.0012345;
  moment_y *= num_boxes;
  moment_y -= global_center_y * global_center_y;
  moment_y *= 1.0000531;
  moment_z *= num_boxes;
  moment_z -= global_center_z * global_center_z;
}
 
template <typename RangeBoxType>
inline void split_boxes(const int &num_boxes, const stk::search::ObjectBoundingBox_T<RangeBoxType> *boxes, float& moment_x, float& moment_y, float& moment_z) {
  float global_center_x(0.0);
  float global_center_y(0.0);
  float global_center_z(0.0);
  //
  //  Calculate the centroid of each component box, the centroid of the current encompassing box, and the
  //  moment of inertia of all particles in the box about (0, 0, 0)
  //
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    float centroid[3];
    centroid[0] = boxes[ibox].GetBox().get_x_min() + boxes[ibox].GetBox().get_x_max();
    centroid[1] = boxes[ibox].GetBox().get_y_min() + boxes[ibox].GetBox().get_y_max();
    centroid[2] = boxes[ibox].GetBox().get_z_min() + boxes[ibox].GetBox().get_z_max();
    global_center_x += centroid[0];
    global_center_y += centroid[1];
    global_center_z += centroid[2];
    moment_x += centroid[0] * centroid[0];
    moment_y += centroid[1] * centroid[1];
    moment_z += centroid[2] * centroid[2];
  }
  adjust_moment_split(num_boxes, global_center_x, global_center_y, global_center_z, moment_x, moment_y, moment_z);
}

#ifdef _OPENMP
template <typename RangeBoxType>
static inline void split_boxes_threaded(const int num_boxes, const stk::search::ObjectBoundingBox_T<RangeBoxType>* const boxes, float& ret_moment_x, float& ret_moment_y, float& ret_moment_z, int numThreadsToUse=-1) {
  float global_center_x(0.0);
  float global_center_y(0.0);
  float global_center_z(0.0);
  float moment_x(0.0);
  float moment_y(0.0);
  float moment_z(0.0);
  //
  //  Calculate the centroid of each component box, the centroid of the current encompassing box, and the
  //  moment of inertia of all particles in the box about (0, 0, 0)
  //
  if(numThreadsToUse < 0) {
    numThreadsToUse = omp_get_max_threads();
  }

#ifdef _OPENMP
#pragma omp parallel for reduction(+: global_center_x, global_center_y, global_center_z, moment_x, moment_y, moment_z) schedule(static) default(none) num_threads(numThreadsToUse) firstprivate(boxes, num_boxes)
#endif
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    float centroid[3];
    centroid[0] = boxes[ibox].GetBox().get_x_min() + boxes[ibox].GetBox().get_x_max();
    centroid[1] = boxes[ibox].GetBox().get_y_min() + boxes[ibox].GetBox().get_y_max();
    centroid[2] = boxes[ibox].GetBox().get_z_min() + boxes[ibox].GetBox().get_z_max();
    global_center_x += centroid[0];
    global_center_y += centroid[1];
    global_center_z += centroid[2];
    moment_x += centroid[0] * centroid[0];
    moment_y += centroid[1] * centroid[1];
    moment_z += centroid[2] * centroid[2];
  }
  adjust_moment_split(num_boxes, global_center_x, global_center_y, global_center_z, moment_x, moment_y, moment_z);
  ret_moment_x = moment_x;
  ret_moment_y = moment_y;
  ret_moment_z = moment_z;
}
#endif

//
//  Store a one noded tree, just set the terminal case
//
template <typename RangeBoxType>
inline void store_1_node_tree(ObjectBoundingBoxHierarchy_T<RangeBoxType> *root_node, const stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes) {
  root_node->set_right_child_offset(-boxes[0].get_object_number());
  root_node->GetBox().set_box(boxes[0].GetBox());
}


template <typename RangeBoxType>
KOKKOS_INLINE_FUNCTION void store_2_node_tree(ObjectBoundingBoxHierarchy_T<RangeBoxType> *root_node, const stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes) {
  //
  //  Special case, only two boxes passed in.  This happens often due to the top to bottom tree structure.
  //  Significant optimizations are avalable for this case
  //  Split the list into two halves.  Put half the boxes into one child and half into the other.  Complete processing
  //  on both halves.
  //
  //  Store self
  //
  root_node->set_right_child_offset(2);
  UnionBoxes(boxes[0].GetBox(), boxes[1].GetBox(), root_node->GetBox());
  //
  //  Store left child
  //
  root_node[1].set_right_child_offset(-boxes[0].get_object_number());
  root_node[1].GetBox().set_box(boxes[0].GetBox());
  //
  //  Store right child
  //
  root_node[2].set_right_child_offset(-boxes[1].get_object_number());
  root_node[2].GetBox().set_box(boxes[1].GetBox());
}

template <typename RangeBoxType>
inline void store_3_node_tree(ObjectBoundingBoxHierarchy_T<RangeBoxType> *root_node, const stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes) {
  //
  //  Special case, exactly three boxes pased in.  This happens resonably often, aditionally, if handle 2 box case, and
  //  3 box case then automatically will handle the one box cases (except for the very special case of exactly one box
  //  in the entire hierarchy, which is handled in the calling routine.)
  //
  //  Determine how to split the boxes.  The closest two boxes will be placed in the right child,
  //  and the third box in the left child.  The index arrays will not be used for this operation
  //  for optimization purposes
  //

  float cent0[3];
  cent0[0] = boxes[0].GetBox().get_x_min() + boxes[0].GetBox().get_x_max();
  cent0[1] = boxes[0].GetBox().get_y_min() + boxes[0].GetBox().get_y_max();
  cent0[2] = boxes[0].GetBox().get_z_min() + boxes[0].GetBox().get_z_max();
  float cent1[3];
  cent1[0] = boxes[1].GetBox().get_x_min() + boxes[1].GetBox().get_x_max();
  cent1[1] = boxes[1].GetBox().get_y_min() + boxes[1].GetBox().get_y_max();
  cent1[2] = boxes[1].GetBox().get_z_min() + boxes[1].GetBox().get_z_max();
  float cent2[3];
  cent2[0] = boxes[2].GetBox().get_x_min() + boxes[2].GetBox().get_x_max();
  cent2[1] = boxes[2].GetBox().get_y_min() + boxes[2].GetBox().get_y_max();
  cent2[2] = boxes[2].GetBox().get_z_min() + boxes[2].GetBox().get_z_max();

  float dist_squared_01 = std::fabs(cent0[0]-cent1[0]) + std::fabs(cent0[1]-cent1[1]) + std::fabs(cent0[2]-cent1[2]);
  float dist_squared_12 = std::fabs(cent1[0]-cent2[0]) + std::fabs(cent1[1]-cent2[1]) + std::fabs(cent1[2]-cent2[2]);
  float dist_squared_20 = std::fabs(cent2[0]-cent0[0]) + std::fabs(cent2[1]-cent0[1]) + std::fabs(cent2[2]-cent0[2]);

  //
  //  Create the left and right hierarchies
  //
  if(dist_squared_01 < dist_squared_12 && dist_squared_01 < dist_squared_20) {
    //
    //  Right child has boxes 0 and 1, left child has box 2
    //
    (root_node+1)->set_right_child_offset(-boxes[2].get_object_number());
    (root_node+1)->GetBox().set_box(boxes[2].GetBox());
    (root_node+3)->set_right_child_offset(-boxes[0].get_object_number());
    (root_node+3)->GetBox().set_box(boxes[0].GetBox());
    (root_node+4)->set_right_child_offset(-boxes[1].get_object_number());
    (root_node+4)->GetBox().set_box(boxes[1].GetBox());
    (root_node+2)->set_right_child_offset(2);
    UnionBoxes(boxes[0].GetBox(),boxes[1].GetBox(), (root_node+2)->GetBox());
  } else {
    if(dist_squared_12 < dist_squared_20) {
      //
      //  Right child has boxes 1 and 2, left child has box 0
      //
      (root_node+1)->set_right_child_offset(-boxes[0].get_object_number());
      (root_node+1)->GetBox().set_box(boxes[0].GetBox());
      (root_node+3)->set_right_child_offset(-boxes[1].get_object_number());
      (root_node+3)->GetBox().set_box(boxes[1].GetBox());
      (root_node+4)->set_right_child_offset(-boxes[2].get_object_number());
      (root_node+4)->GetBox().set_box(boxes[2].GetBox());
      (root_node+2)->set_right_child_offset(2);
      UnionBoxes(boxes[1].GetBox(),boxes[2].GetBox(), (root_node+2)->GetBox());
    } else {
      //
      //  Right child has boxes 0 and 2, left child has box 1
      //
      (root_node+1)->set_right_child_offset(-boxes[1].get_object_number());
      (root_node+1)->GetBox().set_box(boxes[1].GetBox());
      (root_node+3)->set_right_child_offset(-boxes[0].get_object_number());
      (root_node+3)->GetBox().set_box(boxes[0].GetBox());
      (root_node+4)->set_right_child_offset(-boxes[2].get_object_number());
      (root_node+4)->GetBox().set_box(boxes[2].GetBox());
      (root_node+2)->set_right_child_offset(2);
      UnionBoxes(boxes[0].GetBox(), boxes[2].GetBox(), (root_node+2)->GetBox());
    }
  }
  root_node->set_right_child_offset(2);

  UnionBoxes((root_node+1)->GetBox(),(root_node+2)->GetBox(), root_node->GetBox());
}

// ======================================================================================================
//
//  ObjectBoundingBoxHierarchy methods
//
// ======================================================================================================

#ifdef _OPENMP
//
//  Extract 'numLimb' subtrees from the hierarchy_data by drilling down to the lower level trees
//
template <typename RangeBoxType>
inline void update_tree_partial(ObjectBoundingBoxHierarchy_T<RangeBoxType>*  hierarchy_data, 
                         ObjectBoundingBoxHierarchy_T<RangeBoxType>** subTreeHier, 
                         unsigned subdivisionLevel,
                         unsigned mySubTree) {
  //
  //  Terminal leaf case
  //
  if(subdivisionLevel == 1) {
    subTreeHier[mySubTree] = hierarchy_data;
    return;
  }

  const int cur_right_child_offset(hierarchy_data->right_child_offset);
  if(cur_right_child_offset <= 0) {
    //
    //  Early exit small tree terminal leaf case
    //
    subTreeHier[mySubTree] = hierarchy_data;
    return;
  } else {
    ObjectBoundingBoxHierarchy_T<RangeBoxType>* left_child  = hierarchy_data+1;
    update_tree_partial(left_child, subTreeHier, subdivisionLevel/2, mySubTree);
    ObjectBoundingBoxHierarchy_T<RangeBoxType>* right_child = hierarchy_data+cur_right_child_offset;
    update_tree_partial(right_child, subTreeHier, subdivisionLevel/2, mySubTree+subdivisionLevel/2);
  }
}

//
//  Finalize the tree output by computing hierarchal bounding box data on the top few tree levels.  Assume
//  subtrees already complete via threaded processing.
//
template <typename RangeBoxType>
inline void complete_tree_partial(ObjectBoundingBoxHierarchy_T<RangeBoxType>*  hierarchy_data, 
                                                       unsigned subdivisionLevel) {
  if(subdivisionLevel == 1) {
    //  Terminal thread leaf case
    return;
  }

  const int cur_right_child_offset(hierarchy_data->right_child_offset);
  if(cur_right_child_offset <= 0) {
    //
    //  Early exit small tree terminal leaf case
    //
    return;
  }

  ObjectBoundingBoxHierarchy_T<RangeBoxType>* left_child  = hierarchy_data+1;
  complete_tree_partial(left_child, subdivisionLevel/2);
  ObjectBoundingBoxHierarchy_T<RangeBoxType>* right_child = hierarchy_data+cur_right_child_offset;
  complete_tree_partial(right_child, subdivisionLevel/2);

  UnionBoxes(left_child->m_box, right_child->m_box, hierarchy_data->m_box);
}

#endif

template <typename RangeBoxType>
inline void update_tree_internal(ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_data, stk::search::ObjectBoundingBox_T<RangeBoxType> *boxes) {
  const int cur_right_child_offset(hierarchy_data->right_child_offset);
  if(cur_right_child_offset <= 0) {
    hierarchy_data->m_box.set_box(boxes[-cur_right_child_offset].GetBox());
  } else {
    ObjectBoundingBoxHierarchy_T<RangeBoxType>* left_child  = hierarchy_data+1;
    update_tree_internal(left_child, boxes);
    ObjectBoundingBoxHierarchy_T<RangeBoxType>* right_child = hierarchy_data+cur_right_child_offset;
    update_tree_internal(right_child, boxes);

    UnionBoxes(left_child->m_box, right_child->m_box, hierarchy_data->m_box);
  }
}



//
//  Update the bounding boxes that makeup a tree without changing their sorting.  Note,
//  this is only recomennded if box motions are small. If box motions are large the
//  updated tree will still be valid, but will become very suboptimal for searching.
//
template <typename RangeBoxType>
inline void update_tree(ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_data, stk::search::ObjectBoundingBox_T<RangeBoxType> *boxes) {
#ifdef _OPENMP
  unsigned maxNumThread = omp_get_max_threads();
  if(maxNumThread == 1) {
    update_tree_internal(hierarchy_data, boxes);
  }
  //
  //  First generate ~maxNumThread sub trees
  //
  //  The algorithm is optimized for thread counts that are powers of 2, for other counts the algorithm will function
  //  but the final parallel sub-tree update step will be unbalanced and sub-optimal.
  //
  unsigned numLimb = 1;
  while(numLimb < maxNumThread) {
    numLimb*=2;
  }
  std::vector<ObjectBoundingBoxHierarchy_T<RangeBoxType>*> sub_hierarchy_data(numLimb, nullptr);
  update_tree_partial(hierarchy_data, sub_hierarchy_data.data(), numLimb, 0);
  //
  //  In thread-parallel run the serial tree update algorithm on each sub-tree
  //
#pragma omp parallel for default(none)  shared(sub_hierarchy_data, boxes, numLimb)
  for(unsigned i=0; i<numLimb; ++i) {
    if(sub_hierarchy_data[i] == nullptr) continue;
    update_tree_internal(sub_hierarchy_data[i], boxes);
  }
  //
  //  Now do the top level recombination of the subtrees.  Note, assuming this operation is cheap due to number of threads
  //  used generally being much smaller than the number of boxes, thus the operation is just done in serial.
  //
  complete_tree_partial(hierarchy_data, numLimb);
#else
  // Pure serial update
  update_tree_internal(hierarchy_data, boxes);
#endif
}

template <typename RangeBoxType>
struct CreateHierarchyStack {
  ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_loc;
  stk::search::ObjectBoundingBox_T<RangeBoxType> *boxes;
  stk::search::ObjectBoundingBox_T<RangeBoxType> *scratch_boxes;
  int num_boxes;
  int* index_array0_t;
  int* index_array1_t;
  int* index_array2_t;
  int* scratch_index0_array;
  int* scratch_index1_array;
  int* scratch_index2_array;
  int unwind;
};

//
//  Non recursive hierarchy creation routine
//
//  NKC Note, this routine was origonally made non-recursive due to peformance issues with recursion on Janus
//  In ~2004.  Entirely likely new machines have no such issue with recursion and the recursive algorithm may
//  at least as faster or faster than the non-recursive algorithm.  Should probably try it out someday.
//
//  Removing recurse components improves speed.  However, recursion is easier to understand, here is the
//  Recursive puesdo code for this routine:
//
//  ObjectBoundingBoxHierarchy::create_hierarchy(boxes, scratch_boxes, num_boxes, index0, index1, index2) {
//    if(num_boxes == 1) {
//      1) Set the bounding box of this tree node to the bounding box of the object that the tree node represents.
//      2) Set the right child offset of this node to the negative of the object number this node represents
//    }
//    if(num_boxes > 1) {
//      1) Calculate the number of boxes to place into the left and right children
//      2) Determine which direction to split these boxes, this direction is taken as the direction with
//         the largest moment of inertia of box centroids
//      3) Sort the boxes based on the split direction using the index array for that direction.  Place sorted
//         boxes into scratch_boxes.
//      4) Update the sorting lists and create the left and right child sorting lists.
//      5) Recurively call the routine for the left and right child
//      this[1].create_hierachy(left_boxes, left_scratch_boxes, left_num_boxes, left_index0, left_index1, left_index2)
//      this[right_child_offset].create_hierachy(right_boxes, right_scratch_boxes, right_num_boxes, right_index0, right_index1, right_index2)
//      6) Recursive unwinding, set the bounding box of the node as the sum of the bounding boxes of the child nodes
//    }
//  }
//
template <typename RangeBoxType>
inline void create_hierarchy_loop(ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_start_ptr,
                                                       stk::search::ObjectBoundingBox_T<RangeBoxType> *boxes_start_ptr,
                                                       stk::search::ObjectBoundingBox_T<RangeBoxType> *scratch_boxes_start_ptr,
                                                       const int num_boxes_start,
                                                       int* index_array0_t_start_ptr,
                                                       int* index_array1_t_start_ptr,
                                                       int* index_array2_t_start_ptr,
                                                       int* scratch_index0_array_start_ptr,
                                                       int* scratch_index1_array_start_ptr,
                                                       int* scratch_index2_array_start_ptr) {
  //
  //  The tree creation loop handles several special cases, the effect is that the case of one, two and three boxs is not handled.
  //  correctly, if it occurs handle it here and exit tree creation
  //
  if(num_boxes_start == 1) {
    store_1_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  } else if(num_boxes_start ==2) {
    store_2_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  } else if(num_boxes_start ==3) {
    store_3_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  }
  //
  //  Create the object stack variables, this routine is really a recursive routine, however, some performance gain
  //  can be had by rewriting as a loop and managing the recursive stack by hand.
  //
  CreateHierarchyStack<RangeBoxType> object_stack[MAX_TREE_LEVELS];
  CreateHierarchyStack<RangeBoxType> *stack_ptr = object_stack;
  CreateHierarchyStack<RangeBoxType> *current_object = stack_ptr;
  //
  //  Push the head object onto the stack, this object contains all input boxes and will represent the root
  //  node of the object tree.
  //
  current_object->hierarchy_loc = hierarchy_start_ptr;
  current_object->boxes = boxes_start_ptr;
  current_object->scratch_boxes = scratch_boxes_start_ptr;
  current_object->num_boxes = num_boxes_start;
  current_object->index_array0_t = index_array0_t_start_ptr;
  current_object->index_array1_t = index_array1_t_start_ptr;
  current_object->index_array2_t = index_array2_t_start_ptr;
  current_object->scratch_index0_array = scratch_index0_array_start_ptr;
  current_object->scratch_index1_array = scratch_index1_array_start_ptr;
  current_object->scratch_index2_array = scratch_index2_array_start_ptr;
  current_object->unwind = 0;

  do {
    if(current_object->unwind) {
      //
      //  Stack unwinding step.  Compute the current node bounding box from the sum of the two child boxes.  Note, the last step
      //  of tree creation will be the unwinding of the root node, thus check for the terminal condition here.
      //
      ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_loc = current_object->hierarchy_loc;
      RangeBoxType& box1 = (hierarchy_loc+1)->m_box;
      RangeBoxType& box2 = (hierarchy_loc+hierarchy_loc->right_child_offset)->m_box;
      UnionBoxes(box1, box2, hierarchy_loc->m_box);
      if(stack_ptr != object_stack) {
        current_object = --stack_ptr;
        continue;
      } else {
        break;
      }
    } else {
      const int num_boxes = current_object->num_boxes;
      if(num_boxes == 2) {
        store_2_node_tree(current_object->hierarchy_loc, current_object->boxes);
        current_object = --stack_ptr;
        continue;
      } else if(num_boxes == 3) {
        store_3_node_tree(current_object->hierarchy_loc, current_object->boxes);
        current_object = --stack_ptr;
        continue;
      } else if (num_boxes <= 6){
        //
        //  Special case, do not update splitting directions, the unique cases for two and three boxes do not use them
        //  and splitting a N<=6 boxes yields only 3, or 2 length sub trees.
        //
        ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_loc = current_object->hierarchy_loc;
        int *const index_array0_t = current_object->index_array0_t;
        int *const index_array1_t = current_object->index_array1_t;
        int *const index_array2_t = current_object->index_array2_t;
        stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes = current_object->boxes;
        stk::search::ObjectBoundingBox_T<RangeBoxType> *const scratch_boxes = current_object->scratch_boxes;
        //
        //  There are more than 2 boxes, compute an optimal splitting direction for the boxes, and divide them into two halves.
        //  Compute the centroid of each component box, the centroid of the current encompasing box, and the moment of
        //  inertia of all particles in the box about (0,0,0)
        //
        const int right_child_size = num_boxes/2;
        const int left_child_size = num_boxes - right_child_size;
        float moment_x(0.0);
        float moment_y(0.0);
        float moment_z(0.0);
        split_boxes(num_boxes, boxes, moment_x, moment_y, moment_z);
        //
        //  Determine the longest centroid bounding box direction.  This is the direction in which bounding boxes will be
        //  sorted and split.  Reorder the box sorting arrays based on the split direction
        //
        int *primary_index;
        if(moment_x > moment_y && moment_x > moment_z) {
          primary_index = index_array0_t;
        } else {
          if(moment_y > moment_z) {
            primary_index = index_array1_t;
          } else {
            primary_index = index_array2_t;
          }
        }
        //
        //  Sort boxes into array scratch_boxes
        //  Create the new primary index array
        //
        for(int ibox = 0; ibox < num_boxes; ++ibox) {
          scratch_boxes[ibox] = boxes[primary_index[ibox]];
        }
        //
        //  Update the current tree pointer for right child offset
        //
        const int right_child_offset1 = left_child_size*2;
        hierarchy_loc->right_child_offset = right_child_offset1;
        //
        //  Primary computations on this object are complete, thus set the unwinding flag
        //
        current_object->unwind = 1;
        //
        //  Add the right child to the object stack
        //
        current_object = ++stack_ptr;
        current_object->hierarchy_loc = hierarchy_loc+right_child_offset1;
        current_object->boxes = scratch_boxes + left_child_size;
        current_object->num_boxes = right_child_size;
        current_object->unwind = 0;
        //
        //  Add the left child box to the stack,
        //
        current_object = ++stack_ptr;
        //
        //  Update the current object with the rest of the left objects data
        //
        current_object->hierarchy_loc = hierarchy_loc+1;
        current_object->boxes = scratch_boxes;
        current_object->num_boxes = left_child_size;
        current_object->unwind = 0;
        continue;
      } else {
        //
        //  This is the general stack functionallity for more than 6 boxes
        //

        ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_loc = current_object->hierarchy_loc;
        int *const index_array0_t = current_object->index_array0_t;
        int *const index_array1_t = current_object->index_array1_t;
        int *const index_array2_t = current_object->index_array2_t;
        int *const scratch_index0_array = current_object->scratch_index0_array;
        int *const scratch_index1_array = current_object->scratch_index1_array;
        int *const scratch_index2_array = current_object->scratch_index2_array;
        stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes = current_object->boxes;
        stk::search::ObjectBoundingBox_T<RangeBoxType> *const scratch_boxes = current_object->scratch_boxes;

        //
        //  There are more than 6 boxes, compute an optimal splitting direction for the boxes, and divide them into two halves.
        //  Compute the centroid of each component box, the centroid of the current encompasing box, and the moment of
        //  interial of all particles in the box about (0,0,0)
        //
        const int right_child_size = num_boxes/2;
        const int left_child_size = num_boxes - right_child_size;
        float moment_x(0.0);
        float moment_y(0.0);
        float moment_z(0.0);
        split_boxes(num_boxes, boxes, moment_x, moment_y, moment_z);
        //
        //  Determine the longest centroid bounding box direction.  This is the direction in which bounding boxes will be
        //  sorted and split.  Reorder the box sorting arrays based on the split direction
        //
        int *primary_index, *sub_index1, *sub_index2;
        int *primary_sindex, *sub_sindex1, *sub_sindex2;
     
        if(moment_x > moment_y && moment_x > moment_z) {
          primary_index = index_array0_t;
          sub_index1   = index_array1_t;
          sub_index2   = index_array2_t;
          primary_sindex = scratch_index0_array;
          sub_sindex1   = scratch_index1_array;
          sub_sindex2   = scratch_index2_array;
        } else {
          if(moment_y > moment_z) {
            primary_index = index_array1_t;
            sub_index1   = index_array0_t;
            sub_index2   = index_array2_t;
            primary_sindex = scratch_index1_array;
            sub_sindex1   = scratch_index0_array;
            sub_sindex2   = scratch_index2_array;
          } else {
            primary_index = index_array2_t;
            sub_index1   = index_array0_t;
            sub_index2   = index_array1_t;
            primary_sindex = scratch_index2_array;
            sub_sindex1   = scratch_index0_array;
            sub_sindex2   = scratch_index1_array;
          }
        }
        //
        //  Sort boxes into array scratch_boxes
        //  Create the new primary index array
        //
        for(int ibox = 0; ibox < num_boxes; ++ibox) {
          scratch_boxes[ibox] = boxes[primary_index[ibox]];
          primary_sindex[primary_index[ibox]] = ibox;
        }
        //
        //  Reorder secondary arrays to be consistent with the primary array splitting
        //
        int left_pos1 = 0;
        int right_pos1 = left_child_size;
        int left_pos2 = 0;
        int right_pos2 = left_child_size;
        for(int ibox = 0; ibox < num_boxes; ++ibox) {
          const int primary1 = primary_sindex[sub_index1[ibox]];
          const int primary2 = primary_sindex[sub_index2[ibox]];

          if(primary1 < left_child_size) {
            sub_sindex1[left_pos1++] = primary1;
          } else {
            sub_sindex1[right_pos1++] = primary1 - left_child_size;
          }

          if(primary2 < left_child_size) {
            sub_sindex2[left_pos2++] = primary2;
          } else {
            sub_sindex2[right_pos2++] = primary2 - left_child_size;
          }
        }
        //
        //  Update the primary index to be correct with the new box ordering
        //
        for(int ibox = 0; ibox < left_child_size; ++ibox) {
          primary_sindex[ibox] = ibox;
        }
        for(int ibox = left_child_size; ibox < num_boxes; ++ibox) {
          primary_sindex[ibox] = ibox - left_child_size;
        }

        //
        //  Update the current tree pointer for right child offset
        //
        const int right_child_offset1 = left_child_size*2;
        hierarchy_loc->right_child_offset = right_child_offset1;
        //
        //  Primary computations on this object are complete, thus set the unwinding flag
        //
        current_object->unwind = 1;
        //
        //  Add the right child to the object stack
        //
        current_object = ++stack_ptr;
        current_object->hierarchy_loc = hierarchy_loc + right_child_offset1;
        current_object->boxes = scratch_boxes + left_child_size;
        current_object->scratch_boxes = boxes + left_child_size;
        current_object->num_boxes = right_child_size;
        current_object->index_array0_t = scratch_index0_array + left_child_size;
        current_object->index_array1_t = scratch_index1_array + left_child_size;
        current_object->index_array2_t = scratch_index2_array + left_child_size;
        current_object->scratch_index0_array = index_array0_t + left_child_size;
        current_object->scratch_index1_array = index_array1_t + left_child_size;
        current_object->scratch_index2_array = index_array2_t + left_child_size;
        current_object->unwind = 0;
        //
        //  Add the left child box to the stack,
        //
        current_object = ++stack_ptr;
        //
        //  Update the current object with the rest of the left objects data
        //
        current_object->hierarchy_loc = hierarchy_loc+1;
        current_object->boxes = scratch_boxes;
        current_object->scratch_boxes = boxes;
        current_object->num_boxes = left_child_size;
        current_object->index_array0_t = scratch_index0_array;
        current_object->index_array1_t = scratch_index1_array;
        current_object->index_array2_t = scratch_index2_array;
        current_object->scratch_index0_array = index_array0_t;
        current_object->scratch_index1_array = index_array1_t;
        current_object->scratch_index2_array = index_array2_t;
        current_object->unwind = 0;
        continue;
      }
    }
  } while(true);
}

//
// Threaded hierarchy creation routine, basically drill down the tree until have created enough sub-limbs to efficently 
// thread then pass each limb down to a threaded
//
#ifdef _OPENMP
template <typename RangeBoxType>
void create_hierarchy_loop_threaded(ObjectBoundingBoxHierarchy_T<RangeBoxType> *const hierarchy_start_ptr,
                                                                stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes_start_ptr,
                                                                stk::search::ObjectBoundingBox_T<RangeBoxType> *const scratch_boxes_start_ptr,
                                                                const int num_boxes_start,
                                                                int *const index_array0_t_start_ptr,
                                                                int *const index_array1_t_start_ptr,
                                                                int *const index_array2_t_start_ptr,
                                                                int* const scratch_index0_array_start_ptr,
                                                                int* const scratch_index1_array_start_ptr,
                                                                int* const scratch_index2_array_start_ptr) {

  unsigned maxNumThread = omp_get_max_threads();
  if(maxNumThread == 1) {
    create_hierarchy_loop(hierarchy_start_ptr, boxes_start_ptr, scratch_boxes_start_ptr, num_boxes_start,
                          index_array0_t_start_ptr, index_array1_t_start_ptr, index_array2_t_start_ptr,
                          scratch_index0_array_start_ptr, scratch_index1_array_start_ptr, scratch_index2_array_start_ptr);
    return;
  }


  //  The algorithm is optimized for thread counts that are powers of 2, for other counts the algorithm will function
  //  but the final parallel tree creation step will be unbalanced and sub-optimal.
  //
  unsigned numLimb = 1;
  while(numLimb < maxNumThread) {
    numLimb*=2;
  }

  //
  //  Subdivide tree subdivision times
  //

  std::vector<ObjectBoundingBoxHierarchy_T<RangeBoxType>*> limbHier        (numLimb);
  std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType>*>          limbBoxes       (numLimb);
  std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType>*>          limbScratchBoxes(numLimb);
  std::vector<int>                         limbNumBoxes    (numLimb);
  std::vector<int*>                        limbIndex0      (numLimb);
  std::vector<int*>                        limbIndex1      (numLimb);
  std::vector<int*>                        limbIndex2      (numLimb);
  std::vector<int*>                        limbSIndex0      (numLimb);
  std::vector<int*>                        limbSIndex1      (numLimb);
  std::vector<int*>                        limbSIndex2      (numLimb);
  create_hierarchy_partial(hierarchy_start_ptr, boxes_start_ptr, scratch_boxes_start_ptr, num_boxes_start,
                           index_array0_t_start_ptr,
                           index_array1_t_start_ptr,
                           index_array2_t_start_ptr,
                           scratch_index0_array_start_ptr,
                           scratch_index1_array_start_ptr,
                           scratch_index2_array_start_ptr,
                           numLimb,
                           0, limbHier.data(), limbBoxes.data(), limbScratchBoxes.data(), limbNumBoxes.data(),
                           limbIndex0.data(), limbIndex1.data(), limbIndex2.data(),
                           limbSIndex0.data(), limbSIndex1.data(), limbSIndex2.data());
#pragma omp parallel for default(none) shared(numLimb, limbHier, limbBoxes, limbScratchBoxes, limbNumBoxes, limbIndex0, limbIndex1, limbIndex2, limbSIndex0, limbSIndex1, limbSIndex2)
  for(unsigned i=0; i<numLimb; ++i) {
    if(limbHier[i] == nullptr) continue;
    create_hierarchy_loop(limbHier[i], limbBoxes[i], limbScratchBoxes[i], limbNumBoxes[i],
                          limbIndex0[i], limbIndex1[i], limbIndex2[i], 
                          limbSIndex0[i], limbSIndex1[i], limbSIndex2[i]);
  }

  complete_tree_partial(hierarchy_start_ptr,numLimb);
}
#endif

//------------------------------------------------------------------------------------------------------
template <typename RangeBoxType>
inline void create_hierarchy(ObjectBoundingBoxHierarchy_T<RangeBoxType> *hierarchy_data,
                                                  stk::search::ObjectBoundingBox_T<RangeBoxType> *boxes,
                                                  const int num_boxes) {
  //
  //  Generate an initial sorting of the input boxes, sort in each of the three coordinate directions independently
  //  Iterating though the box list in the order given by the 'index' arrays will return the boxes in each of the 
  //  three sorted orders.
  //
  std::vector< std::pair<float,int> > centroid_x(num_boxes);
  std::vector< std::pair<float,int> > centroid_y(num_boxes);
  std::vector< std::pair<float,int> > centroid_z(num_boxes);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static) default(none) shared(centroid_x, centroid_y, centroid_z, boxes) firstprivate(num_boxes)
#endif
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    float centroid[3];
    centroid[0] = boxes[ibox].GetBox().get_x_min() + boxes[ibox].GetBox().get_x_max();
    centroid[1] = boxes[ibox].GetBox().get_y_min() + boxes[ibox].GetBox().get_y_max();
    centroid[2] = boxes[ibox].GetBox().get_z_min() + boxes[ibox].GetBox().get_z_max();
    // you might read the next three lines and want to use std::make_pair instead
    // please do not do this: causes stack use after scope error due to some
    // behavior of std::make_pair (c++11) that I do not understand
    centroid_x[ibox].first = centroid[0]; centroid_x[ibox].second = ibox;
    centroid_y[ibox].first = centroid[1]; centroid_y[ibox].second = ibox;
    centroid_z[ibox].first = centroid[2]; centroid_z[ibox].second = ibox;
  }
  stk::search::threadedSort(centroid_x);
  stk::search::threadedSort(centroid_y);
  stk::search::threadedSort(centroid_z);

  std::vector<int> index_array0(num_boxes);
  std::vector<int> index_array1(num_boxes);
  std::vector<int> index_array2(num_boxes);

#ifdef _OPENMP

#pragma omp parallel for schedule(static) default(none) shared(centroid_x, centroid_y, centroid_z, index_array0, index_array1, index_array2) firstprivate(num_boxes)
#endif
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    index_array0[ibox] = centroid_x[ibox].second;
    index_array1[ibox] = centroid_y[ibox].second;
    index_array2[ibox] = centroid_z[ibox].second;
  }
  //
  //  Scratch workspace for the search tree creation algorithm
  //
  std::vector<int> scratch_index0_array(num_boxes);
  std::vector<int> scratch_index1_array(num_boxes);
  std::vector<int> scratch_index2_array(num_boxes);
  std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType>> scratch_boxes(num_boxes);
  //
  //  Organize the boxes into a searchable hierarchy
  //
#ifdef _OPENMP


  int maxThreadCount= omp_get_max_threads();

  if(maxThreadCount <= 8) {
    //  This method avoids thread nesting and the potential overhead it inccurs.  Overall it is theoretically not 
    //  all that scalable and is targeted to relatively low thread counts.
    create_hierarchy_loop_threaded(hierarchy_data,
                                                   boxes,
                                                   scratch_boxes.data(),
                                                   num_boxes,
                                                   index_array0.data(),
                                                   index_array1.data(),
                                                   index_array2.data(),
                                                   scratch_index0_array.data(),
                                                   scratch_index1_array.data(),
                                                   scratch_index2_array.data());
  } else {
    //  This method attempts to exploit maximum available parallelism.  However, the method requires
    //  nested threading which can have high overhead on some systems, only use if a largeish number
    //  of threads is available.
    int origNesting = omp_get_nested();
    omp_set_nested(1);
    create_hierarchy_fork_threaded(hierarchy_data,
                                                   boxes,
                                                   scratch_boxes.data(),
                                                   num_boxes,
                                                   index_array0.data(),
                                                   index_array1.data(),
                                                   index_array2.data(),
                                                   scratch_index0_array.data(),
                                                   scratch_index1_array.data(),
                                                   scratch_index2_array.data(),
                                                   maxThreadCount);

    omp_set_nested(origNesting);
  }

#else 
  //  Pure serial tree creation
  create_hierarchy_loop(hierarchy_data,
                        boxes,
                        scratch_boxes.data(),
                        num_boxes,
                        index_array0.data(),
                        index_array1.data(),
                        index_array2.data(),
                        scratch_index0_array.data(),
                        scratch_index1_array.data(),
                        scratch_index2_array.data());
#endif
}

#ifdef _OPENMP
//
//  Divide up the first few nodes (subdivisionLevel) of the hierarchy, but do NOT set their bounding boxes.  The remainder
//  of the hierarchy can then be computed in thread separate computations.
//  
template <typename RangeBoxType>
void create_hierarchy_partial(ObjectBoundingBoxHierarchy_T<RangeBoxType> *const hierarchy_start_ptr,
                                                          stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes_start_ptr,
                                                          stk::search::ObjectBoundingBox_T<RangeBoxType> *const scratch_boxes_start_ptr,
                                                          const int num_boxes,
                                                          int *const index_array0_t_start_ptr,
                                                          int *const index_array1_t_start_ptr,
                                                          int *const index_array2_t_start_ptr,
                                                          int *const sindex_array0_t_start_ptr,
                                                          int *const sindex_array1_t_start_ptr,
                                                          int *const sindex_array2_t_start_ptr,
                                                          unsigned subdivisionLevel,
                                                          unsigned mySubTree,
                                                          ObjectBoundingBoxHierarchy_T<RangeBoxType>** subTreeHier,
                                                          stk::search::ObjectBoundingBox_T<RangeBoxType>**          subTreeBoxStart,
                                                          stk::search::ObjectBoundingBox_T<RangeBoxType>**          subTreeBoxScratchStart,
                                                          int*                         subTreeNumBoxes,
                                                          int**                        subTreeIndex0,
                                                          int**                        subTreeIndex1,
                                                          int**                        subTreeIndex2,
                                                          int**                        subTreeSIndex0,
                                                          int**                        subTreeSIndex1,
                                                          int**                        subTreeSIndex2) {



  assert(subdivisionLevel    > 0);
  //
  //  If we have reached the leaf prior to hitting the subdivision level
  //  The current sub-trees will be nullptr and the stack will be unwound
  //
  if(num_boxes == 1) {
    store_1_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    for(unsigned i=0; i<subdivisionLevel; ++i) {
      subTreeHier[mySubTree + i] = nullptr;
    }
    return;
  }
  //
  //  If we have reached the terminal subdivision level, we are done, the tree data will be stored and passed to the next
  //  communication round.
  //
  if(subdivisionLevel == 1) {
    store_1_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    subTreeHier           [mySubTree] = hierarchy_start_ptr;
    subTreeBoxStart       [mySubTree] = boxes_start_ptr;
    subTreeBoxScratchStart[mySubTree] = scratch_boxes_start_ptr;
    subTreeNumBoxes       [mySubTree] = num_boxes;
    subTreeIndex0         [mySubTree] = index_array0_t_start_ptr;
    subTreeIndex1         [mySubTree] = index_array1_t_start_ptr;
    subTreeIndex2         [mySubTree] = index_array2_t_start_ptr;
    subTreeSIndex0        [mySubTree] = sindex_array0_t_start_ptr;
    subTreeSIndex1        [mySubTree] = sindex_array1_t_start_ptr;
    subTreeSIndex2        [mySubTree] = sindex_array2_t_start_ptr;
    return;
  }

  assert(subdivisionLevel%2 == 0);

  //
  //  Continue division of the boxes.
  //  Compute an optimal splitting direction for the boxes, and divide them into two halves.
  //  Compute the centroid of each component box, the centroid of the current encompassing box, and the moment of
  //  inertia of all particles in the box about (0,0,0)
  //
  const int right_child_size = num_boxes/2;
  const int left_child_size = num_boxes - right_child_size;
  //
  //  Determine splitting direction, use moments of inertia
  //

  float moment_x(0.0);
  float moment_y(0.0);
  float moment_z(0.0);
  split_boxes_threaded(num_boxes, boxes_start_ptr, moment_x, moment_y, moment_z);

  //
  //  Determine the longest centroid bounding box direction.  This is the direction in which bounding boxes will be
  //  sorted and split.  Reorder the box sorting arrays based on the split direction
  //
  int *primary_index, *sub_index1, *sub_index2;
  int *primary_sindex;
  int *sub_sindex1, *sub_sindex2;
  if(moment_x > moment_y && moment_x > moment_z) {
    primary_index  = index_array0_t_start_ptr;
    sub_index1    = index_array1_t_start_ptr;
    sub_index2    = index_array2_t_start_ptr;
    primary_sindex = sindex_array0_t_start_ptr;
    sub_sindex1   = sindex_array1_t_start_ptr;
    sub_sindex2   = sindex_array2_t_start_ptr;
  } else {
    if(moment_y > moment_z) {
      primary_index  = index_array1_t_start_ptr;
      sub_index1    = index_array0_t_start_ptr;
      sub_index2    = index_array2_t_start_ptr;
      primary_sindex = sindex_array1_t_start_ptr;
      sub_sindex1   = sindex_array0_t_start_ptr;
      sub_sindex2   = sindex_array2_t_start_ptr;
    } else {
      primary_index  = index_array2_t_start_ptr;
      sub_index1    = index_array0_t_start_ptr;
      sub_index2    = index_array1_t_start_ptr;
      primary_sindex = sindex_array2_t_start_ptr;
      sub_sindex1   = sindex_array0_t_start_ptr;
      sub_sindex2   = sindex_array1_t_start_ptr;
    }
  }
  //
  //  Resort the boxes based off of the primary array.  First use the new primary_sindex array to story a mapping between old box 
  //  positions and new box positions.
  //
#pragma omp parallel for schedule(static) default(none) shared(primary_index, primary_sindex) firstprivate(scratch_boxes_start_ptr, boxes_start_ptr, num_boxes)
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    //  NKC, try memcopy
    const int index = primary_index[ibox];
    scratch_boxes_start_ptr[ibox] = boxes_start_ptr[index];
    primary_sindex[index] = ibox;
  }
  //
  //  Reorder secondary arrays to be consistent with the primary array order
  //  Use at most just four threads here, one starting from the left and one from the right,
  //  can't think of any better way to exploit threading for this particular case
  //
#pragma omp parallel sections default(none) shared(primary_sindex, sub_index1, sub_sindex1, sub_index2, sub_sindex2) firstprivate(left_child_size, num_boxes)
  {
  #pragma omp section
    {
      //  Section 1, process the first half of the sub_index lists, this will fill the output arrays forwards
      int left_pos1 = 0;
      int right_pos1 = left_child_size;
      for(int ibox = 0; ibox < left_child_size; ++ibox) {
        const int primary1 = primary_sindex[sub_index1[ibox]];
        if(primary1 < left_child_size) {
          sub_sindex1[left_pos1++] = primary1;
        } else {
          sub_sindex1[right_pos1++] = primary1 - left_child_size;
        }
      }      
    }
  #pragma omp section
    {
      //  Section 2, process the first half of the sub_index lists, this will fill the output arrays forwards
      int left_pos2 = 0;
      int right_pos2 = left_child_size;
      for(int ibox = 0; ibox < left_child_size; ++ibox) {
        const int primary2 = primary_sindex[sub_index2[ibox]];
        if(primary2 < left_child_size) {
          sub_sindex2[left_pos2++] = primary2;
        } else {
          sub_sindex2[right_pos2++] = primary2 - left_child_size;
        }
      }      
    }
  #pragma omp section
    {
      //  Section 3, process the second half of the sub index lists this will fill the output arrays backwards,
      //  sections 1 and 3 should meet exactly in the middle
      int left_pos1 = left_child_size-1;
      int right_pos1 = num_boxes-1;
      for(int ibox = num_boxes-1; ibox >= left_child_size; --ibox) {
        const int primary1 = primary_sindex[sub_index1[ibox]];
        if(primary1 < left_child_size) {
          sub_sindex1[left_pos1--] = primary1;
        } else {
          sub_sindex1[right_pos1--] = primary1 - left_child_size;
        }
      }      
    }
#pragma omp section
    {
      //  Section 4, process the second half of the sub index lists this will fill the output arrays backwards,
      //  sections 2 and 4 should meet exactly in the middle
      int left_pos2 = left_child_size-1;
      int right_pos2 = num_boxes-1;
      for(int ibox = num_boxes-1; ibox >= left_child_size; --ibox) {
        const int primary2 = primary_sindex[sub_index2[ibox]];
        if(primary2 < left_child_size) {
          sub_sindex2[left_pos2--] = primary2;
        } else {
          sub_sindex2[right_pos2--] = primary2 - left_child_size;
        }
      }      
    }
  }
  //
  //  Update the primary index to be correct with the new box ordering
  //
#pragma omp parallel for schedule(static) default(none) shared(primary_sindex) firstprivate(left_child_size)
  for(int ibox = 0; ibox < left_child_size; ++ibox) {
    primary_sindex[ibox] = ibox;
  }
#pragma omp parallel for schedule(static) default(none) shared(primary_sindex) firstprivate(left_child_size, num_boxes)
  for(int ibox = left_child_size; ibox < num_boxes; ++ibox) {
    primary_sindex[ibox] = ibox - left_child_size;
  }
  //
  //  Update the current tree pointer for right child offset
  //
  const int right_child_offset1 = left_child_size*2;
  hierarchy_start_ptr->right_child_offset = right_child_offset1;

  //
  //  Recusively call on the left and right box sets
  //
  create_hierarchy_partial(hierarchy_start_ptr+1,
                           scratch_boxes_start_ptr,
                           boxes_start_ptr,
                           left_child_size,
                           sindex_array0_t_start_ptr,
                           sindex_array1_t_start_ptr,
                           sindex_array2_t_start_ptr,
                           index_array0_t_start_ptr,
                           index_array1_t_start_ptr,
                           index_array2_t_start_ptr,
                           subdivisionLevel/2,
                           mySubTree,
                           subTreeHier,
                           subTreeBoxStart,
                           subTreeBoxScratchStart,
                           subTreeNumBoxes,
                           subTreeIndex0,
                           subTreeIndex1,
                           subTreeIndex2,
                           subTreeSIndex0,
                           subTreeSIndex1,
                           subTreeSIndex2);

  create_hierarchy_partial(hierarchy_start_ptr      + right_child_offset1,
                           scratch_boxes_start_ptr  + left_child_size,
                           boxes_start_ptr          + left_child_size,
                           right_child_size,
                           sindex_array0_t_start_ptr + left_child_size,
                           sindex_array1_t_start_ptr + left_child_size,
                           sindex_array2_t_start_ptr + left_child_size,
                           index_array0_t_start_ptr + left_child_size,
                           index_array1_t_start_ptr + left_child_size,
                           index_array2_t_start_ptr + left_child_size,
                           subdivisionLevel/2,
                           mySubTree                + subdivisionLevel/2,
                           subTreeHier,
                           subTreeBoxStart,
                           subTreeBoxScratchStart,
                           subTreeNumBoxes,
                           subTreeIndex0,
                           subTreeIndex1,
                           subTreeIndex2,
                           subTreeSIndex0,
                           subTreeSIndex1,
                           subTreeSIndex2);

}
#endif




#ifdef _OPENMP
//
//  Alternative algorithm to threaded tree creation based on forks.  Possibly ideal for large thread counts, also 
//  forking is a natural and easier to understand way to write the recursive algorthm.
//
template <typename RangeBoxType>
inline void create_hierarchy_fork_threaded(ObjectBoundingBoxHierarchy_T<RangeBoxType> *const hierarchy_start_ptr,
                                                                stk::search::ObjectBoundingBox_T<RangeBoxType> *const boxes_start_ptr,
                                                                stk::search::ObjectBoundingBox_T<RangeBoxType> *const scratch_boxes_start_ptr,
                                                                const int num_boxes,
                                                                int *const index_array0_t_start_ptr,
                                                                int *const index_array1_t_start_ptr,
                                                                int *const index_array2_t_start_ptr,
                                                                int* const sindex_array0_t_start_ptr,
                                                                int* const sindex_array1_t_start_ptr,
                                                                int* const sindex_array2_t_start_ptr,
                                                                int numThreadsToUse) { 

  if(numThreadsToUse == 1) {
    //
    //  Terminal case, down to one thread just run serial algorithm on remaining tree branch avoiding any thread
    //  logic overhead.
    //
    create_hierarchy_loop(hierarchy_start_ptr, boxes_start_ptr, scratch_boxes_start_ptr, num_boxes,
                          index_array0_t_start_ptr, index_array1_t_start_ptr, index_array2_t_start_ptr,
                          sindex_array0_t_start_ptr, sindex_array1_t_start_ptr, sindex_array2_t_start_ptr);
    return;
  }

  //
  //  Reached the leaf node store an exit
  //
  if(num_boxes == 1) {
    store_1_node_tree(hierarchy_start_ptr, boxes_start_ptr);
    return;
  }
  //
  //  Continue division of the boxes.
  //  Compute an optimal splitting direction for the boxes, and divide them into two halves.
  //  Division direction is based on approximation of the moments of inertial of the box
  //  sets.
  //
  const int right_child_size = num_boxes/2;
  const int left_child_size = num_boxes - right_child_size;

  float moment_x(0.0);
  float moment_y(0.0);
  float moment_z(0.0);

  split_boxes_threaded(num_boxes, boxes_start_ptr, moment_x, moment_y, moment_z, numThreadsToUse);
  //
  //  Boxes will be split in half in the 'primary' direction.  Sub1 and sub2 are the other two directions.
  //
  int *primary_index, *sub_index1, *sub_index2;
  int *primary_sindex;
  int *sub_sindex1, *sub_sindex2;
  if(moment_x > moment_y && moment_x > moment_z) {
    primary_index  = index_array0_t_start_ptr;
    sub_index1    = index_array1_t_start_ptr;
    sub_index2    = index_array2_t_start_ptr;
    primary_sindex = sindex_array0_t_start_ptr;
    sub_sindex1   = sindex_array1_t_start_ptr;
    sub_sindex2   = sindex_array2_t_start_ptr;
  } else {
    if(moment_y > moment_z) {
      primary_index  = index_array1_t_start_ptr;
      sub_index1    = index_array0_t_start_ptr;
      sub_index2    = index_array2_t_start_ptr;
      primary_sindex = sindex_array1_t_start_ptr;
      sub_sindex1   = sindex_array0_t_start_ptr;
      sub_sindex2   = sindex_array2_t_start_ptr;
    } else {
      primary_index  = index_array2_t_start_ptr;
      sub_index1    = index_array0_t_start_ptr;
      sub_index2    = index_array1_t_start_ptr;
      primary_sindex = sindex_array2_t_start_ptr;
      sub_sindex1   = sindex_array0_t_start_ptr;
      sub_sindex2   = sindex_array1_t_start_ptr;
    }
  }
  //
  //  Resort the boxes based off of the primary array.  First use the new primary_sindex array to story a mapping between old box 
  //  positions and new box positions.  Partion the 'sub' lists to the left and right sub trees.  Note the sub lists still
  //  need to remain in sorted order and are partioned with the below code.
  //

#pragma omp parallel for schedule(static) num_threads(numThreadsToUse) default(none) shared(primary_index, primary_sindex) firstprivate(num_boxes, scratch_boxes_start_ptr, boxes_start_ptr)
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    const int index = primary_index[ibox];
    scratch_boxes_start_ptr[ibox] = boxes_start_ptr[index];
    primary_sindex[index] = ibox;
  }
  //
  //  Reorder secondary arrays to be consistent with the primary array order
  //  Use at most just four threads here, one starting fro the left and one from the right of each of the lists
  //  can't think of any better way to exploit threading any more for this particular operation
  //
#pragma omp parallel sections num_threads(numThreadsToUse) default(none) shared(primary_sindex, sub_index1, sub_sindex1, sub_index2, sub_sindex2) firstprivate(left_child_size, num_boxes)
  {
  #pragma omp section
    {
      //  Section 1, process the first half of the sub_index lists, this will fill the output arrays forwards
      int left_pos1 = 0;
      int right_pos1 = left_child_size;
      for(int ibox = 0; ibox < left_child_size; ++ibox) {
        const int primary1 = primary_sindex[sub_index1[ibox]];
        if(primary1 < left_child_size) {
          sub_sindex1[left_pos1++] = primary1;
        } else {
          sub_sindex1[right_pos1++] = primary1 - left_child_size;
        }
      }      
    }
  #pragma omp section
    {
      //  Section 2, process the first half of the sub_index lists, this will fill the output arrays forwards
      int left_pos2 = 0;
      int right_pos2 = left_child_size;
      for(int ibox = 0; ibox < left_child_size; ++ibox) {
        const int primary2 = primary_sindex[sub_index2[ibox]];
        if(primary2 < left_child_size) {
          sub_sindex2[left_pos2++] = primary2;
        } else {
          sub_sindex2[right_pos2++] = primary2 - left_child_size;
        }
      }      
    }
  #pragma omp section
    {
      //  Section 3, process the second half of the sub index lists this will fill the output arrays backwards,
      //  sections 1 and 3 should meet exactly in the middle
      int left_pos1 = left_child_size-1;
      int right_pos1 = num_boxes-1;
      for(int ibox = num_boxes-1; ibox >= left_child_size; --ibox) {
        const int primary1 = primary_sindex[sub_index1[ibox]];
        if(primary1 < left_child_size) {
          sub_sindex1[left_pos1--] = primary1;
        } else {
          sub_sindex1[right_pos1--] = primary1 - left_child_size;
        }
      }      
    }
  #pragma omp section
    {
      //  Section 4, process the second half of the sub index lists this will fill the output arrays backwards,
      //  sections 2 and 4 should meet exactly in the middle
      int left_pos2 = left_child_size-1;
      int right_pos2 = num_boxes-1;
      for(int ibox = num_boxes-1; ibox >= left_child_size; --ibox) {
        const int primary2 = primary_sindex[sub_index2[ibox]];
        if(primary2 < left_child_size) {
          sub_sindex2[left_pos2--] = primary2;
        } else {
          sub_sindex2[right_pos2--] = primary2 - left_child_size;
        }
      }      
    }
  }
  //
  //  Update the primary index to be correct with the new box ordering
  //
#pragma omp parallel for schedule(static) num_threads(numThreadsToUse) default(none) shared(primary_sindex) firstprivate(left_child_size)
  for(int ibox = 0; ibox < left_child_size; ++ibox) {
    primary_sindex[ibox] = ibox;
  }
#pragma omp parallel for schedule(static) num_threads(numThreadsToUse) default(none) shared(primary_sindex) firstprivate(left_child_size, num_boxes)
  for(int ibox = left_child_size; ibox < num_boxes; ++ibox) {
    primary_sindex[ibox] = ibox - left_child_size;
  }
  //
  //  Update the current tree pointer for right child offset
  //
  const int right_child_offset1 = left_child_size*2;
  hierarchy_start_ptr->right_child_offset = right_child_offset1;

  //
  //  Recusively call on the left and right box sets.  Fork the current thread into two new primary threads.
  //  Each fork will then create new thread teams in the subsequent call.
  //
  int rightNumThreads = numThreadsToUse/2;
  int leftNumThreads = numThreadsToUse-rightNumThreads;

#pragma omp parallel sections num_threads(2) default(none), shared(leftNumThreads, rightNumThreads) firstprivate(hierarchy_start_ptr, scratch_boxes_start_ptr, boxes_start_ptr, left_child_size, sindex_array0_t_start_ptr, sindex_array1_t_start_ptr, sindex_array2_t_start_ptr, index_array0_t_start_ptr, index_array1_t_start_ptr, index_array2_t_start_ptr, right_child_offset1, right_child_size)
  {
#pragma omp section
    {
      create_hierarchy_fork_threaded(hierarchy_start_ptr+1,
                               scratch_boxes_start_ptr,
                               boxes_start_ptr,
                               left_child_size,
                               sindex_array0_t_start_ptr,
                               sindex_array1_t_start_ptr,
                               sindex_array2_t_start_ptr,
                               index_array0_t_start_ptr,
                               index_array1_t_start_ptr,
                               index_array2_t_start_ptr,
                               leftNumThreads);
    }
#pragma omp section
    {
      create_hierarchy_fork_threaded(hierarchy_start_ptr      + right_child_offset1,
                               scratch_boxes_start_ptr  + left_child_size,
                               boxes_start_ptr          + left_child_size,
                               right_child_size,
                               sindex_array0_t_start_ptr + left_child_size,
                               sindex_array1_t_start_ptr + left_child_size,
                               sindex_array2_t_start_ptr + left_child_size,
                               index_array0_t_start_ptr + left_child_size,
                               index_array1_t_start_ptr + left_child_size,
                               index_array2_t_start_ptr + left_child_size,
                               rightNumThreads);
    }
  }

  RangeBoxType& box1 = (hierarchy_start_ptr+1)->m_box;
  RangeBoxType& box2 = (hierarchy_start_ptr+right_child_offset1)->m_box;
  UnionBoxes(box1, box2,  hierarchy_start_ptr->m_box);
}
#endif

//---------------------------------------------------------------------------------------------------------------------

template<typename RangeBoxType>
template<typename DomainBox>
inline void ProximitySearchTree_T<RangeBoxType>::SearchForOverlap(const DomainBox& searchObject, std::vector<int>& returnList) const {
  returnList.clear();
  if(m_tree.empty()) {
    return;
  } else {
    ObjectBoundingBoxHierarchy_T<RangeBoxType> const* current_object  = m_tree.data();
    ObjectBoundingBoxHierarchy_T<RangeBoxType> const* object_stack[MAX_TREE_LEVELS];
    const ObjectBoundingBoxHierarchy_T<RangeBoxType> ** stack_ptr = object_stack;
    do {
      if(intersects(searchObject, current_object->m_box)) {
        const int right_child_offset1 = current_object->right_child_offset;
        if(right_child_offset1 > 0) {
          *(stack_ptr++) = current_object++ + right_child_offset1;
          continue;
        }
        returnList.emplace_back(-right_child_offset1);
      }
      if(stack_ptr != object_stack) {
        current_object = *(--stack_ptr);
        continue;
      }
      break;
    } while(true);
  }
}

 template<typename RangeBoxType>
 template<typename DomainBox>
   inline void ProximitySearchTree_T<RangeBoxType>::SearchForOverlap(const DomainBox&           searchObject,
                                                                                          std::vector<int>&          returnIndexList,
                                                                                          std::vector<RangeBoxType>& returnBoxList) const {
   returnIndexList.clear();
   returnBoxList.clear();
   if(m_tree.empty()) {
     return;
   } else {
     ObjectBoundingBoxHierarchy_T<RangeBoxType> const* current_object  = m_tree.data();
     ObjectBoundingBoxHierarchy_T<RangeBoxType> const* object_stack[MAX_TREE_LEVELS];
     const ObjectBoundingBoxHierarchy_T<RangeBoxType> ** stack_ptr = object_stack;
     do {
       if(intersects(searchObject, current_object->m_box)) {
         const int right_child_offset1 = current_object->right_child_offset;
         if(right_child_offset1 > 0) {
           *(stack_ptr++) = current_object++ + right_child_offset1;
           continue;
         }
         returnIndexList.emplace_back(-right_child_offset1);
         returnBoxList.emplace_back(current_object->m_box);
       }
       if(stack_ptr != object_stack) {
         current_object = *(--stack_ptr);
         continue;
       }
       break;
     } while(true);
   }
 }

template<typename RangeBoxType>
template<typename DomainBox>
inline bool ProximitySearchTree_T<RangeBoxType>::AnyOverlap(const DomainBox& searchObject) const {
  if(m_tree.empty()) {
    return false;
  } else {
    ObjectBoundingBoxHierarchy_T<RangeBoxType> const* current_object  = m_tree.data();
    ObjectBoundingBoxHierarchy_T<RangeBoxType> const* object_stack[MAX_TREE_LEVELS];
    const ObjectBoundingBoxHierarchy_T<RangeBoxType> ** stack_ptr = object_stack;
    do {
      if(intersects(searchObject, current_object->m_box)) {
        const int right_child_offset1 = current_object->right_child_offset;
        if(right_child_offset1 > 0) {
          *(stack_ptr++) = current_object++ + right_child_offset1;
          continue;
        }
        return true;
      }
      if(stack_ptr != object_stack) {
        current_object = *(--stack_ptr);
        continue;
      }
      return false;
    } while(true);

  }
}

template<typename RangeBoxType>
ProximitySearchTree_T<RangeBoxType>::ProximitySearchTree_T(const std::vector<RangeBoxType>& inputBoxes) {
  unsigned numBox = inputBoxes.size();
  std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> > tempBoxes;
  tempBoxes.resize(numBox);
  for(unsigned int ibox = 0; ibox < numBox; ++ibox) {
    tempBoxes[ibox].SetBox(inputBoxes[ibox]);
    tempBoxes[ibox].set_object_number(ibox);
  }
  InitializeSearch(tempBoxes);
}



template<typename RangeBoxType>
ProximitySearchTree_T<RangeBoxType>::ProximitySearchTree_T(std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> >& inputBoxes) {
  InitializeSearch(inputBoxes);
}

//
//  Fill out the tree based on bounding boxes
//
template<typename RangeBoxType>
void ProximitySearchTree_T<RangeBoxType>::InitializeSearch(std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> >& inputBoxes) {
  unsigned numBox = inputBoxes.size();
  if(numBox != 0) {
    m_tree.resize(2*numBox-1);
    create_hierarchy(m_tree.data(), inputBoxes.data(), numBox);
  } else {
    m_tree.clear();
  }
}

//
// Update exiting tree with (subtly) differnt coordinates
//
template<typename RangeBoxType>
int ProximitySearchTree_T<RangeBoxType>::UpdateSearch(std::vector<stk::search::ObjectBoundingBox_T<RangeBoxType> >& updateBoxes) {
  unsigned numTreeBoxes;
  if(m_tree.size() == 0) {
    numTreeBoxes = 0;
  } else {
    numTreeBoxes = (m_tree.size() + 1)/2;
  }
  if(updateBoxes.size() != numTreeBoxes) {
    return -1;
  }
  update_tree(m_tree.data(), updateBoxes.data());
  return 0;
}


  }

}

#endif

