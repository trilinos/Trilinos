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
// 

#ifndef STK_STK_SEARCH_STK_SEARCH_KDTREE_BOUNDINGSPHERES_HPP_
#define STK_STK_SEARCH_STK_SEARCH_KDTREE_BOUNDINGSPHERES_HPP_

#include <cmath>
#include <Kokkos_Core.hpp>
#include <stk_search/Sphere.hpp>
#include <stk_search/kdtree/KDTree.hpp>

namespace stk {
  namespace search {


template <typename NumT>
KOKKOS_FORCEINLINE_FUNCTION
void UnionSpheres(const Sphere<NumT> & sphereA, const Sphere<NumT> & sphereB, Sphere<NumT> & outSphere)
{
  typedef Sphere<NumT>                     sphere_type;
  typedef typename sphere_type::point_type  point_type;
  const unsigned Dim = point_type::Dim;

  const point_type &ptA = sphereA.center();
  const point_type &ptB = sphereB.center();

  NumT vecAB[Dim];
  for (unsigned i = 0; i < Dim; ++i) {
    vecAB[i] = ptB[i] - ptA[i];
  }

  NumT distSq = 0;
  for (unsigned i = 0; i < Dim; ++i) {
    distSq += vecAB[i] * vecAB[i];
  }

  NumT dist = sqrt(distSq);

  NumT radiusA = sphereA.radius();
  NumT radiusB = sphereB.radius();
  NumT alpha = 0.5 * (radiusB - radiusA + dist) / dist;

  if (alpha <= 0) {
    outSphere.center() = sphereA.center();
    outSphere.radius() = sphereA.radius();
  }
  else if (alpha >= 1) {
    outSphere.center() = sphereB.center();
    outSphere.radius() = sphereB.radius();
  }
  else {
    point_type &center = outSphere.center();
    for (unsigned i = 0; i < Dim; ++i) {
      center[i] = ptA[i] + alpha * vecAB[i];
    }
    outSphere.radius() = sphereA.radius() + alpha * dist;
  }
}

///< Set output box to contain both input boxes
template <>
KOKKOS_FORCEINLINE_FUNCTION void UnionBoxes<Sphere<float> >(const Sphere<float> & inBox1, const Sphere<float> & inBox2,
                                                            Sphere<float> & outBox) {
  UnionSpheres<float>(inBox1, inBox2, outBox);
}

template <>
KOKKOS_FORCEINLINE_FUNCTION void UnionBoxes<Sphere<double> >(const Sphere<double> & inBox1, const Sphere<double> & inBox2,
                                                            Sphere<double> & outBox) {
  UnionSpheres<double>(inBox1, inBox2, outBox);
}


template<typename SrcBV_T, typename SphereNum_T>
void BVWrapBV(const SrcBV_T &srcBV, Sphere<SphereNum_T> &wrappingBV) {

  SphereNum_T x_min = srcBV.get_x_min();
  SphereNum_T y_min = srcBV.get_y_min();
  SphereNum_T z_min = srcBV.get_z_min();

  SphereNum_T x_max = srcBV.get_x_max();
  SphereNum_T y_max = srcBV.get_y_max();
  SphereNum_T z_max = srcBV.get_z_max();

  SphereNum_T ctr_x = 0.5 * (x_min + x_max);
  SphereNum_T ctr_y = 0.5 * (y_min + y_max);
  SphereNum_T ctr_z = 0.5 * (z_min + z_max);

  SphereNum_T rad_x = x_max - ctr_x;
  SphereNum_T rad_y = y_max - ctr_y;
  SphereNum_T rad_z = z_max - ctr_z;

  SphereNum_T radiusSq = (rad_x * rad_x) + (rad_y * rad_y) + (rad_z * rad_z);
  SphereNum_T radius = sqrt(radiusSq);

  wrappingBV.set_center(Point<SphereNum_T>(ctr_x, ctr_y, ctr_z));
  wrappingBV.set_radius(radius);
}

template<typename BV_T>
double computeBBoxVolume(const BV_T &bv) {
  typedef typename BV_T::value_type value_type;
  value_type x_diff = bv.get_x_max() - bv.get_x_min();
  value_type y_diff = bv.get_y_max() - bv.get_y_min();
  value_type z_diff = bv.get_z_max() - bv.get_z_min();

  return x_diff * y_diff * z_diff;
}


// Destructively hack a KDTree.
template<typename TightBV_T, typename SloppyBV_T>
bool CheatFromOtherTreesBVs(const ProximitySearchTree_T<TightBV_T> &tight_tree, ProximitySearchTree_T<SloppyBV_T> &edit_tree) {
  const std::vector<ObjectBoundingBoxHierarchy_T<TightBV_T> > &tightTreeData = tight_tree.data();
  std::vector<ObjectBoundingBoxHierarchy_T<SloppyBV_T> >      &editTreeData = edit_tree.data();

  if (tightTreeData.size() != editTreeData.size()) {
    return false;
  }

  const int numNodes = tightTreeData.size();
  for (int i = 0; i < numNodes; ++i) {
    if (tightTreeData[i].right_child_offset != editTreeData[i].right_child_offset) {
      return false;
    }
  }

  for (int i = 0; i < numNodes; ++i) {
    const int right_child_offset = tightTreeData[i].right_child_offset;
    if (right_child_offset <= 0) {
      continue;
    }
    SloppyBV_T candidate;
    BVWrapBV(tightTreeData[i].m_box, candidate);
    if (computeBBoxVolume(candidate) <  computeBBoxVolume(editTreeData[i].m_box)) {
      // std::cout << "Replace " << editTreeData[i].m_box << " with sphere that wraps " << tightTreeData[i].m_box
      //           << ", which is " << candidate << std::endl;
      editTreeData[i].m_box = candidate;
    }
  }
  return true;
}


template<typename Num_T>
void makeSphereNodalKDTree(const std::vector<Sphere<Num_T> > &input, ProximitySearchTree_T<Sphere<Num_T> > &sphereTree) {

  typedef Num_T                             num_type;
  typedef Sphere<num_type>                  sphere_type;

  int numSpheres = input.size();

  typedef ObjectBoundingBox_T<sphere_type> obj_box_type;
  std::vector<obj_box_type> input_boxes;
  input_boxes.reserve(numSpheres);
  for (int i = 0; i < numSpheres; ++i) {
    input_boxes.push_back(obj_box_type(input[i], i));
  }

  sphereTree.InitializeSearch(input_boxes);
}


template<typename Num_T>
void makeAABBTightenedSphereKDTree(const std::vector<Sphere<Num_T> > &input, ProximitySearchTree_T<Sphere<Num_T> > &sphereTree) {

  makeSphereNodalKDTree(input, sphereTree);

  typedef Num_T             num_type;
  typedef Sphere<num_type>  sphere_type;
  typedef Box<num_type>     aabb_type;

  int numSpheres = input.size();

  typedef ObjectBoundingBox_T<aabb_type> obj_box_type;
  std::vector<obj_box_type> obj_aabbs;
  obj_aabbs.reserve(numSpheres);
  for (int i = 0; i < numSpheres; ++i) {
    const sphere_type &s = input[i];
    aabb_type box(s.get_x_min(), s.get_y_min(), s.get_z_min(), s.get_x_max(), s.get_y_max(), s.get_z_max());
    obj_aabbs.push_back(obj_box_type(box, i));
  }

  ProximitySearchTree_T<aabb_type> aabbTree(obj_aabbs);
  CheatFromOtherTreesBVs(aabbTree, sphereTree);
}

}}


#endif /* STK_STK_SEARCH_STK_SEARCH_KDTREE_BOUNDINGSPHERES_HPP_ */
