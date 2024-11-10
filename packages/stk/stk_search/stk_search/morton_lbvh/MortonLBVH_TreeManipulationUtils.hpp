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
//
// The device kernels in this file were ported with minimal changes from
// Sam Mish's (Sandia) code which extends work by Tero Karras.  Mish's
// source file (kernels.cu) includes the following notice:
//
//    Note: these CUDA kernels were written with the input and guidance from
//    Tero Karras' paper, "Maximizing Parallelism in the Construction of BVHs,
//    Octrees and k-d Trees", as well as his explicit sources, which can be
//    found in the NVIDIA_research directory.
//
// NVIDIA CORPORATION provides those sources under the following open-source
// copyright/licence:
//
// -------------------------Begin NVIDIA Notice-----------------------------
// Copyright (c) 2013, NVIDIA CORPORATION. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of NVIDIA CORPORATION nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------End NVIDIA Notice-------------------------------
//

#ifndef MORTONLBVH_TREEMANIPULATIONUTILS_HPP
#define MORTONLBVH_TREEMANIPULATIONUTILS_HPP

#include <stk_search/morton_lbvh/MortonLBVH_CommonTypes.hpp>
#include <stk_search/morton_lbvh/MortonLBVH_BoundingBoxes.hpp>
#include <stk_search/morton_lbvh/MortonLBVH_CollisionList.hpp>
#include <stk_search/morton_lbvh/MortonLBVH_Tree.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Kokkos_Core.hpp>
#include "Kokkos_Sort.hpp"
//#if KOKKOS_VERSION < 40300
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_ROCTHRUST)
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif
//#endif
#include <iostream>
#include <ostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>

//
// Cuda and gcc disagree about whether the argument to clz*(.) is signed or not!
//
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
#define AABB_CLZ(x) (__clz(x))
#define AABB_CLZLL(x) (__clzll(x))
#define AABB_MIN(x, y) (::min(x, y))
#define AABB_MAX(x, y) (::max(x, y))
#else
#define AABB_CLZ(x) (__builtin_clz(x))
#define AABB_CLZLL(x) (__builtin_clzll(x))
#define AABB_MIN(x, y) (::std::min(x, y))
#define AABB_MAX(x, y) (::std::max(x, y))
#endif

namespace stk::search {

template <typename ViewType, typename ExecutionSpace>
struct TotalBoundsFunctor
{
  using size_type = typename ExecutionSpace::size_type;
  using BoxType = typename ViewType::value_type::box_type;
  using RealType = typename BoxType::value_type;

  using value_type = MortonAABox<RealType>;

  TotalBoundsFunctor(const MortonAabbTree<ViewType, ExecutionSpace> &tree);

  KOKKOS_INLINE_FUNCTION
  void init(value_type &update) const;

  static void apply(MortonAabbTree<ViewType, ExecutionSpace> &tree, ExecutionSpace const& execSpace);

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type idx, value_type &update) const;

  KOKKOS_INLINE_FUNCTION
  void join(value_type &update, const value_type &input) const;

  ViewType m_minMaxs;
};

template <typename ViewType, typename ExecutionSpace>
TotalBoundsFunctor<ViewType, ExecutionSpace>::TotalBoundsFunctor(const MortonAabbTree<ViewType, ExecutionSpace> &tree)
  : m_minMaxs(tree.m_minMaxs)
{}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void TotalBoundsFunctor<ViewType, ExecutionSpace>::init(value_type &update) const
{
  update.m_min[0] = FLT_MAX;
  update.m_min[1] = FLT_MAX;
  update.m_min[2] = FLT_MAX;

  update.m_max[0] = -FLT_MAX;
  update.m_max[1] = -FLT_MAX;
  update.m_max[2] = -FLT_MAX;
}

template <typename ViewType, typename ExecutionSpace>
void TotalBoundsFunctor<ViewType, ExecutionSpace>::apply(MortonAabbTree<ViewType, ExecutionSpace> &tree, ExecutionSpace const& execSpace)
{
  value_type retBox;
  retBox.m_min[0] = FLT_MAX;
  retBox.m_min[1] = FLT_MAX;
  retBox.m_min[2] = FLT_MAX;

  retBox.m_max[0] = -FLT_MAX;
  retBox.m_max[1] = -FLT_MAX;
  retBox.m_max[2] = -FLT_MAX;

  if (tree.hm_numLeaves() > 0) {
    const TotalBoundsFunctor tbf(tree);
    const size_t numLeaves = tree.hm_numLeaves();
    auto policy = Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, numLeaves);
    Kokkos::parallel_reduce(policy, tbf, retBox);
  }

  tree.m_globalMinPt[0] = retBox.m_min[0];
  tree.m_globalMinPt[1] = retBox.m_min[1];
  tree.m_globalMinPt[2] = retBox.m_min[2];

  tree.m_globalMaxPt[0] = retBox.m_max[0];
  tree.m_globalMaxPt[1] = retBox.m_max[1];
  tree.m_globalMaxPt[2] = retBox.m_max[2];
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void TotalBoundsFunctor<ViewType, ExecutionSpace>::operator()(size_type idx, value_type &update) const
{
  const BoxType& box = get_box(m_minMaxs(idx));
  update.m_min[0] = fmin(box.get_x_min(), update.m_min[0]);
  update.m_min[1] = fmin(box.get_y_min(), update.m_min[1]);
  update.m_min[2] = fmin(box.get_z_min(), update.m_min[2]);

  update.m_max[0] = fmax(box.get_x_max(), update.m_max[0]);
  update.m_max[1] = fmax(box.get_y_max(), update.m_max[1]);
  update.m_max[2] = fmax(box.get_z_max(), update.m_max[2]);
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void TotalBoundsFunctor<ViewType, ExecutionSpace>::join(value_type &update, const value_type &input) const
{
  update.m_min[0] = fmin(update.m_min[0], input.m_min[0]);
  update.m_min[1] = fmin(update.m_min[1], input.m_min[1]);
  update.m_min[2] = fmin(update.m_min[2], input.m_min[2]);

  update.m_max[0] = fmax(update.m_max[0], input.m_max[0]);
  update.m_max[1] = fmax(update.m_max[1], input.m_max[1]);
  update.m_max[2] = fmax(update.m_max[2], input.m_max[2]);
}


template <typename ViewType, typename ExecutionSpace>
struct MortonEncoder
{
  using value_type = int;
  using BoxType = typename ViewType::value_type::box_type;
  using RealType = typename BoxType::value_type;

  using LBVH_types = MortonLbvhTypes<ExecutionSpace>;

  MortonEncoder(const MortonAabbTree<ViewType, ExecutionSpace> &tree, bool reallyEncode);

  static void apply(const MortonAabbTree<ViewType, ExecutionSpace> &tree, ExecutionSpace const& execSpace, bool reallyEncode = true);

  KOKKOS_INLINE_FUNCTION
  void operator()(unsigned leafIdx) const;

  ViewType m_minMaxs;
  typename LBVH_types::local_ordinals_t m_idsOut;
  typename LBVH_types::aabb_morton_codes_t m_codesOut;
  const LocalOrdinal m_numPts;
  const RealType m_xWidth;
  const RealType m_yWidth;
  const RealType m_zWidth;
  const RealType m_globalXMin;
  const RealType m_globalYMin;
  const RealType m_globalZMin;
  const bool m_reallyDo;
};

template <typename ViewType, typename ExecutionSpace>
MortonEncoder<ViewType, ExecutionSpace>::MortonEncoder(const MortonAabbTree<ViewType, ExecutionSpace> &tree,
                                                       bool reallyEncode)
  : m_minMaxs(tree.m_minMaxs),
    m_idsOut(tree.m_leafIds),
    m_codesOut(tree.m_leafCodes),
    m_numPts(tree.hm_numLeaves()),
    m_xWidth(tree.m_globalMaxPt[0] - tree.m_globalMinPt[0]),
    m_yWidth(tree.m_globalMaxPt[1] - tree.m_globalMinPt[1]),
    m_zWidth(tree.m_globalMaxPt[2] - tree.m_globalMinPt[2]),
    m_globalXMin(tree.m_globalMinPt[0]),
    m_globalYMin(tree.m_globalMinPt[1]),
    m_globalZMin(tree.m_globalMinPt[2]),
    m_reallyDo(reallyEncode)
{
}

template <typename ViewType, typename ExecutionSpace>
void MortonEncoder<ViewType, ExecutionSpace>::apply(const MortonAabbTree<ViewType, ExecutionSpace> &tree,
                                                    ExecutionSpace const& execSpace, bool reallyEncode)
{
  const MortonEncoder op(tree, reallyEncode);
  const size_t numLeaves = tree.hm_numLeaves();
  auto policy = Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, numLeaves);
  Kokkos::parallel_for(policy, op);
}

#ifdef SMALL_MORTON  // 32 bit Morton code

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void MortonEncoder<ExecutionSpace>::operator()(unsigned leafIdx) const
{
  const BoxType& box = get_box(m_minMaxs(leafIdx));
  const RealType ctdX = 0.5 * (box.get_x_min() + box.get_x_max());
  const RealType ctdY = 0.5 * (box.get_y_min() + box.get_y_max());
  const RealType ctdZ = 0.5 * (box.get_z_min() + box.get_z_max());

  // std::cout << "box(" << leafIdx << ") = (" << m_minMax(leafIdx, 0) << " "
  //           <<  m_minMax(leafIdx, 1) << " " <<  m_minMax(leafIdx, 2)
  //           << ") (" << m_minMax(leafIdx, 3) << " " <<  m_minMax(leafIdx, 4)
  //           << " " <<  m_minMax(leafIdx, 5) << ")" <<  std::endl;
  // std::cout << "centroid(" << leafIdx << ") = (" << ctdX << " " << ctdY << " " << ctdZ
  //           << ")" << std::endl;

  //  for a 32-bit morton code

  morton_code_t ux, uy, uz;

  // for a 32-bit morton code, each spatial dimension gets 10 bits.
  // To get the most out of each bit, we normalize the coordinates
  // against 1023.0f (i.e. 2^10 - 1)

  ux = static_cast<morton_code_t>((ctdX - m_globalXMin) * 1023.0f / m_xWidth);
  uy = static_cast<morton_code_t>((ctdY - m_globalYMin) * 1023.0f / m_yWidth);
  uz = static_cast<morton_code_t>((ctdZ - m_globalZMin) * 1023.0f / m_zWidth);

  ux = (ux * 0x00010001u) & 0xFF0000FFu;
  ux = (ux * 0x00000101u) & 0x0F00F00Fu;
  ux = (ux * 0x00000011u) & 0xC30C30C3u;
  ux = (ux * 0x00000005u) & 0x49249249u;

  uy = (uy * 0x00010001u) & 0xFF0000FFu;
  uy = (uy * 0x00000101u) & 0x0F00F00Fu;
  uy = (uy * 0x00000011u) & 0xC30C30C3u;
  uy = (uy * 0x00000005u) & 0x49249249u;

  uz = (uz * 0x00010001u) & 0xFF0000FFu;
  uz = (uz * 0x00000101u) & 0x0F00F00Fu;
  uz = (uz * 0x00000011u) & 0xC30C30C3u;
  uz = (uz * 0x00000005u) & 0x49249249u;

  m_idsOut(leafIdx) = leafIdx;
  m_codesOut(leafIdx) = ux * 4 + uy * 2 + uz;
}

#else  // 64 bit Morton codes

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void MortonEncoder<ViewType, ExecutionSpace>::operator()(unsigned leafIdx) const
{
  m_idsOut(leafIdx) = leafIdx;

  if (m_reallyDo) {
    const BoxType& box = get_box(m_minMaxs(leafIdx));
    const RealType ctdX = 0.5 * (box.get_x_min() + box.get_x_max());
    const RealType ctdY = 0.5 * (box.get_y_min() + box.get_y_max());
    const RealType ctdZ = 0.5 * (box.get_z_min() + box.get_z_max());

    // std::cout << "box(" << leafIdx << ") = (" << m_minMaxs(leafIdx, 0) << " "
    //           <<  m_minMaxs(leafIdx, 1) << " " <<  m_minMaxs(leafIdx, 2)
    //           << ") (" << m_minMaxs(leafIdx, 3) << " " <<  m_minMaxs(leafIdx, 4)
    //           << " " <<  m_minMaxs(leafIdx, 5) << ")" <<  std::endl;
    // std::cout << "centroid(" << leafIdx << ") = (" << ctdX << " " << ctdY << " " << ctdZ
    //           << ")" << std::endl;

    //  for a 64-bit morton code

    morton_code_t ux, uy, uz;

    // for a 64-bit morton code, each spatial dimension gets 21 bits.
    // To get the most out of each bit, we normalize the coordinates
    // against 2097151.0f (i.e. 2^21 - 1)

    ux = static_cast<morton_code_t>((ctdX - m_globalXMin) * 2097151.0f / m_xWidth);
    uy = static_cast<morton_code_t>((ctdY - m_globalYMin) * 2097151.0f / m_yWidth);
    uz = static_cast<morton_code_t>((ctdZ - m_globalZMin) * 2097151.0f / m_zWidth);

    ux = (ux | ux << 32) & 0x001f00000000ffff;
    ux = (ux | ux << 16) & 0x001f0000ff0000ff;
    ux = (ux | ux << 8) & 0x100f00f00f00f00f;
    ux = (ux | ux << 4) & 0x10c30c30c30c30c3;
    ux = (ux | ux << 2) & 0x1249249249249249;

    uy = (uy | uy << 32) & 0x001f00000000ffff;
    uy = (uy | uy << 16) & 0x001f0000ff0000ff;
    uy = (uy | uy << 8) & 0x100f00f00f00f00f;
    uy = (uy | uy << 4) & 0x10c30c30c30c30c3;
    uy = (uy | uy << 2) & 0x1249249249249249;

    uz = (uz | uz << 32) & 0x001f00000000ffff;
    uz = (uz | uz << 16) & 0x001f0000ff0000ff;
    uz = (uz | uz << 8) & 0x100f00f00f00f00f;
    uz = (uz | uz << 4) & 0x10c30c30c30c30c3;
    uz = (uz | uz << 2) & 0x1249249249249249;

    m_codesOut(leafIdx) = ux * 4 + uy * 2 + uz;
  }
}

#endif  // 64 bit Morton code


// Serial sort is the default.
template <typename TreeType, typename ExecutionSpace>
struct SortByCodeIdPair
{
  using LBVH_types = MortonLbvhTypes<ExecutionSpace>;

  SortByCodeIdPair(const TreeType &tree);

  static void apply(const TreeType &tree, bool reallyEncode = true);

  std::vector<morton_code_id_pair> m_buffer;
  typename LBVH_types::local_ordinals_hmt hm_leafIds;
  typename LBVH_types::aabb_morton_codes_hmt hm_leafCodes;
};

template <typename TreeType, typename ExecutionSpace>
SortByCodeIdPair<TreeType, ExecutionSpace>::SortByCodeIdPair(const TreeType &tree)
{
  hm_leafIds = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, tree.m_leafIds);
  hm_leafCodes = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, tree.m_leafCodes);

  Kokkos::deep_copy(hm_leafCodes, tree.m_leafCodes);
  Kokkos::deep_copy(hm_leafIds, tree.m_leafIds);

  LocalOrdinal numLeaves = tree.hm_numLeaves();
  m_buffer.resize(numLeaves);

  for (LocalOrdinal idx = 0; idx < numLeaves; ++idx) {
    m_buffer[idx] = {hm_leafCodes(idx), hm_leafIds(idx)};
  }
}

template <typename TreeType, typename ExecutionSpace>
void SortByCodeIdPair<TreeType, ExecutionSpace>::apply(const TreeType &tree,
                                                       bool reallyEncode)
{
  SortByCodeIdPair tmp(tree);
  std::sort(tmp.m_buffer.begin(), tmp.m_buffer.end());
  LocalOrdinal numLeaves = tree.hm_numLeaves();

  for (LocalOrdinal idx = 0; idx < numLeaves; ++idx) {
    tmp.hm_leafCodes(idx) = tmp.m_buffer[idx].m_code;
    tmp.hm_leafIds(idx) = tmp.m_buffer[idx].m_id;
  }

  Kokkos::deep_copy(tree.m_leafCodes, tmp.hm_leafCodes);
  Kokkos::deep_copy(tree.m_leafIds, tmp.hm_leafIds);
}


template <typename TreeType, typename ExecutionSpace>
struct SortByCode
{
  static void apply(const TreeType &tree, ExecutionSpace const& execSpace)
  {
    if constexpr (Kokkos::SpaceAccessibility<ExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space>::accessible) {
      SortByCodeIdPair<TreeType, ExecutionSpace>::apply(tree);
    }
    else {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_ROCTHRUST)
      const int n = tree.m_leafIds.extent(0);

      morton_code_t *rawLeafCodes = tree.m_leafCodes.data();
      thrust::device_ptr<morton_code_t> rawLeafCodesThr = thrust::device_pointer_cast(rawLeafCodes);
      LocalOrdinal *rawLeafIds = tree.m_leafIds.data();
      thrust::device_ptr<LocalOrdinal> rawLeafIdsThr = thrust::device_pointer_cast(rawLeafIds);
      //thrust::stable_sort_by_key(rawLeafCodesThr, rawLeafCodesThr + n, rawLeafIdsThr);
      thrust::sort_by_key(rawLeafCodesThr, rawLeafCodesThr + n, rawLeafIdsThr);
#else
#if KOKKOS_VERSION >= 40300
      Kokkos::Experimental::sort_by_key(execSpace, tree.m_leafCodes, tree.m_leafIds);
#else
      STK_ThrowErrorMessage("Need at least Kokkos 4.3");
#endif
#endif
    }
  }
};

template <typename RealType, typename TreeType, typename ExecutionSpace>
struct BuildRadixTree
{
  using LBVH_types = MortonLbvhTypes<ExecutionSpace>;

  BuildRadixTree(const TreeType &tree);

  static void apply(const TreeType &tree, ExecutionSpace const& execSpace);

  KOKKOS_INLINE_FUNCTION
  void operator()(unsigned argIdx) const;

  KOKKOS_INLINE_FUNCTION
  int leaves_cpr(LocalOrdinal baseIdx, LocalOrdinal testIdx) const;

  const LocalOrdinal m_numLeaves;
  const LocalOrdinal m_numInternalNodes;
  typename LBVH_types::aabb_morton_codes_tmt tm_leafCodes;
  typename LBVH_types::local_ordinals_tmt tm_leafIds;
  typename LBVH_types::local_ordinal_pairs_t m_nodeChildren;
  typename LBVH_types::local_ordinals_t m_nodeParents;
  typename LBVH_types::local_ordinals_t m_atomicFlags;
};

template <typename RealType, typename TreeType, typename ExecutionSpace>
BuildRadixTree<RealType, TreeType, ExecutionSpace>::BuildRadixTree(const TreeType &tree)
  : m_numLeaves(tree.hm_numLeaves()),
    m_numInternalNodes(tree.hm_numInternalNodes()),
    tm_leafCodes(tree.m_leafCodes),
    tm_leafIds(tree.m_leafIds),
    m_nodeChildren(tree.m_nodeChildren),
    m_nodeParents(tree.m_nodeParents),
    m_atomicFlags(tree.m_atomicFlags)
{}

template <typename RealType, typename TreeType, typename ExecutionSpace>
void BuildRadixTree<RealType, TreeType, ExecutionSpace>::apply(const TreeType &tree, ExecutionSpace const& execSpace)
{
  if (tree.hm_numLeaves() <= 0) {
    return;
  }
  BuildRadixTree<RealType,TreeType,ExecutionSpace> op(tree);
  auto policy = Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, static_cast<unsigned>(tree.hm_numInternalNodes()));
  Kokkos::parallel_for(policy, op);
}

template <typename RealType, typename TreeType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void BuildRadixTree<RealType, TreeType, ExecutionSpace>::operator()(unsigned argIdx) const
{
  LocalOrdinal idx = static_cast<LocalOrdinal>(argIdx);

  if (idx >= m_numInternalNodes) return;

  const uint32_t numLeavesU32 = static_cast<uint32_t>(m_numLeaves);

  // Choose direction
  int prefixPrev = leaves_cpr(idx, idx - 1);
  int prefixNext = leaves_cpr(idx, idx + 1);
  int d = (prefixNext > prefixPrev) ? 1 : -1;
  int prefixMin = AABB_MIN(prefixPrev, prefixNext);

  // Find upper bound for length.
  int lmax = 128 >> 2;  // Shifted back upon entrance in loop.
  uint32_t probe = 0;
  do {
    lmax <<= static_cast<int>(2);
    probe = idx + lmax * d;
  } while ((probe < numLeavesU32) && (leaves_cpr(idx, probe) > prefixMin));

  // Determine length.
  int l = 0;
  for (int t = lmax >> 1; t > 0; t >>= 1) {
    probe = idx + (l + t) * d;
    if ((probe < numLeavesU32) && (leaves_cpr(idx, probe) > prefixMin)) l += t;
  }
  int j = idx + l * d;
  int prefixNode = leaves_cpr(idx, j);

  // Find split point.
  int s = 0;
  int t = l;
  do {
    t = (t + 1) >> 1;
    probe = idx + (s + t) * d;
    if ((probe < numLeavesU32) && (leaves_cpr(idx, probe) > prefixNode)) {
      s += t;
    }
  } while (t > 1);
  int k = idx + s * d + AABB_MIN(d, 0);

  // Output node.
  int lo = AABB_MIN(idx, j);
  int hi = AABB_MAX(idx, j);
  int cx, cy;
  cx = (lo == k) ? k + 0 : (k + 0 + m_numLeaves);
  cy = (hi == k + 1) ? k + 1 : (k + 1 + m_numLeaves);

  // If the children are leaves, their boxes didn't get moved in the sort!
  // In that case, we need to follow the m_leaf_ids map instead of using
  // (k + 0) and (k + 1) directly.
  if (cx < m_numLeaves) {
    cx = tm_leafIds(cx);
  }
  if (cy < m_numLeaves) {
    cy = tm_leafIds(cy);
  }

  m_nodeChildren(idx + m_numLeaves, 0) = cx;
  m_nodeChildren(idx + m_numLeaves, 1) = cy;
  m_nodeParents(cx) = idx + m_numLeaves;
  m_nodeParents(cy) = idx + m_numLeaves;
  m_atomicFlags(idx) = 0;

  if (idx == 0) {
    m_nodeParents(m_numLeaves) = m_numLeaves;
  }
}

#ifdef SMALL_MORTON  // 32 bit Morton

template <typename RealType, typename TreeType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
int BuildRadixTree<RealType, TreeType, ExecutionSpace>::leaves_cpr(LocalOrdinal baseIdx, LocalOrdinal testIdx) const
{
  if (testIdx < 0 || testIdx >= m_numLeaves) {
    return -1;
  }
  typename LBVH_types::aabb_morton_codes_tmt::value_type xorResult =
      tm_leafCodes(baseIdx) ^ tm_leafCodes(testIdx);
  if (xorResult != 0) {
    return AABB_CLZ(xorResult);
  }
  else {
    return 32 + AABB_CLZ(tm_leafIds(baseIdx) ^ tm_leafIds(testIdx));
  }
}

#else  // 64 bit Morton

template <typename RealType, typename TreeType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
int BuildRadixTree<RealType, TreeType, ExecutionSpace>::leaves_cpr(LocalOrdinal baseIdx, LocalOrdinal testIdx) const
{
  if (testIdx < 0 || testIdx >= m_numLeaves) {
    return -1;
  }
  typename LBVH_types::aabb_morton_codes_t::value_type xorResult = tm_leafCodes(baseIdx) ^ tm_leafCodes(testIdx);
  if (xorResult != 0) {
    return AABB_CLZLL(xorResult);
  }
  else {
    return 64 + AABB_CLZ(tm_leafIds(baseIdx) ^ tm_leafIds(testIdx));
  }
}

#endif  // 64 bit Morton


template <typename ViewType, typename ExecutionSpace>
struct UpdateInteriorNodeBVs
{
  using LBVH_types = MortonLbvhTypes<ExecutionSpace>;
  using BoxType = typename ViewType::value_type::box_type;
  using RealType = typename BoxType::value_type;
  using kokkos_aabb_types = MortonAabbTypes<RealType, ExecutionSpace>;
  using bboxes_3d_view_amt = typename kokkos_aabb_types::bboxes_3d_view_amt;

  UpdateInteriorNodeBVs(const MortonAabbTree<ViewType, ExecutionSpace> &tree);

  static void apply(const MortonAabbTree<ViewType, ExecutionSpace> &tree, ExecutionSpace const& execSpace);

  KOKKOS_INLINE_FUNCTION
  void operator()(unsigned argIdx) const;

  KOKKOS_FORCEINLINE_FUNCTION
  void get_box(RealType bvMinMax[6], LocalOrdinal idx, const bboxes_3d_view_amt &boxesMinMax) const;

  KOKKOS_FORCEINLINE_FUNCTION
  void get_stk_box(RealType bvMinMax[6], LocalOrdinal idx, const ViewType &boxesMinMax) const;

  const LocalOrdinal m_numLeaves;
  const LocalOrdinal m_numInternalNodes;
  typename LBVH_types::local_ordinal_pairs_tmt tm_nodeChildren;
  typename LBVH_types::local_ordinals_tmt tm_nodeParents;
  ViewType m_leafMinMaxs;

  // Will write to internal nodes bounding boxes.
  bboxes_3d_view_amt m_nodeMinMaxs;
  typename LBVH_types::local_ordinals_t m_atomicFlags;
};

template <typename ViewType, typename ExecutionSpace>
UpdateInteriorNodeBVs<ViewType, ExecutionSpace>::UpdateInteriorNodeBVs(const MortonAabbTree<ViewType, ExecutionSpace> &tree)
  : m_numLeaves(tree.hm_numLeaves()),
    m_numInternalNodes(tree.hm_numInternalNodes()),
    tm_nodeChildren(tree.m_nodeChildren),
    tm_nodeParents(tree.m_nodeParents),
    m_leafMinMaxs(tree.m_minMaxs),
    m_nodeMinMaxs(tree.m_nodeMinMaxs),
    m_atomicFlags(tree.m_atomicFlags)
{}

template <typename ViewType, typename ExecutionSpace>
void UpdateInteriorNodeBVs<ViewType, ExecutionSpace>::apply(const MortonAabbTree<ViewType, ExecutionSpace> &tree, ExecutionSpace const& execSpace)
{
  const UpdateInteriorNodeBVs<ViewType, ExecutionSpace> op(tree);
  const size_t numLeaves = tree.hm_numLeaves();

  auto policy = Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, numLeaves);
  Kokkos::parallel_for(policy, op);
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void UpdateInteriorNodeBVs<ViewType, ExecutionSpace>::operator()(unsigned argIdx) const
{
  if (m_numLeaves > 1) {
    LocalOrdinal idx = static_cast<LocalOrdinal>(argIdx);

    RealType bvMinMax[6];
    get_stk_box(bvMinMax, idx, m_leafMinMaxs);

    LocalOrdinal parent = tm_nodeParents(idx);
    RealType sibMinMax[6];

    while (Kokkos::atomic_fetch_add(&m_atomicFlags(parent - m_numLeaves), 1) == 1) {
      LocalOrdinal sib = tm_nodeChildren(parent, 0);
      if (sib == idx) {
        sib = tm_nodeChildren(parent, 1);
      }

      if (sib < m_numLeaves) {
        get_stk_box(sibMinMax, sib, m_leafMinMaxs);
      }
      else {
        get_box(sibMinMax, sib-m_numLeaves, m_nodeMinMaxs);
      }

      const LocalOrdinal thisIdx = parent - m_numLeaves;

      bvMinMax[0] = AABB_MIN(bvMinMax[0], sibMinMax[0]);
      bvMinMax[1] = AABB_MIN(bvMinMax[1], sibMinMax[1]);
      bvMinMax[2] = AABB_MIN(bvMinMax[2], sibMinMax[2]);
      bvMinMax[3] = AABB_MAX(bvMinMax[3], sibMinMax[3]);
      bvMinMax[4] = AABB_MAX(bvMinMax[4], sibMinMax[4]);
      bvMinMax[5] = AABB_MAX(bvMinMax[5], sibMinMax[5]);
      m_nodeMinMaxs(thisIdx, 0) = bvMinMax[0];
      m_nodeMinMaxs(thisIdx, 1) = bvMinMax[1];
      m_nodeMinMaxs(thisIdx, 2) = bvMinMax[2];
      m_nodeMinMaxs(thisIdx, 3) = bvMinMax[3];
      m_nodeMinMaxs(thisIdx, 4) = bvMinMax[4];
      m_nodeMinMaxs(thisIdx, 5) = bvMinMax[5];

      idx = parent;
      parent = tm_nodeParents(parent);
      if (idx == parent) {
        return;
      }
    }
  }
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_FORCEINLINE_FUNCTION
void UpdateInteriorNodeBVs<ViewType, ExecutionSpace>::get_box(RealType bvMinMax[6], LocalOrdinal idx,
                                                              const bboxes_3d_view_amt &boxMinMaxs) const
{
  bvMinMax[0] = boxMinMaxs(idx, 0);
  bvMinMax[1] = boxMinMaxs(idx, 1);
  bvMinMax[2] = boxMinMaxs(idx, 2);
  bvMinMax[3] = boxMinMaxs(idx, 3);
  bvMinMax[4] = boxMinMaxs(idx, 4);
  bvMinMax[5] = boxMinMaxs(idx, 5);
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_FORCEINLINE_FUNCTION
void UpdateInteriorNodeBVs<ViewType, ExecutionSpace>::get_stk_box(RealType bvMinMax[6], LocalOrdinal idx,
                                                              const ViewType &boxMinMaxs) const
{
  const BoxType& box = stk::search::get_box(boxMinMaxs(idx));
  bvMinMax[0] = box.get_x_min();
  bvMinMax[1] = box.get_y_min();
  bvMinMax[2] = box.get_z_min();
  bvMinMax[3] = box.get_x_max();
  bvMinMax[4] = box.get_y_max();
  bvMinMax[5] = box.get_z_max();
}


template <typename ExecutionSpace>
class CollisionListCallback
{
  public:
    CollisionListCallback(CollisionList<ExecutionSpace>& collisionList) :
      m_collisionList(collisionList)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(int domainIdx, int rangeIdx) const
    {
      m_collisionList.push_back(domainIdx, rangeIdx);
    }

    bool resize_for_second_pass()
    {
      int numActualCollisions = m_collisionList.get_num_collisions();
      bool needSecondPass = numActualCollisions > m_collisionList.get_capacity();
      if (needSecondPass)
      {
        m_collisionList.reset(numActualCollisions);
      }

      return needSecondPass;
    }

    CollisionList<ExecutionSpace> get_collision_list() const { return m_collisionList; }

  private:
    CollisionList<ExecutionSpace> m_collisionList;
};


template <typename DomainViewType, typename RangeViewType, typename ExecutionSpace, typename Callback>
struct Traverse_MASTB_BVH_Functor
{
  using DomainBoxType = typename DomainViewType::value_type::box_type;
  using RangeBoxType = typename RangeViewType::value_type::box_type;
  using RealType = typename RangeBoxType::value_type;
  using LBVH_types = MortonLbvhTypes<ExecutionSpace>;
  using kokkos_aabb_types      = MortonAabbTypes<RealType, ExecutionSpace>;
  using local_ordinals_tmt     = typename LBVH_types::local_ordinals_tmt;
  using bboxes_3d_view_t       = typename kokkos_aabb_types::bboxes_3d_view_t;
  using bboxes_const_3d_view_t = typename kokkos_aabb_types::bboxes_const_3d_view_t;
  using collision_list_type    = CollisionList<ExecutionSpace>;

  Traverse_MASTB_BVH_Functor(DomainViewType domainMinMaxs,
                             local_ordinals_tmt domainIds,
                             const MortonAabbTree<RangeViewType, ExecutionSpace> &rangeTree,
                             Callback& callback,
                             bool flippedResults = false);

  KOKKOS_INLINE_FUNCTION
  void operator()(unsigned domainIdx) const;

  KOKKOS_FORCEINLINE_FUNCTION
  bool overlaps_range(RealType bvMinMax[6], LocalOrdinal rangeIdx) const;

  KOKKOS_FORCEINLINE_FUNCTION
  bool is_range_leaf(LocalOrdinal rangeIdx) const{ return (rangeIdx < m_rangeRoot); }

  KOKKOS_FORCEINLINE_FUNCTION
  void get_box(RealType bvMinMax[6], LocalOrdinal idx, const DomainViewType &boxMinMaxs) const;

  std::ostream &stream_pair(LocalOrdinal domainIdx, bool overlap, LocalOrdinal rangeIdx, std::ostream &os) const;

  KOKKOS_INLINE_FUNCTION
  void record_result(LocalOrdinal domainIdx, LocalOrdinal rangeIdx, bool flip) const
  {
    LocalOrdinal domainIdxFlipped = flip ? rangeIdx : domainIdx;
    LocalOrdinal rangeIdxFlipped  = flip ? domainIdx : rangeIdx;
    m_callback(domainIdxFlipped, rangeIdxFlipped);
  }

  DomainViewType m_domainMinMaxs;
  typename LBVH_types::local_ordinals_tmt tm_domainIds;

  const LocalOrdinal m_rangeRoot;
  RangeViewType tm_rangeMinMaxs;
  bboxes_const_3d_view_t tm_rangeNodeMinMaxs;
  typename LBVH_types::local_ordinal_pairs_tmt tm_rangeNodeChildren;
  typename LBVH_types::local_ordinals_tmt tm_rangeLeafIds;

  const bool m_flippedResults;
  Callback m_callback;
};

template <typename DomainViewType, typename RangeViewType, typename ExecutionSpace, typename Callback>
Traverse_MASTB_BVH_Functor<DomainViewType, RangeViewType, ExecutionSpace, Callback>::Traverse_MASTB_BVH_Functor(
    DomainViewType domainMinMaxs,
    local_ordinals_tmt domainIds,
    const MortonAabbTree<RangeViewType, ExecutionSpace> &rangeTree,
    Callback& callback,
    bool flippedResults)
  : m_domainMinMaxs(domainMinMaxs),
    tm_domainIds(domainIds),
    m_rangeRoot(rangeTree.hm_numLeaves()),
    tm_rangeMinMaxs(rangeTree.m_minMaxs),
    tm_rangeNodeMinMaxs(rangeTree.m_nodeMinMaxs),
    tm_rangeNodeChildren(rangeTree.m_nodeChildren),
    tm_rangeLeafIds(rangeTree.m_leafIds),
    m_flippedResults(flippedResults),
    m_callback(callback)
{}


template <typename Tree1Type, typename Tree2Type, typename ExecutionSpace, typename Callback>
void search_tree(
    const Tree1Type &domainTree,
    const Tree2Type &rangeTree,
    Callback& callback,
    ExecutionSpace const& execSpace,
    bool flipOutputPairs = false)
{
  Kokkos::Profiling::pushRegion("search_tree");
  if ((domainTree.hm_numLeaves() == 0) || (rangeTree.hm_numLeaves() == 0)) {
    callback.resize_for_second_pass();
    return;
  }

  using DomainViewType = typename Tree1Type::view_type;
  using RangeViewType = typename Tree2Type::view_type;
  Kokkos::Profiling::pushRegion("construct Traverse_MASTB_BVH_Functor");
  const Traverse_MASTB_BVH_Functor<DomainViewType,RangeViewType,ExecutionSpace,Callback> op(domainTree.m_minMaxs, domainTree.m_leafIds, rangeTree,
                                      callback, flipOutputPairs);
  Kokkos::Profiling::popRegion();
  auto policy = Kokkos::RangePolicy<ExecutionSpace>(execSpace, 0, domainTree.hm_numLeaves());
  Kokkos::parallel_for("pll-for Traverse_MASTB_BVH_Functor", policy, op);
  Kokkos::Profiling::pushRegion("search_tree - fence after pll-for");
  execSpace.fence();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("search_tree - resize for second pass");
  if (callback.resize_for_second_pass()) {
    const Traverse_MASTB_BVH_Functor<DomainViewType,RangeViewType,ExecutionSpace,Callback> op2(domainTree.m_minMaxs, domainTree.m_leafIds, rangeTree,
                                         callback, flipOutputPairs);
    Kokkos::parallel_for("Traverse_MASTB_BVH_Functor - pass2", policy, op2);
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("search_tree - last fence");
  execSpace.fence();
  Kokkos::Profiling::popRegion();
  Kokkos::Profiling::popRegion();
}

template <typename DomainViewType, typename RangeViewType, typename ExecutionSpace, typename Callback>
KOKKOS_INLINE_FUNCTION void Traverse_MASTB_BVH_Functor<DomainViewType, RangeViewType, ExecutionSpace, Callback>::operator()(unsigned argDomainIdx) const
{
  LocalOrdinal domainIdx = tm_domainIds(argDomainIdx);

  RealType bvMinMax[6];
  get_box(bvMinMax, domainIdx, m_domainMinMaxs);

  if (m_rangeRoot > 1) {
    int ridxStack[64];
    int* stackPtr = ridxStack;
    *stackPtr++ = -1;

    int nodeIdx = m_rangeRoot;
    do {
      // Check each child node for overlap.
      const int childL = tm_rangeNodeChildren(nodeIdx, 0);
      const int childR = tm_rangeNodeChildren(nodeIdx, 1);
      const bool overlapL = overlaps_range(bvMinMax, childL);
      const bool overlapR = overlaps_range(bvMinMax, childR);

      bool traverseL = false;

      // Query overlaps a leaf node => report collision.
      if (overlapL) {
        if (is_range_leaf(childL)) {
           record_result(domainIdx, childL, m_flippedResults);
        }
        else {
          traverseL = true;
          nodeIdx = childL;
        }
      }

      // Query overlaps and internal node => traverse.
      if (overlapR) {
        if (is_range_leaf(childR)) {
          record_result(domainIdx, childR, m_flippedResults);
          if (!traverseL) {
            nodeIdx = *--stackPtr; // pop
          }
        }
        else {
          if (traverseL) {
            *stackPtr++ = childR;  // push
          }
          else {
            nodeIdx = childR;
          }
        }
      }
      else if (!traverseL) {
        nodeIdx = *--stackPtr; // pop
      }
    } while (nodeIdx >= 0);
  }
  else {
    // Degenerate case of only one leaf node
    bool overlap = overlaps_range(bvMinMax, 0);
    if (overlap) {
      record_result(domainIdx, 0, m_flippedResults);
    }
  }
}

template <typename DomainViewType, typename RangeViewType, typename ExecutionSpace, typename Callback>
KOKKOS_FORCEINLINE_FUNCTION
bool Traverse_MASTB_BVH_Functor<DomainViewType, RangeViewType, ExecutionSpace, Callback>::overlaps_range(RealType bvMinMax[6],
                                                                                                         LocalOrdinal rangeIdx) const
{
  if(rangeIdx < m_rangeRoot) {
    const RangeBoxType& box = stk::search::get_box(tm_rangeMinMaxs(rangeIdx));

    return (bvMinMax[3] < box.get_x_min() ||
            bvMinMax[4] < box.get_y_min() ||
            bvMinMax[5] < box.get_z_min() ||
            bvMinMax[0] > box.get_x_max() ||
            bvMinMax[1] > box.get_y_max() ||
            bvMinMax[2] > box.get_z_max()) ? false : true;
  }

  const LocalOrdinal rangeNodeIdx = rangeIdx - m_rangeRoot;
  return (bvMinMax[3] < tm_rangeNodeMinMaxs(rangeNodeIdx,0) ||
          bvMinMax[4] < tm_rangeNodeMinMaxs(rangeNodeIdx,1) ||
          bvMinMax[5] < tm_rangeNodeMinMaxs(rangeNodeIdx,2) ||
          bvMinMax[0] > tm_rangeNodeMinMaxs(rangeNodeIdx,3) ||
          bvMinMax[1] > tm_rangeNodeMinMaxs(rangeNodeIdx,4) ||
          bvMinMax[2] > tm_rangeNodeMinMaxs(rangeNodeIdx,5)) ? false : true;
}

template <typename DomainViewType, typename RangeViewType, typename ExecutionSpace, typename Callback>
KOKKOS_FORCEINLINE_FUNCTION
void Traverse_MASTB_BVH_Functor<DomainViewType, RangeViewType, ExecutionSpace, Callback>::get_box(RealType bvMinMax[6], LocalOrdinal idx,
                                                                       const DomainViewType &boxMinMaxs) const
{
  const DomainBoxType& box = stk::search::get_box(boxMinMaxs(idx));
  bvMinMax[0] = box.get_x_min();
  bvMinMax[1] = box.get_y_min();
  bvMinMax[2] = box.get_z_min();
  bvMinMax[3] = box.get_x_max();
  bvMinMax[4] = box.get_y_max();
  bvMinMax[5] = box.get_z_max();
}

template <typename DomainViewType, typename RangeViewType, typename ExecutionSpace, typename Callback>
std::ostream &Traverse_MASTB_BVH_Functor<DomainViewType, RangeViewType, ExecutionSpace, Callback>::stream_pair(LocalOrdinal domainIdx, bool overlap,
                                                                                LocalOrdinal rangeIdx, std::ostream &os) const
{
  auto dbox = get_box(m_domainMinMaxs(domainIdx));
  os << " {(" << dbox.get_x_min() << "," << dbox.get_y_min() << "," <<dbox.get_z_min() 
     << ") (" << dbox.get_x_max() << "," << dbox.get_y_max() << "," <<dbox.get_z_max() 
     << ")}";
  os << (overlap ? " overlaps " : " does not overlap ");
  auto rbox = rangeIdx < m_rangeRoot ? get_box(tm_rangeMinMaxs(rangeIdx)) : get_box(tm_rangeNodeMinMaxs(rangeIdx - m_rangeRoot));
  os << " {(" << rbox.get_x_min() << "," << rbox.get_y_min() << "," <<rbox.get_z_min() 
     << ") (" << rbox.get_x_max() << "," << rbox.get_y_max() << "," <<rbox.get_z_max() 
     << ")}";
  os << std::endl;
  return os;
}

}

#endif // MORTONLBVH_TREEMANIPULATIONUTILS_HPP
