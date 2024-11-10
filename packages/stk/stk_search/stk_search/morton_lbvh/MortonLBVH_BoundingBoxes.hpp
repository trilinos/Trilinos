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

#ifndef MORTONLBVH_BOUNDINGBOXES_HPP
#define MORTONLBVH_BOUNDINGBOXES_HPP

#include <stk_search/morton_lbvh/MortonLBVH_CommonTypes.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>
#include <ostream>
#include <string>
#include <algorithm>

namespace stk::search {

template <typename RealType>
struct MortonAABox
{
  RealType m_min[3];
  RealType m_max[3];
};


// Axis Aligned Bounding boxes.
template <typename RealType, typename ExecutionSpace>
struct MortonAabbList
{
  using real_type = RealType;
  using kokkos_aabb_types = MortonAabbTypes<real_type, ExecutionSpace>;

  using aabb_points_t   = typename kokkos_aabb_types::aabb_points_t;
  using aabb_points_hmt = typename kokkos_aabb_types::aabb_points_hmt;
  using device_type     = typename aabb_points_t::device_type;

  MortonAabbList(const std::string &baseName, LocalOrdinal numBoxes = 0);

  KOKKOS_INLINE_FUNCTION
  void host_set_box(LocalOrdinal boxIdx, double minX, double maxX, double minY, double maxY, double minZ, double maxZ);

  void get_bounds(double min[3], double max[3]);

  void reset(LocalOrdinal numBoxes);
  void resize(LocalOrdinal numBoxes);

  void sync_from_device();
  void sync_to_device();

  std::ostream &streamit(std::ostream &os) const;
  std::ostream &streamit(std::ostream &os, size_t boxIdx) const;

  LocalOrdinal m_numBoxes;
  std::string m_minsName;
  std::string m_maxsName;

  aabb_points_t m_mins;
  aabb_points_t m_maxs;
  aabb_points_hmt hm_mins;
  aabb_points_hmt hm_maxs;
};

template <typename RealType, typename ExecutionSpace>
MortonAabbList<RealType, ExecutionSpace>::MortonAabbList(const std::string &baseName, LocalOrdinal numBoxes)
  : m_numBoxes(numBoxes),
    m_minsName(compound_name(baseName, "mins")),
    m_maxsName(compound_name(baseName, "maxs")),
    m_mins(m_minsName, m_numBoxes),
    m_maxs(m_maxsName, m_numBoxes)
{
  std::cout << "Creating a MortonAabbList object" << std::endl;
  hm_mins = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_mins);
  hm_maxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_maxs);

  // Should assert here that the views are not null...
#ifdef BANDC_DEBUG
  std::cout << "Created it." << std::endl;
#endif
}

template <typename RealType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void MortonAabbList<RealType, ExecutionSpace>::host_set_box(LocalOrdinal boxIdx,
                                                            double minX, double maxX,
                                                            double minY, double maxY,
                                                            double minZ, double maxZ)
{
  hm_mins(boxIdx, 0) = minX;
  hm_mins(boxIdx, 1) = minY;
  hm_mins(boxIdx, 2) = minZ;
  hm_maxs(boxIdx, 0) = maxX;
  hm_maxs(boxIdx, 1) = maxY;
  hm_maxs(boxIdx, 2) = maxZ;
}

template <typename RealType, typename ExecutionSpace>
void MortonAabbList<RealType, ExecutionSpace>::get_bounds(double min[3], double max[3])
{
  const LocalOrdinal count = m_numBoxes;
  if (count <= 0) {
    min[0] = min[1] = min[2] = -1;
    max[0] = max[1] = max[2] = 1;
  }
  else {
    for (LocalOrdinal j = 0; j < 3; ++j) {
      min[j] = hm_mins(0, j);
      max[j] = hm_maxs(0, j);
    }
  }
  for (LocalOrdinal i = 1; i < count; ++i) {
    for (LocalOrdinal j = 0; j < 3; ++j) {
      min[j] = std::min(hm_mins(i, j), min[j]);
      max[j] = std::max(hm_maxs(i, j), max[j]);
    }
  }
}

template <typename RealType, typename ExecutionSpace>
void MortonAabbList<RealType, ExecutionSpace>::reset(LocalOrdinal numBoxes)
{
  if ((numBoxes >= 0) && (numBoxes != m_numBoxes)) {
    m_numBoxes = numBoxes;
    m_mins = aabb_points_t(m_minsName, m_numBoxes);
    m_maxs = aabb_points_t(m_maxsName, m_numBoxes);
    hm_mins = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_mins);
    hm_maxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_maxs);
  }
}

template <typename RealType, typename ExecutionSpace>
void MortonAabbList<RealType, ExecutionSpace>::resize(LocalOrdinal numBoxes)
{
  if ((numBoxes >= 0) && (numBoxes != m_numBoxes)) {
    m_numBoxes = numBoxes;
    Kokkos::resize(m_mins, m_numBoxes);
    Kokkos::resize(m_maxs, m_numBoxes);
    hm_mins = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_mins);
    hm_maxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_maxs);
  }
}

template <typename RealType, typename ExecutionSpace>
void MortonAabbList<RealType, ExecutionSpace>::sync_from_device()
{
  Kokkos::deep_copy(hm_mins, m_mins);
  Kokkos::deep_copy(hm_maxs, m_maxs);
}

template <typename RealType, typename ExecutionSpace>
void MortonAabbList<RealType, ExecutionSpace>::sync_to_device()
{
  Kokkos::deep_copy(m_mins, hm_mins);
  Kokkos::deep_copy(m_maxs, hm_maxs);
}

template <typename RealType, typename ExecutionSpace>
std::ostream &MortonAabbList<RealType, ExecutionSpace>::streamit(std::ostream &os) const
{
  LocalOrdinal num_boxes = hm_mins.extent(0);
  os << "{AABBs " << num_boxes << std::endl;

  for (LocalOrdinal idx = 0; idx < num_boxes; ++idx) {
    os << " {(" << hm_mins(idx, 0) << " " << hm_mins(idx, 1) << " " << hm_mins(idx, 2)
       << ") (" << hm_maxs(idx, 0) << " " << hm_maxs(idx, 1) << " " << hm_maxs(idx, 2) << ")}" << std::endl;
  }

  os << "}";
  return os;
}

template <typename RealType, typename ExecutionSpace>
std::ostream &MortonAabbList<RealType, ExecutionSpace>::streamit(std::ostream &os, size_t idx) const
{
  os << " {(" << hm_mins(idx, 0) << " " << hm_mins(idx, 1) << " " << hm_mins(idx, 2)
     << ") (" << hm_maxs(idx, 0) << " " << hm_maxs(idx, 1) << " " << hm_maxs(idx, 2) << ")}";
  return os;
}

}

#endif // MORTONLBVH_BOUNDINGBOXES_HPP
