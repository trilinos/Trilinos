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

#ifndef STK_SEARCH_MORTON_LBVH_MORTONLBVH_TREE_HPP
#define STK_SEARCH_MORTON_LBVH_MORTONLBVH_TREE_HPP

#include <stk_search/morton_lbvh/MortonLBVH_CommonTypes.hpp>
#include <stk_search/Box.hpp>
#include <stk_search/BoxIdent.hpp>
#include <string>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <algorithm>

//#define DEBUG_MORTON_ACCELERATED_SEARCH 1

namespace stk::search {

template<typename BoxType>
KOKKOS_INLINE_FUNCTION
const BoxType& get_box(BoxType& box) { return box; }

template<typename BoxType, typename IdentType>
KOKKOS_INLINE_FUNCTION
const BoxType& get_box(BoxIdent<BoxType,IdentType>& boxIdent) { return boxIdent.box; }

template<typename BoxType, typename IdentProcType>
KOKKOS_INLINE_FUNCTION
const BoxType& get_box(BoxIdentProc<BoxType,IdentProcType>& boxIdentProc) { return boxIdentProc.box; }


template<typename> struct BoxIdentViewTrait {};

template <typename ViewValueType, typename... otherTemplateArgs>
struct BoxIdentViewTrait<Kokkos::View<ViewValueType*, otherTemplateArgs...>>
{
  using InputBoxType = typename ViewValueType::box_type;
  using ValueType = typename InputBoxType::value_type;
  using IdentType = typename ViewValueType::ident_type;
  using ViewType = Kokkos::View<BoxIdent<Box<ValueType>,IdentType>*, otherTemplateArgs...>;
};

template<typename> struct BoxIdentProcViewTrait {};

template <typename ViewValueType, typename... otherTemplateArgs>
struct BoxIdentProcViewTrait<Kokkos::View<ViewValueType*, otherTemplateArgs...>>
{
  using InputBoxType = typename ViewValueType::box_type;
  using ValueType = typename InputBoxType::value_type;
  using IdentProcType = typename ViewValueType::ident_proc_type;
  using ViewType = Kokkos::View<BoxIdentProc<Box<ValueType>,IdentProcType>*, otherTemplateArgs...>;
};


template <typename ViewType, typename ExecutionSpace>
struct MortonAabbTree
{
  using execution_space = ExecutionSpace;
  using view_type = ViewType;
  using BoxViewType     = ViewType;
  using BoxViewType_hmt = typename ViewType::HostMirror;
  using BoxType         = typename BoxViewType::value_type::box_type;
  using real_type = typename BoxType::value_type;
  using LBVH_types = MortonLbvhTypes<ExecutionSpace>;
  using kokkos_aabb_types = MortonAabbTypes<real_type, ExecutionSpace>;

  using local_ordinal_scl_t   = typename LBVH_types::local_ordinal_scl_t;
  using local_ordinal_scl_hmt = typename LBVH_types::local_ordinal_scl_hmt;
  using local_ordinal_scl_tmt = typename LBVH_types::local_ordinal_scl_tmt;
  using local_ordinal_pairs_t = typename LBVH_types::local_ordinal_pairs_t;
  using local_ordinals_t      = typename LBVH_types::local_ordinals_t;

  using aabb_morton_codes_t = typename LBVH_types::aabb_morton_codes_t;

  using bboxes_3d_view_t   = typename kokkos_aabb_types::bboxes_3d_view_t;
  using bboxes_3d_view_hmt = typename kokkos_aabb_types::bboxes_3d_view_hmt;

  MortonAabbTree(const std::string &base_name, LocalOrdinal num_leaves = 0, bool support_host_boxes = true);

  void init(LocalOrdinal num_leaves);

  void reset(LocalOrdinal num_leaves) {
    resize(num_leaves);
  }

  void resize(LocalOrdinal num_leaves);

  KOKKOS_INLINE_FUNCTION
  void host_set_box(LocalOrdinal box_idx, real_type min_x, real_type max_x, real_type min_y, real_type max_y,
                    real_type min_z, real_type max_z);

  KOKKOS_FORCEINLINE_FUNCTION
  void device_set_box(LocalOrdinal box_idx, real_type min_x, real_type min_y, real_type min_z,
                                            real_type max_x, real_type max_y, real_type max_z) const;

  void sync_from_device() const;
  void sync_to_device() const;

  std::ostream &streamit(std::ostream &os) const;
  std::ostream &streamit(std::ostream &os, size_t box_idx) const;

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
  void dump(const std::string& prefix) const;
#endif

  std::string m_baseName;
  real_type m_globalMinPt[3];
  real_type m_globalMaxPt[3];
  bool m_supportHostBoxes;
  bool m_verbose;

  local_ordinal_scl_t m_numLeaves;
  local_ordinal_scl_hmt hm_numLeaves;
  local_ordinal_scl_tmt tm_numLeaves;
  local_ordinal_scl_t m_numInternalNodes;
  local_ordinal_scl_hmt hm_numInternalNodes;
  local_ordinal_scl_tmt tm_numInternalNodes;

  BoxViewType m_minMaxs;
  BoxViewType_hmt hm_minMaxs;
  bboxes_3d_view_t m_nodeMinMaxs;
  bboxes_3d_view_hmt hm_nodeMinMaxs;

  local_ordinal_pairs_t m_nodeChildren;
  local_ordinals_t m_nodeParents;
  local_ordinals_t m_atomicFlags;
  local_ordinals_t m_leafIds;
  aabb_morton_codes_t m_leafCodes;

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
  typename LBVH_types::local_ordinal_pairs_hmt hm_nodeChildren;
  typename LBVH_types::local_ordinals_hmt hm_nodeParents;
  typename LBVH_types::local_ordinals_hmt hm_leafIds;
  typename LBVH_types::aabb_morton_codes_hmt hm_leafCodes;
#endif
};

template <typename ViewType, typename ExecutionSpace>
MortonAabbTree<ViewType, ExecutionSpace>::MortonAabbTree(const std::string &baseName,
                                                         LocalOrdinal numLeaves,
                                                         bool supportHostBoxes)
  : m_baseName(baseName),
    m_supportHostBoxes(supportHostBoxes),
    m_verbose(false)
{
  init(numLeaves);
}

template <typename ViewType, typename ExecutionSpace>
void MortonAabbTree<ViewType, ExecutionSpace>::init(LocalOrdinal numLeaves)
{
  LocalOrdinal numInternalNodes = std::max(numLeaves - 1, 0);
  LocalOrdinal numNodes = numLeaves + numInternalNodes;

  // Allocate views.
  Kokkos::Timer timer0;
  m_numLeaves = no_init<local_ordinal_scl_t>(compound_name(m_baseName, "numLeaves"));
  m_numInternalNodes = no_init<local_ordinal_scl_t>(compound_name(m_baseName, "numInternalNodes"));
  m_minMaxs = no_init<ViewType>(compound_name(m_baseName, "minMaxs"), numLeaves);
  m_nodeMinMaxs = no_init<bboxes_3d_view_t>(compound_name(m_baseName, "nodeMinMaxs"), numInternalNodes);

  m_nodeChildren = no_init<local_ordinal_pairs_t>(compound_name(m_baseName, "children"), numNodes);
  m_nodeParents = no_init<local_ordinals_t>(compound_name(m_baseName, "parents"), numNodes);

  m_atomicFlags = no_init<local_ordinals_t>(compound_name(m_baseName, "atomicFlags"), numInternalNodes);
  m_leafIds = no_init<local_ordinals_t>(compound_name(m_baseName, "ids"), numLeaves);
  m_leafCodes = no_init<aabb_morton_codes_t>(compound_name(m_baseName, "leafCodes"), numLeaves);
  const double time0 = timer0.seconds();

#ifdef MAS_DEBUG
  std::cout << "Points Views layout type is " << print_layout_name<aabb_points_t::traits::array_layout>() << std::endl;

  std::cout << "Points Host Mirror Views layout type is " << print_layout_name<aabb_points_hmt::traits::array_layout>()
            << std::endl;
#endif

  // Want num_leaves, num_intnl_nodes, and num_nodes readable on both host and device.
  Kokkos::Timer timer1;
  hm_numLeaves = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_numLeaves);
  hm_numLeaves() = numLeaves;
  hm_numInternalNodes = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_numInternalNodes);
  hm_numInternalNodes() = numInternalNodes;
  tm_numLeaves = m_numLeaves;
  tm_numInternalNodes = m_numInternalNodes;
  const double time1 = timer1.seconds();

  // Host mirror views really for debugging
  Kokkos::Timer timer2;

  if (m_supportHostBoxes) {
    hm_minMaxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_minMaxs);
    hm_nodeMinMaxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeMinMaxs);
  }

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
  hm_nodeChildren = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeChildren);
  hm_nodeParents = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeParents);
  hm_leafIds = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_leafIds);
  hm_leafCodes = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_leafCodes);
#endif
  const double time2 = timer2.seconds();

  if (m_verbose) {
    std::cout << std::setw(15) << "MortonAabbTree::init " << numLeaves << " leaves " << time0 << " "
              << " " << time1 << " " << time2 << std::endl;
  }
}

template <typename ViewType, typename ExecutionSpace>
void MortonAabbTree<ViewType, ExecutionSpace>::resize(LocalOrdinal numLeaves)
{
  if ((numLeaves < 0) || (numLeaves == hm_numLeaves())) return;
  Kokkos::Profiling::pushRegion("MortonAabbTree::resize");

  LocalOrdinal numInternalNodes = std::max(numLeaves - 1, 0);
  LocalOrdinal numNodes = numLeaves + numInternalNodes;

  // Allocate views.
  Kokkos::resize(m_minMaxs, numLeaves);
  Kokkos::resize(m_nodeMinMaxs, numInternalNodes);
  Kokkos::resize(m_nodeChildren, numNodes);
  Kokkos::resize(m_nodeParents, numNodes);
  Kokkos::resize(m_atomicFlags, numInternalNodes);
  Kokkos::resize(m_leafIds, numLeaves);
  Kokkos::resize(m_leafCodes, numLeaves);

  // Want num_leaves, num_intnl_nodes, and num_nodes readable on both host and device.
  hm_numLeaves() = numLeaves;
  hm_numInternalNodes() = numInternalNodes;
  Kokkos::deep_copy(m_numLeaves, hm_numLeaves);
  Kokkos::deep_copy(m_numInternalNodes, hm_numInternalNodes);
  tm_numLeaves = m_numLeaves;
  tm_numInternalNodes = m_numInternalNodes;

  if (m_supportHostBoxes) {
    // Host mirror views really for debugging
    hm_minMaxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_minMaxs);
    hm_nodeMinMaxs = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeMinMaxs);
  }

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
  hm_nodeChildren = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeChildren);
  hm_nodeParents = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_nodeParents);
  hm_leafIds = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_leafIds);
  hm_leafCodes = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_leafCodes);
#endif
  Kokkos::Profiling::popRegion();
}

template<typename BoxRealType, typename RealType>
KOKKOS_INLINE_FUNCTION
void set_stk_box(Box<BoxRealType>& box, RealType minX, RealType minY, RealType minZ,
                                     RealType maxX, RealType maxY, RealType maxZ)
{
  box.min_corner()[0] = minX;
  box.min_corner()[1] = minY;
  box.min_corner()[2] = minZ;
  box.max_corner()[0] = maxX;
  box.max_corner()[1] = maxY;
  box.max_corner()[2] = maxZ;
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
void MortonAabbTree<ViewType, ExecutionSpace>::host_set_box(LocalOrdinal boxIdx,
                                                            real_type minX, real_type maxX,
                                                            real_type minY, real_type maxY,
                                                            real_type minZ, real_type maxZ)
{
  hm_minMaxs(boxIdx).box.min_corner()[0] = minX;
  hm_minMaxs(boxIdx).box.min_corner()[1] = minY;
  hm_minMaxs(boxIdx).box.min_corner()[2] = minZ;
  hm_minMaxs(boxIdx).box.max_corner()[0] = maxX;
  hm_minMaxs(boxIdx).box.max_corner()[1] = maxY;
  hm_minMaxs(boxIdx).box.max_corner()[2] = maxZ;
}

template <typename ViewType, typename ExecutionSpace>
KOKKOS_FORCEINLINE_FUNCTION
void MortonAabbTree<ViewType, ExecutionSpace>::device_set_box(LocalOrdinal boxIdx,
                                                              real_type minX, real_type minY, real_type minZ,
                                                              real_type maxX, real_type maxY, real_type maxZ) const
{
  m_minMaxs(boxIdx).box.min_corner()[0] = minX;
  m_minMaxs(boxIdx).box.min_corner()[1] = minY;
  m_minMaxs(boxIdx).box.min_corner()[2] = minZ;
  m_minMaxs(boxIdx).box.max_corner()[0] = maxX;
  m_minMaxs(boxIdx).box.max_corner()[1] = maxY;
  m_minMaxs(boxIdx).box.max_corner()[2] = maxZ;
}

template <typename ViewType, typename ExecutionSpace>
void MortonAabbTree<ViewType, ExecutionSpace>::sync_from_device() const
{
  Kokkos::deep_copy(hm_numLeaves, m_numLeaves);
  Kokkos::deep_copy(hm_numInternalNodes, m_numInternalNodes);

  if (m_supportHostBoxes) {
    Kokkos::deep_copy(hm_minMaxs, m_minMaxs);
    Kokkos::deep_copy(hm_nodeMinMaxs, m_nodeMinMaxs);
  }

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
  Kokkos::deep_copy(hm_nodeChildren, m_nodeChildren);
  Kokkos::deep_copy(hm_nodeParents, m_nodeParents);

  Kokkos::deep_copy(hm_leafIds, m_leafIds);
  Kokkos::deep_copy(hm_leafCodes, m_leafCodes);
#endif
}

template <typename ViewType, typename ExecutionSpace>
void MortonAabbTree<ViewType, ExecutionSpace>::sync_to_device() const
{
  Kokkos::deep_copy(m_numLeaves, hm_numLeaves);
  Kokkos::deep_copy(m_numInternalNodes, hm_numInternalNodes);

  if (m_supportHostBoxes) {
    Kokkos::deep_copy(m_minMaxs, hm_minMaxs);
    Kokkos::deep_copy(m_nodeMinMaxs, hm_nodeMinMaxs);
  }

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
  Kokkos::deep_copy(m_nodeChildren, hm_nodeChildren);
  Kokkos::deep_copy(m_nodeParents, hm_nodeParents);

  Kokkos::deep_copy(m_leafIds, hm_leafIds);
  Kokkos::deep_copy(m_leafCodes, hm_leafCodes);
#endif
}

template <typename ViewType, typename ExecutionSpace>
std::ostream &MortonAabbTree<ViewType, ExecutionSpace>::streamit(std::ostream &os) const
{
  os << "{MortonAabbTree " << m_baseName << " " << sizeof(morton_code_t)
     << " (" << m_globalMinPt[0] << " " << m_globalMinPt[1] << " " << m_globalMinPt[2]
     << ") ( " << m_globalMaxPt[0] << " " << m_globalMaxPt[1] << " " << m_globalMaxPt[2]
     << ") " << hm_numLeaves << " " << hm_numInternalNodes << std::endl;

  if (m_supportHostBoxes) {
    for (LocalOrdinal idx = 0; idx < hm_numLeaves(); ++idx) {
      os << " {(" << hm_minMaxs(idx, 0) << " " << hm_minMaxs(idx, 1) << " " << hm_minMaxs(idx, 2)
         << ") (" << hm_minMaxs(idx, 4) << " " << hm_minMaxs(idx, 4) << " " << hm_minMaxs(idx, 5) << ")"
      #ifdef DEBUG_MORTON_ACCELERATED_SEARCH
         << " (" << hm_nodeChildren(idx, 0) << ", " << hm_nodeChildren(idx, 1) << ") " << hm_nodeParents(idx) << " "
         << hm_leafIds(idx) << " 0x" << std::setfill('0') << std::setw(sizeof(morton_code_t) * 2) << std::hex
         << hm_leafCodes(idx) << std::dec
      #endif
         << "}" << std::endl;
    }
    for (LocalOrdinal idx = hm_numLeaves(); idx < hm_numLeaves()+hm_numInternalNodes(); ++idx) {
      LocalOrdinal nodeIdx = idx - hm_numLeaves();
      os << " {(" << hm_nodeMinMaxs(nodeIdx, 0) << " " << hm_nodeMinMaxs(nodeIdx, 1) << " " << hm_nodeMinMaxs(nodeIdx, 2)
         << ") (" << hm_nodeMinMaxs(nodeIdx, 3) << " " << hm_nodeMinMaxs(nodeIdx, 4) << " " << hm_nodeMinMaxs(nodeIdx, 5) << ")"
      #ifdef DEBUG_MORTON_ACCELERATED_SEARCH
         << " (" << hm_nodeChildren(idx, 0) << ", " << hm_nodeChildren(idx, 1) << ") " << hm_nodeParents(idx)
      #endif
         << "}" << std::endl;
    }
  }
  os << "}";
  return os;
}

template <typename ViewType, typename ExecutionSpace>
std::ostream &MortonAabbTree<ViewType, ExecutionSpace>::streamit(std::ostream &os, size_t idx) const
{
  if (m_supportHostBoxes) {
    if (idx < hm_numLeaves()) {
      os << " {(" << hm_minMaxs(idx, 0) << " " << hm_minMaxs(idx, 1) << " " << hm_minMaxs(idx, 2)
         << ") (" << hm_minMaxs(idx, 3) << " " << hm_minMaxs(idx, 4) << " " << hm_minMaxs(idx, 5) << ")}";
    }
    else {
      size_t nodeIdx = idx - hm_numLeaves();
      os << " {(" << hm_nodeMinMaxs(nodeIdx, 0) << " " << hm_nodeMinMaxs(nodeIdx, 1) << " " << hm_nodeMinMaxs(nodeIdx, 2)
         << ") (" << hm_nodeMinMaxs(nodeIdx, 3) << " " << hm_nodeMinMaxs(nodeIdx, 4) << " " << hm_nodeMinMaxs(nodeIdx, 5) << ")}";
    }
  }
  else {
    os << "MortonAabbTree::streamit(.) m_supportHostBoxes = 0. ";
  }
  return os;
}

#ifdef DEBUG_MORTON_ACCELERATED_SEARCH
template<typename ViewType, typename ExecutionSpace>
void MortonAabbTree<ViewType, ExecutionSpace>::dump(const std::string& prefix) const
{
  std::cout << prefix << ": dump, numLeaves=" << hm_numLeaves() << ", numInternalNodes="
            << hm_numInternalNodes() << std::endl;
  sync_from_device();
  {
    std::ofstream of(prefix+"_leafCodes.txt");
    for (size_t i = 0; i < hm_leafCodes.extent(0); ++i) {
      of << i << " " << hm_leafCodes(i) << std::endl;
    }
  }
  {
    std::ofstream of(prefix+"_leafIds.txt");
    for (size_t i = 0; i < hm_leafIds.extent(0); ++i) {
      of << i << " " << hm_leafIds(i) << std::endl;
    }
  }
  {
    std::ofstream of(prefix+"_nodeChildren.txt");
    for (size_t i = 0; i < hm_nodeChildren.extent(0); ++i) {
      of << i << " " << hm_nodeChildren(i,0) << " " << hm_nodeChildren(i,1) << std::endl;
    }
  }
  {
    std::ofstream of(prefix+"_nodeParents.txt");
    for (size_t i = 0; i < hm_nodeParents.extent(0); ++i) {
      of << i << " " << hm_nodeParents(i) << std::endl;
    }
  }
  {
    std::ofstream of(prefix+"_minMaxs.txt");
    for (size_t i = 0; i < hm_minMaxs.extent(0); ++i) {
      of << i;
      for (size_t j = 0; j < 6; ++j) {
        of << " " << std::setprecision(7) << hm_minMaxs(i,j);
      }
      of << std::endl;
    }
  }
  {
    std::ofstream of(prefix+"_nodeMinMaxs.txt");
    for (size_t i = 0; i < hm_nodeMinMaxs.extent(0); ++i) {
      of << i;
      for (size_t j = 0; j < 6; ++j) {
        of << " " << std::setprecision(7) << hm_nodeMinMaxs(i,j);
      }
      of << std::endl;
    }
  }
}
#endif

}  // namespace stk::search

#endif  // STK_SEARCH_MORTON_LBVH_MORTONLBVH_TREE_HPP
