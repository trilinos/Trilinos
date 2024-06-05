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

#ifndef MORTONLBVH_COLLISIONLIST_HPP
#define MORTONLBVH_COLLISIONLIST_HPP

#include <stk_search/morton_lbvh/MortonLBVH_CommonTypes.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>
#include <utility>
#include <set>

namespace stk::search {

template <typename ExecutionSpace>
struct CollisionList
{
  using LVBH_types = MortonLbvhTypes<ExecutionSpace>;
  using memory_space = typename LVBH_types::memory_space;
  using data_val_type      = typename LVBH_types::local_ordinal_pairs_t;
  using host_data_val_type = typename data_val_type::HostMirror;
  using idx_type      = Kokkos::View<int, memory_space>;
  using host_idx_type = typename idx_type::HostMirror;

  struct CorruptorTag {};

  CollisionList(const std::string &dataName, LocalOrdinal capacity = 0);

  void reset(LocalOrdinal newCapacity = -1);

  KOKKOS_INLINE_FUNCTION
  void operator()(const CorruptorTag &, const int &idx) const { m_data(idx, 1) += 10; }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_collision(int idx, int domainIdx, int rangeIdx) const;

  KOKKOS_FORCEINLINE_FUNCTION
  bool push_back(LocalOrdinal domainIdx, LocalOrdinal rangeIdx, bool flipped = false) const;

  KOKKOS_INLINE_FUNCTION
  bool set_num_collisions(int numCollisions) const;

  int get_num_collisions() const;
  int get_capacity() const { return m_capacity; }

  void make_contents_corrupt();

  void sync_from_device();
  void sync_to_device();

  std::string m_dataName;
  std::string m_idxName;
  LocalOrdinal m_capacity;

  data_val_type m_data;
  idx_type m_idx;
  host_data_val_type hm_data;
  host_idx_type hm_idx;
};

template <typename ExecutionSpace>
CollisionList<ExecutionSpace>::CollisionList(const std::string &baseName, LocalOrdinal capacity)
  : m_dataName(compound_name(baseName, "data")),
    m_idxName(compound_name(baseName, "idx")),
    m_capacity(capacity),
    m_data(no_init<data_val_type>(m_dataName, m_capacity)),
    m_idx(m_idxName)
{
#ifdef BANDC_DEBUG
  std::cout << "Creating CollisionList mirror view." << std::endl;
#endif
  // This should probably be done lazily, when sync_from_device() gets called.
  hm_data = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_data);
  hm_idx = Kokkos::create_mirror_view(m_idx);
#ifdef BANDC_DEBUG
  std::cout << "Collision List Mirror View created" << std::endl;
#endif
}

template <typename ExecutionSpace>
void CollisionList<ExecutionSpace>::reset(LocalOrdinal newCapacity)
{
  if ((newCapacity >= 0) && (newCapacity != m_capacity)) {
    m_capacity = newCapacity;
    m_data = no_init<data_val_type>(m_dataName, m_capacity);
    hm_data = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, m_data);
  }
  hm_idx() = 0;
  Kokkos::deep_copy(m_idx, hm_idx);
}

template <typename ExecutionSpace>
KOKKOS_FORCEINLINE_FUNCTION
void CollisionList<ExecutionSpace>::set_collision(int idx, int domainIdx, int rangeIdx) const
{
  m_data(idx, 0) = domainIdx;
  m_data(idx, 1) = rangeIdx;
}

template <typename ExecutionSpace>
KOKKOS_FORCEINLINE_FUNCTION
bool CollisionList<ExecutionSpace>::push_back(LocalOrdinal domainIdx, LocalOrdinal rangeIdx, bool flipped) const
{
  typename idx_type::value_type idx = Kokkos::atomic_fetch_add(&m_idx(), 1);
  if (idx < static_cast<int>(m_data.extent(0))) {
    if (flipped) {
      set_collision(idx, rangeIdx, domainIdx);
    }
    else {
      set_collision(idx, domainIdx, rangeIdx);
    }
    return true;
  }
  return false;
}

template <typename ExecutionSpace>
KOKKOS_INLINE_FUNCTION
bool CollisionList<ExecutionSpace>::set_num_collisions(int numCollisions) const
{
  if (numCollisions <= static_cast<int>(m_data.extent(0))) {
    m_idx() = numCollisions;
    return true;
  }
  else {
    return false;
  }
}

template <typename ExecutionSpace>
int CollisionList<ExecutionSpace>::get_num_collisions() const
{
  if (hm_idx() <= 0) {
    Kokkos::deep_copy(hm_idx, m_idx);
  }
  return hm_idx();
}

template <typename ExecutionSpace>
void CollisionList<ExecutionSpace>::make_contents_corrupt()
{
  Kokkos::parallel_for(Kokkos::RangePolicy<CorruptorTag>(0, m_data.extent(0)), *this);
}

template <typename ExecutionSpace>
void CollisionList<ExecutionSpace>::sync_from_device()
{
  Kokkos::deep_copy(hm_data, m_data);
  Kokkos::deep_copy(hm_idx, m_idx);
}

template <typename ExecutionSpace>
void CollisionList<ExecutionSpace>::sync_to_device()
{
  Kokkos::deep_copy(m_data, hm_data);
  Kokkos::deep_copy(m_idx, hm_idx);
}


template <typename ExecutionSpace>
bool SameContents(CollisionList<ExecutionSpace> &listA, CollisionList<ExecutionSpace> &listB)
{
  listA.sync_from_device();
  listB.sync_from_device();

  int numCollisionsA = listA.get_num_collisions();
  int numCollisionsB = listB.get_num_collisions();

  if (numCollisionsA != numCollisionsB) {
    return false;
  }

  using collision_pair = std::pair<LocalOrdinal, LocalOrdinal>;

  std::set<collision_pair> setA;
  for (LocalOrdinal i = 0; i < numCollisionsA; ++i) {
    setA.insert(collision_pair(listA.hm_data(i, 0), listA.hm_data(i, 1)));
  }
  std::set<collision_pair> setB;
  for (LocalOrdinal i = 0; i < numCollisionsB; ++i) {
    setB.insert(collision_pair(listB.hm_data(i, 0), listB.hm_data(i, 1)));
  }

  return (setA == setB);
}

}

#endif // MORTONLBVH_COLLISIONLIST_HPP
