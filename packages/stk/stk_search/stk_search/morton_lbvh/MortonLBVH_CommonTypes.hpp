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

#ifndef STK_SEARCH_MORTON_LBVH_MORTON_COMMON_TYPES_HPP
#define STK_SEARCH_MORTON_LBVH_MORTON_COMMON_TYPES_HPP

#include <Kokkos_Core.hpp>
#include <string>

namespace stk::search {

// Within an MPI rank ints should suffice.  LocalOrdinal MUST BE SIGNED!
using LocalOrdinal = int;

#ifdef SMALL_MORTON
using morton_code_t = uint32_t;  // 32 bit codes
#else
using morton_code_t = uint64_t;  // 64 bit codes
#endif

template <typename ExecutionSpace>
struct MortonLbvhMemorySpace
{
  using memory_space = typename ExecutionSpace::memory_space;
};

template <class ViewT>
inline ViewT no_init(const std::string &name)
{
  return ViewT(Kokkos::ViewAllocateWithoutInitializing(name));
}

template <class ViewT>
inline ViewT no_init(const std::string &name, unsigned len)
{
  return ViewT(Kokkos::ViewAllocateWithoutInitializing(name), len);
}

template <class ViewT>
inline ViewT with_init(const std::string &name, unsigned len)
{
  return ViewT(name, len);
}

template <typename ExecutionSpace>
struct MortonLbvhTypes
{
  using memory_space = typename ExecutionSpace::memory_space;

  // View of a single LocalOrdinal.
  using local_ordinal_scl_t   = Kokkos::View<LocalOrdinal, memory_space>;
  using local_ordinal_scl_hmt = typename local_ordinal_scl_t::HostMirror;
  using local_ordinal_scl_tmt = Kokkos::View<const LocalOrdinal, memory_space, Kokkos::MemoryRandomAccess>;

  // Will need a view of LocalOrdinal Scalars.
  using local_ordinals_t   = Kokkos::View<LocalOrdinal *, memory_space>;
  using local_ordinals_hmt = typename local_ordinals_t::HostMirror;
  using local_ordinals_tmt = Kokkos::View<LocalOrdinal *, memory_space, Kokkos::MemoryRandomAccess>;

  // Will need a view of LocalOrdinalPairs.
  using local_ordinal_pairs_t   = Kokkos::View<LocalOrdinal * [2], memory_space>;
  using local_ordinal_pairs_hmt = typename local_ordinal_pairs_t::HostMirror;
  using local_ordinal_pairs_tmt = Kokkos::View<LocalOrdinal * [2], memory_space, Kokkos::MemoryRandomAccess>;

  using aabb_morton_codes_t   = Kokkos::View<morton_code_t *, ExecutionSpace>;
  using aabb_morton_codes_hmt = typename aabb_morton_codes_t::HostMirror;
  using aabb_morton_codes_tmt = Kokkos::View<const morton_code_t *, ExecutionSpace, Kokkos::MemoryRandomAccess>;
};

template <typename RealType, typename ExecutionSpace>
struct MortonAabbTypes
{
  using memory_space = typename MortonLbvhTypes<ExecutionSpace>::memory_space;

  // Points
  using aabb_points_t         = Kokkos::View<RealType * [3], memory_space>;
  using aabb_points_hmt       = typename aabb_points_t::HostMirror;
  using aabb_const_points_t   = Kokkos::View<const RealType * [3], memory_space>;
  using aabb_const_points_tmt = Kokkos::View<const RealType * [3], memory_space, Kokkos::MemoryRandomAccess>;

  // We'll use these when convert from using (min_pt, max_pt) pairs.
  using bboxes_3d_view_t       = Kokkos::View<RealType * [6], Kokkos::LayoutRight, memory_space>;
  using bboxes_3d_view_hmt     = typename bboxes_3d_view_t::HostMirror;
  using bboxes_const_3d_view_t = Kokkos::View<const RealType * [6], Kokkos::LayoutRight, memory_space>;
  using bboxes_3d_view_amt     = Kokkos::View<RealType * [6], Kokkos::LayoutRight, memory_space,
                                              Kokkos::MemoryTraits<Kokkos::Atomic>>;
};

struct morton_code_id_pair
{
  morton_code_t m_code;
  LocalOrdinal m_id;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator<(const morton_code_id_pair &rhs) const
  {
    return (m_code < rhs.m_code) || ((m_code == rhs.m_code) && (m_id < rhs.m_id));
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator>(const morton_code_id_pair &rhs) const
  {
    return (m_code > rhs.m_code) || ((m_code == rhs.m_code) && (m_id > rhs.m_id));
  }
};

inline std::string compound_name(const std::string &baseName, const std::string &componentName)
{
  std::string sep("_");
  return baseName + sep + componentName;
}

}  // namespace stk::search

#endif  // STK_SEARCH_MORTON_LBVH_MORTON_COMMON_TYPES_HPP
