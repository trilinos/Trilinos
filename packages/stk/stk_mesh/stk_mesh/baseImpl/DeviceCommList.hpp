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

#ifndef STK_MESH_DEVICECOMMLIST_HPP
#define STK_MESH_DEVICECOMMLIST_HPP

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/NgpTypes.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include "Kokkos_Core.hpp"
#include <cstddef>

namespace stk::mesh::impl {

// Re-gather each comm map slot's FastMeshIndex from the entity recorded in it.  A device modification
// that does not change who communicates what renumbers buckets but keeps every slot's entity, so this
// keeps the map's layout, which must stay matched position-by-position with every peer's.
template <typename NgpMemSpace>
void refresh_device_comm_map_meshIndices(const FastSharedCommMapViewType<NgpMemSpace>& rankCommMap, const MeshIndexType<NgpMemSpace>& meshIndices, const Kokkos::View<Entity*, NgpMemSpace>& rankCommMapEntities) {
  using MeshExecSpace = typename NgpMemSpace::execution_space;
  Kokkos::parallel_for("refresh_device_comm_map_meshIndices", Kokkos::RangePolicy<MeshExecSpace>(0, rankCommMap.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      rankCommMap(i) = meshIndices(rankCommMapEntities(i).local_offset());
    }
  );
}

// Number of entries in rankCommMapEntities whose entity is in sortedEntities, which must be sorted by
// local offset and unique (as get_destroyable_entities returns them).
template <typename NgpMemSpace>
unsigned count_comm_map_entries_in_entities_list(const Kokkos::View<Entity*, NgpMemSpace>& rankCommMapEntities,
                                                 const Kokkos::View<Entity*, NgpMemSpace>& sortedEntities) {
  using MeshExecSpace = typename NgpMemSpace::execution_space;
  const size_t numSorted = sortedEntities.extent(0);
  unsigned numFound = 0;
  Kokkos::parallel_reduce("count_comm_map_entries_for",
    Kokkos::RangePolicy<MeshExecSpace>(0, rankCommMapEntities.extent(0)),
    KOKKOS_LAMBDA(const size_t i, unsigned& localCount) {
      const auto offset = rankCommMapEntities(i).local_offset();
      size_t first = 0;
      size_t last = numSorted;
      while (first < last) {
        const size_t mid = first + (last - first) / 2;
        if (sortedEntities(mid).local_offset() < offset) { first = mid + 1; }
        else { last = mid; }
      }
      if (first < numSorted && sortedEntities(first).local_offset() == offset) { ++localCount; }
    }, Kokkos::Sum<unsigned>(numFound)
  );
  return numFound;
}

template <typename OffsetViewType, typename NumSharedViewType>
KOKKOS_INLINE_FUNCTION
Kokkos::pair<size_t, size_t> comm_map_range(const OffsetViewType& offset, const NumSharedViewType& numShared,
                                            int proc, bool includeGhosts)
{
  const size_t dataBegin = offset(proc);
  const size_t dataEnd = includeGhosts ? offset(proc+1) : dataBegin + numShared(proc);
  return Kokkos::pair<size_t, size_t>(dataBegin, dataEnd);
}

template <typename DstArrayType, typename SrcArrayType>
void copy_per_rank_comm_views(DstArrayType& dst, const SrcArrayType& src, EntityRank endRank) {
  for (EntityRank rank = stk::topology::BEGIN_RANK; rank < endRank; ++rank) {
    if (dst[rank].extent(0) != src[rank].extent(0)) {
      Kokkos::realloc(Kokkos::WithoutInitializing, dst[rank], src[rank].extent(0));
    }
    Kokkos::deep_copy(dst[rank], src[rank]);
  }
}

}  // namespace stk::mesh::impl

#endif
