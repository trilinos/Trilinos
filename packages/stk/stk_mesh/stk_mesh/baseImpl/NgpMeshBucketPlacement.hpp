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

#ifndef STK_MESH_NGPMESHBUCKETPLACEMENT_HPP
#define STK_MESH_NGPMESHBUCKETPLACEMENT_HPP

#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/NgpMeshTypes.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"

namespace stk::mesh::impl {

struct NumNewBucketsToAddPerPartition {
  EntityRank rank = topology::INVALID_RANK;
  unsigned partitionId = INVALID_INDEX;
  unsigned numBucketsToAdd = 0;

  KOKKOS_INLINE_FUNCTION
  bool operator==(NumNewBucketsToAddPerPartition const& rhs) const
  {
    return (rank == rhs.rank) && (partitionId == rhs.partitionId) && (numBucketsToAdd == rhs.numBucketsToAdd);
  }
};

struct GrowLastBucketInPartition {
  EntityRank rank = topology::INVALID_RANK;
  unsigned partitionId = INVALID_INDEX;
  bool growLastBucket = false;

  KOKKOS_INLINE_FUNCTION
  bool operator==(GrowLastBucketInPartition const& rhs) const
  {
    return (rank == rhs.rank) && (partitionId == rhs.partitionId) && (growLastBucket == rhs.growLastBucket);
  }
};

template <typename EntitySrcDestType>
KOKKOS_INLINE_FUNCTION
bool is_same_dest_partition_group(EntitySrcDestType const& lhs, EntitySrcDestType const& rhs)
{
  return (lhs.rank == rhs.rank) && (lhs.destPartitionId == rhs.destPartitionId);
}

template <typename MeshType, typename EntitySrcDestView, typename NewBucketsToAddViewType,
          typename GrowLastBucketViewType>
void assign_dest_bucket_id_and_ordinal(MeshType const& ngpMesh, EntitySrcDestView& entitySrcDestView,
                                       NewBucketsToAddViewType& numNewBucketsToAddInPartitions,
                                       GrowLastBucketViewType& growLastBucketInPartition)
{
  using ExecSpace = typename MeshType::MeshExecSpace;
  using MemSpace = typename EntitySrcDestView::memory_space;
  using UnsignedViewType = Kokkos::View<unsigned*, MemSpace>;

  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();
  Kokkos::Profiling::pushRegion("check_mesh_consistency");
  deviceBucketRepo.check_mesh_consistency();
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("remove_unmoved_entities");
  auto begin = Kokkos::Experimental::begin(entitySrcDestView);
  auto end = Kokkos::Experimental::remove_if(
      ExecSpace{}, entitySrcDestView,
      KOKKOS_LAMBDA(typename EntitySrcDestView::value_type const& srcDest) {
        return srcDest.srcPartitionId == srcDest.destPartitionId;
      });

  auto newLength = Kokkos::Experimental::distance(begin, end);
  Kokkos::resize(entitySrcDestView, newLength); 
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("sort_entity_src_dest");
  if (!Kokkos::Experimental::is_sorted(ExecSpace{}, entitySrcDestView)) {
    Kokkos::sort(entitySrcDestView);
  }
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("assign_dest_buckets_per_partition_group");
  auto numEntities = entitySrcDestView.extent(0);
  UnsignedViewType destPartitionGroupId(Kokkos::view_alloc(Kokkos::WithoutInitializing, "destPartitionGroupId"), numEntities);
  UnsignedViewType destPartitionGroupStartIdx(Kokkos::view_alloc(Kokkos::WithoutInitializing, "destPartitionGroupStartIdx"), numEntities);

  // For numbering destination partition groups (id and starting idx)
  Kokkos::parallel_scan("number_dest_partition_groups", Kokkos::RangePolicy<ExecSpace>(0, numEntities),
    KOKKOS_LAMBDA(const unsigned i, unsigned& update, const bool isFinal) {
      const bool isGroupStart = (i == 0) || !is_same_dest_partition_group(entitySrcDestView(i), entitySrcDestView(i-1));
      if (isFinal && isGroupStart) {
        destPartitionGroupStartIdx(update) = i;
      }
      if (isGroupStart) { ++update; }
      if (isFinal) {
        destPartitionGroupId(i) = update - 1;
      }
    }
  );

  Kokkos::parallel_for("assign_dest_bucket_id_and_ordinal", Kokkos::RangePolicy<ExecSpace>(0, numEntities),
    KOKKOS_LAMBDA(const unsigned i) {
      auto& srcDest = entitySrcDestView(i);
      auto rank = srcDest.rank;
      auto destPartitionId = srcDest.destPartitionId;
      auto groupId = destPartitionGroupId(i);
      auto posInGroup = i - destPartitionGroupStartIdx(groupId);

      auto maxBucketCapacity = deviceBucketRepo.get_maximum_bucket_capacity();
      auto partition = deviceBucketRepo.get_partition(rank, destPartitionId);
      auto lastBucketIdx = partition->get_last_avail_bucket_index(maxBucketCapacity);
      auto appendToLastBucket = (lastBucketIdx != INVALID_INDEX);

      unsigned startBucketIdxInPartition = partition->m_buckets.size();
      unsigned startBucketOrd = 0;
      unsigned startBucketCapacity = maxBucketCapacity;

      if (appendToLastBucket) {
        auto destBucket = partition->m_buckets[lastBucketIdx];
        startBucketIdxInPartition = lastBucketIdx;
        startBucketOrd = destBucket->size();
        startBucketCapacity = destBucket->capacity();
      }

      auto posFromStartBucket = startBucketOrd + posInGroup;
      srcDest.destBucketId = INVALID_BUCKET_ID;  // to be determined later; can't be done now
      srcDest.destBucketIndexInPartition = startBucketIdxInPartition + posFromStartBucket / maxBucketCapacity;
      srcDest.destBucketOrd = posFromStartBucket % maxBucketCapacity;

      auto isGroupEnd = (i == numEntities - 1) || !is_same_dest_partition_group(entitySrcDestView(i), entitySrcDestView(i+1));
      if (isGroupEnd) {
        STK_NGP_ThrowAssertMsg(!appendToLastBucket ||
                               partition->m_buckets[lastBucketIdx]->size() == partition->m_buckets[lastBucketIdx]->get_active_entity_span(),
                               "Appending to a bucket in the middle of a partition. sync_from_partitions was not called.");

        auto numBucketsToAdd = (appendToLastBucket ? 0 : 1) + posFromStartBucket / maxBucketCapacity;
        numNewBucketsToAddInPartitions(groupId) = {rank, destPartitionId, numBucketsToAdd};

        if (appendToLastBucket && startBucketCapacity < maxBucketCapacity &&
            posFromStartBucket >= startBucketCapacity) {
          growLastBucketInPartition(groupId) = {rank, destPartitionId, true};
        }
      }
    }
  );
  Kokkos::fence();
  Kokkos::Profiling::popRegion();
}

template <typename MeshType, typename EntitySrcDestView>
void set_dest_bucket_ids(MeshType const& ngpMesh, EntitySrcDestView& entitySrcDestView)
{
  auto& deviceBucketRepo = ngpMesh.get_device_bucket_repository();

  Kokkos::parallel_for(entitySrcDestView.extent(0),
    KOKKOS_LAMBDA(const int i) {
      auto& entitySrcDest = entitySrcDestView(i);
      auto rank = entitySrcDest.rank;
      auto destPartitionId = entitySrcDest.destPartitionId;
      auto destPartition = deviceBucketRepo.get_partition(rank, destPartitionId);
      auto destBucketId = destPartition->get_bucket(entitySrcDest.destBucketIndexInPartition)->bucket_id();
      entitySrcDest.destBucketId = destBucketId;
    }
  );
}

} // namespace stk::mesh::impl

#endif
