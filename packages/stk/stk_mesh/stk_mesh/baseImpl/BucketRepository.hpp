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

#ifndef stk_mesh_BucketRepository_hpp
#define stk_mesh_BucketRepository_hpp

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/Types.hpp>      // for EntityRank, OrdinalVector, etc
#include <vector>                       // for vector
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class EntitySorterBase; } }
namespace stk { namespace mesh { namespace impl { class Partition; } } }

namespace stk {
namespace mesh {
namespace impl {


class BucketRepository
{
public:
  //------------------------------------
  /** \brief  Query the upper bound on the number of mesh entities
    *         that may be associated with a single bucket.
    */
  BucketRepository(
      BulkData & mesh,
      unsigned entity_rank_count,
      unsigned initialBucketCapacity = get_default_initial_bucket_capacity(),
      unsigned maximumBucketCapacity = get_default_maximum_bucket_capacity()
      );

  ~BucketRepository();

  /** \brief  Query to get all buckets of a given entity rank.
   *
   *  Don't call inside BucketRepository member functions!
   */
  const BucketVector & buckets( EntityRank rank ) const;

  BulkData& mesh() const { return m_mesh; }

  void set_needs_to_be_sorted(stk::mesh::Bucket &bucket, bool needsSorting);
  void internal_default_sort_bucket_entities(bool mustSortFacesByNodeIds=false);
  void internal_custom_sort_bucket_entities(const EntitySorterBase& sorter);

  Bucket *get_bucket(EntityRank entity_rank, int bucket_id) const;

  friend class Partition;
  template<typename NgpMemSpace> friend class stk::mesh::DeviceMeshT;

  void add_entity_with_part_memberships(const Entity entity,
                                        const EntityRank arg_entity_rank,
                                        const OrdinalVector &parts);

  void change_entity_part_membership(const MeshIndex &meshIndex, const OrdinalVector &parts);

  void remove_entity(const MeshIndex &meshIndex);

  Partition *get_or_create_partition(const EntityRank arg_entity_rank ,
                                     const OrdinalVector &parts);

  Partition *get_partition(const EntityRank arg_entity_rank ,
                           const OrdinalVector &parts,
                           std::vector<Partition*>::iterator& ik);

  Partition *create_partition(const EntityRank arg_entity_rank ,
                              const OrdinalVector &parts,
                              std::vector<Partition*>::iterator& ik);

  // For use by BulkData::internal_modification_end().
  void internal_modification_end();

  // Update m_buckets from the partitions.
  void sync_from_partitions();
  void sync_from_partitions(EntityRank rank);

  const std::vector<Partition *>& get_partitions(EntityRank rank) const;

  Partition* get_partition(const EntityRank arg_entity_rank, const OrdinalVector &parts);

  bool being_destroyed() const { return m_being_destroyed; }

  unsigned get_bucket_capacity() const { return m_maximumBucketCapacity; }
  unsigned get_initial_bucket_capacity() const { return m_initialBucketCapacity; }
  unsigned get_maximum_bucket_capacity() const { return m_maximumBucketCapacity; }

  void delete_bucket(Bucket * bucket);

  void set_need_sync_from_partitions(EntityRank entityRank) { m_need_sync_from_partitions[entityRank] = true; }

  void set_remove_mode_tracking();
  void set_remove_mode_fill_and_sort();

private:
  BucketRepository();

  Bucket *allocate_bucket(EntityRank entityRank,
                          const std::vector<unsigned> & key,
                          unsigned initialCapacity,
                          unsigned maximumCapacity);

  void deallocate_bucket(Bucket *bucket);

  void sync_bucket_ids(EntityRank entity_rank);

  void ensure_data_structures_sized();

  BulkData & m_mesh ;

  // Vector of bucket pointers for each rank.
  std::vector< BucketVector >   m_buckets ;

  std::vector<std::vector<Partition *> > m_partitions;
  std::vector<bool> m_need_sync_from_partitions;

  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;
  bool m_being_destroyed;
};

inline
Bucket *BucketRepository::get_bucket(EntityRank entity_rank, int bucket_id) const
{
  const BucketVector & all_buckets_for_rank = m_buckets[entity_rank];
  STK_ThrowAssert(static_cast<size_t>(bucket_id) < all_buckets_for_rank.size());
  return all_buckets_for_rank[bucket_id];
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_BucketRepository_hpp
