// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include "stk_mesh/base/ConnectivityMap.hpp"  // for ConnectivityMap
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { namespace impl { class Partition; } } }
namespace stk { namespace mesh { namespace utest { struct SyncToPartitions; } } }


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
  static const unsigned max_bucket_capacity = 1024;
  static const unsigned default_bucket_capacity = 512;

  BucketRepository(
      BulkData & mesh,
      unsigned entity_rank_count,
      const ConnectivityMap & connectivity_map,
      unsigned bucket_capacity = default_bucket_capacity
      );

  ~BucketRepository();

  /** \brief  Query to get all buckets of a given entity rank.
   *
   *  Don't call inside BucketRepository member functions!
   */
  const BucketVector & buckets( EntityRank rank ) const
  {
    static const BucketVector emptyBucketVector;

    if( rank < m_buckets.size() )
    {
      if (m_need_sync_from_partitions[rank])
      {
        const_cast<BucketRepository *>(this)->sync_from_partitions(rank);
      }

      return m_buckets[ rank ];
    }
    else
    {
      return emptyBucketVector;
    }
  }

  BulkData& mesh() const { return m_mesh; }

  //------------------------------------
  size_t total_field_data_footprint(const FieldBase &f, EntityRank rank) const;

  void set_needs_to_be_sorted(stk::mesh::Bucket &bucket, bool needsSorting);
  void internal_default_sort_bucket_entities();
  void internal_custom_sort_bucket_entities(const EntitySorterBase& sorter);

  Bucket *get_bucket(EntityRank entity_rank, int bucket_id) const;

  template <class RankType>
  inline
  Bucket *get_bucket(RankType rank_id) const;

  ////
  //// Partitions are now the primary location of buckets.
  ////

  friend class Partition;
  friend struct stk::mesh::utest::SyncToPartitions;

  void add_entity_with_part_memberships(const Entity entity,
                                        const EntityRank arg_entity_rank,
                                        const OrdinalVector &parts);

  void change_entity_part_membership(const MeshIndex &meshIndex, const OrdinalVector &parts);

  void remove_entity(const MeshIndex &meshIndex);

  Partition *get_or_create_partition(const EntityRank arg_entity_rank ,
                                     const OrdinalVector &parts);

  // For use by BulkData::internal_modification_end().
  void internal_modification_end();

  // Update m_buckets from the partitions.
  void sync_from_partitions();
  void sync_from_partitions(EntityRank rank);

  // Used in unit tests.  Returns the current partitions.
  std::vector<Partition *> get_partitions(EntityRank rank) const;

  const ConnectivityMap& connectivity_map() const { return m_connectivity_map; }

  bool being_destroyed() const { return m_being_destroyed; }

  unsigned get_bucket_capacity() const { return m_bucket_capacity; }

  void delete_bucket(Bucket * bucket);

private:
  BucketRepository();

  Bucket *allocate_bucket(EntityRank arg_entity_rank,
                          const std::vector<unsigned> & arg_key,
                          size_t arg_capacity);

  void deallocate_bucket(Bucket *bucket);

  void sync_bucket_ids(EntityRank entity_rank);

  void ensure_data_structures_sized();

  BulkData & m_mesh ; // Associated Bulk Data Aggregate

  // Vector of bucket pointers by rank.  This is now a cache and no longer the primary
  // location of Buckets when USE_STK_MESH_IMPL_PARTITION is #defined.
  std::vector< BucketVector >   m_buckets ;

  std::vector<std::vector<Partition *> > m_partitions;
  std::vector<bool> m_need_sync_from_partitions;

  ConnectivityMap m_connectivity_map;

  unsigned m_bucket_capacity;

  bool m_being_destroyed;
};

inline
Bucket *BucketRepository::get_bucket(EntityRank entity_rank, int bucket_id) const
{
  const BucketVector & all_buckets_for_rank = m_buckets[entity_rank];
  ThrowAssert(static_cast<size_t>(bucket_id) < all_buckets_for_rank.size());
  return all_buckets_for_rank[bucket_id];
}

#undef RANK_DEPENDENT_GET_BUCKET_FN_DEF

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_BucketRepository_hpp
