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

#include <stk_mesh/baseImpl/BucketRepository.hpp>
#include <algorithm>                    // for copy, remove_if
#include <new>                          // for operator new
#include <sstream>                      // for operator<<, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, raw_part_equal
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition, lower_bound
#include <stk_mesh/baseImpl/ForEachEntityLoopAbstractions.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/SortAndUnique.hpp"

namespace stk {
namespace mesh {
namespace impl {



BucketRepository::BucketRepository(BulkData & mesh,
                                   unsigned entity_rank_count,
                                   unsigned initialBucketCapacity,
                                   unsigned maximumBucketCapacity)
  : m_mesh(mesh),
    m_buckets(entity_rank_count),
    m_partitions(entity_rank_count),
    m_need_sync_from_partitions(entity_rank_count, false),
    m_initialBucketCapacity(initialBucketCapacity),
    m_maximumBucketCapacity(maximumBucketCapacity),
    m_being_destroyed(false)
{
  STK_ThrowRequireMsg(initialBucketCapacity > 0,
                      "The initial Bucket capacity (" << initialBucketCapacity << ") must be positive");
  STK_ThrowRequireMsg(initialBucketCapacity <= maximumBucketCapacity,
                      "The initial Bucket capacity (" << initialBucketCapacity
                      << ") cannot be larger than the maximum Bucket capacity (" << maximumBucketCapacity << ")");
}

BucketRepository::~BucketRepository()
{
  // Destroy buckets, which were *not* allocated by the set.

  m_being_destroyed = true;

  try {

    for ( std::vector<std::vector<Partition *> >::iterator pv_i = m_partitions.begin();
          pv_i != m_partitions.end(); ++pv_i)
    {
      for (std::vector<Partition *>::iterator p_j = pv_i->begin();
           p_j != pv_i->end(); ++p_j)
      {
        Partition * tmp = *p_j;
        delete tmp;
      }
      pv_i->clear();
    }
    m_partitions.clear();
    m_buckets.clear();

    const FieldVector& fields = m_mesh.mesh_meta_data().get_fields();
    for(size_t i=0; i<fields.size(); ++i) {
      fields[i]->get_meta_data_for_field().clear();
    }

  } catch(...) {}
}

void BucketRepository::set_needs_to_be_sorted(stk::mesh::Bucket &bucket, bool needsSorting)
{
    bucket.getPartition()->set_flag_needs_to_be_sorted(needsSorting);
}

void BucketRepository::internal_default_sort_bucket_entities(bool mustSortFacesByNodeIds)
{
    for(std::vector<Partition*>& partitionVector : m_partitions)
        for(Partition* partition : partitionVector)
            partition->default_sort_if_needed(mustSortFacesByNodeIds);
}

void BucketRepository::internal_custom_sort_bucket_entities(const EntitySorterBase& sorter)
{
    for(std::vector<Partition*>& partitionVector : m_partitions)
        for(Partition* partition : partitionVector)
            partition->sort(sorter);
}

void BucketRepository::add_entity_with_part_memberships(const stk::mesh::Entity entity,
                                                        const EntityRank arg_entity_rank,
                                                        const OrdinalVector &parts)
{
    Partition *partition = get_or_create_partition(arg_entity_rank, parts);
    partition->add(entity);
}

void BucketRepository::change_entity_part_membership(const MeshIndex &meshIndex, const OrdinalVector &parts)
{
    Bucket *bucket = meshIndex.bucket;
    Partition *destinationPartition = get_or_create_partition(bucket->entity_rank(), parts);
    Entity entity = get_entity(meshIndex);
    Partition *sourcePartition = bucket->getPartition();
    sourcePartition->move_to(entity, *destinationPartition);
}

void BucketRepository::remove_entity(const MeshIndex &meshIndex)
{
    Bucket *bucket = meshIndex.bucket;
    Partition *partition = bucket->getPartition();
    Entity entity = get_entity(meshIndex);
    partition->remove(entity);
}

void BucketRepository::ensure_data_structures_sized()
{
    if(m_buckets.empty())
    {
        size_t entity_rank_count = m_mesh.mesh_meta_data().entity_rank_count();
        STK_ThrowRequireMsg(entity_rank_count > 0, "MetaData doesn't have any entity-ranks! Did you forget to "
                            "initialize MetaData before creating BulkData?");
        m_buckets.resize(entity_rank_count);
        m_partitions.resize(entity_rank_count);
        m_need_sync_from_partitions.resize(entity_rank_count, false);
    }
}

////
//// Note that we need to construct a key vector that the particular
//// format so we can use the lower_bound(..) function to lookup the
//// partition.  Because we are using partitions now instead of
//// buckets, it should be possible to do without that vector and
//// instead do the lookup directly from the OrdinalVector.
////

Partition *BucketRepository::get_or_create_partition(
  const EntityRank arg_entity_rank ,
  const OrdinalVector &parts)
{
  const unsigned maxKeyTmpBufferSize = 64;
  PartOrdinal keyTmpBuffer[maxKeyTmpBufferSize];
  OrdinalVector keyTmpVec;

  PartOrdinal* keyPtr = nullptr;
  PartOrdinal* keyEnd = nullptr;

  fill_key_ptr(parts, &keyPtr, &keyEnd, maxKeyTmpBufferSize, keyTmpBuffer, keyTmpVec);

  std::vector<Partition *>::iterator ik;

  Partition* partition = get_partition(arg_entity_rank, parts, ik, keyPtr, keyEnd);

  if(partition == nullptr) {
    partition = create_partition(arg_entity_rank, parts, ik, keyPtr, keyEnd);
  }
  return partition;
}

void BucketRepository::fill_key_ptr(const OrdinalVector& parts, PartOrdinal** keyPtr, PartOrdinal** keyEnd,
                                    const unsigned maxKeyTmpBufferSize, PartOrdinal* keyTmpBuffer, OrdinalVector& keyTmpVec)
{
  const size_t part_count = parts.size();

  const size_t keyLen = 2 + part_count;

  *keyPtr = keyTmpBuffer;
  *keyEnd = *keyPtr+keyLen;

  if (keyLen >= maxKeyTmpBufferSize) {
    keyTmpVec.resize(keyLen);
    *keyPtr = keyTmpVec.data();
    *keyEnd = *keyPtr+keyLen;
  }

  //----------------------------------
  // Key layout:
  // { part_count + 1 , { part_ordinals } , partition_count }
  // Thus partition_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key for a bucket in the partition.
  const unsigned max = static_cast<unsigned>(-1);
  (*keyPtr)[0] = part_count+1;
  (*keyPtr)[ (*keyPtr)[0] ] = max ;

  {
    for ( unsigned i = 0 ; i < part_count ; ++i ) { (*keyPtr)[i+1] = parts[i] ; }
  }
}

Partition *BucketRepository::get_partition(const EntityRank arg_entity_rank, const OrdinalVector &parts)
{
  PartOrdinal* keyPtr = nullptr;
  PartOrdinal* keyEnd = nullptr;
  std::vector<impl::Partition*>::iterator ik;

  const unsigned maxKeyTmpBufferSize = 64;
  PartOrdinal keyTmpBuffer[maxKeyTmpBufferSize];
  OrdinalVector keyTmpVec;

  fill_key_ptr(parts, &keyPtr, &keyEnd, maxKeyTmpBufferSize, keyTmpBuffer, keyTmpVec);

  return get_partition(arg_entity_rank, parts, ik, keyPtr, keyEnd);
}

Partition *BucketRepository::get_partition(
  const EntityRank arg_entity_rank ,
  const OrdinalVector &parts,
  std::vector<Partition*>::iterator& ik,
  PartOrdinal* keyPtr,
  PartOrdinal* keyEnd)
{
  STK_ThrowRequireMsg(m_mesh.mesh_meta_data().check_rank(arg_entity_rank), "Entity rank " << arg_entity_rank
                      << " is invalid");

  ensure_data_structures_sized();

  std::vector<Partition *> & partitions = m_partitions[ arg_entity_rank ];

  // If the partition is found, the iterator will be right after it, thanks to the
  // trickiness above.
  ik = lower_bound( partitions , keyPtr );
  const bool partition_exists =
    (ik != partitions.begin()) && raw_part_equal( ik[-1]->key() , keyPtr );

  if (partition_exists)
  {
    return ik[-1];
  }

  return nullptr;
}

Partition* BucketRepository::create_partition(
  const EntityRank arg_entity_rank,
  const OrdinalVector& parts,
  std::vector<Partition*>::iterator& ik,
  PartOrdinal* keyPtr,
  PartOrdinal* keyEnd)
{
  keyPtr[keyPtr[0]] = 0;

  Partition *partition = new Partition(m_mesh, this, arg_entity_rank, keyPtr, keyEnd);
  STK_ThrowRequire(partition != nullptr);

  m_need_sync_from_partitions[arg_entity_rank] = true;
  m_partitions[arg_entity_rank].insert( ik , partition );

  return partition;
}

void BucketRepository::internal_modification_end()
{
    sync_from_partitions();

    for(EntityRank from_rank = stk::topology::NODE_RANK; from_rank < stk::topology::NUM_RANKS; ++from_rank)
    {
        const BucketVector &buckets = this->buckets(from_rank);
        unsigned num_buckets = buckets.size();
        for(unsigned j = 0; j < num_buckets; ++j)
        {
            STK_ThrowAssert(buckets[j] != nullptr);
            Bucket &bucket = *buckets[j];

            for(EntityRank to_rank = stk::topology::NODE_RANK; to_rank < stk::topology::NUM_RANKS; ++to_rank)
            {
                if (from_rank == to_rank) {
                    continue;
                }

                switch(to_rank)
                {
                    case stk::topology::NODE_RANK:
                        bucket.m_fixed_node_connectivity.end_modification(&bucket.m_mesh);
                        break;
                    case stk::topology::EDGE_RANK:
                        bucket.m_dynamic_edge_connectivity.compress_connectivity();
                        break;
                    case stk::topology::FACE_RANK:
                        bucket.m_dynamic_face_connectivity.compress_connectivity();
                        break;
                    case stk::topology::ELEMENT_RANK:
                        bucket.m_dynamic_element_connectivity.compress_connectivity();
                        break;
                    case stk::topology::INVALID_RANK:
                        break;
                    default:
                        bucket.m_dynamic_other_connectivity.compress_connectivity();
                        break;
                }
            }
        }
    }
}

void BucketRepository::sync_from_partitions()
{
  for (EntityRank rank = stk::topology::NODE_RANK; rank < static_cast<EntityRank>(m_partitions.size()); ++rank)
  {
    sync_from_partitions(rank);
  }
}

namespace {

inline bool is_null(stk::mesh::impl::Partition *p) { return (p ? false : true);}

struct bucket_less_by_first_entity_identifier
{
    bool operator()(const Bucket* first, const Bucket* second) const
    {
        bool result = false;
        if ((first->size() > 0) && (second->size() > 0))
        {
            const stk::mesh::BulkData& mesh = first->mesh();
            result = EntityLess(mesh)((*first)[0], (*second)[0]);
        }
        return (first->size() == 0) || result;
    }
};

}

void BucketRepository::sync_from_partitions(EntityRank rank)
{
  if (m_need_sync_from_partitions[rank])
  {
      std::vector<Partition *> &partitions = m_partitions[rank];

      size_t num_partitions = partitions.size();
      size_t num_buckets = 0;
      for (size_t p_i = 0; p_i < num_partitions; ++p_i)
      {
        if (!partitions[p_i]->empty())
        {
          num_buckets += partitions[p_i]->num_buckets();
        }
      }

      m_buckets[rank].resize(num_buckets);

      bool has_hole = false;
      BucketVector::iterator bkts_i = m_buckets[rank].begin();
      for (size_t p_i = 0; p_i < num_partitions; ++p_i)
      {
        Partition &partition = *partitions[p_i];

        if (partition.empty())
        {
          delete partitions[p_i];
          partitions[p_i] = 0;
          has_hole = true;
          continue;
        }
        size_t num_bkts_in_partition = partition.num_buckets();
        std::copy(partition.begin(), partition.end(), bkts_i);
        bkts_i += num_bkts_in_partition;
      }

      if (has_hole)
      {
        std::vector<Partition *>::iterator new_end;
        new_end = std::remove_if(partitions.begin(), partitions.end(), is_null);
        size_t new_size = new_end - partitions.begin();  // OK because has_hole is true.
        partitions.resize(new_size);
      }
  }

  if(m_mesh.should_sort_buckets_by_first_entity_identifier())
  {
      std::sort(m_buckets[rank].begin(), m_buckets[rank].end(), bucket_less_by_first_entity_identifier());
  }

  if (m_need_sync_from_partitions[rank] == true || m_mesh.should_sort_buckets_by_first_entity_identifier())
  {
      sync_bucket_ids(rank);
  }

  m_need_sync_from_partitions[rank] = false;
}

Bucket *BucketRepository::allocate_bucket(EntityRank entityRank,
                                          const std::vector<unsigned> & key,
                                          unsigned initialCapacity,
                                          unsigned maximumCapacity)
{
  STK_ThrowAssertMsg(stk::util::is_sorted_and_unique(std::vector<unsigned>(key.begin()+1,key.end()-1),std::less<unsigned>()),
                     "bucket created with 'key' vector that's not sorted and unique");
  BucketVector &bucket_vec = m_buckets[entityRank];
  const unsigned bucket_id = bucket_vec.size();
  Bucket * new_bucket = new Bucket(m_mesh, entityRank, key, initialCapacity, maximumCapacity, bucket_id);
  STK_ThrowRequire(new_bucket != nullptr);

  bucket_vec.push_back(new_bucket);
  m_need_sync_from_partitions[entityRank] = true;

  return new_bucket;
}

void BucketRepository::deallocate_bucket(Bucket *b)
{
  STK_ThrowAssertMsg(b != nullptr, "BucketRepository::deallocate_bucket(.) m_buckets invariant broken.");

  const unsigned bucket_id = b->bucket_id();
  const EntityRank bucket_rank = b->entity_rank();

  STK_ThrowAssertMsg(b == m_buckets[bucket_rank][bucket_id],
                     "BucketRepository::deallocate_bucket(.) m_buckets invariant broken.");

  m_buckets[bucket_rank][bucket_id] = nullptr; // space will be reclaimed by sync_from_partitions
  m_need_sync_from_partitions[bucket_rank] = true;
  delete b;
}

void BucketRepository::sync_bucket_ids(EntityRank entity_rank)
{
  BucketVector &buckets = m_buckets[entity_rank];
  unsigned num_buckets = buckets.size();
  std::vector<unsigned> id_map(num_buckets);

  for (unsigned i = 0; i < num_buckets; ++i)
  {
    STK_ThrowAssertMsg(buckets[i] != nullptr, "BucketRepository::sync_bucket_ids() called when m_buckets["
                       << entity_rank << "] is not dense.");
    id_map[i] = buckets[i]->bucket_id();
    buckets[i]->m_bucket_id = i;
  }

  m_mesh.reorder_buckets_callback(entity_rank, id_map);
}

std::vector<Partition *> BucketRepository::get_partitions(EntityRank rank) const
{
  std::vector<Partition *> retval;
  std::vector<Partition *> const& bf_vec = m_partitions[rank];
  for (size_t i = 0; i < bf_vec.size(); ++i)
  {
    retval.push_back(bf_vec[i]);
  }
  return retval;
}

void BucketRepository::delete_bucket(Bucket * bucket)
{
    bucket->getPartition()->delete_bucket(bucket);
}

void BucketRepository::set_remove_mode_tracking()
{
  for(std::vector<Partition*>& partitionVector : m_partitions) {
    for(Partition* partition : partitionVector) {
      partition->set_remove_mode(Partition::TRACK_THEN_SLIDE);
    }
  }
}

void BucketRepository::set_remove_mode_fill_and_sort()
{
  for(std::vector<Partition*>& partitionVector : m_partitions) {
    for(Partition* partition : partitionVector) {
      partition->set_remove_mode(Partition::FILL_HOLE_THEN_SORT);
    }
  }
}

} // namespace impl
} // namespace mesh
} // namespace stk
