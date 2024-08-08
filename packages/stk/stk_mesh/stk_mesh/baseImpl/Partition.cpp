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

#include <stk_mesh/baseImpl/Partition.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_topology/topology.hpp>    // for topology, operator<<, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartOrdinal, etc
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/GlobalIdEntitySorter.hpp>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc

namespace stk { namespace mesh { class FieldBase; } }


// Testing this!
#define PARTITION_HOLD_EMPTY_BUCKET_OPTIMIZATION

using namespace stk::mesh::impl;

Partition::Partition(BulkData& mesh, BucketRepository *repo, EntityRank rank,
                     const PartOrdinal* keyBegin, const PartOrdinal* keyEnd)
  : m_mesh(mesh),
    m_repository(repo)
  , m_rank(rank)
  , m_extPartitionKey(keyBegin, keyEnd)
  , m_size(0)
  , m_updated_since_sort(false)
  , m_removeMode(FILL_HOLE_THEN_SORT)
  , m_removedEntities()
{
}

Partition::~Partition()
{
  size_t num_bkts = m_buckets.size();
  for (size_t i = 0; i < num_bkts; ++i)
  {
    delete m_buckets[i];
  }
}

void Partition::remove_internal(Bucket& bucket, unsigned bucketOrd)
{
  overwrite_from_end(bucket, bucketOrd);
  remove_impl();
}

bool Partition::remove(Entity entity)
{
  STK_ThrowAssert(belongs(m_mesh.bucket(entity)));

  Bucket &bucket   = m_mesh.bucket(entity);
  unsigned ordinal = m_mesh.bucket_ordinal(entity);

  if (m_removeMode == TRACK_THEN_SLIDE) {
    bool foundBucket = false;
    for(unsigned myBktIdx=0; myBktIdx<m_buckets.size(); ++myBktIdx) {
      if (m_buckets[myBktIdx]->bucket_id() == bucket.bucket_id()) {
        const bool lastBucketLastEntity =
          (myBktIdx==m_buckets.size()-1) && (ordinal==m_buckets[myBktIdx]->size()-1);
        if (lastBucketLastEntity) {
          remove_impl();
        }
        else {
          m_removedEntities.push_back(FastMeshIndex{myBktIdx, ordinal});
          m_buckets[myBktIdx]->m_entities[ordinal] = Entity(); //Partition is friend of Bucket!
          m_updated_since_sort = true;
        }
        foundBucket = true;
        break;
      }
    }
    STK_ThrowRequireMsg(foundBucket, "Failed to find bucket in partition for entity that is being removed.");
  }
  else {
    remove_internal(bucket, ordinal);
  }
  return true;
}

bool Partition::add(Entity entity)
{
  if (m_mesh.bucket_ptr(entity))
  {
    // If an entity already belongs to a partition, it cannot be added to one.
    return false;
  }

  // If the last bucket is full, automatically create a new one.
  Bucket *bucket = get_bucket_for_adds();
  bucket->add_entity(entity);
  ++m_size;

  m_removeMode = FILL_HOLE_THEN_SORT;
  m_updated_since_sort = true;

  internal_check_invariants();

  return true;
}

bool Partition::move_to(Entity entity, Partition &dst_partition)
{
  STK_ThrowAssert(belongs(m_mesh.bucket(entity)));

  Bucket *src_bucket   = m_mesh.bucket_ptr(entity);
  unsigned src_ordinal = m_mesh.bucket_ordinal(entity);
  if (src_bucket && (src_bucket->getPartition() == &dst_partition))
  {
    return false;
  }

  STK_ThrowRequireMsg(src_bucket && (src_bucket->getPartition() == this),
                  "Partition::move_to cannot move an entity that does not belong to it.");

  // If the last bucket is full, automatically create a new one.
  Bucket *dst_bucket = nullptr;
  try {
      dst_bucket = dst_partition.get_bucket_for_adds();
  }
  catch(std::bad_alloc& ex) {
      std::ostringstream strout;
      strout << "Error adding "<<m_mesh.entity_key(entity)<< " to mesh:"
             << "ERROR: The current simulation ran out of memory." << std::endl
             << "Consider requesting less processors per node and resubmit your job." << std::endl;
      throw std::bad_alloc();
    }
  catch(std::exception& e) {
      std::ostringstream os;
      os << "Error adding "<<m_mesh.entity_key(entity)<< " to mesh: " << e.what();
      throw std::runtime_error(os.str());
  }

  STK_ThrowErrorMsgIf(src_bucket && src_bucket->topology().is_valid() && (src_bucket->topology() != dst_bucket->topology()),
                  "Error: cannot change topology of entity (rank: "
                  << static_cast<stk::topology::rank_t>(m_mesh.entity_rank(entity))
                  << ", global_id: " << m_mesh.identifier(entity) << ") from "
                  << src_bucket->topology() << "to " << dst_bucket->topology() << "."
                  );

  // Copy the entity's data to the new bucket before removing the entity from its old bucket.
  dst_bucket->copy_entity(entity);

  if (m_removeMode == TRACK_THEN_SLIDE) {
    bool foundBucket = false;
    for(unsigned myBktIdx=0; myBktIdx<m_buckets.size(); ++myBktIdx) {
      if (m_buckets[myBktIdx]->bucket_id() == src_bucket->bucket_id()) {
        const bool lastBucketLastEntity =
          (myBktIdx==m_buckets.size()-1) && (src_ordinal==m_buckets[myBktIdx]->size()-1);
        if (lastBucketLastEntity) {
          remove_impl();
        }
        else {
          m_removedEntities.push_back(FastMeshIndex{myBktIdx, src_ordinal});
          m_buckets[myBktIdx]->m_entities[src_ordinal] = Entity(); //Partition is friend of Bucket!
          m_updated_since_sort = true;
        }
        foundBucket = true;
        break;
      }
    } 
    STK_ThrowRequireMsg(foundBucket, "Failed to find bucket in partition for entity that is being removed.");
  }
  else {
    overwrite_from_end(*src_bucket, src_ordinal);
    remove_impl();
  }

  dst_partition.m_updated_since_sort = true;
  dst_partition.m_size++;
  dst_partition.internal_check_invariants();

  return true;
}

void Partition::overwrite_from_end(Bucket& bucket, unsigned ordinal)
{
  Bucket *last = *(end() - 1);

  const bool NOT_last_entity_in_last_bucket =
    (last != &bucket) || (bucket.size() != ordinal + 1);
  if ( NOT_last_entity_in_last_bucket )
  {
    // Copy last entity to spot being vacated.
    Entity e_swap = (*last)[ last->size() - 1 ];
    bucket.overwrite_entity(ordinal, e_swap );
    m_updated_since_sort = true;
  }
}

void Partition::delete_bucket(Bucket* bucket)
{
  m_size -= bucket->size();

  auto iter = std::find(m_buckets.begin(), m_buckets.end(), bucket);
  m_repository->deallocate_bucket(bucket);
  m_buckets.erase(iter, iter+1);
}

void Partition::remove_bucket(Bucket* bucket)
{
  m_size -= bucket->size();

  auto iter = std::find(m_buckets.begin(), m_buckets.end(), bucket);
  m_buckets.erase(iter, iter+1);
}

void Partition::add_bucket(Bucket* bucket)
{
  clear_pending_removes_by_filling_from_end();
  m_size += bucket->size();
  bucket->m_partition = this;
  stk::util::insert_keep_sorted_and_unique(bucket, m_buckets);
}

void Partition::remove_impl()
{
  STK_ThrowAssert(!empty());

  Bucket &last_bucket   = **(end() - 1);

  last_bucket.remove_entity();

  if ( 0 == last_bucket.size() )
  {
    size_t num_buckets = m_buckets.size();

    // Don't delete the last bucket now --- might want it later in this modification cycle.
    if (num_buckets > 1)
    {
      m_repository->deallocate_bucket( m_buckets.back() );
      m_buckets.pop_back();
    }
    else
    {
#ifndef PARTITION_HOLD_EMPTY_BUCKET_OPTIMIZATION
      m_repository->deallocate_bucket(m_buckets.back()());
      m_buckets.pop_back();
#else
      m_repository->m_need_sync_from_partitions[m_rank] = true;
#endif
    }
  }

  --m_size;

  internal_check_invariants();
}

void Partition::reset_partition_key(const std::vector<unsigned>& newKey)
{
  m_extPartitionKey = newKey;
}

void Partition::default_sort_if_needed(bool mustSortFacesByNodeIds)
{
  if (!empty() && m_updated_since_sort)
  {
    if (m_removeMode == FILL_HOLE_THEN_SORT) {
      sort(GlobalIdEntitySorter(mustSortFacesByNodeIds));
    }
    else {
      finalize_pending_removes_by_sliding_memory();
    }
  }
}

stk::mesh::FieldVector get_fields_for_bucket(const stk::mesh::BulkData& mesh,
                                             const stk::mesh::Bucket& bkt)
{
  stk::mesh::FieldVector fields;

  if(mesh.is_field_updating_active()) {
    if(bkt.entity_rank() < stk::topology::NUM_RANKS) {
      const stk::mesh::FieldVector& rankFields = mesh.mesh_meta_data().get_fields(bkt.entity_rank());
      fields.reserve(rankFields.size());

      for(unsigned ifield=0; ifield<rankFields.size(); ++ifield) {
        if(field_bytes_per_entity(*rankFields[ifield], bkt)) {
          fields.push_back(rankFields[ifield]);
        }
      }
    }
  }

  return fields;
}

void Partition::sort(const EntitySorterBase& sorter)
{
  std::vector<unsigned> partition_key = get_legacy_partition_id();
  //index of bucket in partition
  partition_key[ partition_key[0] ] = 0;

  std::vector<Entity> entities(m_size);

  BucketVector::iterator buckets_begin, buckets_end;
  buckets_begin = begin();
  buckets_end = end();

  // Copy all the entities in the Partition into a vector for sorting.
  size_t new_i = 0;
  for (BucketVector::iterator b_i = buckets_begin; b_i != buckets_end; ++b_i)
  {
    Bucket &b = **b_i;
    size_t b_size = b.size();
    std::copy(&b.m_entities[0], &b.m_entities[0] + b_size, entities.data() + new_i);
    new_i += b_size;
  }

  sorter.sort(m_mesh, entities);

  // Make sure that there is a vacancy somewhere.
  //
  stk::mesh::Bucket *vacancy_bucket = *(buckets_end - 1);
  unsigned vacancy_ordinal = vacancy_bucket->size();
  stk::mesh::Bucket *tmp_bucket = 0;

  if (vacancy_ordinal >= vacancy_bucket->capacity())
  {
    // If we need a temporary bucket, it only needs to hold one entity and
    // the corresponding field data.
    tmp_bucket = m_repository->allocate_bucket(m_rank, partition_key, 1, 1);
    vacancy_bucket = tmp_bucket;
    vacancy_ordinal = 0;
  }

  // Allocate space to copy in to
  vacancy_bucket->add_entity();

  // Now that we have the entities sorted, we need to put them and their data
  // in the right order in the buckets.

  std::vector<Entity>::iterator sorted_ent_vector_itr = entities.begin();

  Bucket* orig_vacancy_bucket = vacancy_bucket;


  FieldVector reduced_fields = get_fields_for_bucket(m_mesh, *vacancy_bucket);

  for (BucketVector::iterator bucket_itr = begin(); bucket_itr != buckets_end; ++bucket_itr)
  {
      Bucket &curr_bucket = **bucket_itr;
      const unsigned n = *bucket_itr == orig_vacancy_bucket ? curr_bucket.size() -1 : curr_bucket.size(); // skip very last entity in partition

      for ( unsigned curr_bucket_ord = 0; curr_bucket_ord < n ; ++curr_bucket_ord , ++sorted_ent_vector_itr ) {
          STK_ThrowAssert(sorted_ent_vector_itr != entities.end());

          Entity curr_entity = curr_bucket[curr_bucket_ord];
          STK_ThrowAssert(m_mesh.is_valid(curr_entity));

          if ( curr_entity != *sorted_ent_vector_itr ) // check if we need to move
          {
              // Move current entity to the vacant spot
              if (vacancy_bucket != &curr_bucket || vacancy_ordinal != curr_bucket_ord) {
                  vacancy_bucket->overwrite_entity( vacancy_ordinal, curr_entity, &reduced_fields );
              }

              // Set the vacant spot to where the required entity is now.
              MeshIndex& meshIndex = m_mesh.mesh_index(*sorted_ent_vector_itr);
              vacancy_bucket  = meshIndex.bucket;
              vacancy_ordinal = meshIndex.bucket_ordinal;

              // Move required entity to the required spot
              curr_bucket.overwrite_entity( curr_bucket_ord, *sorted_ent_vector_itr, &reduced_fields );
          }
      }
  }

  m_updated_since_sort = false;

  orig_vacancy_bucket->remove_entity();

  if (tmp_bucket) {
    m_repository->deallocate_bucket(tmp_bucket);
  }

  if (not m_buckets.empty()) {
    m_buckets.back()->reset_empty_space(reduced_fields);
  }

  internal_check_invariants();
}

bool out_of_range(const stk::mesh::FastMeshIndex& fastMeshIndex, const stk::mesh::BucketVector& buckets)
{
  return fastMeshIndex.bucket_id >= buckets.size() ||
         fastMeshIndex.bucket_ord >= buckets[fastMeshIndex.bucket_id]->size();
}

stk::mesh::Entity get_last_entity(const stk::mesh::BucketVector& buckets)
{
  const stk::mesh::Bucket* lastBucket = buckets.back();
  const size_t size = lastBucket->size();
  if (size == 0) return stk::mesh::Entity();
  return (*lastBucket)[lastBucket->size()-1];
}

void Partition::clear_pending_removes_by_filling_from_end()
{
  if (m_removeMode == TRACK_THEN_SLIDE && !m_removedEntities.empty()) {
    stk::util::sort_and_unique(m_removedEntities);
    for(const FastMeshIndex& fastMeshIndex : m_removedEntities) {
      while(!empty() && !m_mesh.is_valid(get_last_entity(m_buckets))) {
        remove_impl();
      }

      if (empty()) { break; }

      if (out_of_range(fastMeshIndex, m_buckets)) { continue; }

      Bucket& bucket = *m_buckets[fastMeshIndex.bucket_id];
      unsigned bucketOrd = fastMeshIndex.bucket_ord;
      remove_internal(bucket, bucketOrd);
    }

    m_removedEntities.clear();
  }
  m_removeMode = FILL_HOLE_THEN_SORT;
}

void Partition::check_sorted(const std::string& prefixMsg)
{
  EntityLess entityLess(m_mesh);
  Entity prevEnt = Entity();
  unsigned prevBktId = 0;
  unsigned prevBktOffset = 0;
  for(const Bucket* bptr : m_buckets) {
    unsigned offset = 0;
    for(Entity ent : *bptr) {
      if(m_mesh.is_valid(ent)){
      if (m_mesh.is_valid(prevEnt)) {
        STK_ThrowRequireMsg(entityLess(prevEnt, ent), "check_sorted "<<prefixMsg<<", m_buckets.size()="<<m_buckets.size()<<", found ent ID "<<m_mesh.identifier(ent)<<" at bktId "<<bptr->bucket_id()<<" offset "<<offset<<", after prevEnt ID "<<m_mesh.identifier(prevEnt)<<" at bktId "<<prevBktId<<" offset "<<prevBktOffset);
      }
      prevEnt = ent;
      prevBktId = bptr->bucket_id();
      prevBktOffset = offset++;
      }
    }
  }
}

void Partition::finalize_pending_removes_by_sliding_memory()
{
  if (m_removeMode == TRACK_THEN_SLIDE && !m_removedEntities.empty()) {
    stk::util::sort_and_unique(m_removedEntities);

    std::vector<FastMeshIndex>::iterator rmEntIter = m_removedEntities.begin();
    FastMeshIndex slotToFill = *rmEntIter;
    unsigned bktIdx = slotToFill.bucket_id;
    const bool lastBucket = bktIdx == m_buckets.size()-1;
    unsigned bktOrd = slotToFill.bucket_ord+1 ;
    bool lastEntity = bktOrd >= m_buckets[bktIdx]->size();
    if (!lastBucket && lastEntity) {
      ++bktIdx;
      bktOrd = 0;
      lastEntity = bktOrd >= m_buckets[bktIdx]->size();
    }

    FieldVector reduced_fields = get_fields_for_bucket(m_mesh, *m_buckets[bktIdx]);

    while(bktIdx < m_buckets.size()) {
      Bucket& bkt = *m_buckets[bktIdx];

      while(bktOrd < bkt.size()) {
        Entity entity = bkt[bktOrd];

        if (m_mesh.is_valid(entity)) {
          m_buckets[slotToFill.bucket_id]->overwrite_entity(slotToFill.bucket_ord, entity, &reduced_fields);

          ++slotToFill.bucket_ord;
          if (slotToFill.bucket_ord >= m_buckets[slotToFill.bucket_id]->size()) {
            ++slotToFill.bucket_id;
            slotToFill.bucket_ord = 0;
          }
        }
        ++bktOrd;
      }
      ++bktIdx;
      bktOrd = 0;
    }

    unsigned numRemoved = m_removedEntities.size();
    for(unsigned i=0; i<numRemoved; ++i) {
      remove_impl();
    }
    m_removedEntities.clear();
  }

  m_updated_since_sort = false;
}

void Partition::set_remove_mode(RemoveMode removeMode)
{
  if(m_removeMode!=removeMode && m_updated_since_sort) {
    default_sort_if_needed();
  }

  m_removeMode = removeMode;
}

stk::mesh::Bucket *Partition::get_bucket_for_adds()
{
  clear_pending_removes_by_filling_from_end();

  if (no_buckets()) {
    std::vector<unsigned> partition_key = get_legacy_partition_id();
    partition_key[ partition_key[0] ] = 0;
    Bucket *bucket = m_repository->allocate_bucket(m_rank, partition_key,
                                                   m_repository->get_initial_bucket_capacity(),
                                                   m_repository->get_maximum_bucket_capacity());
    bucket->m_partition = this;
    m_buckets.push_back(bucket);

    return bucket;
  }

  Bucket *bucket = *(end() - 1);  // Last bucket of the partition.

  if (bucket->size() == bucket->capacity()) {
    if (bucket->size() == m_repository->get_maximum_bucket_capacity()) {
      std::vector<unsigned> partition_key = get_legacy_partition_id();
      partition_key[ partition_key[0] ] = m_buckets.size();
      bucket = m_repository->allocate_bucket(m_rank, partition_key,
                                             m_repository->get_initial_bucket_capacity(),
                                             m_repository->get_maximum_bucket_capacity());
      bucket->m_partition = this;
      m_buckets.push_back(bucket);
    }
    else {
      bucket->grow_capacity();
    }
  }

  return bucket;
}

