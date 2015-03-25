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

#include <stk_mesh/baseImpl/Partition.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for EntityLess, BulkData
#include <stk_topology/topology.hpp>    // for topology, operator<<, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartOrdinal, etc
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssert, etc
#include "stk_util/util/TrackingAllocator.hpp"  // for tracking_allocator
namespace stk { namespace mesh { class FieldBase; } }


// Testing this!
#define PARTITION_HOLD_EMPTY_BUCKET_OPTIMIZATION

namespace stk {
namespace mesh {
namespace impl {

std::ostream &operator<<(std::ostream &os, const stk::mesh::impl::Partition &bf)
{
  return bf.streamit(os);
}

} // impl
} // mesh
} // stk

using namespace stk::mesh::impl;

std::ostream &Partition::streamit(std::ostream &os) const
{
  const MetaData & mesh_meta_data = m_mesh.mesh_meta_data();

  os << "{Partition m_rank = " << m_rank << ", m_size = " << m_size;

  os << "  legacy partition id : {";
  const std::vector<unsigned> &family_key = get_legacy_partition_id();
  for (size_t i = 0; i < family_key.size(); ++i)
  {
    os << " " << family_key[i];
    if ((i > 0) && (i < family_key.size() - 1))
    {
      const Part & part = mesh_meta_data.get_part( family_key[i] );
      os << " " << part.name();
    }
  }
  os << " }}";

  return os;
}

std::ostream &Partition::dumpit(std::ostream &os) const
{
  os << "{ Partition (rank = " << m_rank << ")  \n";
  for (BucketVector::const_iterator b_i = begin(); b_i != end(); ++b_i)
  {
    Bucket &b = **b_i;
    print(os, "  ", b );
  }
  os << "}\n";
  return os;
}

std::string Partition::dumpit() const
{
  std::ostringstream output;
  dumpit(output);

  return output.str();
}

Partition::Partition(BulkData& mesh, BucketRepository *repo, EntityRank rank,
                     const std::vector<PartOrdinal> &key)
  : m_mesh(mesh),
    m_repository(repo)
  , m_rank(rank)
  , m_extPartitionKey(key)
  , m_size(0)
  , m_updated_since_compress(false)
  , m_updated_since_sort(false)
{
  // Nada.
}

// Only the BucketRepository will delete a Partition.
Partition::~Partition()
{
  typedef tracking_allocator<Bucket, BucketTag> bucket_allocator;
  size_t num_bkts = m_buckets.size();
  for (size_t i = 0; i < num_bkts; ++i)
  {
    m_buckets[i]->~Bucket();
    bucket_allocator().deallocate(m_buckets[i],1);
  }
}

BucketRepository &Partition::getRepository(stk::mesh::BulkData &mesh)
{
  return mesh.m_bucket_repository;
}

bool Partition::remove(Entity entity)
{
  ThrowAssert(belongs(m_mesh.bucket(entity)));

  Bucket &bucket   = m_mesh.bucket(entity);
  unsigned ordinal = m_mesh.bucket_ordinal(entity);
  overwrite_from_end(bucket, ordinal);

  m_mesh.set_mesh_index(entity, 0, 0);

  remove_impl();
  internal_check_invariants();
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
  bucket->mesh().mark_entity_and_upward_related_entities_as_modified(entity);
  ++m_size;

  m_updated_since_compress = m_updated_since_sort = true;

  internal_check_invariants();

  return true;
}

bool Partition::move_to(Entity entity, Partition &dst_partition)
{
  ThrowAssert(belongs(m_mesh.bucket(entity)));

  Bucket *src_bucket   = m_mesh.bucket_ptr(entity);
  unsigned src_ordinal = m_mesh.bucket_ordinal(entity);
  if (src_bucket && (src_bucket->getPartition() == &dst_partition))
  {
    return false;
  }

  ThrowRequireMsg(src_bucket && (src_bucket->getPartition() == this),
                  "Partition::move_to cannot move an entity that does not belong to it.");

  // If the last bucket is full, automatically create a new one.
  Bucket *dst_bucket = dst_partition.get_bucket_for_adds();

  ThrowErrorMsgIf(src_bucket && src_bucket->topology().is_valid() && (src_bucket->topology() != dst_bucket->topology()),
                  "Error: cannot change topology of entity (rank: "
                  << static_cast<stk::topology::rank_t>(m_mesh.entity_rank(entity))
                  << ", global_id: " << m_mesh.identifier(entity) << ") from "
                  << src_bucket->topology() << "to " << dst_bucket->topology() << "."
                  );

  dst_bucket->mesh().mark_entity_and_upward_related_entities_as_modified(entity);

  // Copy the entity's data to the new bucket before removing the entity from its old bucket.
  dst_bucket->copy_entity(entity);

  overwrite_from_end(*src_bucket, src_ordinal);

  dst_partition.m_updated_since_compress = dst_partition.m_updated_since_sort = true;
  dst_partition.m_size++;

  remove_impl();

  m_updated_since_compress = m_updated_since_sort = true;

  internal_check_invariants();
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

    // Entity field data has relocated.
  }
}

void Partition::remove_impl()
{
  ThrowAssert(!empty());

  Bucket &last_bucket   = **(end() - 1);

  last_bucket.remove_entity();

  if ( 0 == last_bucket.size() )
  {
    size_t num_buckets = m_buckets.size();

    // Don't delete the last bucket now --- might want it later in this modification cycle.
    if (num_buckets > 1)
    {
      // The current 'last' bucket in a partition is to be deleted.
      // The previous 'last' bucket becomes the new 'last' bucket in the partition.
      Bucket *new_last = m_buckets[num_buckets - 2];
      m_buckets[0]->set_last_bucket_in_partition(new_last);
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

  m_updated_since_compress = m_updated_since_sort = true;
  --m_size;

  internal_check_invariants();
}

void Partition::compress(bool force)
{
  // An easy optimization.
  if (!force && (empty() || !m_updated_since_compress))
  {
    return;
  }

  std::vector<unsigned> partition_key = get_legacy_partition_id();
  //index of bucket in partition
  partition_key[ partition_key[0] ] = 0;

  const size_t num_entities = m_size; // m_size will change when we add new buckets
                                      // so need to store it here

  std::vector<Entity> entities(num_entities);

  // Copy the entities (but not their data) into a vector, where they will be sorted.
  //
  size_t new_i = 0;
  for (BucketVector::iterator b_i = begin(); b_i != end(); ++b_i)
  {
    Bucket &b = **b_i;
    size_t b_size = b.size();
    std::copy(&b.m_entities[0], &b.m_entities[b_size], &entities[new_i]);
    new_i += b_size;
  }

  //sort entities
  std::sort( entities.begin(), entities.end(), EntityLess(m_mesh) );
  m_updated_since_sort = false;

  // ceiling
  const size_t num_new_buckets = (num_entities + (m_repository->max_bucket_capacity - 1u)) / m_repository->max_bucket_capacity;

  BucketVector tmp_buckets(num_new_buckets);

  size_t curr_entity = 0;
  for (size_t bi=0; bi<num_new_buckets; ++bi) {

    const size_t bucket_size = ((num_entities - curr_entity) < m_repository->max_bucket_capacity) ? num_entities - curr_entity : m_repository->max_bucket_capacity;
    tmp_buckets[bi] = m_repository->allocate_bucket( m_rank, partition_key, bucket_size );
    Bucket & new_bucket = *tmp_buckets[bi];
    new_bucket.set_first_bucket_in_partition(tmp_buckets[0]);
    new_bucket.m_partition = this;

    for (size_t i=0; i<bucket_size; ++i, ++curr_entity) {
      Entity entity = entities[curr_entity];
      new_bucket.copy_entity( entity );
      // Don't worry about deleting from old buckets, they are being deallocated right after this loop
    }
  }


  //remove old buckets
  for (BucketVector::iterator b_i = begin(); b_i != end(); ++b_i) {
    m_repository->deallocate_bucket(*b_i);
  }

  m_size = num_entities;
  m_buckets.swap(tmp_buckets);

  m_updated_since_compress = false;
  m_updated_since_sort = false;

  internal_check_invariants();
}

void Partition::sort(bool force)
{
  if (!force && (empty() || !m_updated_since_sort))
  {
    return;
  }

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
    std::copy(&b.m_entities[0], &b.m_entities[0] + b_size, &entities[0] + new_i);
    new_i += b_size;
  }

  std::sort( entities.begin(), entities.end(), EntityLess(m_mesh) );

  // Make sure that there is a vacancy somewhere.
  //
  stk::mesh::Bucket *vacancy_bucket = *(buckets_end - 1);
  size_t vacancy_ordinal = vacancy_bucket->size();
  stk::mesh::Bucket *tmp_bucket = 0;

  if (vacancy_ordinal >= vacancy_bucket->capacity())
  {
    // If we need a temporary bucket, it only needs to hold one entity and
    // the corresponding field data.
    tmp_bucket = m_repository->allocate_bucket(m_rank, partition_key, 1 /* capacity */);
    vacancy_bucket = tmp_bucket;
    vacancy_ordinal = 0;
  }

  // Allocate space to copy in to
  vacancy_bucket->add_entity();

  // Now that we have the entities sorted, we need to put them and their data
  // in the right order in the buckets.

  std::vector<Entity>::iterator sorted_ent_vector_itr = entities.begin();

  Bucket* orig_vacancy_bucket = vacancy_bucket;


  FieldVector reduced_fields;
  if(m_mesh.is_field_updating_active()) {
    if(buckets_begin != buckets_end && (stk::topology::rank_t)(*buckets_begin)->entity_rank() < stk::topology::NUM_RANKS) {
      const FieldVector& rankFields = m_mesh.mesh_meta_data().get_fields((stk::topology::rank_t)(*buckets_begin)->entity_rank());
      reduced_fields.reserve(rankFields.size());
      for(unsigned ifield=0; ifield<rankFields.size(); ++ifield) {
	if(field_bytes_per_entity(*rankFields[ifield], *(*buckets_begin))) {
	  reduced_fields.push_back(rankFields[ifield]);
	}
      }
    }
  }

  for (BucketVector::iterator bucket_itr = begin(); bucket_itr != buckets_end; ++bucket_itr)
  {
    Bucket &curr_bucket = **bucket_itr;
    const unsigned n = *bucket_itr == orig_vacancy_bucket ? curr_bucket.size() -1 : curr_bucket.size(); // skip very last entity in partition

    //if(m_mesh.is_field_updating_active()) { 
      for ( unsigned curr_bucket_ord = 0; curr_bucket_ord < n ; ++curr_bucket_ord , ++sorted_ent_vector_itr ) {
	ThrowAssert(sorted_ent_vector_itr != entities.end());

	Entity curr_entity = curr_bucket[curr_bucket_ord];
	ThrowAssert(m_mesh.is_valid(curr_entity));

	if ( curr_entity != *sorted_ent_vector_itr ) // check if we need to move
	  {
	    // Move current entity to the vacant spot
	    if (vacancy_bucket != &curr_bucket || vacancy_ordinal != curr_bucket_ord) {
	      vacancy_bucket->overwrite_entity( vacancy_ordinal, curr_entity, &reduced_fields );
	    }

	    // Set the vacant spot to where the required entity is now.
	    vacancy_bucket  = & (m_mesh.bucket(*sorted_ent_vector_itr));
	    vacancy_ordinal = m_mesh.bucket_ordinal(*sorted_ent_vector_itr);

	    // Move required entity to the required spot
	    curr_bucket.overwrite_entity( curr_bucket_ord, *sorted_ent_vector_itr, &reduced_fields );
	  }
      }
      //} else {
      /*
      for ( unsigned curr_bucket_ord = 0; curr_bucket_ord < n ; ++curr_bucket_ord , ++sorted_ent_vector_itr ) {
	ThrowAssert(sorted_ent_vector_itr != entities.end());

	Entity curr_entity = curr_bucket[curr_bucket_ord];
	ThrowAssert(m_mesh.is_valid(curr_entity));

	if ( curr_entity != *sorted_ent_vector_itr ) // check if we need to move
	  {
	    // Move current entity to the vacant spot
	    if (vacancy_bucket != &curr_bucket || vacancy_ordinal != curr_bucket_ord) {
	      vacancy_bucket->overwrite_entity( vacancy_ordinal, curr_entity );
	    }

	    // Set the vacant spot to where the required entity is now.
	    vacancy_bucket  = & (m_mesh.bucket(*sorted_ent_vector_itr));
	    vacancy_ordinal = m_mesh.bucket_ordinal(*sorted_ent_vector_itr);

	    // Move required entity to the required spot
	    curr_bucket.overwrite_entity( curr_bucket_ord, *sorted_ent_vector_itr);
	  }
      }
      */
      //}


  }

  m_updated_since_sort = false;

  orig_vacancy_bucket->remove_entity();

  if (tmp_bucket) {
    m_repository->deallocate_bucket(tmp_bucket);
  }

  internal_check_invariants();
}

stk::mesh::Bucket *Partition::get_bucket_for_adds()
{
  if (no_buckets())
  {
    std::vector<unsigned> partition_key = get_legacy_partition_id();
    partition_key[ partition_key[0] ] = 0;
    Bucket *bucket = m_repository->allocate_bucket(m_rank, partition_key,
                                                   m_repository->get_bucket_capacity());
    bucket->m_partition = this;
    bucket->set_first_bucket_in_partition(bucket);
    m_buckets.push_back(bucket);

    return bucket;
  }

  Bucket *bucket = *(end() - 1);  // Last bucket of the partition.

  if (bucket->size() == bucket->capacity())
  {
    std::vector<unsigned> partition_key = get_legacy_partition_id();
    partition_key[ partition_key[0] ] = m_buckets.size();
    bucket = m_repository->allocate_bucket(m_rank, partition_key,
                                           m_repository->get_bucket_capacity());
    bucket->m_partition = this;
    bucket->set_first_bucket_in_partition(m_buckets[0]);
    m_buckets[0]->set_last_bucket_in_partition(bucket);
    m_buckets.push_back(bucket);
  }

  return bucket;
}

size_t Partition::field_data_footprint(const FieldBase& f) const
{
  size_t retval = 0;

  size_t num_bkts = m_buckets.size();
  for (size_t i = 0; i < num_bkts; ++i)
  {
    Bucket *b_ptr = m_buckets[i];
    if (b_ptr)
    {
      retval += b_ptr->capacity() * field_bytes_per_entity(f, *b_ptr);
    }
  }

  return retval;
}
