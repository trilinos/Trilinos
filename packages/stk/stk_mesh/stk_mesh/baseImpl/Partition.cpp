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
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
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
  , m_updated_since_sort(false)
{
  // Nada.
}

// Only the BucketRepository will delete a Partition.
Partition::~Partition()
{
  size_t num_bkts = m_buckets.size();
  for (size_t i = 0; i < num_bkts; ++i)
  {
    delete m_buckets[i];
  }
}

bool Partition::remove(Entity entity)
{
  ThrowAssert(belongs(m_mesh.bucket(entity)));

  Bucket &bucket   = m_mesh.bucket(entity);
  unsigned ordinal = m_mesh.bucket_ordinal(entity);
  overwrite_from_end(bucket, ordinal);

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
  ++m_size;

  m_updated_since_sort = true;

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
  Bucket *dst_bucket = nullptr;
  try {
      dst_bucket = dst_partition.get_bucket_for_adds();
  }
  catch(std::exception& e) {
      std::ostringstream os;
      os << "Error adding "<<m_mesh.entity_key(entity)<< " to mesh: " << e.what();
      throw std::runtime_error(os.str());
  }

  ThrowErrorMsgIf(src_bucket && src_bucket->topology().is_valid() && (src_bucket->topology() != dst_bucket->topology()),
                  "Error: cannot change topology of entity (rank: "
                  << static_cast<stk::topology::rank_t>(m_mesh.entity_rank(entity))
                  << ", global_id: " << m_mesh.identifier(entity) << ") from "
                  << src_bucket->topology() << "to " << dst_bucket->topology() << "."
                  );

  // Copy the entity's data to the new bucket before removing the entity from its old bucket.
  dst_bucket->copy_entity(entity);

  overwrite_from_end(*src_bucket, src_ordinal);

  dst_partition.m_updated_since_sort = true;
  dst_partition.m_size++;

  remove_impl();

  m_updated_since_sort = true;

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

void Partition::delete_bucket(Bucket * bucket)
{
    m_size -= bucket->size();

    auto iter = std::find(m_buckets.begin(), m_buckets.end(), bucket);
    m_repository->deallocate_bucket(bucket);
    m_buckets.erase(iter, iter+1);
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

  m_updated_since_sort = true;
  --m_size;

  internal_check_invariants();
}

void Partition::default_sort_if_needed()
{
  if (!empty() && m_updated_since_sort)
  {
      sort(GlobalIdEntitySorter());
  }
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
