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

#include <stk_mesh/baseImpl/BucketRepository.hpp>
#include <algorithm>                    // for copy, remove_if
#include <new>                          // for operator new
#include <sstream>                      // for operator<<, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, raw_part_equal
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition, lower_bound
#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc

namespace stk {
namespace mesh {
namespace impl {

BucketRepository::BucketRepository(BulkData & mesh,
                                   unsigned entity_rank_count,
                                   const ConnectivityMap& connectivity_map,
                                   unsigned bucket_capacity)
  : m_mesh(mesh),
    m_buckets(entity_rank_count),
    m_partitions(entity_rank_count),
    m_need_sync_from_partitions(entity_rank_count, false),
    m_connectivity_map(connectivity_map),
    m_bucket_capacity(bucket_capacity),
    m_being_destroyed(false)
{
  // Nada.
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

size_t BucketRepository::total_field_data_footprint(const FieldBase& f, EntityRank rank) const
{
  if (rank > m_partitions.size() || static_cast<unsigned>(f.entity_rank()) != rank)
  {
    return 0;
  }

  size_t retval = 0;
  const std::vector<Partition *> &r_partitions = m_partitions[rank];
  size_t num_partitions = r_partitions.size();
  for (size_t i = 0; i < num_partitions; ++i)
  {
    retval += r_partitions[i]->field_data_footprint(f);
  }
  return retval;
}

void BucketRepository::internal_sort_bucket_entities()
{
  for (std::vector<std::vector<Partition *> >::const_iterator
         i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
  {
    const std::vector<Partition *> & pset = *i ;
    for ( std::vector<Partition*>::const_iterator
            ip = pset.begin() ; ip != pset.end() ; ++ip )
    {
      (*ip)->sort();
    }
  }
}

void BucketRepository::optimize_buckets()
{
  for (std::vector<std::vector<Partition *> >::const_iterator
         i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
  {
    const std::vector<Partition *> & pset = *i ;
    for ( std::vector<Partition*>::const_iterator
            ip = pset.begin() ; ip != pset.end() ; ++ip )
    {
      (*ip)->compress();
    }
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
  enum { KEY_TMP_BUFFER_SIZE = 64 };

  ThrowRequireMsg(MetaData::get(m_mesh).check_rank(arg_entity_rank),
                  "Entity rank " << arg_entity_rank << " is invalid");

  if (m_buckets.empty()) {
    size_t entity_rank_count = m_mesh.mesh_meta_data().entity_rank_count();
    ThrowRequireMsg( entity_rank_count > 0,
                   "MetaData doesn't have any entity-ranks! Did you forget to initialize MetaData before creating BulkData?");
    m_buckets.resize(entity_rank_count);
    m_partitions.resize(entity_rank_count);
    m_need_sync_from_partitions.resize(entity_rank_count, false);
  }

  std::vector<Partition *> & partitions = m_partitions[ arg_entity_rank ];

  const size_t part_count = parts.size();
  std::vector<unsigned> key(2 + part_count) ;

  //----------------------------------
  // Key layout:
  // { part_count + 1 , { part_ordinals } , partition_count }
  // Thus partition_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key for a bucket in the partition.
  const unsigned max = static_cast<unsigned>(-1);
  key[0] = part_count+1;
  key[ key[0] ] = max ;

  {
    for ( unsigned i = 0 ; i < part_count ; ++i ) { key[i+1] = parts[i] ; }
  }

  // If the partition is found, the iterator will be right after it, thanks to the
  // trickiness above.
  const std::vector<Partition *>::iterator ik = lower_bound( partitions , &key[0] );
  const bool partition_exists =
    (ik != partitions.begin()) && raw_part_equal( ik[-1]->key() , &key[0] );

  if (partition_exists)
  {
    return ik[-1];
  }

  key[key[0]] = 0;

  Partition *partition = new Partition(m_mesh, this, arg_entity_rank, key);
  ThrowRequire(partition != NULL);

  m_need_sync_from_partitions[arg_entity_rank] = true;
  partitions.insert( ik , partition );

  return partition ;
}

void BucketRepository::internal_modification_end()
{
  sync_from_partitions();

  // What needs to be done depends on the connectivity map.
  for (EntityRank from_rank = stk::topology::NODE_RANK;
        from_rank < m_connectivity_map.m_map.size();
        ++from_rank)
  {
    const BucketVector &buckets = m_buckets[from_rank];
    unsigned num_buckets = buckets.size();
    for (unsigned j = 0; j < num_buckets; ++j)
    {
      ThrowAssert(buckets[j] != NULL);
      Bucket &bucket = *buckets[j];

      // Update the hop-saving connectivity data on this bucket.
      //
      for (EntityRank to_rank = stk::topology::NODE_RANK;
          to_rank < m_connectivity_map.m_map[from_rank].size();
          ++to_rank)
      {
        switch (m_connectivity_map.m_map[from_rank][to_rank])
        {
        case FIXED_CONNECTIVITY:
          switch (to_rank)
          {
          case stk::topology::NODE_RANK:
            bucket.m_fixed_node_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::EDGE_RANK:
            bucket.m_fixed_edge_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::FACE_RANK:
            bucket.m_fixed_face_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::ELEMENT_RANK:
            bucket.m_fixed_element_connectivity.end_modification(&bucket.m_mesh);
            break;
          default:
            break;
          }
          break;
        case DYNAMIC_CONNECTIVITY:
          switch (to_rank)
          {
          case stk::topology::NODE_RANK:
            bucket.m_dynamic_node_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::EDGE_RANK:
            bucket.m_dynamic_edge_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::FACE_RANK:
            bucket.m_dynamic_face_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::ELEMENT_RANK:
            bucket.m_dynamic_element_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::INVALID_RANK:
            break;
          default:
            bucket.m_dynamic_other_connectivity.end_modification(&bucket.m_mesh);
            break;
          }
          break;
        case INVALID_CONNECTIVITY_TYPE:
        default:
          break;
        }
      }
    }
  }
}

void BucketRepository::sync_from_partitions()
{
  for (EntityRank rank = stk::topology::NODE_RANK; rank < m_partitions.size(); ++rank)
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
        if (first->size() == 0)
        {
            return true;
        }
        else if (second->size() == 0)
        {
            return false;
        }
        else
        {
            const stk::mesh::BulkData& mesh = first->mesh();
            stk::mesh::EntityId firstId = mesh.identifier((*first)[0]);
            stk::mesh::EntityId secondId = mesh.identifier((*second)[0]);
            return firstId < secondId;
       }
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

Bucket *BucketRepository::allocate_bucket(EntityRank arg_entity_rank,
                                          const std::vector<unsigned> & arg_key,
                                          size_t arg_capacity )
{
  BucketVector &bucket_vec = m_buckets[arg_entity_rank];
  const unsigned bucket_id = bucket_vec.size();

  Bucket * new_bucket = new Bucket(m_mesh, arg_entity_rank, arg_key, arg_capacity, m_connectivity_map, bucket_id);
  ThrowRequire(new_bucket != NULL);

  bucket_vec.push_back(new_bucket);
  m_need_sync_from_partitions[arg_entity_rank] = true;

  return new_bucket;
}

void BucketRepository::deallocate_bucket(Bucket *b)
{
  ThrowAssertMsg(b != NULL,
                 "BucketRepository::deallocate_bucket(.) m_buckets invariant broken.");

  const unsigned bucket_id = b->bucket_id();
  const EntityRank bucket_rank = b->entity_rank();

  ThrowAssertMsg(b == m_buckets[bucket_rank][bucket_id],
                 "BucketRepository::deallocate_bucket(.) m_buckets invariant broken.");

  m_buckets[bucket_rank][bucket_id] = NULL; // space will be reclaimed by sync_from_partitions
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
    ThrowAssertMsg(buckets[i] != NULL,
                   "BucketRepository::sync_bucket_ids() called when m_buckets["
                   << entity_rank << "] is not dense.");
    id_map[i] = buckets[i]->bucket_id();
    buckets[i]->m_bucket_id = i;
  }

  m_mesh.reorder_buckets_callback(entity_rank, id_map);
}

std::vector<Partition *> BucketRepository::get_partitions(EntityRank rank) const
{
  if (!m_mesh.in_synchronized_state())
  {
    std::vector<Partition *>();
  }
  std::vector<Partition *> retval;
  std::vector<Partition *> const& bf_vec = m_partitions[rank];
  for (size_t i = 0; i < bf_vec.size(); ++i)
  {
    retval.push_back(bf_vec[i]);
  }
  return retval;
}

} // namespace impl
} // namespace mesh
} // namespace stk
