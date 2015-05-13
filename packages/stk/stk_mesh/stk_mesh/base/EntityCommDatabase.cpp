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

#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <sstream>                      // for operator<<, basic_ostream
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Relation.hpp>   // for Relation
#include <string>                       // for operator<<
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssertMsg
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer


namespace stk {
namespace mesh {

//----------------------------------------------------------------------------

namespace {

unsigned count_parallel_consistent_parts(const MetaData & meta, const unsigned* first, const unsigned* last) {
    unsigned count = 0;
    for (unsigned part_index=0; part_index < last-first; ++part_index) {
        const unsigned part_ordinal = first[part_index];
        if ( (part_ordinal != meta.locally_owned_part().mesh_meta_data_ordinal()) &&
                (part_ordinal != meta.globally_shared_part().mesh_meta_data_ordinal()) &&
                (meta.get_parts()[part_ordinal]->entity_membership_is_parallel_consistent() )) {
            ++count;
        }
    }
    return count;
}

void pack_bucket_part_list(const Bucket & bucket, CommBuffer & buf ) {
    const MetaData & meta = bucket.mesh().mesh_meta_data();
    const std::pair<const unsigned *, const unsigned *>
      part_ordinals = bucket.superset_part_ordinals();
    buf.pack<unsigned>( count_parallel_consistent_parts(meta, part_ordinals.first, part_ordinals.second) );
    unsigned nparts = part_ordinals.second - part_ordinals.first;
    for (unsigned part_index=0; part_index < nparts; ++part_index) {
        const unsigned part_ordinal = part_ordinals.first[part_index];
        if ( (part_ordinal != meta.locally_owned_part().mesh_meta_data_ordinal()) &&
             (part_ordinal != meta.globally_shared_part().mesh_meta_data_ordinal()) &&
             (meta.get_parts()[part_ordinal]->entity_membership_is_parallel_consistent() )) {
            buf.pack<unsigned>(part_ordinal);
        }
    }
}

}

void pack_entity_info(const BulkData& mesh, CommBuffer & buf , const Entity entity )
{
  const EntityKey & key   = mesh.entity_key(entity);
  const unsigned    owner = mesh.parallel_owner_rank(entity);

  buf.pack<EntityKey>( key );
  buf.pack<unsigned>( owner );
  pack_bucket_part_list(mesh.bucket(entity), buf);

  const unsigned tot_rel = mesh.count_relations(entity);
  buf.pack<unsigned>( tot_rel );

  Bucket& bucket = mesh.bucket(entity);
  unsigned ebo   = mesh.bucket_ordinal(entity);

  ThrowAssertMsg(mesh.is_valid(entity), "BulkData at " << &mesh << " does not know Entity " << entity.local_offset());
  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    const unsigned nrel = bucket.num_connectivity(ebo, irank);
    if (nrel > 0) {
      Entity const *rel_entities = bucket.begin(ebo, irank);
      ConnectivityOrdinal const *rel_ordinals = bucket.begin_ordinals(ebo, irank);
      Permutation const *rel_permutations = bucket.begin_permutations(ebo, irank);
      ThrowAssert(rel_ordinals);
      ThrowAssert(rel_permutations);

      for ( unsigned i = 0 ; i < nrel ; ++i ) {
        if (mesh.is_valid(rel_entities[i])) {
          buf.pack<EntityKey>( mesh.entity_key(rel_entities[i]) );
          buf.pack<unsigned>( rel_ordinals[i] );
          if (bucket.has_permutation(irank)) {
            buf.pack<unsigned>( rel_permutations[i] );
          } else {
            buf.pack<unsigned>(0u);
          }
        }
      }
    }
  }
}

void unpack_entity_info(
  CommBuffer       & buf,
  const BulkData   & mesh ,
  EntityKey        & key ,
  int              & owner ,
  PartVector       & parts ,
  std::vector<Relation> & relations )
{
  unsigned nparts = 0 ;
  unsigned nrel = 0 ;

  buf.unpack<EntityKey>( key );
  buf.unpack<int>( owner );
  buf.unpack<unsigned>( nparts );

  parts.resize( nparts );

  for ( unsigned i = 0 ; i < nparts ; ++i ) {
    unsigned part_ordinal = ~0u ;
    buf.unpack<unsigned>( part_ordinal );
    parts[i] = & MetaData::get(mesh).get_part( part_ordinal );
  }

  buf.unpack( nrel );

  relations.clear();
  relations.reserve( nrel );

  for ( unsigned i = 0 ; i < nrel ; ++i ) {
    EntityKey rel_key ;
    unsigned rel_id = 0 ;
    unsigned rel_attr = 0 ;
    buf.unpack<EntityKey>( rel_key );
    buf.unpack<unsigned>( rel_id );
    buf.unpack<unsigned>( rel_attr );
    Entity const entity =
      mesh.get_entity( rel_key.rank(), rel_key.id() );
    if ( mesh.is_valid(entity) ) {
      Relation rel(entity, mesh.entity_rank(entity), rel_id );
      rel.set_attribute(rel_attr);
      relations.push_back( rel );
    }
  }
}


//----------------------------------------------------------------------

void pack_field_values(const BulkData& mesh, CommBuffer & buf , Entity entity )
{
    if (!mesh.is_field_updating_active()) {
        return;
    }
    const Bucket   & bucket = mesh.bucket(entity);
    const MetaData & mesh_meta_data = MetaData::get(mesh);
    const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields(bucket.entity_rank());
    for ( std::vector< FieldBase * >::const_iterator
            i = fields.begin() ; i != fields.end() ; ++i ) {
        const FieldBase & f = **i ;
        if ( f.data_traits().is_pod ) {
            const unsigned size = field_bytes_per_entity( f, bucket );
#ifndef NDEBUG
            buf.pack<unsigned>( size );
#endif
            if ( size ) {
                unsigned char * const ptr =
                        reinterpret_cast<unsigned char *>( stk::mesh::field_data( f , entity ) );
                buf.pack<unsigned char>( ptr , size );
            }
        }
    }
}

bool unpack_field_values(const BulkData& mesh,
                         CommBuffer & buf , Entity entity , std::ostream & error_msg )
{
    if (!mesh.is_field_updating_active()) {
        return true;
    }
    const Bucket   & bucket = mesh.bucket(entity);
    const MetaData & mesh_meta_data = MetaData::get(mesh);
    const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields(bucket.entity_rank());
    const std::vector< FieldBase * >::const_iterator i_end = fields.end();
    const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();
    std::vector< FieldBase * >::const_iterator i ;
    bool ok = true ;
    for ( i = i_beg ; i_end != i ; ) {
        const FieldBase & f = **i ; ++i ;
        if ( f.data_traits().is_pod ) {
            const unsigned size = field_bytes_per_entity( f, bucket );
#ifndef NDEBUG
            unsigned recv_data_size = 0 ;
            buf.unpack<unsigned>( recv_data_size );
            if ( size != recv_data_size ) {
                if ( ok ) {
                    ok = false ;
                    error_msg << mesh.identifier(entity);
                }
                error_msg << " " << f.name();
                error_msg << " " << size ;
                error_msg << " != " << recv_data_size ;
                buf.skip<unsigned char>( recv_data_size );
            }
#endif
            if ( size )
            { // Non-zero and equal
                unsigned char * ptr =
                        reinterpret_cast<unsigned char *>( stk::mesh::field_data( f , entity ) );
                buf.unpack<unsigned char>( ptr , size );
            }
        }
    }
    return ok ;
}

//----------------------------------------------------------------------

// A cached find function that stores the result in m_last_lookup if successful and returns true.
// Otherwise, the find failed and it returns false.
bool EntityCommDatabase::cached_find(const EntityKey& key) const
{
  if (m_last_lookup != m_comm_map.end() && key == m_last_lookup->first) {
    return true;
  }

  map_type::iterator find_it = const_cast<map_type&>(m_comm_map).find(key);
  if (find_it == m_comm_map.end()) {
    return false;
  }
  else {
    m_last_lookup = find_it;
    return true;
  }
}


void EntityCommDatabase::insert(const EntityKey& key)
{
  if (!cached_find(key)) {
    m_last_lookup = m_comm_map.insert(std::make_pair(key, EntityComm())).first;
  }
}


int EntityCommDatabase::owner_rank( const EntityKey & key ) const
{
  if (!cached_find(key)) return InvalidProcessRank;

  return m_last_lookup->second.owner_rank;
}

PairIterEntityComm EntityCommDatabase::shared_comm_info( const EntityKey & key ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  return shared_comm_info_range(m_last_lookup->second.comm_map);
}


PairIterEntityComm EntityCommDatabase::comm( const EntityKey & key ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  const EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;
  return PairIterEntityComm(comm_map);
}

const EntityComm* EntityCommDatabase::entity_comm( const EntityKey & key ) const
{
  if (!cached_find(key)) return NULL;

  const EntityComm& entity_comm = m_last_lookup->second;
  return &entity_comm;
}


PairIterEntityComm EntityCommDatabase::comm( const EntityKey & key, const Ghosting & sub ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  const EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  const EntityCommInfo s_begin( sub.ordinal() ,     0 );
  const EntityCommInfo s_end(   sub.ordinal() + 1 , 0 );

  EntityCommInfoVector::const_iterator i = comm_map.begin();
  EntityCommInfoVector::const_iterator e = comm_map.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  return PairIterEntityComm( i , e );
}


bool EntityCommDatabase::insert( const EntityKey & key, const EntityCommInfo & val, int owner )
{
  insert(key);
  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;
  m_last_lookup->second.owner_rank = owner;

  EntityCommInfoVector::iterator i =
    std::lower_bound( comm_map.begin() , comm_map.end() , val );

  const bool result = ((i == comm_map.end()) || (val != *i));

  if ( result ) {
    comm_map.insert( i , val );
  }

  return result;
}


bool EntityCommDatabase::erase( const EntityKey & key, const EntityCommInfo & val )
{
  if (!cached_find(key)) return false;

  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  EntityCommInfoVector::iterator i =
    std::lower_bound( comm_map.begin() , comm_map.end() , val );

  const bool result = ( (i != comm_map.end()) && (val == *i) ) ;

  if ( result ) {
    comm_map.erase( i );
    if (comm_map.empty()) {
      m_last_lookup = m_comm_map.erase(m_last_lookup);
    }
  }

  return result ;
}


bool EntityCommDatabase::erase( const EntityKey & key, const Ghosting & ghost )
{
  if (!cached_find(key)) return false;

  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  const EntityCommInfo s_begin( ghost.ordinal() ,     0 );
  const EntityCommInfo s_end(   ghost.ordinal() + 1 , 0 );

  EntityCommInfoVector::iterator i = comm_map.begin();
  EntityCommInfoVector::iterator e = comm_map.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  const bool result = i != e ;

  if ( result ) {
    comm_map.erase( i , e );
    if (comm_map.empty()) {
      m_last_lookup = m_comm_map.erase(m_last_lookup);
    }
  }

  // if there is no more comm info, just remove it from the map?

  return result ;
}


void EntityCommDatabase::comm_clear_ghosting(const EntityKey & key)
{
  if (!cached_find(key)) return;

  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  EntityCommInfoVector::iterator j = comm_map.begin();
  while ( j != comm_map.end() && j->ghost_id == 0 ) { ++j ; }
  comm_map.erase( j , comm_map.end() );

  if (comm_map.empty()) {
    m_last_lookup = m_comm_map.erase(m_last_lookup);
  }
}


void EntityCommDatabase::comm_clear(const EntityKey & key)
{
  if (!cached_find(key)) return;

  m_last_lookup = m_comm_map.erase(m_last_lookup);
}


bool EntityCommDatabase::change_owner_rank(const EntityKey& key, int owner)
{
  // Do not add key to map, only update rank if it was already in the map
  if (cached_find(key)) {
    const int orig_owner = m_last_lookup->second.owner_rank;
    m_last_lookup->second.owner_rank = owner;
    return orig_owner != owner;
  }
  return false;
}



} // namespace mesh
}

