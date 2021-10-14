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

#include <stk_mesh/base/Bucket.hpp>
#include <stdlib.h>                     // for nullptr
#include <algorithm>                    // for copy, swap, lower_bound, etc
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, field_data, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/ConnectedTopologyNodes.hpp>
#include <stk_mesh/base/ConnectedSparseNodes.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology, etc
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include <stk_mesh/base/FieldTraits.hpp>
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector, etc
#include "stk_topology/topology.hpp"    // for topology::num_nodes
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include "stk_util/util/SortAndUnique.hpp"
namespace stk { namespace mesh { namespace impl { template <EntityRank TargetRank, stk::mesh::ConnectivityType> class BucketConnectivity; } } }


//----------------------------------------------------------------------
namespace stk {
namespace mesh {

namespace {

void setup_connectivity(stk::topology bucket_topology,
                        EntityRank from_rank,
                        EntityRank to_rank,
                        ConnectivityType& conn_type)
{
  if (bucket_topology != stk::topology::END_TOPOLOGY &&
      bucket_topology.num_sub_topology(to_rank) > 0 &&
      to_rank == stk::topology::NODE_RANK) {
    conn_type = FIXED_CONNECTIVITY;
  }
  else if (from_rank > stk::topology::ELEMENT_RANK ||
           to_rank > stk::topology::ELEMENT_RANK ||
           from_rank != to_rank) {
    conn_type = DYNAMIC_CONNECTIVITY;
  }
}

} //namespace anonymous

//----------------------------------------------------------------------

bool raw_part_equal( const unsigned * lhs , const unsigned * rhs )
{
  bool result = true ;
  {
    const unsigned * const end_lhs = lhs + *lhs ;
    while ( result && end_lhs != lhs ) {
      result = *lhs == *rhs ;
      ++lhs ; ++rhs ;
    }
  }
  return result ;
}

inline
bool bucket_key_less( const unsigned * lhs , const unsigned * rhs )
{
  const unsigned * const last_lhs = lhs + ( *lhs < *rhs ? *lhs : *rhs );
  while ( last_lhs != lhs && *lhs == *rhs ) { ++lhs ; ++rhs ; }
  return *lhs < *rhs ;
}

// The part count and part ordinals are less
bool BucketLess::operator()( const Bucket * lhs_bucket ,
                             const unsigned * rhs ) const
{ return bucket_key_less( lhs_bucket->key() , rhs ); }

bool BucketLess::operator()( const unsigned * lhs ,
                             const Bucket * rhs_bucket ) const
{ return bucket_key_less( lhs , rhs_bucket->key() ); }

Bucket::Bucket(BulkData & arg_mesh,
               EntityRank arg_entity_rank,
               const std::vector<unsigned> & arg_key,
               size_t arg_capacity,
               unsigned bucket_id)
  : m_mesh(arg_mesh),
    m_sparse_connectivity(arg_mesh.get_sparse_connectivity()),
    m_entity_rank(arg_entity_rank),
    m_key(arg_key),
    m_partOrdsBeginEnd(m_key.data()+1,m_key.data()+m_key[0]),
    m_topology(get_topology(m_mesh.mesh_meta_data(), arg_entity_rank, superset_part_ordinals())),
    m_capacity(arg_capacity),
    m_size(0),
    m_bucket_id(bucket_id),
    m_ngp_bucket_id(INVALID_BUCKET_ID),
    m_is_modified(true),
    m_entities(arg_capacity),
    m_partition(nullptr),
    m_node_kind(INVALID_CONNECTIVITY_TYPE),
    m_edge_kind(INVALID_CONNECTIVITY_TYPE),
    m_face_kind(INVALID_CONNECTIVITY_TYPE),
    m_topoNodes(m_topology, m_entities),
    m_sparseNodes(m_sparse_connectivity, m_entity_rank, m_entities),
    m_owned(has_superset(*this, m_mesh.mesh_meta_data().locally_owned_part())),
    m_shared(has_superset(*this, m_mesh.mesh_meta_data().globally_shared_part())),
    m_aura(has_superset(*this, m_mesh.mesh_meta_data().aura_part()))
{
  ThrowAssertMsg(arg_capacity != 0, "Buckets should never have zero capacity");
  setup_connectivity(m_topology, arg_entity_rank, stk::topology::NODE_RANK, m_node_kind);
  setup_connectivity(m_topology, arg_entity_rank, stk::topology::EDGE_RANK, m_edge_kind);
  setup_connectivity(m_topology, arg_entity_rank, stk::topology::FACE_RANK, m_face_kind);

  m_parts.reserve(m_key.size());
  supersets(m_parts);
  m_mesh.new_bucket_callback(m_entity_rank, m_parts, m_capacity, this);

  initialize_ngp_field_bucket_ids();
}

Bucket::~Bucket()
{
  m_mesh.destroy_bucket_callback(m_entity_rank, *this, m_capacity);
}

void Bucket::change_connected_nodes(unsigned bucket_ordinal, stk::mesh::Entity* new_nodes)
{
  const unsigned numNodes = m_topology.num_nodes();
  for (unsigned i=0; i<numNodes; ++i) {
    CONN_TYPE(m_node_kind, m_topoNodes.set_connected_node(bucket_ordinal, new_nodes[i], i), m_sparseNodes.reset_connected_node(bucket_ordinal, new_nodes[i], i));
  }
}

void Bucket::change_existing_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity, ConnectivityOrdinal elem_connectivity_ordinal, Permutation permut)
{
  m_sparse_connectivity.replace_permutation(entity_rank(),
                                            m_entities[bucket_ordinal_of_lower_ranked_entity],
                                            stk::topology::ELEM_RANK,
                                            elem_connectivity_ordinal, permut);
}

void Bucket::change_existing_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, ConnectivityOrdinal face_connectivity_ordinal, Permutation permut)
{
  m_sparse_connectivity.replace_permutation(entity_rank(),
                                            m_entities[bucket_ordinal_of_higher_ranked_entity],
                                            stk::topology::FACE_RANK,
                                            face_connectivity_ordinal, permut);
}

void Bucket::change_existing_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, ConnectivityOrdinal edge_connectivity_ordinal, Permutation permut)
{
  m_sparse_connectivity.replace_permutation(entity_rank(),
                                            m_entities[bucket_ordinal_of_higher_ranked_entity],
                                            stk::topology::EDGE_RANK,
                                            edge_connectivity_ordinal, permut);
}

bool Bucket::member_all( const PartVector & parts ) const
{
  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;

  bool result_all = true ;

  for ( ; result_all && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    result_all = member(ord);
  }
  return result_all ;
}

bool Bucket::member_any( const PartVector & parts ) const
{
  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    result_none = !member(ord);
  }
  return ! result_none ;
}

bool Bucket::member_any( const OrdinalVector & parts ) const
{
  const OrdinalVector::const_iterator ip_end = parts.end();
        OrdinalVector::const_iterator ip     = parts.begin() ;

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = *ip;
    result_none = !member(ord);
  }
  return ! result_none ;
}

unsigned char* Bucket::field_data_location(const FieldBase& field) const{
  return reinterpret_cast<unsigned char*>(stk::mesh::field_data(field, *this, 0));
}

//----------------------------------------------------------------------

bool has_superset( const Bucket & bucket , const PartVector & ps )
{
  const std::pair<const unsigned *, const unsigned *>
    part_ord = bucket.superset_part_ordinals();

  bool result = ! ps.empty();

  for ( PartVector::const_iterator
        i = ps.begin() ; result && i != ps.end() ; ++i ) {

    const unsigned ordinal = (*i)->mesh_meta_data_ordinal();

    const unsigned * iter =
      std::lower_bound( part_ord.first , part_ord.second , ordinal );

    result = iter < part_ord.second && ordinal == *iter ;
  }
  return result ;
}

void Bucket::supersets( PartVector & ps ) const
{
  const MetaData & mesh_meta_data = m_mesh.mesh_meta_data();

  std::pair<const unsigned *, const unsigned *>
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = & mesh_meta_data.get_part( * part_ord.first );
  }
}

void Bucket::supersets( OrdinalVector & ps ) const
{
  std::pair<const unsigned *, const unsigned *>
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = *part_ord.first;
  }
}

//----------------------------------------------------------------------

bool Bucket::field_data_is_allocated(const FieldBase& field) const
{
    return field_is_allocated_for_bucket(field, *this);
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Bucket & k )
{
  const MetaData & mesh_meta_data = k.mesh().mesh_meta_data();
  const std::string & entity_rank_name =
    mesh_meta_data.entity_rank_names()[ k.entity_rank() ];

  const PartVector& parts = k.supersets();

  s << "Bucket( " << entity_rank_name << " : " ;
  for (auto p : parts) {
    s << p->name() << " " ;
  }
  s << ")" ;

  return s ;
}

std::ostream &
print( std::ostream & os , const std::string & indent , const Bucket & bucket )
{
  const MetaData & mesh_meta_data = bucket.mesh().mesh_meta_data();
  const BulkData & mesh = bucket.mesh();
  const std::string & entity_rank_name =
    mesh_meta_data.entity_rank_names()[ bucket.entity_rank() ];

  const std::pair<const unsigned *, const unsigned *>
    part_ids = bucket.superset_part_ordinals();

  os << "Bucket(size = " << bucket.size() << std::endl << indent << "Part intersection {" ;

  for ( const unsigned * i = part_ids.first ; i < part_ids.second ; ++i ) {
    const Part & part = mesh_meta_data.get_part( *i );
    os << " " << part.name();
  }

  os << " }" << std::endl << indent << entity_rank_name << " members {" ;

  for ( unsigned j = 0 ; j < bucket.size() ; ++j ) {
    const EntityId id = mesh.identifier(bucket[j]);
    os << " " << id ;
  }
  os << " }" << std::endl ;

  return os ;
}

struct EntityRankLess
{
  inline bool operator()(Entity entity, EntityRank rank) const {
    return m_mesh->entity_rank(entity) < rank;
  }

  inline bool operator()(EntityRank rank, Entity entity) const {
    return rank < m_mesh->entity_rank(entity);
  }

  const BulkData *m_mesh;
};

//----------------------------------------------------------------------
void Bucket::initialize_ngp_field_bucket_ids()
{
  const MetaData& meta = mesh().mesh_meta_data();
  const FieldVector& allFields = meta.get_fields();
  m_ngp_field_bucket_id.resize(allFields.size());
  m_ngp_field_is_modified.resize(allFields.size());

  for(FieldBase* field : allFields) {
    m_ngp_field_bucket_id[field->mesh_meta_data_ordinal()] = INVALID_BUCKET_ID;
    m_ngp_field_is_modified[field->mesh_meta_data_ordinal()] = false;
  }
}

void Bucket::set_ngp_field_bucket_id(unsigned fieldOrdinal, unsigned ngpFieldBucketId)
{
  ThrowRequire(fieldOrdinal < m_ngp_field_bucket_id.size());
  m_ngp_field_bucket_id[fieldOrdinal] = ngpFieldBucketId;
  m_ngp_field_is_modified[fieldOrdinal] = false;
}

unsigned Bucket::get_ngp_field_bucket_id(unsigned fieldOrdinal) const
{
  ThrowRequire(fieldOrdinal < m_ngp_field_bucket_id.size());
  return m_ngp_field_bucket_id[fieldOrdinal];
}

unsigned Bucket::get_ngp_field_bucket_is_modified(unsigned fieldOrdinal) const
{
  return m_ngp_field_is_modified[fieldOrdinal];
}

void Bucket::reset_part_ord_begin_end()
{
  m_partOrdsBeginEnd.first = m_key.data()+1;
  m_partOrdsBeginEnd.second = m_key.data()+m_key[0];
}

void Bucket::reset_bucket_key(const OrdinalVector& newPartOrdinals)
{
  unsigned partitionCount = m_key[m_key.size() - 1];
  unsigned newPartCount = newPartOrdinals.size();

  m_key.resize(newPartCount + 2);
  m_key[0] = newPartCount + 1;
  m_key[newPartCount+1] = partitionCount;

  for(unsigned i = 0; i < newPartCount; i++) {
    m_key[i+1] = newPartOrdinals[i];
  }
}

void Bucket::reset_bucket_parts(const OrdinalVector& newPartOrdinals)
{
  reset_bucket_key(newPartOrdinals);
  reset_part_ord_begin_end();
  supersets(m_parts);
}

void Bucket::initialize_slot(unsigned bucketOrdinal, Entity entity)
{
  mark_for_modification();
  m_entities[bucketOrdinal]    = entity;
  if (mesh().is_valid(entity)) {
    mesh().set_state(entity, Created);
  }

  if (m_topology != stk::topology::INVALID_TOPOLOGY &&
      entity_rank() > stk::topology::NODE_RANK &&
      entity_rank() <= stk::topology::ELEM_RANK) {
    const unsigned numNodes = m_topology.num_nodes();
    ThrowRequireMsg(m_node_kind == FIXED_CONNECTIVITY, "Bucket topo="<<m_topology<<" has wrong m_node_kind.");
    for (unsigned ord = 0; ord < numNodes; ++ord) {
      m_topoNodes.set_connected_node(bucketOrdinal, Entity(), ord);
    }
  }
}

int Bucket::parallel_owner_rank(unsigned ordinal) const
{
  return m_mesh.parallel_owner_rank(m_entities[ordinal]);
}

void Bucket::reset_entity_location(Entity entity, unsigned to_ordinal, const FieldVector* fields)
{
  mark_for_modification();
  MeshIndex& meshIndex = mesh().mesh_index(entity);
  Bucket & from_bucket = *meshIndex.bucket;
  const unsigned from_ordinal = meshIndex.bucket_ordinal;

  m_entities[to_ordinal]    = entity;

  const unsigned fromNumNodes = from_bucket.num_nodes(from_ordinal);
  const unsigned numNodes = m_topology != stk::topology::INVALID_TOPOLOGY ? m_topology.num_nodes() : fromNumNodes;
  if (numNodes > 0 && fromNumNodes == numNodes) {
    const Entity* entityNodes = from_bucket.begin_nodes(from_ordinal);
    if (m_topology == stk::topology::INVALID_TOPOLOGY || from_bucket.topology() == stk::topology::INVALID_TOPOLOGY) {
      const ConnectivityOrdinal* ords = from_bucket.begin_node_ordinals(from_ordinal);
      for(unsigned i=0; i<fromNumNodes; ++i) {
        CONN_TYPE(m_node_kind, m_topoNodes.set_connected_node(to_ordinal, entityNodes[i], ords[i]), m_sparseNodes.reset_connected_node(to_ordinal, entityNodes[i], ords[i]));
      }
    }
    else {
      ThrowAssertMsg(entityNodes != nullptr, "from_bucket.begin_nodes(from_ordinal) == nullptr!");
      CONN_TYPE(m_node_kind, m_topoNodes.set_connected_nodes(to_ordinal, numNodes, entityNodes), m_sparseNodes.reset_connected_nodes(to_ordinal, numNodes, entityNodes));
    }
  }
  else {
    ThrowRequireMsg(fromNumNodes == 0, "Skipping node copy even though from-entity has "<<fromNumNodes<<" nodes");
  }

  meshIndex.bucket = this;
  meshIndex.bucket_ordinal = to_ordinal;

  m_mesh.copy_entity_fields_callback(m_entity_rank, m_bucket_id, to_ordinal,
                                     from_bucket.m_bucket_id, from_ordinal, fields);
}

void Bucket::add_entity(Entity entity)
{
  ThrowAssert(m_size < m_capacity);
  ThrowAssert(!mesh().is_valid(m_entities[m_size]));
  ThrowAssert(!mesh().is_valid(entity) || mesh().bucket_ptr(entity) == nullptr);
  ThrowAssert(!mesh().is_valid(entity) || mesh().entity_rank(entity) == m_entity_rank);

  initialize_slot(m_size, entity);

  if (mesh().is_valid(entity)) {
    mesh().set_mesh_index(entity, this, m_size);
  }

  this->mesh().add_entity_callback(entity_rank(), bucket_id(), m_size);
  ++m_size;
}

bool Bucket::destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id )
{
  EntityRank toRank = mesh().entity_rank(e_to);
  if (toRank != stk::topology::NODE_RANK) {
    EntityRank fromRank = mesh().entity_rank(e_from);
    return m_sparse_connectivity.remove_connectivity(fromRank, e_from, toRank, e_to, local_id);
  }

  Bucket& fromBucket = mesh().bucket(e_from);
  if (fromBucket.bucket_id() != bucket_id()) {
    return fromBucket.destroy_relation(e_from, e_to, local_id);
  }

  const unsigned from_bucket_ordinal = mesh().bucket_ordinal(e_from);
  
  CONN_TYPE(m_node_kind, m_topoNodes.set_connected_node(from_bucket_ordinal, Entity(), local_id), m_sparseNodes.set_connected_node(from_bucket_ordinal, Entity(), local_id));
  return true;
}

bool Bucket::declare_relation(unsigned bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation )
{
  EntityRank toRank = mesh().entity_rank(e_to);
  if (toRank == stk::topology::NODE_RANK) {
    return (m_node_kind==FIXED_CONNECTIVITY ? m_topoNodes.set_connected_node(bucket_ordinal, e_to, ordinal) : m_sparseNodes.reset_connected_node(bucket_ordinal, e_to, ordinal));
  }

  return m_sparse_connectivity.add_connectivity(entity_rank(), m_entities[bucket_ordinal],
                                                           toRank, e_to, ordinal, permutation);
}

void Bucket::remove_entity()
{
  ThrowAssert(m_size > 0);

  mark_for_modification();
  mesh().remove_entity_field_data_callback(entity_rank(), bucket_id(), m_size-1);
  --m_size;

  initialize_slot(m_size, Entity());
}

void Bucket::copy_entity(Entity entity)
{
  ThrowAssert(m_size < m_capacity);
  ThrowAssert(mesh().is_valid(entity));
  ThrowAssert(!mesh().is_valid(m_entities[m_size]));
  ThrowAssert(mesh().bucket_ptr(entity) != nullptr);
  ThrowAssert(mesh().bucket_ptr(entity) != this);
  ThrowAssert(mesh().entity_rank(entity) == m_entity_rank);

  mark_for_modification();

  this->mesh().add_entity_callback(this->entity_rank(), this->bucket_id(), m_size);
  reset_entity_location(entity, m_size);

  ++m_size;
}

void Bucket::overwrite_entity(unsigned to_ordinal, Entity entity, const FieldVector* fields)
{
  ThrowAssert(to_ordinal < m_capacity);
  ThrowAssert(mesh().is_valid(entity));
  ThrowAssert(mesh().bucket_ptr(entity) != nullptr);
  ThrowAssert(mesh().entity_rank(entity) == m_entity_rank);

  reset_entity_location(entity, to_ordinal, fields);
}


void Bucket::parent_topology( EntityRank parent_rank, std::vector<stk::topology> & parent_topologies) const
{
  PartVector parts;
  supersets(parts);

  parent_topologies.clear();

  for (auto& part : parts ) {
    Part const & p = *part;
    if ( (p.primary_entity_rank() == parent_rank) && p.topology().is_valid()) {
      parent_topologies.push_back(p.topology());
    }
  }

  stk::util::sort_and_unique(parent_topologies);
}

} // namespace mesh
} // namespace stk
