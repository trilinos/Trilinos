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
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology, etc
#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector, etc
#include "stk_topology/topology.hpp"    // for topology::num_nodes
#include "stk_util/util/ReportHandler.hpp"  // for STK_ThrowAssert, etc
#include "stk_util/util/SortAndUnique.hpp"
namespace stk { namespace mesh { namespace impl { template <EntityRank TargetRank, stk::mesh::ConnectivityType> class BucketConnectivity; } } }



//----------------------------------------------------------------------
namespace stk {
namespace mesh {

namespace {

struct AddEntityFunctor
{
  template <EntityRank Rank, ConnectivityType OtherType>
  void operator()(Bucket&, impl::BucketConnectivity<Rank,OtherType>& connectivity, Bucket*)
  { connectivity.add_entity(); }

  template <EntityRank Rank, ConnectivityType OtherType>
  void operator()(Bucket& thisBucket, impl::BucketConnDynamic& connectivity, Bucket*)
  {connectivity.grow_if_necessary(thisBucket.size());}

  bool is_modifying() const { return true; }
};

struct RemoveEntityFunctor
{
  RemoveEntityFunctor(unsigned bktOrdinal)
  : m_bucketOrdinal(bktOrdinal)
  {}

  template <EntityRank Rank, ConnectivityType OtherType>
  void operator()(Bucket&, impl::BucketConnectivity<Rank,OtherType>& connectivity, Bucket*)
  { connectivity.remove_entity(); }

  template <EntityRank Rank, ConnectivityType OtherType>
  void operator()(Bucket&, impl::BucketConnDynamic& connectivity, Bucket*)
  { connectivity.remove_connectivity(m_bucketOrdinal); }

  bool is_modifying() const { return true; }

  unsigned m_bucketOrdinal;
};

struct ClearEntityFunctor
{
  ClearEntityFunctor(unsigned bktOrdinal)
  : m_bucketOrdinal(bktOrdinal)
  {}

  template<typename ConnectivityType>
  void operator()(Bucket&, ConnectivityType& connectivity)
  {}

  void operator()(Bucket&, impl::BucketConnDynamic& connectivity)
  { connectivity.remove_connectivity(m_bucketOrdinal); }

  bool is_modifying() const { return true; }

  unsigned m_bucketOrdinal;
};

struct ReplaceEntityFunctor
{
  ReplaceEntityFunctor(unsigned srcOrdinal, unsigned destOrdinal)
  : m_srcOrdinal(srcOrdinal), m_destOrdinal(destOrdinal)
  {}

  template <EntityRank Rank, ConnectivityType OtherType>
  void operator()(Bucket&, impl::BucketConnectivity<Rank,OtherType>& connectivity, Bucket*)
  {
    connectivity.replace_connectivity(m_destOrdinal,
                      connectivity.num_connectivity(m_srcOrdinal),
                      connectivity.begin(m_srcOrdinal),
                      connectivity.begin_ordinals(m_srcOrdinal),
                      connectivity.begin_permutations(m_srcOrdinal));
  }

  template <EntityRank Rank, ConnectivityType OtherType>
  void operator()(Bucket&, impl::BucketConnDynamic& connectivity, Bucket*)
  { connectivity.swap_connectivity(m_srcOrdinal, m_destOrdinal); }

  bool is_modifying() const { return true; }

  unsigned m_srcOrdinal;
  unsigned m_destOrdinal;
};

struct DeclareRelationFunctor
{
  DeclareRelationFunctor(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal,
                         Permutation permutation)
    : m_bucket_ordinal(bucket_ordinal),
      m_to(to),
      m_ordinal(ordinal),
      m_permutation(permutation),
      m_modified(false)
  {}

  template <typename Connectivity>
  void operator()(Bucket& bucket, Connectivity& connectivity)
  {
    STK_ThrowAssert(!m_modified);
    m_modified = connectivity.add_connectivity(m_bucket_ordinal, m_to, m_ordinal, m_permutation);
  }

  unsigned m_bucket_ordinal;
  Entity m_to;
  ConnectivityOrdinal m_ordinal;
  Permutation m_permutation;
  bool m_modified;
};

struct DestroyRelationFunctor
{
  DestroyRelationFunctor(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal)
    : m_bucket_ordinal(bucket_ordinal),
      m_to(to),
      m_ordinal(ordinal),
      m_modified(false)
  {}

  template <typename Connectivity>
  void operator()(Bucket& bucket, Connectivity& connectivity)
  {
    STK_ThrowAssert(!m_modified);
    m_modified = connectivity.remove_connectivity(m_bucket_ordinal, m_to, m_ordinal);
  }

  unsigned m_bucket_ordinal;
  Entity m_to;
  ConnectivityOrdinal m_ordinal;
  bool m_modified;
};

struct ReplaceRelationFunctor
{
  ReplaceRelationFunctor(unsigned bucket_ordinal,
                         unsigned numConnectivity,
                         const Entity* connectivity,
                         const ConnectivityOrdinal* ordinals,
                         const Permutation* permutations)
    : m_bucket_ordinal(bucket_ordinal),
      m_numConnectivity(numConnectivity),
      m_connectivity(connectivity),
      m_ordinals(ordinals),
      m_permutations(permutations),
      m_modified(false)
  {}

  template <typename Connectivity>
  void operator()(Bucket& bucket, Connectivity& connectivity)
  {
    STK_ThrowAssert(!m_modified);
    m_modified = connectivity.replace_connectivity(m_bucket_ordinal, m_numConnectivity,
                                  m_connectivity, m_ordinals, m_permutations);
  }

  unsigned m_bucket_ordinal;
  unsigned m_numConnectivity;
  const Entity* m_connectivity;
  const ConnectivityOrdinal* m_ordinals;
  const Permutation* m_permutations;
  bool m_modified;
};

template <typename FixedConnectivity>
void setup_connectivity(stk::topology bucket_topology,
                        EntityRank from_rank,
                        EntityRank to_rank,
                        ConnectivityType& conn_type,
                        FixedConnectivity& fixed_conn)
{
  if (bucket_topology != stk::topology::END_TOPOLOGY &&
      bucket_topology.num_sub_topology(to_rank) > 0 &&
      to_rank == stk::topology::NODE_RANK) {
    fixed_conn.set_num_connectivity(bucket_topology.num_sub_topology(to_rank));
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

namespace impl {

static const unsigned default_initial_bucket_capacity = 16;
static const unsigned default_maximum_bucket_capacity = 512;

}

unsigned get_default_bucket_capacity() { return impl::default_maximum_bucket_capacity; }
unsigned get_default_initial_bucket_capacity() { return impl::default_initial_bucket_capacity; }
unsigned get_default_maximum_bucket_capacity() { return impl::default_maximum_bucket_capacity; }

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

//----------------------------------------------------------------------

Bucket::Bucket(BulkData & mesh,
               EntityRank entityRank,
               const std::vector<unsigned> & key,
               unsigned initialCapacity,
               unsigned maximumCapacity,
               unsigned bucketId)
  : m_mesh(mesh),
    m_entity_rank(entityRank),
    m_topology(),
    m_key(key),
    m_partOrdsBeginEnd(m_key.data()+1,m_key.data()+m_key[0]),
    m_capacity(initialCapacity),
    m_maxCapacity(maximumCapacity),
    m_size(0),
    m_bucket_id(bucketId),
    m_ngp_bucket_id(INVALID_BUCKET_ID),
    m_is_modified(true),
    m_entities(maximumCapacity),
    m_partition(nullptr),
    m_node_kind(INVALID_CONNECTIVITY_TYPE),
    m_edge_kind(INVALID_CONNECTIVITY_TYPE),
    m_face_kind(INVALID_CONNECTIVITY_TYPE),
    m_element_kind(INVALID_CONNECTIVITY_TYPE),
    m_fixed_node_connectivity(),
    m_fixed_edge_connectivity(),
    m_fixed_face_connectivity(),
    m_fixed_element_connectivity(),
    m_dynamic_node_connectivity(initialCapacity),
    m_dynamic_edge_connectivity(initialCapacity, should_store_permutations(entityRank, stk::topology::EDGE_RANK)),
    m_dynamic_face_connectivity(initialCapacity, should_store_permutations(entityRank, stk::topology::FACE_RANK)),
    m_dynamic_element_connectivity(initialCapacity, should_store_permutations(entityRank, stk::topology::ELEM_RANK)),
    m_dynamic_other_connectivity(initialCapacity),
    m_owned(has_superset(*this, m_mesh.mesh_meta_data().locally_owned_part())),
    m_shared(has_superset(*this, m_mesh.mesh_meta_data().globally_shared_part())),
    m_aura(has_superset(*this, m_mesh.mesh_meta_data().aura_part()))
{
  STK_ThrowAssertMsg(initialCapacity != 0, "Buckets should never have zero capacity");
  STK_ThrowAssertMsg(initialCapacity <= maximumCapacity, "Initial capacity cannot exceed the maximum capacity");

  m_topology = get_topology(m_mesh.mesh_meta_data(), entityRank, superset_part_ordinals());

  setup_connectivity(m_topology, entityRank, stk::topology::NODE_RANK, m_node_kind, m_fixed_node_connectivity);
  setup_connectivity(m_topology, entityRank, stk::topology::EDGE_RANK, m_edge_kind, m_fixed_edge_connectivity);
  setup_connectivity(m_topology, entityRank, stk::topology::FACE_RANK, m_face_kind, m_fixed_face_connectivity);
  setup_connectivity(m_topology, entityRank, stk::topology::ELEMENT_RANK, m_element_kind, m_fixed_element_connectivity);

  m_parts.reserve(m_key.size());
  supersets(m_parts);
  m_mesh.new_bucket_callback(m_entity_rank, m_parts, m_capacity, this);

  initialize_ngp_field_bucket_ids();
}

Bucket::~Bucket()
{
  m_mesh.destroy_bucket_callback(m_entity_rank, *this, m_capacity);
}

size_t Bucket::memory_size_in_bytes() const
{
  size_t bytes = sizeof(Bucket);
  bytes += impl::capacity_in_bytes(m_entities);
  bytes += m_fixed_node_connectivity.heap_memory_in_bytes();
  bytes += m_fixed_edge_connectivity.heap_memory_in_bytes();
  bytes += m_fixed_face_connectivity.heap_memory_in_bytes();
  bytes += m_fixed_element_connectivity.heap_memory_in_bytes();
  bytes += m_dynamic_node_connectivity.heap_memory_in_bytes();
  bytes += m_dynamic_edge_connectivity.heap_memory_in_bytes();
  bytes += m_dynamic_face_connectivity.heap_memory_in_bytes();
  bytes += m_dynamic_element_connectivity.heap_memory_in_bytes();
  bytes += m_dynamic_other_connectivity.heap_memory_in_bytes();
  return bytes;
}

void Bucket::change_existing_connectivity(unsigned bucket_ordinal, stk::mesh::Entity* new_nodes)
{
    unsigned num_nodes = this->num_nodes(bucket_ordinal);
    Entity *nodes=0;
    if (m_node_kind == FIXED_CONNECTIVITY)
    {
        nodes = m_fixed_node_connectivity.begin(bucket_ordinal);
    }
    else
    {
        nodes = m_dynamic_node_connectivity.begin(bucket_ordinal);
    }

    for (unsigned i=0;i<num_nodes;++i)
    {
        nodes[i] = new_nodes[i];
    }
}

void Bucket::change_existing_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity, unsigned elem_connectivity_ordinal, stk::mesh::Permutation permut)
{
    stk::mesh::Permutation *perms = m_element_kind == FIXED_CONNECTIVITY ?
        m_fixed_element_connectivity.begin_permutations(bucket_ordinal_of_lower_ranked_entity)
        :
        m_dynamic_element_connectivity.begin_permutations(bucket_ordinal_of_lower_ranked_entity);

    if (perms)
    {
        perms[elem_connectivity_ordinal] = permut;
    }
}

void Bucket::change_existing_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned face_connectivity_ordinal, stk::mesh::Permutation permut)
{
    stk::mesh::Permutation *perms = (m_face_kind == FIXED_CONNECTIVITY) ?
        m_fixed_face_connectivity.begin_permutations(bucket_ordinal_of_higher_ranked_entity)
        :
        m_dynamic_face_connectivity.begin_permutations(bucket_ordinal_of_higher_ranked_entity);

    if (perms)
    {
        perms[face_connectivity_ordinal] = permut;
    }
}

void Bucket::change_existing_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned edge_connectivity_ordinal, stk::mesh::Permutation permut)
{
    stk::mesh::Permutation *perms=nullptr;
    if (m_edge_kind == FIXED_CONNECTIVITY)
    {
        perms = m_fixed_edge_connectivity.begin_permutations(bucket_ordinal_of_higher_ranked_entity);
    }
    else
    {
        perms = const_cast<Permutation*>(m_dynamic_edge_connectivity.begin_permutations(bucket_ordinal_of_higher_ranked_entity));
    }

    if (perms)
    {
        perms[edge_connectivity_ordinal] = permut;
    }
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

  bucket.debug_dump(os);
  os << std::endl;

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

unsigned Bucket::get_others_begin_index(unsigned bucket_ordinal, EntityRank rank) const
{
  Entity const * const ents_begin = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const ents_end = m_dynamic_other_connectivity.end(bucket_ordinal);

  EntityRankLess cmp = {&m_mesh};
  Entity const *probe = ents_begin;
  probe = std::lower_bound(probe, ents_end, rank, cmp);

  return probe - ents_begin;
}

unsigned Bucket::get_others_end_index(unsigned bucket_ordinal, EntityRank rank) const
{
  Entity const * const ents_begin = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const ents_end = m_dynamic_other_connectivity.end(bucket_ordinal);

  EntityRankLess cmp = {&m_mesh};
  Entity const *probe = ents_begin;
  probe = std::upper_bound(probe, ents_end, rank, cmp);

  return probe - ents_begin;
}

unsigned Bucket::get_others_index_count(unsigned bucket_ordinal, EntityRank rank) const
{
  Entity const * const ents_begin = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const ents_end = m_dynamic_other_connectivity.end(bucket_ordinal);

  EntityRankLess cmp = {&m_mesh};
  Entity const *probe = ents_begin;
  Entity const *const saved_lower = std::lower_bound(probe, ents_end, rank, cmp);
  probe = std::upper_bound(probe, ents_end, rank, cmp);

  return probe - saved_lower;
}

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

void Bucket::grow_ngp_field_bucket_ids()
{
  const MetaData& meta = mesh().mesh_meta_data();
  const FieldVector& allFields = meta.get_fields();

  const unsigned oldNumFields = m_ngp_field_bucket_id.size();
  const unsigned newNumFields = allFields.size();
  STK_ThrowRequire(newNumFields >= oldNumFields);

  const unsigned numNewFields = newNumFields - oldNumFields;
  for (unsigned i = 0; i < numNewFields; ++i) {
    m_ngp_field_bucket_id.push_back(INVALID_BUCKET_ID);
    m_ngp_field_is_modified.push_back(false);
  }
}

void Bucket::set_ngp_field_bucket_id(unsigned fieldOrdinal, unsigned ngpFieldBucketId)
{
  STK_ThrowRequire(fieldOrdinal < m_ngp_field_bucket_id.size());
  m_ngp_field_bucket_id[fieldOrdinal] = ngpFieldBucketId;
  m_ngp_field_is_modified[fieldOrdinal] = false;
}

unsigned Bucket::get_ngp_field_bucket_id(unsigned fieldOrdinal) const
{
  STK_ThrowRequire(fieldOrdinal < m_ngp_field_bucket_id.size());
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

void Bucket::initialize_slot(unsigned ordinal, Entity entity)
{
  mark_for_modification();
  m_entities[ordinal]    = entity;
  if (mesh().is_valid(entity) && mesh().state(entity) == Unchanged) {
    mesh().set_state(entity, Modified);
  }
}

int Bucket::parallel_owner_rank(unsigned ordinal) const
{
  return m_mesh.parallel_owner_rank(m_entities[ordinal]);
}

void Bucket::reset_entity_location(Entity entity, unsigned to_ordinal, const FieldVector* fields)
{
  mark_for_modification();
  Bucket & from_bucket = mesh().bucket(entity);
  const unsigned from_ordinal = mesh().bucket_ordinal(entity);

  m_entities[to_ordinal]    = entity;

  mesh().set_mesh_index(entity, this, to_ordinal);

  m_mesh.copy_entity_fields_callback(m_entity_rank, m_bucket_id, to_ordinal,
                                     from_bucket.m_bucket_id, from_ordinal, fields);
}

void Bucket::reset_empty_space(const FieldVector & fields)
{
  m_mesh.reset_empty_field_data_callback(m_entity_rank, m_bucket_id, m_size, m_capacity, fields);
}

void Bucket::add_entity(Entity entity)
{
  STK_ThrowAssert(m_size < m_capacity);
  STK_ThrowAssert(!mesh().is_valid(m_entities[m_size]));
  STK_ThrowAssert(!mesh().is_valid(entity) || mesh().bucket_ptr(entity) == nullptr);
  STK_ThrowAssert(!mesh().is_valid(entity) || mesh().entity_rank(entity) == m_entity_rank);

  initialize_slot(m_size, entity);

  if (mesh().is_valid(entity)) {
    mesh().set_mesh_index(entity, this, m_size);
  }

  this->mesh().add_entity_callback(entity_rank(), bucket_id(), capacity(), m_size);
  ++m_size;

  AddEntityFunctor functor;
  process_all_connectivity(functor);
}

void
Bucket::grow_capacity()
{
  STK_ThrowAssert(m_capacity < std::numeric_limits<unsigned>::max()/2);
  m_capacity = std::min(2 * m_capacity, m_maxCapacity);
  m_dynamic_node_connectivity.increase_bucket_capacity(m_capacity);
  m_dynamic_edge_connectivity.increase_bucket_capacity(m_capacity);
  m_dynamic_face_connectivity.increase_bucket_capacity(m_capacity);
  m_dynamic_element_connectivity.increase_bucket_capacity(m_capacity);
  m_dynamic_other_connectivity.increase_bucket_capacity(m_capacity);
}

bool Bucket::destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id )
{
  const unsigned from_bucket_ordinal = mesh().bucket_ordinal(e_from);
  DestroyRelationFunctor functor(from_bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id));
  modify_connectivity(functor, m_mesh.entity_rank(e_to));

  return functor.m_modified;
}

bool Bucket::declare_relation(unsigned bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation )
{
  DeclareRelationFunctor functor(bucket_ordinal, e_to, ordinal, permutation);
  modify_connectivity(functor, m_mesh.entity_rank(e_to));

  return functor.m_modified;
}

bool Bucket::replace_relations(unsigned bucketOrdinal,
                               EntityRank rank,
                               unsigned numConnectivity,
                               const Entity* connectivity,
                               const ConnectivityOrdinal* ordinals,
                               const Permutation* permutations)
{
  if (numConnectivity > 0) {
    ReplaceRelationFunctor functor(bucketOrdinal, numConnectivity, connectivity,
                                   ordinals, permutations);
    modify_connectivity(functor, rank);
    return functor.m_modified;
  }
  return false;
}

void Bucket::remove_entity()
{
  STK_ThrowAssert(m_size > 0);

  mark_for_modification();
  mesh().remove_entity_field_data_callback(entity_rank(), bucket_id(), m_size-1);
  const unsigned bktOrdinal = m_size-1;
  --m_size;

  initialize_slot(m_size, Entity());

  RemoveEntityFunctor functor(bktOrdinal);
  process_all_connectivity(functor);
}

void Bucket::copy_entity(Entity entity)
{
  STK_ThrowAssert(m_size < m_capacity);
  STK_ThrowAssert(mesh().is_valid(entity));
  STK_ThrowAssert(!mesh().is_valid(m_entities[m_size]));
  STK_ThrowAssert(mesh().bucket_ptr(entity) != nullptr);
  STK_ThrowAssert(mesh().bucket_ptr(entity) != this);
  STK_ThrowAssert(mesh().entity_rank(entity) == m_entity_rank);

  mark_for_modification();
  Bucket* old_bucket = mesh().bucket_ptr(entity);
  const unsigned old_ordinal = mesh().bucket_ordinal(entity);

  this->mesh().add_entity_callback(entity_rank(), bucket_id(), capacity(), m_size);
  const unsigned newOrdinal = m_size;
  reset_entity_location(entity, newOrdinal);

  ++m_size;

  AddEntityFunctor functor;
  process_all_connectivity(functor);

  EntityRank endRank = static_cast<EntityRank>(mesh().mesh_meta_data().entity_rank_count());

  for(EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
    const unsigned numConn = old_bucket->num_connectivity(old_ordinal, rank);
    if (numConn > 0) {
      const Entity* conn = old_bucket->begin(old_ordinal, rank);
      const ConnectivityOrdinal* ordinals = old_bucket->begin_ordinals(old_ordinal, rank);
      const Permutation* perms = old_bucket->begin_permutations(old_ordinal, rank);
      replace_relations(newOrdinal, rank, numConn, conn, ordinals, perms);
    }
  }
}

void Bucket::overwrite_entity(unsigned to_ordinal, Entity entity, const FieldVector* fields)
{
  STK_ThrowAssert(to_ordinal < m_capacity);
  STK_ThrowAssert(mesh().is_valid(entity));
  STK_ThrowAssert(mesh().bucket_ptr(entity) != nullptr);
  STK_ThrowAssert(mesh().entity_rank(entity) == m_entity_rank);

  const MeshIndex from_index = m_mesh.mesh_index(entity);
  const Bucket* old_bucket = from_index.bucket;
  const unsigned old_ordinal = from_index.bucket_ordinal;
  reset_entity_location(entity, to_ordinal, fields);

  if (bucket_id() == old_bucket->bucket_id()) {
    ReplaceEntityFunctor functor(old_ordinal, to_ordinal);
    process_all_connectivity(functor);
  }
  else {
    ClearEntityFunctor functor(to_ordinal);
    EntityRank endRank = static_cast<EntityRank>(mesh().mesh_meta_data().entity_rank_count());

    for(EntityRank rank = stk::topology::NODE_RANK; rank<endRank; ++rank) {
      const unsigned numConn = old_bucket->num_connectivity(old_ordinal, rank);
      if (numConn > 0) {
        const Entity* conn = old_bucket->begin(old_ordinal, rank);
        const ConnectivityOrdinal* ordinals = old_bucket->begin_ordinals(old_ordinal, rank);
        const Permutation* perms = old_bucket->begin_permutations(old_ordinal, rank);
        replace_relations(to_ordinal, rank, numConn, conn, ordinals, perms);
      }
      else {
        modify_connectivity(functor, rank);
      }
    }
  }
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

void Bucket::check_size_invariant() const
{
}

void Bucket::debug_dump(std::ostream& out, unsigned ordinal) const
{
}

void Bucket::debug_check_for_invalid_connectivity_request(ConnectivityType const* type) const
{
#ifndef NDEBUG
    EntityRank rank = stk::topology::END_RANK;
    if (type == &m_node_kind) {
      rank = stk::topology::NODE_RANK;
    }
    else if (type == &m_edge_kind) {
      rank = stk::topology::EDGE_RANK;
    }
    else if (type == &m_face_kind) {
      rank = stk::topology::FACE_RANK;
    }
    else if (type == &m_element_kind) {
      rank = stk::topology::ELEMENT_RANK;
    }
    else {
      STK_ThrowAssert(false);
    }
    // Asking for connectivity between entities of equal rank is always invalid and ok to ask for
    // Asking for connectivity between for FACE_RANK in 2d is always invalid and ok to ask for
    bool isThisEntityAskingForConnectivityToItsOwnRank = entity_rank() == rank;
    bool isThisEntityAskingForFaceConnectivityOnTwoDimensionalMesh = rank == stk::topology::FACE_RANK && mesh().mesh_meta_data().spatial_dimension() == 2;
    STK_ThrowAssert( isThisEntityAskingForConnectivityToItsOwnRank || isThisEntityAskingForFaceConnectivityOnTwoDimensionalMesh);
#endif
}

} // namespace mesh
} // namespace stk
