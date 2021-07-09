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

// NOTE: This is not a stand alone header file. It is meant to help clean up BulkData.hpp
// This is currently just housing the mutli-lined inlined methods.

#ifndef BULKDATA_INLINED_METHODS_HPP
#define BULKDATA_INLINED_METHODS_HPP

// IWYU pragma: private, include "stk_mesh/base/BulkData.hpp"

namespace stk {
namespace mesh {

inline
unsigned BulkData::num_connectivity(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
}

inline
unsigned BulkData::find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  unsigned num_rels = mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
  ConnectivityOrdinal const *ords = mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);
  ThrowAssert(ords);

  unsigned i = 0;
  for (; i < num_rels; ++i)
  {
    if (ords[i] == ordinal)
      break;
  }
  return i;
}

inline
Entity const* BulkData::begin(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin(mesh_idx.bucket_ordinal, rank);
}

inline
Entity const* BulkData::begin_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_nodes(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edges(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_faces(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_elements(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_ordinals(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline
ConnectivityOrdinal const* BulkData::begin_node_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_edge_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_face_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_element_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_ordinals(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_permutations(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_permutations(mesh_idx.bucket_ordinal, rank);
}

inline
Permutation const* BulkData::begin_node_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_edge_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_face_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_element_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_permutations(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_nodes(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_edges(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_faces(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_elements(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end(mesh_idx.bucket_ordinal, rank);
}

inline
Entity const* BulkData::end_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_nodes(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edges(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_faces(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_elements(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_ordinals(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline
ConnectivityOrdinal const* BulkData::end_node_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_edge_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_face_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_element_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_ordinals(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_permutations(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_permutations(mesh_idx.bucket_ordinal, rank);
}

inline
Permutation const* BulkData::end_node_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_edge_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_face_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_element_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_permutations(mesh_idx.bucket_ordinal);
}

inline
bool BulkData::has_permutation(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  return bucket(entity).has_permutation(rank);
}

inline
bool BulkData::in_shared(Entity entity) const
{
  if (m_entitycomm[entity.local_offset()] != nullptr) {
    return m_entitycomm[entity.local_offset()]->isShared;
  }
  return false;
}

inline
bool BulkData::in_receive_ghost( EntityKey key ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>=AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], key) )
          return true;
  }
  return false;
}

inline
bool BulkData::in_receive_ghost( Entity entity ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>=AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], entity) )
          return true;
  }
  return false;
}

inline
bool BulkData::in_receive_custom_ghost( EntityKey key ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], key) )
          return true;
  }
  return false;
}

inline
bool BulkData::in_receive_ghost( const Ghosting & ghost , EntityKey key ) const
{
  const int owner_rank = parallel_owner_rank(get_entity(key));
  return in_ghost( ghost , key , owner_rank );
}

inline
bool BulkData::in_receive_ghost( const Ghosting & ghost , Entity entity ) const
{
  if (m_entitycomm[entity.local_offset()] == nullptr) {
    return false;
  }

  const int owner_rank = parallel_owner_rank(entity);
  if (owner_rank == parallel_rank()) {
    return false;
  }

  const EntityCommInfoVector& vec = m_entitycomm[entity.local_offset()]->comm_map;
  EntityCommInfoVector::const_iterator i = vec.begin();
  EntityCommInfoVector::const_iterator end = vec.end();
  for(; i!=end; ++i) {
    if (i->ghost_id == ghost.ordinal()) {
      return true;
    }
  }

  return false;
}

inline
bool BulkData::in_send_ghost( EntityKey key) const
{
    const int owner_rank = parallel_owner_rank(get_entity(key));
    for ( PairIterEntityComm ec = internal_entity_comm_map(key); ! ec.empty() ; ++ec )
    {
      if ( ec->ghost_id != 0 &&
           ec->proc     != owner_rank)
      {
        return true;
      }
    }
    return false;
}

inline
bool BulkData::in_send_ghost( Entity entity) const
{
    const int owner_rank = parallel_owner_rank(entity);
    for ( PairIterEntityComm ec = internal_entity_comm_map(entity); ! ec.empty() ; ++ec )
    {
      if ( ec->ghost_id != 0 &&
           ec->proc     != owner_rank)
      {
        return true;
      }
    }
    return false;
}

inline
void BulkData::internal_check_unpopulated_relations(Entity entity, EntityRank rank) const
{
#ifndef NDEBUG
  if (m_check_invalid_rels) {
    const MeshIndex &mesh_idx = mesh_index(entity);
    const Bucket &b = *mesh_idx.bucket;
    Bucket::size_type bucket_ord = mesh_idx.bucket_ordinal;
    ThrowAssertMsg(count_valid_connectivity(entity, rank) == b.num_connectivity(bucket_ord, rank),
                   count_valid_connectivity(entity,rank) << " = count_valid_connectivity("<<entity_key(entity)<<","<<rank<<") != b.num_connectivity("<<bucket_ord<<","<<rank<<") = " << b.num_connectivity(bucket_ord,rank);
                  );

  }
#endif
}

////////////////

inline void BulkData::copy_entity_fields(Entity src, Entity dst)
{
  if (src == dst) return;

  //TODO fix const correctness for src
  MeshIndex & src_mesh_idx = mesh_index(src);
  MeshIndex & dst_mesh_idx = mesh_index(dst);

  copy_entity_fields_callback(dst_mesh_idx.bucket->entity_rank(),
                              dst_mesh_idx.bucket->bucket_id(),
                              dst_mesh_idx.bucket_ordinal,
                              src_mesh_idx.bucket->bucket_id(),
                              src_mesh_idx.bucket_ordinal);
}

inline bool BulkData::relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id )
{
  bool found = false;
  Entity const * rel_entity_it = bucket(entity).begin(bucket_ordinal(entity),subcell_rank);
  const unsigned num_rel = bucket(entity).num_connectivity(bucket_ordinal(entity),subcell_rank);
  ConnectivityOrdinal const * rel_ord_it = bucket(entity).begin_ordinals(bucket_ordinal(entity),subcell_rank);

  for (unsigned i=0 ; i < num_rel ; ++i) {
    if (rel_ord_it[i] == static_cast<ConnectivityOrdinal>(subcell_id) && is_valid(rel_entity_it[i])) {
      found = true;
      break;
    }
  }

  return found;
}

inline VolatileFastSharedCommMapOneRank const& BulkData::volatile_fast_shared_comm_map(EntityRank rank) const
{
  ThrowAssert(this->in_synchronized_state());
  ThrowAssertMsg(rank < stk::topology::ELEMENT_RANK, "Cannot shared entities of rank: " << rank);
  if (m_volatile_fast_shared_comm_map_sync_count < synchronized_count()) {
    internal_update_fast_comm_maps();
  }
  return m_volatile_fast_shared_comm_map[rank];
}

inline Part& BulkData::ghosting_part(const Ghosting& ghosting) const
{
  ThrowRequireMsg(ghosting.ordinal() < m_ghost_parts.size(), "BulkData::ghosting_part ERROR, no part corresponds to ghosting with name="<<ghosting.name()<<" and ordinal="<<ghosting.ordinal());
  return *m_ghost_parts[ghosting.ordinal()];
}

inline bool BulkData::in_index_range(Entity entity) const
{
  return entity.local_offset() < m_entity_keys.size();
}

inline bool BulkData::is_valid(Entity entity) const
{
  ThrowAssertMsg(in_index_range(entity),
                 "Error in stk::mesh::BulkData::is_valid, entity not in index range. "
                 " entity.local_offset()="<<entity.local_offset()<<", valid range is < "
                 << get_size_of_entity_index_space());
  return !m_meshModification.is_entity_deleted(entity.local_offset());
}

inline const MeshIndex& BulkData::mesh_index(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_mesh_indexes[entity.local_offset()];
}

inline MeshIndex& BulkData::mesh_index(Entity entity)
{
  ThrowAssert(entity.local_offset() > 0); // setter check due to non-const

  return m_mesh_indexes[entity.local_offset()];
}

inline EntityId BulkData::identifier(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_entity_keys[entity.local_offset()].id();
}

inline EntityRank BulkData::entity_rank(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_entity_keys[entity.local_offset()].rank();
}

inline EntityKey BulkData::entity_key(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_entity_keys[entity.local_offset()];
}

inline EntityState BulkData::state(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif
  return m_meshModification.get_entity_state(entity.local_offset());
}

inline void BulkData::internal_mark_entity(Entity entity, EntitySharing sharedType)
{
    m_mark_entity[entity.local_offset()] = sharedType;
}

inline BulkData::EntitySharing BulkData::internal_is_entity_marked(Entity entity) const
{
    return m_mark_entity[entity.local_offset()];
}

inline bool BulkData::internal_add_node_sharing_called() const
{
  return m_add_node_sharing_called;
}

inline Bucket & BulkData::bucket(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return *mesh_index(entity).bucket;
}

inline Bucket * BulkData::bucket_ptr(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return mesh_index(entity).bucket;
}

inline Bucket::size_type BulkData::bucket_ordinal(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return mesh_index(entity).bucket_ordinal;
}

inline int BulkData::parallel_owner_rank(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_owner[entity.local_offset()];
}

inline unsigned BulkData::local_id(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_local_ids[entity.local_offset()];
}

#ifdef SIERRA_MIGRATION
inline BulkData::FmwkId BulkData::global_id(stk::mesh::Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_fmwk_global_ids[entity.local_offset()];
}

inline const RelationVector& BulkData::aux_relations(Entity entity) const
{
  ThrowAssert(m_add_fmwk_data);
  ThrowAssert(entity.local_offset() > 0);

  if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
    m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
  }
  return *m_fmwk_aux_relations[entity.local_offset()];
}

inline RelationVector& BulkData::aux_relations(Entity entity)
{
  ThrowAssert(m_add_fmwk_data);
  ThrowAssert(entity.local_offset() > 0);

  if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
    m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
  }
  return *m_fmwk_aux_relations[entity.local_offset()];
}

inline void BulkData::set_global_id(stk::mesh::Entity entity, BulkData::FmwkId id)
{
  ThrowAssert(entity.local_offset() > 0);

  m_modSummary.track_set_global_id(entity, id);

  m_fmwk_global_ids[entity.local_offset()] = id;
}

inline RelationIterator BulkData::internal_begin_relation(Entity entity, const Relation::RelationType relation_type) const
{
  ThrowAssert(m_add_fmwk_data);
  if (impl::internal_is_handled_generically(relation_type)) {
    ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
    return RelationIterator();
  }
  else {
    return aux_relations(entity).begin();
  }
}

inline RelationIterator BulkData::internal_end_relation(Entity entity, const Relation::RelationType relation_type) const
{
  ThrowAssert(m_add_fmwk_data);
  if (impl::internal_is_handled_generically(relation_type)) {
    ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
    return RelationIterator();
  }
  else {
    return aux_relations(entity).end();
  }
}

inline void BulkData::compress_relation_capacity(Entity entity)
{
  RelationVector &rels = aux_relations(entity);
  RelationVector tmp(rels);
  tmp.swap(rels);
}
#endif

inline void BulkData::set_mesh_index(Entity entity, Bucket * in_bucket, Bucket::size_type ordinal )
{
  ThrowAssert(entity.local_offset() > 0);

  if (in_bucket != NULL) {
    ThrowAssertMsg(in_bucket->size() >= ordinal, "Detected bad bucket/ordinal.");
  }
  MeshIndex &mesh_idx = mesh_index(entity);
  mesh_idx.bucket = in_bucket;
  mesh_idx.bucket_ordinal = ordinal;
}

inline void BulkData::set_entity_key(Entity entity, EntityKey key)
{
  ThrowAssert(entity.local_offset() > 0);

  m_entity_keys[entity.local_offset()] = key;
}

inline void BulkData::set_state(Entity entity, EntityState entity_state)
{
  ThrowAssert(entity.local_offset() > 0);

  m_meshModification.set_entity_state(entity.local_offset(), entity_state);
  m_mark_entity[entity.local_offset()] = NOT_MARKED;
}

inline void BulkData::set_local_id(Entity entity, unsigned id)
{
  ThrowAssert(entity.local_offset() > 0);

  m_local_ids[entity.local_offset()] = id;
}

inline void BulkData::log_created_parallel_copy(Entity entity)
{
  if (state(entity) == Created) {
    set_state(entity, Modified);
  }
}

inline bool BulkData::is_valid_connectivity(Entity entity, EntityRank rank) const
{
  if (!is_valid(entity)) return false;
  if (bucket_ptr(entity) == NULL) return false;
  internal_check_unpopulated_relations(entity, rank);
  return true;
}

}
}

#endif
