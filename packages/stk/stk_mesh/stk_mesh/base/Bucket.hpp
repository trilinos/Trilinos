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

#ifndef stk_mesh_base_Bucket_hpp
#define stk_mesh_base_Bucket_hpp

#include <stddef.h>                     // for size_t, NULL
#include <algorithm>                    // for lower_bound
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Part.hpp>       // for contains_ordinal, Part
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/SparseConnectivity.hpp>
#include <stk_mesh/base/ConnectedTopologyNodes.hpp>
#include <stk_mesh/base/ConnectedSparseNodes.hpp>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowAssert, etc
#include <stk_util/util/StridedArray.hpp>
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class DeviceMesh; } }
namespace stk { namespace mesh { namespace impl { class BucketRepository; } } }
namespace stk { namespace mesh { namespace impl { class Partition; } } }

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Print the \ref stk::mesh::Part "part" names
 *          for which this bucket is a subset.
 */
std::ostream & operator << ( std::ostream & , const Bucket & );

/** \brief  Print the parts and entities of this bucket */
std::ostream &
print( std::ostream & , const std::string & indent , const Bucket & );

// The part count and parts are equal
bool raw_part_equal( const unsigned * lhs , const unsigned * rhs );

#define CONN_TYPE(conn_type, fixed_conn, sparse_conn) \
  ((conn_type == FIXED_CONNECTIVITY) ? fixed_conn : sparse_conn)

#define RANK_SWITCH(rank, begin_or_end, item_type, bucket_ordinal)    \
                                                        \
  switch(rank) {                                        \
  case stk::topology::NODE_RANK:                        \
       switch(m_node_kind) {                                 \
       case FIXED_CONNECTIVITY: return m_topoNodes.begin_or_end##item_type(bucket_ordinal);     \
       case DYNAMIC_CONNECTIVITY: return m_sparseNodes.begin_or_end##item_type(bucket_ordinal); \
       default: return nullptr;                           \
       }                                                   \
  default:                                                  \
    return m_sparse_connectivity.begin_or_end##item_type(m_entities[bucket_ordinal], rank); \
}

//----------------------------------------------------------------------
/** \brief  A container for a homogeneous collection of
 *          \ref stk::mesh::Entity "entities".
 *
 *  The entities are homogeneous in that they are of the same rank and topology
 *  and are members of the same parts.
 */
class Bucket
{

public:
  typedef size_t size_type;

  stk::topology topology() const { return m_topology; }

  void parent_topology( EntityRank parent_rank, std::vector<stk::topology> & parent_topologies) const;

  bool owned() const { return m_owned; }
  bool shared() const { return m_shared; }
  bool in_aura() const { return m_aura; }

  //--------------------------------
  // Container-like types and methods:

  typedef const Entity * iterator;

  /** \brief Beginning of the bucket */
  inline iterator begin() const { return &m_entities[0]; }

  /** \brief End of the bucket */
  inline iterator end() const { return &m_entities[0] + m_size; }

  /** \brief  Number of entities associated with this bucket */
  size_type size() const { return m_size ; }

  /** \brief  Capacity of this bucket */
  size_t capacity() const { return m_capacity ; }

  /** \brief  Query the i^th entity */
  Entity operator[] ( size_t i ) const {
    ThrowAssertMsg( i < m_entities.size(), "Index " << i << " is out of bounds");
    return m_entities[i];
  }

  ConnectivityType connectivity_type(EntityRank rank) const;

  //--------------------------------
  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & mesh() const { return m_mesh ; }

  /** \brief  Type of entities in this bucket */
  EntityRank entity_rank() const { return m_entity_rank ; }

  unsigned bucket_id() const { return m_bucket_id; }

  /** \brief  This bucket is a subset of these \ref stk::mesh::Part "parts"
   * WARNING: if the bucket is deleted, this reference will no longer
   * point to valid memory.
   */
  const PartVector& supersets() const { return m_parts; }
  void supersets( OrdinalVector & ) const ;

  //--------------------------------
  /** \brief  Bucket is a subset of the given part */
  bool member( const Part & part) const
  {
    return member(part.mesh_meta_data_ordinal());
  }

  bool member( PartOrdinal partOrdinal ) const
  {
    return std::binary_search(m_partOrdsBeginEnd.first, m_partOrdsBeginEnd.second, partOrdinal);
  }

  /** \brief  Bucket is a subset of all of the given parts */
  bool member_all( const PartVector & ) const ;
  bool member_all( const OrdinalVector & ) const ;

  /** \brief  Bucket is a subset of any of the given parts */
  bool member_any( const PartVector & ) const ;
  bool member_any( const OrdinalVector & ) const ;

  //--------------------------------
  /** Query bucket's supersets' ordinals. */

  std::pair<const unsigned *, const unsigned *>
  superset_part_ordinals() const { return m_partOrdsBeginEnd; }

#ifndef DOXYGEN_COMPILE
  const unsigned * key() const { return m_key.data() ; }
#endif /* DOXYGEN_COMPILE */

  /** \brief  The allocation size, in bytes, of this bucket */
  unsigned allocation_size() const { return 0 ; }

  bool is_empty() const { return size() == 0; }

  impl::Partition *getPartition() const { return m_partition; }

  unsigned char* field_data_location(const FieldBase& field) const;

  bool field_data_is_allocated(const FieldBase& field) const;

  ///
  /// Entity member functions are moved here:
  ///

  int parallel_owner_rank(unsigned ordinal) const;

  //generic rank connectivity calls
  Entity const* begin(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, _connectivity, bucket_ordinal) }
  ConnectivityOrdinal const* begin_ordinals(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, _ordinals, bucket_ordinal) }
  Permutation const* begin_permutations(unsigned bucket_ordinal, EntityRank rank) const
  { return m_sparse_connectivity.begin_permutations(m_entities[bucket_ordinal], rank); }

  Entity const* end(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, _connectivity, bucket_ordinal) }
  ConnectivityOrdinal const* end_ordinals(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, _ordinals, bucket_ordinal) }
  Permutation const* end_permutations(unsigned bucket_ordinal, EntityRank rank) const
  { return m_sparse_connectivity.end_permutations(m_entities[bucket_ordinal], rank); }

  unsigned num_connectivity(unsigned bucket_ordinal, EntityRank rank) const;

  Entity const* begin_nodes(unsigned bucket_ordinal) const
  { return CONN_TYPE(m_node_kind, m_topoNodes.begin_connectivity(bucket_ordinal), m_sparseNodes.begin_connectivity(bucket_ordinal)); }
  Entity const* begin_edges(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_connectivity(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  Entity const* begin_faces(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_connectivity(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  Entity const* begin_elements(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_connectivity(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  ConnectivityOrdinal const* begin_node_ordinals(unsigned bucket_ordinal) const
  { return CONN_TYPE(m_node_kind, m_topoNodes.begin_ordinals(bucket_ordinal), m_sparseNodes.begin_ordinals(bucket_ordinal)); }
  ConnectivityOrdinal const* begin_edge_ordinals(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_ordinals(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  ConnectivityOrdinal const* begin_face_ordinals(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_ordinals(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  ConnectivityOrdinal const* begin_element_ordinals(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_ordinals(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  Permutation const* begin_node_permutations(unsigned bucket_ordinal) const { return nullptr; }
  Permutation const* begin_edge_permutations(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_permutations(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  Permutation const* begin_face_permutations(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_permutations(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  Permutation const* begin_element_permutations(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.begin_permutations(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  unsigned num_nodes(unsigned bucket_ordinal) const
  { return CONN_TYPE(m_node_kind, m_topoNodes.num_nodes_per_entity(bucket_ordinal), m_sparseNodes.num_nodes_per_entity(bucket_ordinal)); }
  unsigned num_edges(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.num_connectivity(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  unsigned num_faces(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.num_connectivity(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  unsigned num_elements(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.num_connectivity(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  Entity const* end_nodes(unsigned bucket_ordinal) const
  { return CONN_TYPE(m_node_kind, m_topoNodes.end_connectivity(bucket_ordinal), m_sparseNodes.end_connectivity(bucket_ordinal)); }
  Entity const* end_edges(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_connectivity(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  Entity const* end_faces(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_connectivity(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  Entity const* end_elements(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_connectivity(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  ConnectivityOrdinal const* end_node_ordinals(unsigned bucket_ordinal) const
  { return CONN_TYPE(m_node_kind, m_topoNodes.end_ordinals(bucket_ordinal), m_sparseNodes.end_ordinals(bucket_ordinal)); }
  ConnectivityOrdinal const* end_edge_ordinals(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_ordinals(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  ConnectivityOrdinal const* end_face_ordinals(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_ordinals(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  ConnectivityOrdinal const* end_element_ordinals(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_ordinals(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  Permutation const* end_node_permutations(unsigned bucket_ordinal) const { return nullptr; }
  Permutation const* end_edge_permutations(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_permutations(m_entities[bucket_ordinal], stk::topology::EDGE_RANK); }
  Permutation const* end_face_permutations(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_permutations(m_entities[bucket_ordinal], stk::topology::FACE_RANK); }
  Permutation const* end_element_permutations(unsigned bucket_ordinal) const
  { return m_sparse_connectivity.end_permutations(m_entities[bucket_ordinal], stk::topology::ELEM_RANK); }

  bool has_permutation(EntityRank rank) const;

  void debug_dump(std::ostream& out, unsigned ordinal = -1u) const;

  /* NGP Bucket methods */
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after November 2021
  STK_DEPRECATED unsigned get_num_nodes_per_entity() const { return topology().num_nodes(); }
#endif

  ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
    return ConnectedEntities(begin(offsetIntoBucket, connectedRank),
                             num_connectivity(offsetIntoBucket, connectedRank));
  }

  ConnectedOrdinals get_connected_ordinals(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const {
    return ConnectedOrdinals(begin_ordinals(offsetIntoBucket, connectedRank),
                             num_connectivity(offsetIntoBucket, connectedRank));
  }

  ConnectedNodes get_nodes(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::NODE_RANK);
  }

  ConnectedEntities get_edges(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::EDGE_RANK);
  }

  ConnectedEntities get_faces(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::FACE_RANK);
  }

  ConnectedEntities get_elements(unsigned offsetIntoBucket) const {
    return get_connected_entities(offsetIntoBucket, stk::topology::ELEM_RANK);
  }

  stk::mesh::Entity host_get_entity(unsigned offsetIntoBucket) const {
    return (*this)[offsetIntoBucket];
  }

  void set_ngp_field_bucket_id(unsigned fieldOrdinal, unsigned ngpFieldBucketId);
  unsigned get_ngp_field_bucket_id(unsigned fieldOrdinal) const;
  unsigned get_ngp_field_bucket_is_modified(unsigned fieldOrdinal) const;

  void reset_part_ord_begin_end();

  void reset_bucket_key(const OrdinalVector& newPartOrdinals);

  void reset_bucket_parts(const OrdinalVector& newPartOrdinals);

protected:
  void change_connected_nodes(unsigned bucket_ordinal, Entity* new_nodes);
  void change_existing_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity, ConnectivityOrdinal elem_connectivity_ordinal, Permutation permut);
  void change_existing_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, ConnectivityOrdinal edge_connectivity_ordinal, Permutation permut);
  void change_existing_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, ConnectivityOrdinal face_connectivity_ordinal, Permutation permut);
  virtual ~Bucket();

protected:
  bool destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id );

  bool declare_relation(unsigned bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation);

  void supersets( PartVector & ) const ;

  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & bulk_data() const { return mesh(); }

  bool is_modified() const { return m_is_modified; }
  unsigned ngp_bucket_id() const { return m_ngp_bucket_id; }
  void set_ngp_bucket_id(unsigned ngpBucketId) {
    m_ngp_bucket_id = ngpBucketId;
    m_is_modified = false;
  }

  void mark_for_modification()
  {
  #ifdef STK_USE_DEVICE_MESH
    m_is_modified = true;
    std::fill(m_ngp_field_is_modified.begin(), m_ngp_field_is_modified.end(), true);
  #endif
  }

  void initialize_ngp_field_bucket_ids();

  Bucket();

  Bucket( const Bucket & );
  Bucket & operator = ( const Bucket & );

  Bucket( BulkData & arg_mesh ,
          EntityRank arg_entity_rank,
          const std::vector<unsigned> & arg_key,
          size_t arg_capacity,
          unsigned bucket_id);

  const std::vector<unsigned> & key_vector() const { return m_key; }

  // Add a new entity to end of bucket
  void add_entity(Entity entity = Entity());

  // Remove an entity from the end of this bucket
  void remove_entity();

  // Copy an existing entity to the end of this bucket
  void copy_entity(Entity entity);

  // overwrites existing entity at ordinal with entity
  // bucket[to_ordinal] = entity;
  // whatever was there before is lost
  //  With optional fields argument only copy listed fields
  void overwrite_entity(unsigned to_ordinal, Entity entity, const std::vector<FieldBase*>* fields = nullptr);

  void initialize_slot(unsigned ordinal, Entity entity);
  //  Optional fields argument, only copy listed fields
  void reset_entity_location(Entity entity, unsigned to_ordinal, const std::vector<FieldBase*>* fields = nullptr);

  friend class BulkData;
  friend class impl::Partition;
  friend class impl::BucketRepository;
  friend class DeviceMesh;

  BulkData             & m_mesh ;        // Where this bucket resides
  SparseConnectivity   & m_sparse_connectivity;
  const EntityRank       m_entity_rank ; // Type of entities for this bucket
  std::vector<unsigned>  m_key ;         // REFACTOR
  std::pair<const unsigned*,const unsigned*> m_partOrdsBeginEnd;
  stk::topology          m_topology ;    // The topology of this bucket
  const size_t           m_capacity ;    // Capacity for entities
  size_type              m_size ;        // Number of entities
  unsigned               m_bucket_id;
  unsigned               m_ngp_bucket_id;
  bool                   m_is_modified;
  PartVector             m_parts;

  // Entity data
  std::vector<Entity>    m_entities;

  impl::Partition    *m_partition;

  ConnectivityType m_node_kind;
  ConnectivityType m_edge_kind;
  ConnectivityType m_face_kind;
  ConnectedTopologyNodes m_topoNodes;
  ConnectedSparseNodes   m_sparseNodes;

  bool m_owned;
  bool m_shared;
  bool m_aura;

  std::vector<unsigned> m_ngp_field_bucket_id;
  std::vector<bool> m_ngp_field_is_modified;
};

/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part" by partID
 */
inline
bool has_superset( const Bucket & bucket,  const unsigned & ordinal )
{
  return bucket.member(ordinal);
}

//----------------------------------------------------------------------
/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part"
 */

inline
bool has_superset( const Bucket & bucket ,  const Part & p )
{
  return bucket.member(p.mesh_meta_data_ordinal());
}

/** \brief  Is this bucket a subset of all of the given
 *          \ref stk::mesh::Part "parts"
 */
bool has_superset( const Bucket & bucket , const PartVector & parts );


struct BucketLess {
  bool operator()( const Bucket * lhs_bucket , const unsigned * rhs ) const ;
  bool operator()( const unsigned * lhs , const Bucket * rhs_bucket ) const ;
};

inline
BucketVector::iterator
lower_bound( BucketVector & v , const unsigned * key )
{ return std::lower_bound( v.begin() , v.end() , key , BucketLess() ); }

struct BucketIdComparator
{
  bool operator()(Bucket const* lhs, Bucket const* rhs) const
  {
    ThrowAssertMsg(lhs->entity_rank() == rhs->entity_rank(), "Cannot compare buckets of different rank");
    return lhs->bucket_id() < rhs->bucket_id();
  }

  bool operator()(unsigned bucket_id, Bucket const* rhs) const
  {
    return bucket_id < rhs->bucket_id();
  }

  bool operator()(Bucket const* lhs, unsigned bucket_id) const
  {
    return lhs->bucket_id() < bucket_id;
  }
};

/** \} */

inline
bool Bucket::member_all( const OrdinalVector& parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const OrdinalVector::const_iterator ip_end = parts.end();
        OrdinalVector::const_iterator ip     = parts.begin() ;

  bool result_all = true ;

  for ( ; result_all && ip_end != ip ; ++ip ) {
    const unsigned ord = *ip;
    result_all = contains_ordinal(i_beg, i_end, ord);
  }
  return result_all ;
}

inline
unsigned Bucket::num_connectivity(unsigned bucket_ordinal, EntityRank rank) const
{
  switch(rank) {
  case stk::topology::NODE_RANK:  return num_nodes(bucket_ordinal);
  default:
    return m_sparse_connectivity.num_connectivity(m_entities[bucket_ordinal], rank);
  }
}

inline
bool Bucket::has_permutation(EntityRank rank) const
{
  return SparseConnectivity::has_permutation(entity_rank(), rank);
}

inline
ConnectivityType Bucket::connectivity_type(EntityRank rank) const
{
  switch(rank) {
  case stk::topology::NODE_RANK:
    return m_node_kind;
  case stk::topology::EDGE_RANK:
    return m_edge_kind;
  case stk::topology::FACE_RANK:
    return m_face_kind;
  default:
    return DYNAMIC_CONNECTIVITY;
  }
}

typedef Bucket::iterator BucketIterator;

} // namespace mesh
} // namespace stk

#endif 

