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
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/BucketConnDynamic.hpp>
#include <stk_mesh/base/BucketConnectivity.hpp>  // for BucketConnectivity
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Part.hpp>       // for contains_ordinal, Part
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for STK_ThrowAssert, etc
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class DeviceMesh; } }
namespace stk { namespace mesh { namespace impl { class BucketRepository; } } }
namespace stk { namespace mesh { namespace impl { class Partition; } } }
namespace stk { namespace mesh { namespace impl { struct OverwriteEntityFunctor; } } }

namespace stk {
namespace mesh {

unsigned get_default_bucket_capacity();
unsigned get_default_initial_bucket_capacity();
unsigned get_default_maximum_bucket_capacity();

constexpr
inline
bool does_rank_have_valid_permutations(stk::mesh::EntityRank rank)
{
    return rank > stk::topology::NODE_RANK && rank < stk::topology::CONSTRAINT_RANK;
}

constexpr
inline
bool should_store_permutations(EntityRank fromRank, EntityRank toRank)
{
    return does_rank_have_valid_permutations(fromRank)
        && does_rank_have_valid_permutations(toRank);
}

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

#define CONNECTIVITY_TYPE_SWITCH(entity_kind, fixed_func_sig, dynamic_func_sig, check_invalid) \
  switch(entity_kind) {                                                 \
  case FIXED_CONNECTIVITY:                                              \
    return fixed_func_sig;                                              \
  case DYNAMIC_CONNECTIVITY:                                            \
    return dynamic_func_sig;                                            \
  default:                                                              \
    if (check_invalid) {                                                \
      check_for_invalid_connectivity_request(&entity_kind);             \
    }                                                                   \
    return 0;                                                           \
  }

#define RANK_SWITCH(rank, begin_or_end, postfix, bucket_ordinal)  \
                                                        \
  switch(rank) {                                                          \
  case stk::topology::NODE_RANK:    return begin_or_end##_node##postfix(bucket_ordinal); \
  case stk::topology::EDGE_RANK:    return begin_or_end##_edge##postfix(bucket_ordinal); \
  case stk::topology::FACE_RANK:    return begin_or_end##_face##postfix(bucket_ordinal); \
  case stk::topology::ELEMENT_RANK: return begin_or_end##_element##postfix(bucket_ordinal); \
  default:                                                              \
    return begin_other##postfix(bucket_ordinal) + get_others_##begin_or_end##_index(bucket_ordinal, rank);  \
}

//----------------------------------------------------------------------
/** \brief  A container for the connectivity for a homogeneous collection of
 *          \ref stk::mesh::Entity "entities".
 *
 *  The entities are homogeneous in that they are of the same entity type
 *  and are members of the same of parts.
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

  size_t memory_size_in_bytes() const;

  /** \brief  Capacity of this bucket */
  size_t capacity() const { return m_capacity; }

  /** \brief  Query the i^th entity */
  Entity operator[] ( size_t i ) const {
    STK_ThrowAssertMsg( i < m_entities.size(), "Index " << i << " is out of bounds");
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
    return contains_ordinal(m_partOrdsBeginEnd.first, m_partOrdsBeginEnd.second, partOrdinal);
  }

  /** \brief  Bucket is a subset of all of the given parts */
  bool member_all( const PartVector & ) const ;
  bool member_all( const OrdinalVector & ) const ;

  /** \brief  Bucket is a subset of any of the given parts */
  template<typename PARTVECTOR>
  bool member_any(const PARTVECTOR & parts) const
  {
    for(const Part* part : parts) {
      if (member(*part)) {
        return true;
      }
    }
    return false ;
  }

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

  void check_size_invariant() const;

  //generic rank connectivity calls
  Entity const* begin(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, s, bucket_ordinal) }
  ConnectivityOrdinal const* begin_ordinals(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, _ordinals, bucket_ordinal) }
  Permutation const* begin_permutations(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, _permutations, bucket_ordinal) }

  Entity const* end(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, s, bucket_ordinal) }
  ConnectivityOrdinal const* end_ordinals(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, _ordinals, bucket_ordinal) }
  Permutation const* end_permutations(unsigned bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, _permutations, bucket_ordinal) }

  unsigned num_connectivity(unsigned bucket_ordinal, EntityRank rank) const;

  Entity const* begin_nodes(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.begin(bucket_ordinal), m_dynamic_node_connectivity.begin(bucket_ordinal), true) }
  Entity const* begin_edges(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.begin(bucket_ordinal), m_dynamic_edge_connectivity.begin(bucket_ordinal), true) }
  Entity const* begin_faces(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.begin(bucket_ordinal), m_dynamic_face_connectivity.begin(bucket_ordinal), true) }
  Entity const* begin_elements(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.begin(bucket_ordinal), m_dynamic_element_connectivity.begin(bucket_ordinal), true) }

  ConnectivityOrdinal const* begin_node_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_node_connectivity.begin_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* begin_edge_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_edge_connectivity.begin_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* begin_face_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_face_connectivity.begin_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* begin_element_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_element_connectivity.begin_ordinals(bucket_ordinal), true) }

  Permutation const* begin_node_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.begin_permutations(bucket_ordinal), m_dynamic_node_connectivity.begin_permutations(bucket_ordinal), true) }
  Permutation const* begin_edge_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.begin_permutations(bucket_ordinal), m_dynamic_edge_connectivity.begin_permutations(bucket_ordinal), true) }
  Permutation const* begin_face_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.begin_permutations(bucket_ordinal), m_dynamic_face_connectivity.begin_permutations(bucket_ordinal), true) }
  Permutation const* begin_element_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.begin_permutations(bucket_ordinal), m_dynamic_element_connectivity.begin_permutations(bucket_ordinal), true) }

  unsigned num_nodes(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.num_connectivity(bucket_ordinal), m_dynamic_node_connectivity.num_connectivity(bucket_ordinal), false) }
  unsigned num_edges(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.num_connectivity(bucket_ordinal), m_dynamic_edge_connectivity.num_connectivity(bucket_ordinal), false) }
  unsigned num_faces(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.num_connectivity(bucket_ordinal), m_dynamic_face_connectivity.num_connectivity(bucket_ordinal), false) }
  unsigned num_elements(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.num_connectivity(bucket_ordinal), m_dynamic_element_connectivity.num_connectivity(bucket_ordinal), false) }

  Entity const* end_nodes(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.end(bucket_ordinal), m_dynamic_node_connectivity.end(bucket_ordinal), true) }
  Entity const* end_edges(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.end(bucket_ordinal), m_dynamic_edge_connectivity.end(bucket_ordinal), true) }
  Entity const* end_faces(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.end(bucket_ordinal), m_dynamic_face_connectivity.end(bucket_ordinal), true) }
  Entity const* end_elements(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.end(bucket_ordinal), m_dynamic_element_connectivity.end(bucket_ordinal), true) }

  ConnectivityOrdinal const* end_node_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.end_ordinals(bucket_ordinal), m_dynamic_node_connectivity.end_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* end_edge_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.end_ordinals(bucket_ordinal), m_dynamic_edge_connectivity.end_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* end_face_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.end_ordinals(bucket_ordinal), m_dynamic_face_connectivity.end_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* end_element_ordinals(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.end_ordinals(bucket_ordinal), m_dynamic_element_connectivity.end_ordinals(bucket_ordinal), true) }

  Permutation const* end_node_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.end_permutations(bucket_ordinal), m_dynamic_node_connectivity.end_permutations(bucket_ordinal), true) }
  Permutation const* end_edge_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.end_permutations(bucket_ordinal), m_dynamic_edge_connectivity.end_permutations(bucket_ordinal), true) }
  Permutation const* end_face_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.end_permutations(bucket_ordinal), m_dynamic_face_connectivity.end_permutations(bucket_ordinal), true) }
  Permutation const* end_element_permutations(unsigned bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.end_permutations(bucket_ordinal), m_dynamic_element_connectivity.end_permutations(bucket_ordinal), true) }

  bool has_permutation(EntityRank rank) const;

  void debug_dump(std::ostream& out, unsigned ordinal = -1u) const;

  /* NGP Bucket methods */

  ConnectedEntities get_connected_entities(unsigned offsetIntoBucket, stk::mesh::EntityRank connectedRank) const
  {
    STK_ThrowAssertMsg(offsetIntoBucket < size(),"Bucket::get_connected_entities offsetIntoBucket="<<offsetIntoBucket<<", size()="<<size());
    ConnectivityType connType = connectivity_type(connectedRank);
    return connType == FIXED_CONNECTIVITY ? get_fixed_connectivity(offsetIntoBucket, connectedRank)
                                          : get_dynamic_connectivity(offsetIntoBucket, connectedRank);
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

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after 2024/06/26
  STK_DEPRECATED
  stk::mesh::Entity host_get_entity(unsigned offsetIntoBucket) const {
    return (*this)[offsetIntoBucket];
  }
#endif

  void set_ngp_field_bucket_id(unsigned fieldOrdinal, unsigned ngpFieldBucketId);
  unsigned get_ngp_field_bucket_id(unsigned fieldOrdinal) const;
  unsigned get_ngp_field_bucket_is_modified(unsigned fieldOrdinal) const;

  void reset_part_ord_begin_end();

  void reset_bucket_key(const OrdinalVector& newPartOrdinals);

  void reset_bucket_parts(const OrdinalVector& newPartOrdinals);

protected:
  void change_existing_connectivity(unsigned bucket_ordinal, stk::mesh::Entity* new_nodes);
  void change_existing_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity, unsigned elem_connectivity_ordinal, stk::mesh::Permutation permut);
  void change_existing_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned edge_connectivity_ordinal, stk::mesh::Permutation permut);
  void change_existing_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned face_connectivity_ordinal, stk::mesh::Permutation permut);
  virtual ~Bucket();

private:

  void grow_capacity();

  bool destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id );

  bool declare_relation(unsigned bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation);

  bool replace_relations(unsigned bucketOrdinal,
                         EntityRank rank,
                         unsigned numConnectivity,
                         const Entity* connectivity,
                         const ConnectivityOrdinal* ordinals,
                         const Permutation* permutations);

  // The following *_other* functions should not be made available externally, in
  // order to avoid external confusion with "constraint" and "other" connectivities.
  // They are currently used within this class to provide connectivities
  // externally through another interface.
  Entity const* begin_others(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.begin(bucket_ordinal);
  }
  ConnectivityOrdinal const* begin_other_ordinals(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.begin_ordinals(bucket_ordinal);
  }
  Permutation const* begin_other_permutations(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.begin_permutations(bucket_ordinal);
  }
  unsigned num_other(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.num_connectivity(bucket_ordinal);
  }
  Entity const* end_others(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.end(bucket_ordinal);
  }
  ConnectivityOrdinal const* end_other_ordinals(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.end_ordinals(bucket_ordinal);
  }
  Permutation const* end_other_permutations(unsigned bucket_ordinal) const {
    return m_dynamic_other_connectivity.end_permutations(bucket_ordinal);
  }

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
  void grow_ngp_field_bucket_ids();

  Bucket();

  Bucket( const Bucket & );
  Bucket & operator = ( const Bucket & );

  Bucket(BulkData & mesh,
         EntityRank entityRank,
         const std::vector<unsigned> & key,
         unsigned initialCapacity,
         unsigned maximumCapacity,
         unsigned bucketId);

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

  unsigned get_others_begin_index(unsigned bucket_ordinal, EntityRank rank) const;
  unsigned get_others_end_index(unsigned bucket_ordinal, EntityRank rank) const;
  unsigned get_others_index_count(unsigned bucket_ordinal, EntityRank rank) const;

  template <typename T>
  void modify_connectivity(T& callable, EntityRank rank);

  template <typename T>
  void process_all_connectivity(T& callable, Bucket* other_bucket = nullptr);

  void check_for_invalid_connectivity_request(ConnectivityType const* type) const
  {
#ifndef NDEBUG
    debug_check_for_invalid_connectivity_request(type);
#endif
  }

  void debug_check_for_invalid_connectivity_request(ConnectivityType const* type) const;

  void reset_empty_space(const FieldVector & fields);

  friend class impl::BucketRepository;
  friend class impl::Partition;
  friend struct impl::OverwriteEntityFunctor;
  friend class BulkData;
  friend class DeviceMesh;

  BulkData             & m_mesh;
  const EntityRank       m_entity_rank;
  stk::topology          m_topology;
  std::vector<unsigned>  m_key;
  std::pair<const unsigned*,const unsigned*> m_partOrdsBeginEnd;
  unsigned               m_capacity;
  unsigned               m_maxCapacity;
  size_type              m_size;
  unsigned               m_bucket_id;
  unsigned               m_ngp_bucket_id;
  bool                   m_is_modified;
  PartVector             m_parts;

  std::vector<Entity> m_entities;

  impl::Partition *m_partition;

  ConnectivityType m_node_kind;
  ConnectivityType m_edge_kind;
  ConnectivityType m_face_kind;
  ConnectivityType m_element_kind;

  ConnectedEntities get_fixed_connectivity(unsigned bucket_ordinal, EntityRank rank) const
  {
    switch(rank) {
    case stk::topology::NODE_RANK: return m_fixed_node_connectivity.get_connected_entities(bucket_ordinal);
    case stk::topology::EDGE_RANK: return m_fixed_edge_connectivity.get_connected_entities(bucket_ordinal);
    case stk::topology::FACE_RANK: return m_fixed_face_connectivity.get_connected_entities(bucket_ordinal);
    case stk::topology::ELEM_RANK: return m_fixed_element_connectivity.get_connected_entities(bucket_ordinal);
    default: return m_dynamic_other_connectivity.get_connected_entities(bucket_ordinal);
    }
  }

  ConnectedEntities get_dynamic_connectivity(unsigned bucket_ordinal, EntityRank rank) const
  {
    switch(rank) {
    case stk::topology::NODE_RANK: return m_dynamic_node_connectivity.get_connected_entities(bucket_ordinal);
    case stk::topology::EDGE_RANK: return m_dynamic_edge_connectivity.get_connected_entities(bucket_ordinal);
    case stk::topology::FACE_RANK: return m_dynamic_face_connectivity.get_connected_entities(bucket_ordinal);
    case stk::topology::ELEM_RANK: return m_dynamic_element_connectivity.get_connected_entities(bucket_ordinal);
    default: return m_dynamic_other_connectivity.get_connected_entities(bucket_ordinal);
    }
  }

  impl::BucketConnectivity<stk::topology::NODE_RANK,    FIXED_CONNECTIVITY> m_fixed_node_connectivity;
  impl::BucketConnectivity<stk::topology::EDGE_RANK,    FIXED_CONNECTIVITY> m_fixed_edge_connectivity;
  impl::BucketConnectivity<stk::topology::FACE_RANK,    FIXED_CONNECTIVITY> m_fixed_face_connectivity;
  impl::BucketConnectivity<stk::topology::ELEMENT_RANK, FIXED_CONNECTIVITY> m_fixed_element_connectivity;

  impl::BucketConnDynamic m_dynamic_node_connectivity;
  impl::BucketConnDynamic m_dynamic_edge_connectivity;
  impl::BucketConnDynamic m_dynamic_face_connectivity;
  impl::BucketConnDynamic m_dynamic_element_connectivity;
  impl::BucketConnDynamic m_dynamic_other_connectivity;

  bool m_owned;
  bool m_shared;
  bool m_aura;

  std::vector<unsigned> m_ngp_field_bucket_id;
  std::vector<bool> m_ngp_field_is_modified;
};
#undef CONNECTIVITY_TYPE_SWITCH
#undef RANK_SWITCH

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
    STK_ThrowAssertMsg(lhs->entity_rank() == rhs->entity_rank(), "Cannot compare buckets of different rank");
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
  const unsigned* beg = m_partOrdsBeginEnd.first;
  const unsigned* end = m_partOrdsBeginEnd.second;
  for (unsigned ord : parts) {
    if (!contains_ordinal(beg, end, ord)) {
      return false;
    }
  }
  return true ;
}

inline
unsigned Bucket::num_connectivity(unsigned bucket_ordinal, EntityRank rank) const
{
  switch(rank) {
  case stk::topology::NODE_RANK:    return num_nodes(bucket_ordinal);
  case stk::topology::EDGE_RANK:    return num_edges(bucket_ordinal);
  case stk::topology::FACE_RANK:    return num_faces(bucket_ordinal);
  case stk::topology::ELEMENT_RANK: return num_elements(bucket_ordinal);
  default:
    return get_others_index_count(bucket_ordinal, rank);
  }
}

inline
bool Bucket::has_permutation(EntityRank rank) const
{
  switch(rank) {
  case stk::topology::NODE_RANK:
    return m_node_kind == FIXED_CONNECTIVITY ? m_fixed_node_connectivity.has_permutation() : m_dynamic_node_connectivity.has_permutation();
  case stk::topology::EDGE_RANK:
    return m_edge_kind == FIXED_CONNECTIVITY ? m_fixed_edge_connectivity.has_permutation() : m_dynamic_edge_connectivity.has_permutation();
  case stk::topology::FACE_RANK:
    return m_face_kind == FIXED_CONNECTIVITY ? m_fixed_face_connectivity.has_permutation() : m_dynamic_face_connectivity.has_permutation();
  case stk::topology::ELEMENT_RANK:
    return m_element_kind == FIXED_CONNECTIVITY ? m_fixed_element_connectivity.has_permutation() : m_dynamic_element_connectivity.has_permutation();
  case stk::topology::CONSTRAINT_RANK:
  default:
    return false;
  }
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
  case stk::topology::ELEMENT_RANK:
    return m_element_kind;
  default:
    return DYNAMIC_CONNECTIVITY;
  }
}

template <typename T>
inline
void Bucket::process_all_connectivity(T& callable, Bucket* other_bucket)
{
  if (callable.is_modifying()) {
    mark_for_modification();
  }

  switch(m_node_kind) {
  case FIXED_CONNECTIVITY:
    callable.template operator()<stk::topology::NODE_RANK, FIXED_CONNECTIVITY>(*this, m_fixed_node_connectivity, other_bucket); break;
  case DYNAMIC_CONNECTIVITY:
    callable.template operator()<stk::topology::NODE_RANK, DYNAMIC_CONNECTIVITY>(*this, m_dynamic_node_connectivity, other_bucket); break;
  default: break;
  }

  switch(m_edge_kind) {
  case FIXED_CONNECTIVITY:
    callable.template operator()<stk::topology::EDGE_RANK, FIXED_CONNECTIVITY>(*this, m_fixed_edge_connectivity, other_bucket); break;
  case DYNAMIC_CONNECTIVITY:
    callable.template operator()<stk::topology::EDGE_RANK, DYNAMIC_CONNECTIVITY>(*this, m_dynamic_edge_connectivity, other_bucket); break;
  default: break;
  }

  switch(m_face_kind) {
  case FIXED_CONNECTIVITY:
    callable.template operator()<stk::topology::FACE_RANK, FIXED_CONNECTIVITY>(*this, m_fixed_face_connectivity, other_bucket); break;
  case DYNAMIC_CONNECTIVITY:
    callable.template operator()<stk::topology::FACE_RANK, DYNAMIC_CONNECTIVITY>(*this, m_dynamic_face_connectivity, other_bucket); break;
  default: break;
  }

  switch(m_element_kind) {
  case FIXED_CONNECTIVITY:
    callable.template operator()<stk::topology::ELEM_RANK, FIXED_CONNECTIVITY>(*this, m_fixed_element_connectivity, other_bucket); break;
  case DYNAMIC_CONNECTIVITY:
    callable.template operator()<stk::topology::ELEM_RANK, DYNAMIC_CONNECTIVITY>(*this, m_dynamic_element_connectivity, other_bucket); break;
  default: break;
  }

  callable.template operator()<stk::topology::INVALID_RANK, DYNAMIC_CONNECTIVITY>(*this, m_dynamic_other_connectivity, other_bucket);
}

template <typename T>
inline
void Bucket::modify_connectivity(T& callable, EntityRank rank)
{
  switch(rank) {
  case stk::topology::NODE_RANK:
    mark_for_modification();

    switch(m_node_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_node_connectivity);   break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_node_connectivity); break;
    default: break;
    }
    break;
  case stk::topology::EDGE_RANK:
    switch(m_edge_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_edge_connectivity);   break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_edge_connectivity); break;
    default: break;
    }
    break;
  case stk::topology::FACE_RANK:
    switch(m_face_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_face_connectivity); break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_face_connectivity); break;
    default: break;
    }
    break;
  case stk::topology::ELEMENT_RANK:
    switch(m_element_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_element_connectivity);   break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_element_connectivity); break;
    default: break;
    }
    break;
  default:
    callable(*this, m_dynamic_other_connectivity);
    break;
  }
}

typedef Bucket::iterator BucketIterator;

} // namespace mesh
} // namespace stk



#endif 
