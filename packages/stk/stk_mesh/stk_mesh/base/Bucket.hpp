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
#include <stk_mesh/base/BucketConnectivity.hpp>  // for BucketConnectivity
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Part.hpp>       // for contains_ordinal, Part
#include <stk_mesh/base/Types.hpp>
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
namespace stk { namespace mesh { namespace impl { struct OverwriteEntityFunctor; } } }
namespace stk { namespace mesh { namespace utest { struct ReversePartition; } } }
namespace stk { namespace mesh { namespace utest { struct SyncToPartitions; } } }
namespace stk { namespace mesh { struct ConnectivityMap; } }

namespace stk {
namespace mesh {

namespace impl {
class Partition;
class BucketRepository;
struct OverwriteEntityFunctor;
} // namespace impl

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
  bool member( const Part & ) const ;

  /** \brief  Bucket is a subset of all of the given parts */
  bool member_all( const PartVector & ) const ;
  bool member_all( const OrdinalVector & ) const ;

  /** \brief  Bucket is a subset of any of the given parts */
  bool member_any( const PartVector & ) const ;
  bool member_any( const OrdinalVector & ) const ;

  //--------------------------------
  /** Query bucket's supersets' ordinals. */

  std::pair<const unsigned *, const unsigned *>
  superset_part_ordinals() const
  {
    return std::pair<const unsigned *, const unsigned *>
      ( key() + 1 , key() + key()[0] );
  }

#ifndef DOXYGEN_COMPILE
  const unsigned * key() const { return m_key.data() ; }
#endif /* DOXYGEN_COMPILE */

  /** \brief  The allocation size, in bytes, of this bucket */
  unsigned allocation_size() const { return 0 ; }

  /** \brief  A method to assist in unit testing - accesses private data as necessary. */
  bool assert_correct() const;

  bool is_empty() const { return size() == 0; }

  impl::Partition *getPartition() const { return m_partition; }

  unsigned char* field_data_location(const FieldBase& field) const;

  bool field_data_is_allocated(const FieldBase& field) const;

  ///
  /// Entity member functions are moved here:
  ///

  int parallel_owner_rank(size_type ordinal) const;

  void check_size_invariant() const;

  //generic rank connectivity calls
  Entity const* begin(size_type bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, s, bucket_ordinal) }
  ConnectivityOrdinal const* begin_ordinals(size_type bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, _ordinals, bucket_ordinal) }
  Permutation const* begin_permutations(size_type bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, begin, _permutations, bucket_ordinal) }

  Entity const* end(size_type bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, s, bucket_ordinal) }
  ConnectivityOrdinal const* end_ordinals(size_type bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, _ordinals, bucket_ordinal) }
  Permutation const* end_permutations(size_type bucket_ordinal, EntityRank rank) const
  { RANK_SWITCH(rank, end, _permutations, bucket_ordinal) }

  unsigned num_connectivity(size_type bucket_ordinal, EntityRank rank) const;

  Entity const* begin_nodes(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.begin(bucket_ordinal), m_dynamic_node_connectivity.begin(bucket_ordinal), true) }
  Entity const* begin_edges(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.begin(bucket_ordinal), m_dynamic_edge_connectivity.begin(bucket_ordinal), true) }
  Entity const* begin_faces(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.begin(bucket_ordinal), m_dynamic_face_connectivity.begin(bucket_ordinal), true) }
  Entity const* begin_elements(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.begin(bucket_ordinal), m_dynamic_element_connectivity.begin(bucket_ordinal), true) }

  ConnectivityOrdinal const* begin_node_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_node_connectivity.begin_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* begin_edge_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_edge_connectivity.begin_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* begin_face_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_face_connectivity.begin_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* begin_element_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.begin_ordinals(bucket_ordinal), m_dynamic_element_connectivity.begin_ordinals(bucket_ordinal), true) }

  Permutation const* begin_node_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.begin_permutations(bucket_ordinal), m_dynamic_node_connectivity.begin_permutations(bucket_ordinal), true) }
  Permutation const* begin_edge_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.begin_permutations(bucket_ordinal), m_dynamic_edge_connectivity.begin_permutations(bucket_ordinal), true) }
  Permutation const* begin_face_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.begin_permutations(bucket_ordinal), m_dynamic_face_connectivity.begin_permutations(bucket_ordinal), true) }
  Permutation const* begin_element_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.begin_permutations(bucket_ordinal), m_dynamic_element_connectivity.begin_permutations(bucket_ordinal), true) }

  unsigned num_nodes(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.num_connectivity(bucket_ordinal), m_dynamic_node_connectivity.num_connectivity(bucket_ordinal), false) }
  unsigned num_edges(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.num_connectivity(bucket_ordinal), m_dynamic_edge_connectivity.num_connectivity(bucket_ordinal), false) }
  unsigned num_faces(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.num_connectivity(bucket_ordinal), m_dynamic_face_connectivity.num_connectivity(bucket_ordinal), false) }
  unsigned num_elements(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.num_connectivity(bucket_ordinal), m_dynamic_element_connectivity.num_connectivity(bucket_ordinal), false) }

  Entity const* end_nodes(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.end(bucket_ordinal), m_dynamic_node_connectivity.end(bucket_ordinal), true) }
  Entity const* end_edges(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.end(bucket_ordinal), m_dynamic_edge_connectivity.end(bucket_ordinal), true) }
  Entity const* end_faces(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.end(bucket_ordinal), m_dynamic_face_connectivity.end(bucket_ordinal), true) }
  Entity const* end_elements(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.end(bucket_ordinal), m_dynamic_element_connectivity.end(bucket_ordinal), true) }

  ConnectivityOrdinal const* end_node_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.end_ordinals(bucket_ordinal), m_dynamic_node_connectivity.end_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* end_edge_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.end_ordinals(bucket_ordinal), m_dynamic_edge_connectivity.end_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* end_face_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.end_ordinals(bucket_ordinal), m_dynamic_face_connectivity.end_ordinals(bucket_ordinal), true) }
  ConnectivityOrdinal const* end_element_ordinals(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.end_ordinals(bucket_ordinal), m_dynamic_element_connectivity.end_ordinals(bucket_ordinal), true) }

  Permutation const* end_node_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_node_kind, m_fixed_node_connectivity.end_permutations(bucket_ordinal), m_dynamic_node_connectivity.end_permutations(bucket_ordinal), true) }
  Permutation const* end_edge_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_edge_kind, m_fixed_edge_connectivity.end_permutations(bucket_ordinal), m_dynamic_edge_connectivity.end_permutations(bucket_ordinal), true) }
  Permutation const* end_face_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_face_kind, m_fixed_face_connectivity.end_permutations(bucket_ordinal), m_dynamic_face_connectivity.end_permutations(bucket_ordinal), true) }
  Permutation const* end_element_permutations(size_type bucket_ordinal) const
  { CONNECTIVITY_TYPE_SWITCH(m_element_kind, m_fixed_element_connectivity.end_permutations(bucket_ordinal), m_dynamic_element_connectivity.end_permutations(bucket_ordinal), true) }

  bool has_permutation(EntityRank rank) const;

  void debug_dump(std::ostream& out, unsigned ordinal = -1u) const;

  /* NGP Bucket methods */
  using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
  using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
  using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

  unsigned get_num_nodes_per_entity() const { return topology().num_nodes(); }

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

  bool member(stk::mesh::PartOrdinal partOrdinal) const;

  void set_ngp_field_bucket_id(unsigned fieldOrdinal, unsigned ngpFieldBucketId);
  unsigned get_ngp_field_bucket_id(unsigned fieldOrdinal) const;
  unsigned get_ngp_field_bucket_is_modified(unsigned fieldOrdinal) const;

protected:
  void change_existing_connectivity(unsigned bucket_ordinal, stk::mesh::Entity* new_nodes);
  void change_existing_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity, unsigned elem_connectivity_ordinal, stk::mesh::Permutation permut);
  void change_existing_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned edge_connectivity_ordinal, stk::mesh::Permutation permut);
  void change_existing_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned face_connectivity_ordinal, stk::mesh::Permutation permut);
  virtual ~Bucket();

private:

  bool destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id );

  bool declare_relation(size_type bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation);

  // The following *_other* functions should not be made available externally, in
  // order to avoid external confusion with "constraint" and "other" connectivities.
  // They are currently used within this class to provide connectivities
  // externally through another interface.
  Entity const* begin_others(size_type bucket_ordinal) const {
    return m_dynamic_other_connectivity.begin(bucket_ordinal);
  }
  ConnectivityOrdinal const* begin_other_ordinals(size_type bucket_ordinal) const {
    return m_dynamic_other_connectivity.begin_ordinals(bucket_ordinal);
  }
  Permutation const* begin_other_permutations(size_type bucket_ordinal) const {
    return m_dynamic_other_connectivity.begin_permutations(bucket_ordinal);
  }
  unsigned num_other(size_type bucket_ordinal) const {
    return m_dynamic_other_connectivity.num_connectivity(bucket_ordinal);
  }
  Entity const* end_others(size_type bucket_ordinal) const {
    return m_dynamic_other_connectivity.end(bucket_ordinal);
  }
  ConnectivityOrdinal const* end_other_ordinals(size_type bucket_ordinal) const {
    return m_dynamic_other_connectivity.end_ordinals(bucket_ordinal);
  }
  Permutation const* end_other_permutations(size_type bucket_ordinal) const {
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

  void mark_for_modification();

  void initialize_ngp_field_bucket_ids();

  Bucket();

  Bucket( const Bucket & );
  Bucket & operator = ( const Bucket & );

  Bucket( BulkData & arg_mesh ,
          EntityRank arg_entity_rank,
          const std::vector<unsigned> & arg_key,
          size_t arg_capacity,
          const ConnectivityMap& connectivity_map,
          unsigned bucket_id
        );

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
  void overwrite_entity(size_type to_ordinal, Entity entity, const std::vector<FieldBase*>* fields=NULL);

  void initialize_slot(size_type ordinal, Entity entity);
  //  Optional fields argument, only copy listed fields
  void reset_entity_location(Entity entity, size_type to_ordinal, const std::vector<FieldBase*>* fields = NULL);

  size_type get_others_begin_index(size_type bucket_ordinal, EntityRank rank) const;
  size_type get_others_end_index(size_type bucket_ordinal, EntityRank rank) const;
  size_type get_others_index_count(size_type bucket_ordinal, EntityRank rank) const;

  template <typename T>
  void modify_connectivity(T& callable, EntityRank rank);

  template <typename T>
  void modify_all_connectivity(T& callable, Bucket* other_bucket=NULL);

  void check_for_invalid_connectivity_request(ConnectivityType const* type) const
  {
#ifndef NDEBUG
    debug_check_for_invalid_connectivity_request(type);
#endif
  }

  void debug_check_for_invalid_connectivity_request(ConnectivityType const* type) const;

  friend class impl::BucketRepository;
  friend class impl::Partition;
  friend struct impl::OverwriteEntityFunctor;
  friend class BulkData;                // Replacement friend.
  friend struct Entity;
  friend struct utest::ReversePartition;
  friend struct utest::SyncToPartitions;
  friend class stk::mesh::DeviceMesh;

  BulkData             & m_mesh ;        // Where this bucket resides
  const EntityRank       m_entity_rank ; // Type of entities for this bucket
  stk::topology          m_topology ;    // The topology of this bucket
  std::vector<unsigned>  m_key ;         // REFACTOR
  const size_t           m_capacity ;    // Capacity for entities
  size_type              m_size ;        // Number of entities
  unsigned               m_bucket_id;    // Index into its BucketRepository's m_bucket[entity_rank()], these are NOT unique
  unsigned               m_ngp_bucket_id;
  bool                   m_is_modified;
  PartVector             m_parts;

  // Entity data
  std::vector<Entity>    m_entities;    // Array of entity handles; will be removed soon

  impl::Partition    *m_partition;

  ConnectivityType m_node_kind;
  ConnectivityType m_edge_kind;
  ConnectivityType m_face_kind;
  ConnectivityType m_element_kind;

  impl::BucketConnectivity<stk::topology::NODE_RANK,    FIXED_CONNECTIVITY> m_fixed_node_connectivity; // fixed connectivity to nodes
  impl::BucketConnectivity<stk::topology::EDGE_RANK,    FIXED_CONNECTIVITY> m_fixed_edge_connectivity; // fixed connectivity to edges
  impl::BucketConnectivity<stk::topology::FACE_RANK,    FIXED_CONNECTIVITY> m_fixed_face_connectivity; // fixed connectivity to faces
  impl::BucketConnectivity<stk::topology::ELEMENT_RANK, FIXED_CONNECTIVITY> m_fixed_element_connectivity; // fixed connectivity to elements

  impl::BucketConnectivity<stk::topology::NODE_RANK,    DYNAMIC_CONNECTIVITY> m_dynamic_node_connectivity; // dynamic connectivity to nodes
  impl::BucketConnectivity<stk::topology::EDGE_RANK,    DYNAMIC_CONNECTIVITY> m_dynamic_edge_connectivity; // dynamic connectivity to edges
  impl::BucketConnectivity<stk::topology::FACE_RANK,    DYNAMIC_CONNECTIVITY> m_dynamic_face_connectivity; // dynamic connectivity to faces
  impl::BucketConnectivity<stk::topology::ELEMENT_RANK, DYNAMIC_CONNECTIVITY> m_dynamic_element_connectivity; // dynamic connectivity to elements

  impl::BucketConnectivity<stk::topology::INVALID_RANK, DYNAMIC_CONNECTIVITY> m_dynamic_other_connectivity; // dynamic connectivity to everything else

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
  std::pair<const unsigned *, const unsigned *>
    part_ord = bucket.superset_part_ordinals();

  part_ord.first =
    std::lower_bound( part_ord.first , part_ord.second , ordinal );

  return part_ord.first < part_ord.second && ordinal == *part_ord.first ;
}

//----------------------------------------------------------------------
/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part"
 */

inline
bool has_superset( const Bucket & bucket ,  const Part & p )
{
  return has_superset(bucket,p.mesh_meta_data_ordinal());
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
unsigned Bucket::num_connectivity(size_type bucket_ordinal, EntityRank rank) const
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
void Bucket::modify_all_connectivity(T& callable, Bucket* other_bucket)
{
  mark_for_modification();

  switch(m_node_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_node_connectivity,   T::template generate_args<stk::topology::NODE_RANK, FIXED_CONNECTIVITY>(other_bucket)); break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_node_connectivity, T::template generate_args<stk::topology::NODE_RANK, DYNAMIC_CONNECTIVITY>(other_bucket)); break;
  default: break;
  }

  switch(m_edge_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_edge_connectivity,   T::template generate_args<stk::topology::EDGE_RANK, FIXED_CONNECTIVITY>(other_bucket)); break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_edge_connectivity, T::template generate_args<stk::topology::EDGE_RANK, DYNAMIC_CONNECTIVITY>(other_bucket)); break;
  default: break;
  }

  switch(m_face_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_face_connectivity,   T::template generate_args<stk::topology::FACE_RANK, FIXED_CONNECTIVITY>(other_bucket)); break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_face_connectivity, T::template generate_args<stk::topology::FACE_RANK, DYNAMIC_CONNECTIVITY>(other_bucket)); break;
  default: break;
  }

  switch(m_element_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_element_connectivity,   T::template generate_args<stk::topology::ELEMENT_RANK, FIXED_CONNECTIVITY>(other_bucket)); break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_element_connectivity, T::template generate_args<stk::topology::ELEMENT_RANK, DYNAMIC_CONNECTIVITY>(other_bucket)); break;
  default: break;
  }

  callable(*this, m_dynamic_other_connectivity, T::template generate_args<stk::topology::INVALID_RANK, DYNAMIC_CONNECTIVITY>(other_bucket));
}

template <typename T>
inline
void Bucket::modify_connectivity(T& callable, EntityRank rank)
{
  switch(rank) {
  case stk::topology::NODE_RANK:
    ThrowAssert(m_node_kind != INVALID_CONNECTIVITY_TYPE);
    mark_for_modification();

    switch(m_node_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_node_connectivity);   break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_node_connectivity); break;
    default: break;
    }
    break;
  case stk::topology::EDGE_RANK:
    ThrowAssert(m_edge_kind != INVALID_CONNECTIVITY_TYPE);
    switch(m_edge_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_edge_connectivity);   break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_edge_connectivity); break;
    default: break;
    }
    break;
  case stk::topology::FACE_RANK:
    ThrowAssert(m_face_kind != INVALID_CONNECTIVITY_TYPE);
    switch(m_face_kind) {
    case FIXED_CONNECTIVITY:   callable(*this, m_fixed_face_connectivity); break;
    case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_face_connectivity); break;
    default: break;
    }
    break;
  case stk::topology::ELEMENT_RANK:
    ThrowAssert(m_element_kind != INVALID_CONNECTIVITY_TYPE);
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
