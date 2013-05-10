/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/MetaData.hpp>

#ifdef SIERRA_MIGRATION

namespace sierra { namespace Fmwk {
class MeshBulkData;
}}

#endif // SIERRA_MIGRATION

#define STK_MESH_BACKWARDS_COMPATIBILITY

namespace stk {
namespace mesh {

namespace utest {
struct ReversePartition;
struct SyncToPartitions;
}

namespace impl {
class Partition;
class BucketRepository;
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

//----------------------------------------------------------------------
/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part"
 */
bool has_superset( const Bucket & ,  const Part & p );

/** \brief  Is this bucket a subset of the given
 *          \ref stk::mesh::Part "part" by partID
 */
bool has_superset( const Bucket & ,  const unsigned & ordinal );

/** \brief  Is this bucket a subset of all of the given
 *          \ref stk::mesh::Part "parts"
 */
bool has_superset( const Bucket & , const PartVector & );


//----------------------------------------------------------------------
/** \brief  A container for the connectivity for a homogeneous collection of
 *          \ref stk::mesh::Entity "entities".
 *
 *  The entities are homogeneous in that they are of the same entity type
 *  and are members of the same of parts.
 */
class Bucket
{
private:
  friend class impl::BucketRepository;
  friend class impl::Partition;
  friend class BulkData;                // Replacement friend.
  friend union Entity;
  friend class utest::ReversePartition;
  friend struct utest::SyncToPartitions;

  BulkData             & m_mesh ;        // Where this bucket resides
  const EntityRank       m_entity_rank ; // Type of entities for this bucket
  stk::topology          m_topology ;    // The topology of this bucket
  std::vector<unsigned>  m_key ;         // REFACTOR
  const size_t           m_capacity ;    // Capacity for entities
  size_t                 m_size ;        // Number of entities
  Bucket               * m_bucket ;      // Pointer to head of partition, but head points to tail
  unsigned               m_bucket_id;    // Index into its BucketRepository's m_bucket[entity_rank()], these are NOT unique

  //m_nodes_per_entity is how many connected-nodes (downward relations) exist
  //for each entity in this bucket.
  //If this is a node bucket or a bucket without a valid topology, then
  //m_nodes_per_entity will be 0.
  //TODO: This can go away once BucketConnectivity is fully live
  unsigned               m_nodes_per_entity;

  // Entity data
  std::vector<Entity>    m_entities;    // Array of entity handles; will be removed soon
  std::vector<int>       m_owner_ranks;

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

#ifdef SIERRA_MIGRATION
  const void*            m_fmwk_mesh_bulk_data;
#endif

public:

  stk::topology topology() const { return m_topology; }

  void parent_topology( EntityRank parent_rank, std::vector<stk::topology> & parent_topologies) const;

  bool owned() const { return m_owned; }
  bool shared() const { return m_shared; }

  //--------------------------------
  // Container-like types and methods:

  typedef const Entity * iterator;

  /** \brief Beginning of the bucket */
  inline iterator begin() const { return &m_entities[0]; }

  /** \brief End of the bucket */
  inline iterator end() const { return &m_entities[0] + m_size; }

  /** \brief  Number of entities associated with this bucket */
  size_t size() const { return m_size ; }

  /** \brief  Capacity of this bucket */
  size_t capacity() const { return m_capacity ; }

  /** \brief  Query the i^th entity */
  Entity operator[] ( size_t i ) const { return m_entities[i]; }

  ConnectivityType connectivity_type(EntityRank rank) const;

  //--------------------------------
  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & mesh() const { return m_mesh ; }

  /** \brief  Type of entities in this bucket */
  unsigned entity_rank() const { return m_entity_rank ; }

  unsigned bucket_id() const { return m_bucket_id; }

  /** \brief  This bucket is a subset of these \ref stk::mesh::Part "parts" */
  void supersets( PartVector & ) const ;
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

  /** \brief Equivalent buckets have the same parts
   */
  bool in_same_partition( const Bucket& b ) const {
    return first_bucket_in_partition() == b.first_bucket_in_partition();
  }


#ifndef DOXYGEN_COMPILE
  const unsigned * key() const { return &m_key[0] ; }
#endif /* DOXYGEN_COMPILE */

  /** \brief  The allocation size, in bytes, of this bucket */
  unsigned allocation_size() const { return 0 ; }

  /** \brief  A method to assist in unit testing - accesses private data as necessary. */
  bool assert_correct() const;

#ifdef SIERRA_MIGRATION
  typedef std::pair<iterator, iterator> EntityRange;

  bool is_empty() const { return size() == 0; }

  const sierra::Fmwk::MeshBulkData* get_bulk_data() const
  {
    return static_cast<const sierra::Fmwk::MeshBulkData*>(m_fmwk_mesh_bulk_data);
  }

  template <class T>
  void set_bulk_data(const T* bulk_ptr) { m_fmwk_mesh_bulk_data = bulk_ptr; }

#endif

  impl::Partition *getPartition() const { return m_partition; }

#ifdef STK_MESH_BACKWARDS_COMPATIBILITY
  unsigned char* field_data_location(const FieldBase& field) const;
#endif

  ///
  /// Entity member functions are moved here:
  ///

  int parallel_owner_rank(unsigned ordinal) const
  {
    return m_owner_ranks[ordinal];
  }

  void check_size_invariant() const;

  // Samba-like interface...
  // TODO: support beyond-element rank (e.g. constaint) connectivity generally

  Node const* begin_nodes(unsigned bucket_ordinal) const;
  Node const* end_nodes  (unsigned bucket_ordinal) const;

  Edge const* begin_edges(unsigned bucket_ordinal) const;
  Edge const* end_edges  (unsigned bucket_ordinal) const;

  Face const* begin_faces(unsigned bucket_ordinal) const;
  Face const* end_faces  (unsigned bucket_ordinal) const;

  Element const* begin_elements(unsigned bucket_ordinal) const;
  Element const* end_elements  (unsigned bucket_ordinal) const;

  Entity const* begin_entities(unsigned bucket_ordinal, EntityRank rank) const;
  Entity const* begin_node_entities(unsigned bucket_ordinal) const;
  Entity const* begin_edge_entities(unsigned bucket_ordinal) const;
  Entity const* begin_face_entities(unsigned bucket_ordinal) const;
  Entity const* begin_element_entities(unsigned bucket_ordinal) const;
  Entity const* begin_other_entities(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.begin_entities(bucket_ordinal); }

  bool other_entities_have_single_rank(unsigned bucket_ordinal, EntityRank rank) const;

  ConnectivityOrdinal const* begin_ordinals(unsigned bucket_ordinal, EntityRank rank) const;
  ConnectivityOrdinal const* begin_node_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* begin_edge_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* begin_face_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* begin_element_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* begin_other_ordinals(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.begin_ordinals(bucket_ordinal); }

  Permutation const* begin_permutations(unsigned bucket_ordinal, EntityRank rank) const;
  Permutation const* begin_node_permutations(unsigned bucket_ordinal) const;
  Permutation const* begin_edge_permutations(unsigned bucket_ordinal) const;
  Permutation const* begin_face_permutations(unsigned bucket_ordinal) const;
  Permutation const* begin_element_permutations(unsigned bucket_ordinal) const;
  Permutation const* begin_other_permutations(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.begin_permutations(bucket_ordinal); }

  unsigned num_connectivity(unsigned bucket_ordinal, EntityRank rank) const;
  unsigned num_nodes(unsigned bucket_ordinal) const;
  unsigned num_edges(unsigned bucket_ordinal) const;
  unsigned num_faces(unsigned bucket_ordinal) const;
  unsigned num_elements(unsigned bucket_ordinal) const;
  unsigned num_other(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.num_connectivity(bucket_ordinal); }

  Entity const* end_entities(unsigned bucket_ordinal, EntityRank rank) const;
  Entity const* end_node_entities(unsigned bucket_ordinal) const;
  Entity const* end_edge_entities(unsigned bucket_ordinal) const;
  Entity const* end_face_entities(unsigned bucket_ordinal) const;
  Entity const* end_element_entities(unsigned bucket_ordinal) const;
  Entity const* end_other_entities(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.end_entities(bucket_ordinal); }

  ConnectivityOrdinal const* end_ordinals(unsigned bucket_ordinal, EntityRank rank) const;
  ConnectivityOrdinal const* end_node_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* end_edge_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* end_face_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* end_element_ordinals(unsigned bucket_ordinal) const;
  ConnectivityOrdinal const* end_other_ordinals(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.end_ordinals(bucket_ordinal); }

  Permutation const* end_permutations(unsigned bucket_ordinal, EntityRank rank) const;
  Permutation const* end_node_permutations(unsigned bucket_ordinal) const;
  Permutation const* end_edge_permutations(unsigned bucket_ordinal) const;
  Permutation const* end_face_permutations(unsigned bucket_ordinal) const;
  Permutation const* end_element_permutations(unsigned bucket_ordinal) const;
  Permutation const* end_other_permutations(unsigned bucket_ordinal) const
  { return m_dynamic_other_connectivity.end_permutations(bucket_ordinal); }

  bool has_permutation(EntityRank rank) const;

  bool destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id );

  bool declare_relation(unsigned bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation);

private:
  /** \brief  The \ref stk::mesh::BulkData "bulk data manager"
   *          that owns this bucket.
   */
  BulkData & bulk_data() const { return mesh(); }

  ~Bucket();

  Bucket();
  Bucket( const Bucket & );
  Bucket & operator = ( const Bucket & );

  Bucket( BulkData & arg_mesh ,
          EntityRank arg_entity_rank,
          const std::vector<unsigned> & arg_key,
          size_t arg_capacity,
          const ConnectivityMap& connectivity_map
        );

  const std::vector<unsigned> & key_vector() const { return m_key; }

  void add_entity(Entity entity);
  void add_entity(); // do not use except when sorting
  void remove_entity(bool due_to_move=false);
  void move_entity(Entity entity);
  void replace_entity(unsigned ordinal, Entity entity); // overwrites existing entity

  void initialize_slot(unsigned ordinal, Entity entity);
  void internal_move_entity(Entity entity, unsigned to_ordinal);

  // BucketKey key = ( part-count , { part-ordinals } , counter )
  //  key[ key[0] ] == counter
  unsigned bucket_counter() const { return m_key[ m_key[0] ]; }

  size_t get_others_begin_index(unsigned bucket_ordinal, EntityRank rank) const;
  size_t get_others_end_index(unsigned bucket_ordinal, EntityRank rank) const;
  size_t get_others_index_count(unsigned bucket_ordinal, EntityRank rank) const;

  Bucket * last_bucket_in_partition() const;
  Bucket * first_bucket_in_partition() const;
  void set_last_bucket_in_partition( Bucket * last_bucket );
  void set_first_bucket_in_partition( Bucket * first_bucket );
  void set_partition_pointer( Bucket * bucket ) { m_bucket = bucket; }
  const Bucket * get_partition_pointer() const { return m_bucket; }

  Bucket * last_bucket_in_partition_impl() const;

  template <typename T>
  void modify_connectivity(T& callable, Bucket* other_bucket);

};

struct BucketLess {
  bool operator()( const Bucket * lhs_bucket , const unsigned * rhs ) const ;
  bool operator()( const unsigned * lhs , const Bucket * rhs_bucket ) const ;
};

inline
std::vector<Bucket*>::iterator
lower_bound( std::vector<Bucket*> & v , const unsigned * key )
{ return std::lower_bound( v.begin() , v.end() , key , BucketLess() ); }

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

//
// Define a begin/end pair
//
#define BEGIN_END_PAIR(rank_name, return_type, data_type)       \
                                                                \
inline                                                          \
return_type const* Bucket::begin_##rank_name##_##data_type(unsigned bucket_ordinal) const \
{                                                                       \
  switch(m_##rank_name##_kind) {                                        \
  case FIXED_CONNECTIVITY:                                              \
    return m_fixed_##rank_name##_connectivity.begin_##data_type(bucket_ordinal); \
  case DYNAMIC_CONNECTIVITY:                                            \
    return m_dynamic_##rank_name##_connectivity.begin_##data_type(bucket_ordinal); \
  default:                                                                \
    return NULL;                                                          \
  }                                                                     \
}                                                                       \
                                                                        \
inline                                                          \
return_type const* Bucket::end_##rank_name##_##data_type(unsigned bucket_ordinal) const \
{                                                                       \
  switch(m_##rank_name##_kind) {                                        \
  case FIXED_CONNECTIVITY:                                              \
    return m_fixed_##rank_name##_connectivity.end_##data_type(bucket_ordinal); \
  case DYNAMIC_CONNECTIVITY:                                            \
    return m_dynamic_##rank_name##_connectivity.end_##data_type(bucket_ordinal); \
  default:                                                                \
    return NULL;                                                          \
  }                                                                     \
}

#define FAST_BEGIN_END_PAIR(rank_name, type_name)                      \
inline                                                                 \
type_name const* Bucket::begin_##rank_name##s(unsigned bucket_ordinal) const    \
{                                                                      \
  switch(m_##rank_name##_kind) {                                        \
  case FIXED_CONNECTIVITY:                                              \
    return m_fixed_##rank_name##_connectivity.begin(bucket_ordinal); \
  case DYNAMIC_CONNECTIVITY:                                            \
    return m_dynamic_##rank_name##_connectivity.begin(bucket_ordinal); \
  default:                                                                \
    return NULL;                                                          \
  }                                                                     \
}                                                                      \
                                                                       \
inline                                                                 \
type_name const* Bucket::end_##rank_name##s(unsigned bucket_ordinal) const      \
{                                                                      \
  switch(m_##rank_name##_kind) {                                        \
  case FIXED_CONNECTIVITY:                                              \
    return m_fixed_##rank_name##_connectivity.end(bucket_ordinal); \
  case DYNAMIC_CONNECTIVITY:                                            \
    return m_dynamic_##rank_name##_connectivity.end(bucket_ordinal); \
  default:                                                                \
    return NULL;                                                          \
  }                                                                     \
}

//
// Define all methods for a rank
//
#define RANK_FUNCTION_DEFS(rank_name, type_name)                       \
                                                                       \
FAST_BEGIN_END_PAIR(rank_name, type_name)                              \
                                                                       \
BEGIN_END_PAIR(rank_name, Entity, entities)                            \
                                                                       \
BEGIN_END_PAIR(rank_name, ConnectivityOrdinal, ordinals)               \
                                                                       \
BEGIN_END_PAIR(rank_name, Permutation, permutations)                   \
                                                                       \
inline                                                                 \
unsigned Bucket::num_##rank_name##s(unsigned bucket_ordinal) const                \
{                                                                       \
  switch(m_##rank_name##_kind) {                                        \
  case FIXED_CONNECTIVITY:                                              \
    return m_fixed_##rank_name##_connectivity.num_connectivity(bucket_ordinal); \
  case DYNAMIC_CONNECTIVITY:                                            \
    return m_dynamic_##rank_name##_connectivity.num_connectivity(bucket_ordinal); \
  default:                                                                \
    return 0;                                                          \
  }                                                                    \
}

//
// Define method for runtime rank
//
#define FUNCTION_DEF(begin_str, end_str, return_type)   \
                                                        \
inline                                                  \
return_type const* Bucket::begin_str##_##end_str(unsigned bucket_ordinal, EntityRank rank) const \
{                                                                       \
  switch(rank) {                                                          \
  case stk::topology::NODE_RANK:    return begin_str##_node_##end_str(bucket_ordinal); \
  case stk::topology::EDGE_RANK:    return begin_str##_edge_##end_str(bucket_ordinal); \
  case stk::topology::FACE_RANK:    return begin_str##_face_##end_str(bucket_ordinal); \
  case stk::topology::ELEMENT_RANK: return begin_str##_element_##end_str(bucket_ordinal); \
  default:                                                              \
    return begin##_other_##end_str(bucket_ordinal) + get_others_##begin_str##_index(bucket_ordinal, rank);  \
  }                                                                     \
}

//
// Methods defined here!
//

RANK_FUNCTION_DEFS(node, Node);
RANK_FUNCTION_DEFS(edge, Edge);
RANK_FUNCTION_DEFS(face, Face);
RANK_FUNCTION_DEFS(element, Element);

FUNCTION_DEF(begin, entities, Entity);
FUNCTION_DEF(end,   entities, Entity);
FUNCTION_DEF(begin, ordinals, ConnectivityOrdinal);
FUNCTION_DEF(end,   ordinals, ConnectivityOrdinal);
FUNCTION_DEF(begin, permutations, Permutation);
FUNCTION_DEF(end,   permutations, Permutation);

#undef BEGIN_END_PAIR
#undef FAST_BEGIN_END_PAIR
#undef RANK_FUNCTION_DEFS
#undef FUNCTION_DEF

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
    return m_edge_kind == FIXED_CONNECTIVITY ? m_fixed_node_connectivity.has_permutation() : m_dynamic_node_connectivity.has_permutation();
  case stk::topology::FACE_RANK:
    return m_face_kind == FIXED_CONNECTIVITY ? m_fixed_node_connectivity.has_permutation() : m_dynamic_node_connectivity.has_permutation();
  case stk::topology::ELEMENT_RANK:
    return m_element_kind == FIXED_CONNECTIVITY ? m_fixed_node_connectivity.has_permutation() : m_dynamic_node_connectivity.has_permutation();
  default:
    return m_dynamic_other_connectivity.has_permutation();
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
void Bucket::modify_connectivity(T& callable, Bucket* other_bucket = NULL)
{
  switch(m_node_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_node_connectivity,   other_bucket == NULL ? NULL : &other_bucket->m_fixed_node_connectivity);   break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_node_connectivity, other_bucket == NULL ? NULL : &other_bucket->m_dynamic_node_connectivity); break;
  default: break;
  }

  switch(m_edge_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_edge_connectivity,   other_bucket == NULL ? NULL : &other_bucket->m_fixed_edge_connectivity);   break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_edge_connectivity, other_bucket == NULL ? NULL : &other_bucket->m_dynamic_edge_connectivity); break;
  default: break;
  }

  switch(m_face_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_face_connectivity,   other_bucket == NULL ? NULL : &other_bucket->m_fixed_face_connectivity);   break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_face_connectivity, other_bucket == NULL ? NULL : &other_bucket->m_dynamic_face_connectivity); break;
  default: break;
  }

  switch(m_element_kind) {
  case FIXED_CONNECTIVITY:   callable(*this, m_fixed_element_connectivity,   other_bucket == NULL ? NULL : &other_bucket->m_fixed_element_connectivity);   break;
  case DYNAMIC_CONNECTIVITY: callable(*this, m_dynamic_element_connectivity, other_bucket == NULL ? NULL : &other_bucket->m_dynamic_element_connectivity); break;
  default: break;
  }

  callable(*this, m_dynamic_other_connectivity, other_bucket == NULL ? NULL : &other_bucket->m_dynamic_other_connectivity);
}

typedef Bucket::iterator BucketIterator;

} // namespace mesh
} // namespace stk
