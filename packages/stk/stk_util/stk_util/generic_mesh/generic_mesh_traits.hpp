#ifndef STK_UTIL_GENERIC_MESH_TRAITS_HPP
#define STK_UTIL_GENERIC_MESH_TRAITS_HPP


namespace stk {

template <typename Mesh>
struct generic_mesh_traits {
  // persistant through modifications (up to change owner):
  //typedef typename Mesh::entity_global_id                entity_global_id;

  // persistant through modifications (up to change owner):
  typedef typename Mesh::entity_local_id                entity_local_id;
  // not-persistant through modifications (up to change owner):
  typedef typename Mesh::entity_descriptor              entity_descriptor;
  // entity_descriptor_iterator de-references to an entity_descriptor:
  typedef typename Mesh::entity_descriptor_iterator     entity_descriptor_iterator;
  typedef typename Mesh::entity_value                   entity_value;
  typedef typename Mesh::entity_size_type               entity_size_type; // not used yet


  typedef typename Mesh::part_descriptor                part_descriptor;
  // part_descriptor_iterator de-references to a part_descriptor:
  typedef typename Mesh::part_descriptor_iterator       part_descriptor_iterator;
  typedef typename Mesh::part_value                     part_value;
  typedef typename Mesh::part_size_type                 part_size_type; // not used yet

  typedef typename Mesh::bucket_descriptor              bucket_descriptor;
  // bucket_descriptor_iterator de-references to a bucket_descriptor:
  typedef typename Mesh::bucket_descriptor_iterator     bucket_descriptor_iterator;
  typedef typename Mesh::bucket_value                   bucket_value; // not used yet
  typedef typename Mesh::bucket_size_type               bucket_size_type; // not used yet
  // potentially heavy bucket_descriptor_iterator:
  typedef typename Mesh::selected_bucket_descriptor_iterator selected_bucket_descriptor_iterator;

  // potentially heavy bucket_descriptor_iterator:
  typedef typename Mesh::part_bucket_descriptor_iterator part_bucket_descriptor_iterator;
  typedef typename Mesh::part_bucket_size_type          part_bucket_size_type; // not used yet

  // potentially heavy part_descriptor_iterator:
  typedef typename Mesh::bucket_part_descriptor_iterator bucket_part_descriptor_iterator;
  typedef typename Mesh::bucket_part_size_type          bucket_part_size_type; // not used yet

  // potentially heavy entity_descriptor_iterator:
  typedef typename Mesh::bucket_entity_descriptor_iterator bucket_entity_descriptor_iterator;
  typedef typename Mesh::bucket_entity_size_type        bucket_entity_size_type; // not used yet

  inline part_descriptor                                universal_part();

  static inline part_descriptor                         null_part();
  static inline bucket_descriptor                       null_bucket();
  static inline entity_local_id                         null_entity();
  static inline relation_descriptor                     null_relation();

  typedef typename Mesh::relation_descriptor            relation_descriptor;
  typedef typename Mesh::relation_descriptor_iterator   relation_descriptor_iterator;
  typedef typename Mesh::selected_relation_descriptor_iterator selected_relation_descriptor_iterator
  typedef typename Mesh::relation_value                 relation_value; // not used yet
  typedef typename Mesh::relation_size_type             relation_size_type; // not used yet
};


// Note:  A Selector is a unary-predicate that accepts the appropriate type
// for the function being used and returns a boolean.

template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::part_descriptor
generic_mesh_traits<Mesh>::universal_part()
{ return Mesh::universal_part(); }

template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::part_descriptor
generic_mesh_traits<Mesh>::null_part()
{ return Mesh::null_part(); }

template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::bucket_descriptor
generic_mesh_traits<Mesh>::null_bucket()
{ return Mesh::null_bucket(); }

template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_local_id
generic_mesh_traits<Mesh>::null_entity()
{ return Mesh::null_entity(); }

template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::relation_descriptor
generic_mesh_traits<Mesh>::null_relation()
{ return Mesh::null_relation(); }


} // namespace stk


#endif // STK_UTIL_GENERIC_MESH_TRAITS_HPP
