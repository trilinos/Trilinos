#ifndef STK_UTIL_GENERIC_MESH_FUNCTIONS_HPP
#define STK_UTIL_GENERIC_MESH_FUNCTIONS_HPP

namespace stk {

// Generic API:  add_entity that uses default constructable entity_value type.
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_local_id
add_entity( Mesh & mesh );

// Generic API:  Generic add_entity that takes a pre-constructed entity_value type.
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_local_id
add_entity( const typename generic_mesh_traits<Mesh>::entity_value & entity_value,
             Mesh & mesh );


// Convenience function that also adds the entity to the range of parts.
template <typename Mesh, typename PartInputIterator>
inline typename generic_mesh_traits<Mesh>::entity_local_id
add_entity(  PartInputIterator first, PartInputIterator last,
              Mesh & mesh );


// Convenience function that also adds the entity to the range of parts.
template <typename Mesh, typename PartInputIterator>
inline typename generic_mesh_traits<Mesh>::entity_local_id
add_entity( const typename generic_mesh_traits<Mesh>::entity_value & entity_value,
              PartInputIterator first, PartInputIterator last,
              Mesh & mesh );



// Generic API:  remove this entity from the Mesh.
template <typename Mesh>
inline void
remove_entity( typename generic_mesh_traits<Mesh>::entity_local_id entity_lid,
                Mesh & mesh );

// Generic API:  add this relation to the Mesh with default-constructed relation_value type.
template <typename Mesh>
inline generic_mesh_traits<Mesh>::relation_descriptor add_relation(
    generic_mesh_traits<Mesh>::entity_local_id entity_from,
    generic_mesh_traits<Mesh>::entity_local_id entity_to,
    Mesh & mesh
    );

// Generic API:  add this relation to the Mesh with pre-constructed relation_value type.
template <typename Mesh>
inline generic_mesh_traits<Mesh>::relation_descriptor add_relation(
    generic_mesh_traits<Mesh>::entity_local_id entity_from,
    generic_mesh_traits<Mesh>::entity_local_id entity_to,
    const generic_mesh_traits<Mesh>::relation_value & relation,
    Mesh & mesh
    );


// Generic API:  Remove this relation from the mesh.
template <typename Mesh>
inline void remove_relation( generic_mesh_traits<Mesh>::relation_descriptor relation_d, Mesh & mesh );


// Generic API:  Get a const reference to the Entity from the entity_local_id.
template <typename Mesh>
inline const generic_mesh_traits<Mesh>::entity_value & get_entity(
    generic_mesh_traits<Mesh>::entity_local_id entity_lid,
    const Mesh & Mesh
    );


// Generic API:  Get a const reference to the Entity from the entity_descriptor.
template <typename Mesh>
inline const generic_mesh_traits<Mesh>::entity_value & get_entity(
    generic_mesh_traits<Mesh>::entity_descriptor entity_d,
    const Mesh & Mesh
    );


// Generic API:  Get an entity_descriptor from an entity_local_id
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_local_id entity_descriptor_to_local_id(
    typename generic_mesh_traits<Mesh>::entity_descriptor entity_d,
    const Mesh & mesh
    );


// Generic API:  Get an entity_local_id from an entity_descriptor.
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_descriptor entity_local_id_to_descriptor(
    typename generic_mesh_traits<Mesh>::entity_local_id entity_lid,
    const Mesh & mesh
    );



// Generic API:  Get a range to all entities in the mesh.
template <typename Mesh>
inline std::pair<
                  typename generic_mesh_traits<Mesh>::entity_descriptor_iterator,
                  typename generic_mesh_traits<Mesh>::entity_descriptor_iterator
                >
  get_entities(Mesh & mesh);


// Generic API:
template <typename Mesh>
inline
std::pair<
          typename generic_mesh_traits<Mesh>::bucket_entity_descriptor_iterator,
          typename generic_mesh_traits<Mesh>::bucket_entity_descriptor_iterator
         >
  get_entities( typename generic_mesh_traits<Mesh>::bucket_descriptor bucket_descriptor,
                Mesh & mesh
              );



// Generic API:  Get a range to all the relations for this entity_local_id
template <typename Mesh>
inline std::pair<
                  typename generic_mesh_traits<Mesh>::relation_descriptor_iterator,
                  typename generic_mesh_traits<Mesh>::relation_descriptor_iterator
                >
  get_relations( generic_mesh_traits<Mesh>::entity_local_id entity_lid, Mesh & mesh );


// Generic API:  Get a range to all the relations for this entity_descriptor
template <typename Mesh>
inline std::pair<generic_mesh_traits<Mesh>::relation_descriptor_iterator, generic_mesh_traits<Mesh>::relation_descriptor_iterator>
  get_relations( generic_mesh_traits<Mesh>::entity_descriptor entity_d, Mesh & mesh );



// Generic API:  Get a range to selected relations for this entity_local_id
// Selector is a unary-predicate that takes a relation_descriptor and returns true/false
template <typename Mesh, typename Selector>
inline std::pair<generic_mesh_traits<Mesh>::selected_relation_descriptor_iterator,generic_mesh_traits<Mesh>::selected_relation_descriptor_iterator>
  get_relations(
      generic_mesh_traits<Mesh>::entity_local_id entity_lid,
      Selector & selector,
      Mesh & mesh
      );


// Generic API:  Get a range to selected relations for this entity_descriptor
// Selector is a unary-predicate that takes a relation_descriptor and returns true/false
template <typename Mesh, typename Selector>
inline std::pair<generic_mesh_traits<Mesh>::selected_relation_descriptor_iterator,generic_mesh_traits<Mesh>::selected_relation_descriptor_iterator>
  get_relations(
      generic_mesh_traits<Mesh>::entity_descriptor entity_d,
      Selector & selector,
      Mesh & mesh
      );


// Generic API:  Get a bucket for an entity_descriptor
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::bucket_descriptor
get_bucket( typename generic_mesh_traits<Mesh>::entity_descriptor entity,
               Mesh & mesh );


// Generic API:  Get all buckets for the mesh
template <typename Mesh>
inline
std::pair<
          typename generic_mesh_traits<Mesh>::bucket_descriptor_iterator,
          typename generic_mesh_traits<Mesh>::bucket_descriptor_iterator
         >
get_buckets( const Mesh & mesh );


// Generic API:  Get buckets associated with a Selector.
// Selector is a unary-predicate that takes a bucket_descriptor and returns true/false
template <typename Mesh, typename Selector >
inline
std::pair<
          typename generic_mesh_traits<Mesh>::selected_bucket_descriptor_iterator,
          typename generic_mesh_traits<Mesh>::selected_bucket_descriptor_iterator
         >
get_buckets( const Selector & selector, Mesh & mesh );


// Generic API:  Get buckets for a particular part_descriptor.
template <typename Mesh>
inline
std::pair<
          typename generic_mesh_traits<Mesh>::part_bucket_descriptor_iterator,
          typename generic_mesh_traits<Mesh>::part_bucket_descriptor_iterator
         >
get_buckets( typename generic_mesh_traits<Mesh>::part_descriptor part_descriptor,
               Mesh & mesh );



// Generic API:  add this part to the Mesh.
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::part_descriptor
add_part( const typename generic_mesh_traits<Mesh>::part_value & part_value,
         Mesh & mesh );


// Generic API:  remove this part from the Mesh.
template <typename Mesh>
inline void
remove_part( typename generic_mesh_traits<Mesh>::part_descriptor part_descriptor,
            Mesh & mesh );



// Generic API:  Move entity so it
// sits in Parts defined by AddPartInputIterator and
// so it does not sit in Parts defined by RemovePartInputIterator
template <typename Mesh, typename AddPartInputIterator, typename RemovePartInputIterator>
inline typename generic_mesh_traits<Mesh>::bucket_descriptor
move_entity( typename generic_mesh_traits<Mesh>::entity_descriptor entity_descriptor,
              AddPartInputIterator add_first, AddPartInputIterator add_last,
              RemovePartInputIterator remove_first, RemovePartInputIterator remove_last,
              Mesh & mesh );



// Generic API:  Get all parts on the mesh.
template <typename Mesh>
inline
std::pair<
          typename generic_mesh_traits<Mesh>::part_descriptor_iterator,
          typename generic_mesh_traits<Mesh>::part_descriptor_iterator
         >
get_parts( const Mesh & mesh );


// Generic API:  Get all Parts associated with a bucket.
template <typename Mesh>
inline
std::pair<
          typename generic_mesh_traits<Mesh>::bucket_part_descriptor_iterator,
          typename generic_mesh_traits<Mesh>::bucket_part_descriptor_iterator
         >
get_parts(
    typename generic_mesh_traits<Mesh>::bucket_descriptor bucket_descriptor,
    Mesh & mesh
    );


} // namespace stk

#endif // STK_UTIL_GENERIC_MESH_FUNCTIONS_HPP

