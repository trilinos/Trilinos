#ifndef SIERRA_MESH_MODIFIABLE_MESH_API_HPP
#define SIERRA_MESH_MODIFIABLE_MESH_API_HPP

#include <sierra/mesh/modifiable_mesh_traits.hpp>
#include <sierra/mesh/modifiable/modifiable_mesh.hpp>

#include <stk_util/environment/ReportHandler.hpp>

namespace sierra {
namespace mesh {

// Add entity to mesh
inline
mesh_traits<modifiable_mesh>::entity_key
add_entity( modifiable_mesh & mesh )
{
  return mesh.add_entity();
}

// Add entity with property to mesh
inline
mesh_traits<modifiable_mesh>::entity_key
add_entity( const mesh_traits<modifiable_mesh>::entity_property & prop , modifiable_mesh & mesh )
{
  return mesh.add_entity(prop);
}

// Destroy_entity
inline
bool
remove_entity( const mesh_traits<modifiable_mesh>::entity_key & entity_key,
               modifiable_mesh & mesh )
{
  return mesh.remove_entity(entity_key);
}

//   // Add relation from entity_from to entity_to without position concept
//   inline
//   mesh_traits<modifiable_mesh>::relation_descriptor
//   add_relation(const mesh_traits<modifiable_mesh>::entity_key& entity_from,
//                const mesh_traits<modifiable_mesh>::entity_key& entity_to,
//                modifiable_mesh & mesh)
//   {
//     ThrowRequireMsg(false, "Error!  modifiable_mesh only supports add_relation with a relation position!");
//     return mesh_traits<modifiable_mesh>::relation_descriptor();
//   }

// Add relation from entity_from to entity_to with position
inline
mesh_traits<modifiable_mesh>::relation_descriptor
add_relation(const mesh_traits<modifiable_mesh>::entity_key& entity_from,
             const mesh_traits<modifiable_mesh>::entity_key& entity_to,
             const mesh_traits<modifiable_mesh>::relation_position & relation_position,
             modifiable_mesh & mesh)
{
  return mesh.add_relation(entity_from, entity_to, relation_position);
}

// Remove this relation from the mesh.
inline
bool
remove_relation( const mesh_traits<modifiable_mesh>::relation_descriptor& relation_d,
                 modifiable_mesh & mesh )
{
  return mesh.remove_relation(relation_d);
}

// Given a relation, get the targetted entity
inline
mesh_traits<modifiable_mesh>::entity_descriptor
target_entity( const mesh_traits<modifiable_mesh>::relation_descriptor& relation,
        const modifiable_mesh & mesh )
{
  return mesh.target_entity(relation);
}

// Generic API:  Get an entity_descriptor from an entity_key
inline
mesh_traits<modifiable_mesh>::entity_descriptor
entity_key_to_entity_descriptor(const mesh_traits<modifiable_mesh>::entity_key& entity_key,
                                const modifiable_mesh & mesh)
{
  return mesh.get_entity_descriptor(entity_key);
}

// Generic API:  Get an entity_key from an entity_descriptor
inline
mesh_traits<modifiable_mesh>::entity_key
entity_descriptor_to_entity_key(const mesh_traits<modifiable_mesh>::entity_descriptor& entity_d,
                                const modifiable_mesh & mesh)
{
  return mesh.get_entity_key(entity_d);
}

// Generic API:  Get a range to all entities in the mesh.
// Concept: EntityListmodifiable_mesh (similar to VertexListGraph concept for a graph).
inline
mesh_traits<modifiable_mesh>::entity_descriptor_range
get_entities(const modifiable_mesh & mesh)
{
  return mesh.get_entities();
}

// Generic API: Get all entities in a bucket
// Concept: HomogeneousSubsetmodifiable_mesh
inline
mesh_traits<modifiable_mesh>::bucket_entity_range
get_entities(const mesh_traits<modifiable_mesh>::bucket_key& bucket_key,
             const modifiable_mesh & mesh)
{
  return mesh.get_entities(bucket_key);
}

//   // Generic API:  Get a range to all the relations for this entity_local_id
//   // Concept: Connectivitymodifiable_mesh (similar to AdjacencyGraph)
//   inline
//   mesh_traits<modifiable_mesh>::relation_range
//   get_relations(const mesh_traits<modifiable_mesh>::entity_descriptor& entity,
//                 const modifiable_mesh & mesh )
//   {
//     return mesh.get_out_relations(entity);
//   }
//   
//   // Generic API:  Get a range to selected relations for this entity_local_id
//   // selector is a unary-predicate that takes a relation_descriptor and returns true/false
//   // Concept: ??? TODO
//   template <class RelationSelector>
//   inline
//   mesh_traits<modifiable_mesh>::relation_range
//   get_relations(const mesh_traits<modifiable_mesh>::entity_descriptor& entity,
//                 const RelationSelector & select,
//                 const modifiable_mesh & mesh)
//   {
//     // TODO
//     ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
//   }

// Generic API:  Get a bucket for an entity_descriptor
// Concept: HomogeneousSubsetmodifiable_mesh
inline
mesh_traits<modifiable_mesh>::bucket_key
get_bucket(const mesh_traits<modifiable_mesh>::entity_descriptor& entity,
           const modifiable_mesh & mesh)
{
  return mesh.get_bucket_location(entity).bucket();
}

// Generic API:  Get a bucket for an entity_descriptor
// Concept: HomogeneousSubsetmodifiable_mesh
inline
mesh_traits<modifiable_mesh>::bucket_location
get_bucket_location(const mesh_traits<modifiable_mesh>::entity_descriptor& entity,
                    const modifiable_mesh & mesh)
{
  return mesh.get_bucket_location(entity);
}

// Get all buckets in a mesh
// Concept: HomogeneousSubsetmodifiable_mesh
inline
mesh_traits<modifiable_mesh>::bucket_range
get_buckets( const modifiable_mesh & mesh )
{
  return mesh.get_buckets();
}

// Generic API:  Get buckets associated with a selector.
// selector is a unary-predicate that takes a bucket_descriptor and returns true/false
// Concept: HomogeneousSubsetmodifiable_mesh + selector
template <class selector>
inline
mesh_traits<modifiable_mesh>::selected_bucket_range
get_buckets( const selector & select,
             const modifiable_mesh & mesh )
{
  return get_selected_buckets(select, mesh);
}

// Generic API:  Get buckets for a particular part
// Concept: Bidirectional HomogeneousSubsetmodifiable_mesh (similar to bidirectional graph)
inline
mesh_traits<modifiable_mesh>::part_bucket_range
get_buckets( const mesh_traits<modifiable_mesh>::part_key& part,
             const modifiable_mesh & mesh )
{
  const modifiable_mesh::bucket_set& buckets = mesh.get_buckets(part);
  return std::make_pair(buckets.begin(), buckets.end());
}

// Generic API:  add part with property to the modifiable_mesh.
inline
mesh_traits<modifiable_mesh>::part_key
add_part( const mesh_traits<modifiable_mesh>::part_identifier& part_name,
          const mesh_traits<modifiable_mesh>::part_property& property,
          modifiable_mesh & mesh )
{
  return mesh.declare_part(part_name, property);
}

// Generic API:  remove this part from the modifiable_mesh.
inline
bool
remove_part( const mesh_traits<modifiable_mesh>::part_key& part,
             modifiable_mesh & mesh )
{
  // TODO - Dynamicmodifiable_mesh does not yet have a remove-part method
  return false;
}

// Move entity so it
// sits in Parts defined by AddPartInputIterator and
// so it does not sit in Parts defined by RemovePartInputIterator
template <class AddPartInputIterator, class RemovePartInputIterator>
inline
mesh_traits<modifiable_mesh>::bucket_key
move_entity( const mesh_traits<modifiable_mesh>::entity_key& entity_key,
             AddPartInputIterator add_first, AddPartInputIterator add_last,
             RemovePartInputIterator remove_first, RemovePartInputIterator remove_last,
             modifiable_mesh & mesh )
{
  mesh.change_entity_parts(entity_key, add_first, add_last, remove_first, remove_last);

  return mesh.get_bucket_location(entity_key).bucket();
}

// Generic API:  Get all parts on the mesh.
// Concept: Subsets
inline
mesh_traits<modifiable_mesh>::part_range
get_parts( const modifiable_mesh & mesh )
{
  return mesh.get_parts();
}

// Generic API:  Get all Parts associated with a bucket.
// Concept: HomogeneousSubsetmodifiable_mesh
inline
mesh_traits<modifiable_mesh>::bucket_part_range
get_parts(const mesh_traits<modifiable_mesh>::bucket_key& bucket,
          const modifiable_mesh & mesh)
{
  const modifiable_mesh::part_set& parts = mesh.get_parts(bucket);
  return std::make_pair(parts.begin(), parts.end());
}

// Generic API:  Begin modification cycle. Returns true if cycle was not already in progress
inline
bool
modification_begin( modifiable_mesh & mesh )
{
  // TODO - How to switch from CSR to modifiable with this interface? Use adaptor?
  return true;
}

// Generic API:  End modification cycle. Returns true if cycle was in progress.
inline
bool
modification_end( modifiable_mesh & mesh )
{
  // TODO - Same issue as with begin
  return true;
}

// Generic API:  Query if we are in a modification cycle
inline
bool
is_modifiable( const modifiable_mesh & mesh )
{
  return true;
}

inline
void
set_number_entity_ranks( int num_ranks, modifiable_mesh & mesh)
{
  mesh.set_num_entity_ranks(num_ranks);
}

} // namespace mesh
} // namespace sierra

#include <sierra/mesh/generic_mesh_utilities.hpp>

#endif // SIERRA_MESH_MODIFIABLE_MESH_API_HPP

