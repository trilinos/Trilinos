#ifndef SIERRA_MESH_STKMESH_MESH_FUNCTIONS_HPP
#define SIERRA_MESH_STKMESH_MESH_FUNCTIONS_HPP

#include <sierra/mesh/stkmesh_mesh_traits.hpp>
#include <sierra/mesh/stkmesh_field_traits.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stdexcept>

namespace sierra {
namespace mesh {

// Need a way to ask if descriptors are valid
inline
bool
is_valid(mesh_traits<stk::mesh::BulkData>::relation_descriptor relation)
{
  return relation != mesh_traits<stk::mesh::BulkData>::relation_descriptor();
}

// Add entity with property to mesh
inline
mesh_traits<stk::mesh::BulkData>::entity_key
add_entity( mesh_traits<stk::mesh::BulkData>::entity_rank entity_rank,
            mesh_traits<stk::mesh::BulkData>::entity_property entity_id,
            stk::mesh::BulkData & mesh )
{
  stk::mesh::PartVector empty_parts;
  stk::mesh::Entity& new_entity = mesh.declare_entity(entity_rank, entity_id, empty_parts);
  return &new_entity;
}

// Add entity with property and parts to mesh. Optimization for STK_Mesh
inline
mesh_traits<stk::mesh::BulkData>::entity_descriptor
add_entity( mesh_traits<stk::mesh::BulkData>::entity_rank entity_rank,
            mesh_traits<stk::mesh::BulkData>::entity_property entity_id,
            stk::mesh::PartVector& parts,
            stk::mesh::BulkData & mesh )
{
  return &mesh.declare_entity(entity_rank, entity_id, parts);
}

// Destroy_entity
inline
bool
remove_entity( mesh_traits<stk::mesh::BulkData>::entity_key entity,
               stk::mesh::BulkData & mesh )
{
  return mesh.destroy_entity(*entity);
}

// Add relation from entity_from to entity_to with position
inline
mesh_traits<stk::mesh::BulkData>::relation_descriptor
add_relation(mesh_traits<stk::mesh::BulkData>::entity_key entity_from,
             mesh_traits<stk::mesh::BulkData>::entity_key entity_to,
             mesh_traits<stk::mesh::BulkData>::relation_position relation_position,
             stk::mesh::BulkData & mesh)
{
  try {
    mesh.declare_relation(*entity_from, *entity_to, relation_position);
  }
  catch (std::runtime_error&) {
    return stk::mesh::Relation(); // invalid relation descriptor
  }

  stk::mesh::PairIterRelation rels = entity_from->relations(entity_to->entity_rank());
  for ( ; !rels.empty(); ++rels) {
    if (rels->identifier() == relation_position && rels->entity() == entity_to) {
      return *rels;
    }
  }
  ThrowRequireMsg(false, "Failed to find added entity");
  return stk::mesh::Relation();
}

// Get relation position from relation-descriptor. Does this need to be in the generic API?
inline
mesh_traits<stk::mesh::BulkData>::relation_position
get_relation_position(const mesh_traits<stk::mesh::BulkData>::relation_descriptor& relation,
                      stk::mesh::BulkData& mesh)
{
  return relation.identifier();
}

// Remove this relation from the mesh.
inline
bool
remove_relation( const mesh_traits<stk::mesh::BulkData>::relation_descriptor& relation,
                 stk::mesh::BulkData & mesh )
{
  // TODO: How do we get the from ptr from Relation? Using an adaptor is an option, but will
  // cause performance problems since ranges of Relation will have to be converted to
  // ranges of the adaptor class, potentially requiring a "heavy" iterator.

  //return mesh.destroy_relation(*relation.m_entity_from,
  //                             *relation.m_entity_to,
  //                             relation.m_relation_id);
  return false;
}

// Given a relation, get the targetted entity
inline
mesh_traits<stk::mesh::BulkData>::entity_descriptor
target_entity( const mesh_traits<stk::mesh::BulkData>::relation_descriptor& relation,
               const stk::mesh::BulkData & mesh )
{
  return relation.entity();
}

// Generic API:  Get an entity_descriptor from an entity_key
inline
mesh_traits<stk::mesh::BulkData>::entity_descriptor
entity_key_to_entity_descriptor(mesh_traits<stk::mesh::BulkData>::entity_key entity,
                                const stk::mesh::BulkData & mesh)
{
  return entity;
}

// Generic API:  Get an entity_key from an entity_descriptor
inline
mesh_traits<stk::mesh::BulkData>::entity_key
entity_descriptor_to_entity_key(mesh_traits<stk::mesh::BulkData>::entity_descriptor entity,
                                const stk::mesh::BulkData & mesh)
{
  return entity;
}

// TEMPORARY workaround until we figure out how to handle relation descriptors. See TODO above.
inline
bool
remove_relation( mesh_traits<stk::mesh::BulkData>::entity_key entity_from,
                 mesh_traits<stk::mesh::BulkData>::entity_key entity_to,
                 mesh_traits<stk::mesh::BulkData>::relation_position position,
                 stk::mesh::BulkData & mesh)
{
  return mesh.destroy_relation(*entity_from, *entity_to, position);
}

// Generic API:  Get a range to all entities in the mesh.
// Concept: EntityListMesh (similar to VertexListGraph concept for a graph).
inline
mesh_traits<stk::mesh::BulkData>::entity_descriptor_range
get_entities(const mesh_traits<stk::mesh::BulkData>::entity_rank entity_rank, const stk::mesh::BulkData & mesh)
{
  return stk::mesh::get_entities(entity_rank, mesh);
}

// Generic API: Get all entities in a bucket
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<stk::mesh::BulkData>::bucket_entity_range
get_entities(mesh_traits<stk::mesh::BulkData>::bucket_key bucket,
             const stk::mesh::BulkData & mesh)
{
  return std::make_pair(bucket->begin(), bucket->end());
}

// Generic API:  Get a range to all the relations for this entity_local_id
// Concept: ConnectivityMesh (similar to AdjacencyGraph)
inline
mesh_traits<stk::mesh::BulkData>::relation_range
get_relations(mesh_traits<stk::mesh::BulkData>::entity_descriptor entity,
              const stk::mesh::BulkData & mesh )
{
  return entity->relations();
}

// Generic API:  Get a range to selected relations for this entity_local_id
// selector is a unary-predicate that takes a relation_descriptor and returns true/false
template <class RelationSelector>
inline
mesh_traits<stk::mesh::BulkData>::relation_range
get_relations(mesh_traits<stk::mesh::BulkData>::entity_descriptor entity,
              const RelationSelector & select,
              const stk::mesh::BulkData & mesh)
{
  return entity->relations(select);
}

// Generic API:  Get a bucket for an entity_descriptor
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<stk::mesh::BulkData>::bucket_key
get_bucket(mesh_traits<stk::mesh::BulkData>::entity_descriptor entity,
           const stk::mesh::BulkData & mesh)
{
  return & entity->bucket();
}

// Generic API:  Get a bucket for an entity_descriptor
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<stk::mesh::BulkData>::bucket_location
get_bucket_location(mesh_traits<stk::mesh::BulkData>::entity_descriptor entity,
                    const stk::mesh::BulkData & mesh)
{
  return entity->bucket_ordinal();
}

// Get all buckets in a mesh
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<stk::mesh::BulkData>::bucket_range
get_buckets( mesh_traits<stk::mesh::BulkData>::entity_rank ent_rank, const stk::mesh::BulkData & mesh )
{
  return std::make_pair(mesh.buckets(ent_rank).begin(), mesh.buckets(ent_rank).end());
}

// Get rank of an entity. Does this belong in the generic API?
inline
mesh_traits<stk::mesh::BulkData>::entity_rank
get_entity_rank( mesh_traits<stk::mesh::BulkData>::entity_key entity,
                 const stk::mesh::BulkData & mesh )
{
  return entity->entity_rank();
}

// Generic API:  Get buckets associated with a selector.
// selector is a unary-predicate that takes a bucket_descriptor and returns true/false
// Concept: HomogeneousSubsetMesh + selector
template <class Selector>
inline
mesh_traits<stk::mesh::BulkData>::selected_bucket_range
get_buckets( const Selector & selector,
             const stk::mesh::BulkData & mesh )
{
  return stk::mesh::get_buckets(selector, mesh);
}

// Generic API:  add part to the mesh.
inline
mesh_traits<stk::mesh::BulkData>::part_key
add_part( stk::mesh::BulkData & mesh )
{
  // TODO: STK_Mesh have both properties, names, and potentially EntityRanks; how are we going to support that?
  //mesh.meta_data().declare_part();
  return mesh_traits<stk::mesh::BulkData>::part_key();
}

// Generic API:  add part with property to the mesh.
inline
mesh_traits<stk::mesh::BulkData>::part_key
add_part( const mesh_traits<stk::mesh::BulkData>::part_property& property, stk::mesh::BulkData & mesh )
{
  // TODO: STK_Mesh have both properties, names, and potentially EntityRanks; how are we going to support that?
  //mesh.meta_data().declare_part();
  return mesh_traits<stk::mesh::BulkData>::part_key();
}

// Move entity so it
// sits in Parts defined by AddPartInputIterator and
// so it does not sit in Parts defined by RemovePartInputIterator
//
//Special sierra-framework-migration argument: the optional 'always_propagate_internal_changes'
//argument is true except when called from the sierra framework. The fmwk does its own part-change
//propagation redundantly, and it is a significant performance optimization to avoid having
//stk-mesh do it.
//We definitely don't want this flag in this function, but to remove it we need to remove more
//of the fmwk/stk duplication in sierra.
template <typename AddPartInputIterator, typename RemovePartInputIterator>
inline
mesh_traits<stk::mesh::BulkData>::bucket_key
move_entity( mesh_traits<stk::mesh::BulkData>::entity_key entity,
             AddPartInputIterator add_first, AddPartInputIterator add_last,
             RemovePartInputIterator remove_first, RemovePartInputIterator remove_last,
             stk::mesh::BulkData & mesh,
             bool always_propagate_internal_changes = true )
{
  mesh.change_entity_parts( *entity, add_first, add_last, remove_first, remove_last, always_propagate_internal_changes );
  stk::mesh::Bucket & bucket = entity->bucket();
  return & bucket;
}

// Generic API:  Get all parts on the mesh.
// Concept: Subsets
inline
mesh_traits<stk::mesh::BulkData>::part_range
get_parts( const stk::mesh::BulkData & mesh )
{
  const stk::mesh::PartVector& part_vector = mesh.mesh_meta_data().get_parts();
  return std::make_pair(boost::make_transform_iterator(part_vector.begin(), To_Ordinal()),
                        boost::make_transform_iterator(part_vector.end(), To_Ordinal()));
}

// Generic API:  Get all Parts associated with a bucket.
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<stk::mesh::BulkData>::bucket_part_range
get_parts(mesh_traits<stk::mesh::BulkData>::bucket_key bucket,
          const stk::mesh::BulkData & mesh)
{
  return bucket->superset_part_ordinals();
}

// Generic API:  Begin modification cycle. Returns true if cycle was not already in progress
inline
bool
modification_begin( stk::mesh::BulkData & mesh )
{
  return mesh.modification_begin();
}

// Generic API:  End modification cycle. Returns true if cycle was in progress.
inline
bool
modification_end( stk::mesh::BulkData & mesh )
{
  return mesh.modification_end();
}

// Generic API:  Query if we are in a modification cycle
inline
bool
is_modifiable( const stk::mesh::BulkData & mesh )
{
  return mesh.synchronized_state() == stk::mesh::BulkData::MODIFIABLE;
}

inline
void
set_number_entity_ranks( int num_ranks, stk::mesh::BulkData & mesh)
{
  //mesh.set_number_entity_ranks(num_ranks);
}

// Access field data
template<class Field>
inline
typename field_traits<Field>::data_type*
get_field_data( const Field& field,
                mesh_traits<stk::mesh::BulkData>::entity_descriptor entity,
                const stk::mesh::BulkData& mesh )
{
  return field_data(field, *entity);
}

} // namespace mesh
} // namespace sierra

#include <sierra/mesh/generic_mesh_utilities.hpp>

#endif // SIERRA_MESH_STKMESH_MESH_FUNCTIONS_HPP
