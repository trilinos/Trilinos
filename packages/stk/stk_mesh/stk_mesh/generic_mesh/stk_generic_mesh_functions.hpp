#ifndef STK_MESH_GENERIC_MESH_FUNCTIONS_HPP
#define STK_MESH_GENERIC_MESH_FUNCTIONS_HPP

#include "stk_mesh/generic_mesh/stk_generic_mesh_traits.hpp"

namespace stk {
namespace mesh {

// Not supported by stk_mesh.
inline STKGenericMesh::entity_local_id
add_entity( STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  STK_Mesh only supports add_entity with a stk::mesh::EntityKey!");
  return NULL;
}

// An API that is almost the same as the Generic Mesh API:
inline STKGenericMesh::entity_local_id
add_entity( const EntityKey & entity_key, STKGenericMesh & mesh )
{
  PartVector empty_part_vector;
  Entity & entity = mesh.m_bulk_data.declare_entity( entity_rank(entity_key), entity_id(entity_key), empty_part_vector );
  return & entity;
}


// Not supported by stk_mesh.
template <typename PartInputIterator>
inline STKGenericMesh::entity_local_id
add_entity(  PartInputIterator first, PartInputIterator last,
              STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  STK_Mesh only supports add_entity with a stk::mesh::EntityKey!");
  return NULL;
}

// An API that is almost the same as the Generic Mesh API:
template <typename PartInputIterator>
inline STKGenericMesh::entity_local_id
add_entity( const EntityKey & entity_key,
              PartInputIterator first, PartInputIterator last,
              STKGenericMesh & mesh )
{
  PartVector part_vector;
  std::copy(first,last,std::back_inserter(part_vector));
  Entity & entity = mesh.m_bulk_data.declare_entity( entity_rank(entity_key), entity_id(entity_key), part_vector );
  STKGenericMesh::entity_descriptor entity_lid = & entity;
  return entity_lid;
}

// destroy_entity:
inline void remove_entity( STKGenericMesh::entity_local_id entity_lid, STKGenericMesh & mesh )
{
  mesh.m_bulk_data.destroy_entity(entity_lid);
}


// Not supported by stk_mesh.
inline STKGenericMesh::relation_descriptor add_relation(
    STKGenericMesh::entity_local_id entity_from,
    STKGenericMesh::entity_local_id entity_to,
    STKGenericMesh & mesh
    )
{
  ThrowRequireMsg(false, "Error!  STK_Mesh only supports add_relation with a stk::mesh::RelationIdentifier!");
  return NULL;
}


// An API that is almost the same as the Generic Mesh API:
inline STKGenericMesh::relation_descriptor add_relation(
    STKGenericMesh::entity_local_id x_entity_from,
    STKGenericMesh::entity_local_id x_entity_to,
    const RelationIdentifier & relation_id,
    STKGenericMesh & mesh
    )
{
  const Entity & entity_from = get_entity(x_entity_from,mesh);
  const Entity & entity_to   = get_entity(x_entity_to,mesh);
  Entity & nonconst_entity_from = const_cast<Entity&>(entity_from);
  Entity & nonconst_entity_to = const_cast<Entity&>(entity_to);
  mesh.m_bulk_data.declare_relation(nonconst_entity_from,nonconst_entity_to,relation_id);
  Internal_STK_Relation_Descriptor_Adapter relation_d(x_entity_from, x_entity_to, relation_id);
  return relation_d;
}

// Remove this relation from the mesh.
inline void remove_relation( STKGenericMesh::relation_descriptor relation_d, STKGenericMesh & mesh )
{
  mesh.m_bulk_data.destroy_relation(*relation_d.m_entity_from, *relation_d.m_entity_to, relation_d.m_relation_id );
}


// Generic API:  Get an entity_descriptor from an entity_local_id
inline STKGenericMesh::entity_local_id entity_descriptor_to_local_id(
    STKGenericMesh::entity_descriptor entity_d,
    const STKGenericMesh & mesh
    )
{
  return entity_d;
}


// Generic API:  Get an entity_local_id from an entity_descriptor.
inline STKGenericMesh::entity_descriptor entity_local_id_to_descriptor(
    STKGenericMesh::entity_local_id entity_lid,
    const STKGenericMesh & mesh
    )
{
  return entity_lid;
}



// Generic API:  Get a range to all entities in the mesh.
inline std::pair<
                  STKGenericMesh::entity_descriptor_iterator,
                  STKGenericMesh::entity_descriptor_iterator
                >
  get_entities(STKGenericMesh & mesh)
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::entity_descriptor_iterator,STKGenericMesh::entity_descriptor_iterator>(NULL,NULL);
}


// Generic API:
inline
std::pair<
          STKGenericMesh::bucket_entity_descriptor_iterator,
          STKGenericMesh::bucket_entity_descriptor_iterator
         >
  get_entities( STKGenericMesh::bucket_descriptor bucket_descriptor,
                STKGenericMesh & mesh
              )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::bucket_entity_descriptor_iterator,STKGenericMesh::bucket_entity_descriptor_iterator>(NULL,NULL);
}



// Generic API:  Get a range to all the relations for this entity_local_id
inline std::pair<
                  STKGenericMesh::relation_descriptor_iterator,
                  STKGenericMesh::relation_descriptor_iterator
                >
  get_relations( STKGenericMesh::entity_local_id entity_lid, STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::relation_descriptor_iterator,STKGenericMesh::relation_descriptor_iterator>(NULL,NULL);
}


// Generic API:  Get a range to all the relations for this entity_descriptor
//inline std::pair<
//                  STKGenericMesh::relation_descriptor_iterator,
//                  STKGenericMesh::relation_descriptor_iterator
//                >
//  get_relations( STKGenericMesh::entity_descriptor entity_d, STKGenericMesh & mesh )
//{
//  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
//  return std::pair<STKGenericMesh::relation_descriptor_iterator,STKGenericMesh::relation_descriptor_iterator>(NULL,NULL);
//}



// Generic API:  Get a range to selected relations for this entity_local_id
// Selector is a unary-predicate that takes a relation_descriptor and returns true/false
template <typename Selector>
inline std::pair<
                  STKGenericMesh::selected_relation_descriptor_iterator,
                  STKGenericMesh::selected_relation_descriptor_iterator
                >
  get_relations(
      STKGenericMesh::entity_local_id entity_lid,
      Selector & selector,
      STKGenericMesh & mesh
      )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::selected_relation_descriptor_iterator,STKGenericMesh::selected_relation_descriptor_iterator>(NULL,NULL);
}


// Generic API:  Get a range to selected relations for this entity_descriptor
// Selector is a unary-predicate that takes a relation_descriptor and returns true/false
//template <typename Selector>
//inline std::pair<
//                  STKGenericMesh::selected_relation_descriptor_iterator,
//                  STKGenericMesh::selected_relation_descriptor_iterator
//                >
//  get_relations(
//      STKGenericMesh::entity_descriptor entity_d,
//      Selector & selector,
//      STKGenericMesh & mesh
//      )
//{
//  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
//  return std::pair<STKGenericMesh::selected_relation_descriptor_iterator,STKGenericMesh::selected_relation_descriptor_iterator>(NULL,NULL);
//}


// Generic API:  Get a bucket for an entity_descriptor
inline STKGenericMesh::bucket_descriptor
get_bucket( STKGenericMesh::entity_descriptor entity, STKGenericMesh & mesh )
{
  return &entity->bucket();
}


// Not supported by stk_mesh:
inline
std::pair<
          STKGenericMesh::bucket_descriptor_iterator,
          STKGenericMesh::bucket_descriptor_iterator
         >
get_buckets( const STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  STK_Mesh only supports get_buckets with an EntityRank specification!");
  return std::pair<STKGenericMesh::bucket_descriptor_iterator,STKGenericMesh::bucket_descriptor_iterator>(NULL,NULL);
}

// An API that is almost the same as the Generic Mesh API:
inline
std::pair<
          STKGenericMesh::bucket_descriptor_iterator,
          STKGenericMesh::bucket_descriptor_iterator
         >
get_buckets( EntityRank entity_rank, const STKGenericMesh & mesh )
{
  //const BucketVector & buckets = mesh.m_bulk_data.buckets(entity_rank);
  //return std::make_pair(buckets.begin(),buckets.end());
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::bucket_descriptor_iterator,STKGenericMesh::bucket_descriptor_iterator>(NULL,NULL);
}

// Generic API:  Get buckets associated with a Selector.
// Selector is a unary-predicate that takes a bucket_descriptor and returns true/false
template <typename Selector >
inline
std::pair<
          STKGenericMesh::selected_bucket_descriptor_iterator,
          STKGenericMesh::selected_bucket_descriptor_iterator
         >
get_buckets( const Selector & selector, STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::selected_relation_descriptor_iterator,STKGenericMesh::selected_relation_descriptor_iterator>(NULL,NULL);
}


// Generic API:  Get buckets for a particular part_descriptor.
inline
std::pair<
          STKGenericMesh::part_bucket_descriptor_iterator,
          STKGenericMesh::part_bucket_descriptor_iterator
         >
get_buckets( STKGenericMesh::part_descriptor part_descriptor, STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::part_bucket_descriptor_iterator,STKGenericMesh::part_bucket_descriptor_iterator>(NULL,NULL);
}


// Generic API:  add this part to the Mesh.
inline STKGenericMesh::part_descriptor
add_part( const STKGenericMesh::part_descriptor & part, STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return NULL;
}


// Generic API:  remove this part from the Mesh.
inline void
remove_part( STKGenericMesh::part_descriptor part_descriptor, STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
}


// Move entity so it
// sits in Parts defined by AddPartInputIterator and
// so it does not sit in Parts defined by RemovePartInputIterator
template <typename AddPartInputIterator, typename RemovePartInputIterator>
inline STKGenericMesh::bucket_descriptor
move_entity( STKGenericMesh::entity_descriptor entity_descriptor,
              AddPartInputIterator add_first, AddPartInputIterator add_last,
              RemovePartInputIterator remove_first, RemovePartInputIterator remove_last,
              STKGenericMesh & mesh )
{
  PartVector add_parts;
  PartVector remove_parts;
  std::copy(add_first,add_last,std::back_inserter(add_parts));
  std::copy(remove_first,remove_last,std::back_inserter(remove_parts));
  Entity & entity = get_entity(entity_descriptor, mesh);
  mesh.m_bulk_data.change_entity_parts( entity, add_parts, remove_parts );
  Bucket & bucket = entity.bucket();
  STKGenericMesh::bucket_descriptor bucket_d = &bucket;
  return bucket_d;
}


// Generic API:  Get all parts on the mesh.
inline
std::pair<
          STKGenericMesh::part_descriptor_iterator,
          STKGenericMesh::part_descriptor_iterator
         >
get_parts( const STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  This function is not implemented yet!");
  return std::pair<STKGenericMesh::part_descriptor_iterator,STKGenericMesh::part_descriptor_iterator>(NULL,NULL);
}


// Generic API:  Get all Parts associated with a bucket.
inline
std::pair<
          STKGenericMesh::bucket_part_descriptor_iterator,
          STKGenericMesh::bucket_part_descriptor_iterator
         >
get_parts(
    STKGenericMesh::bucket_descriptor bucket_descriptor, STKGenericMesh & mesh
    )
{
  return bucket_descriptor->superset_part_ordinals();
}

inline
STKGenericMesh::part_descriptor
get_part_descriptor(STKGenericMesh::bucket_part_descriptor descriptor, STKGenericMesh & mesh)
{
  return &(mesh.m_bulk_data.mesh_meta_data().get_part(descriptor));
}

// Generic API:  Begin modification cycle. Returns true if cycle was not already in progress
inline
bool
modification_begin( STKGenericMesh & mesh )
{
  return mesh.m_bulk_data.modification_begin();
}

// Generic API:  End modification cycle. Returns true if cycle was in progress.
inline
bool
modification_end( STKGenericMesh & mesh )
{
  return mesh.m_bulk_data.modification_end();
}

// Generic API:  Query if we are in a modification cycle
inline
bool
is_modifiable( const STKGenericMesh & mesh )
{
  return mesh.m_bulk_data.synchronized_state() == stk::mesh::BulkData::MODIFIABLE;
}

// Generic API:  Rotate the field data of multistate fields.
inline
void
rotate_multistate_fields( STKGenericMesh & mesh )
{
  mesh.m_bulk_data.update_field_data_states();
}

} // namespace mesh
} // namespace stk

#endif // STK_MESH_GENERIC_MESH_FUNCTIONS_HPP
