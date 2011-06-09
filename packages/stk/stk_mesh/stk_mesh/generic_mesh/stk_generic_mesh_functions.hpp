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

// Not supported by stk_mesh.
inline STKGenericMesh::entity_local_id
add_entity( const STKGenericMesh::entity_value & entity_value, STKGenericMesh & mesh )
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


// Not supported by stk_mesh.
template <typename PartInputIterator>
inline STKGenericMesh::entity_local_id
add_entity( const STKGenericMesh::entity_value & entity_value,
              PartInputIterator first, PartInputIterator last,
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
  STKGenericMesh::entity_value & entity = mesh.m_bulk_data.declare_entity( entity_rank(entity_key), entity_id(entity_key), part_vector );
  STKGenericMesh::entity_descriptor entity_lid = & entity;
  return entity_lid;
}



// destroy_entity:
inline void remove_entity( typename STKGenericMesh::entity_local_id entity_lid, STKGenericMesh & mesh )
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


// Not supported by stk_mesh.
inline STKGenericMesh::relation_descriptor add_relation(
    STKGenericMesh::entity_local_id entity_from,
    STKGenericMesh::entity_local_id entity_to,
    const STKGenericMesh::relation_value & relation,
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
  STKGenericMesh::entity_value & entity_from = get_entity(x_entity_from,mesh);
  STKGenericMesh::entity_value & entity_to   = get_entity(x_entity_to,mesh);
  mesh.declare_relation(entity_from,entity_to,relation_id);
  Internal_STK_Relation_Descriptor_Adapter relation_d(x_entity_from, x_entity_to, relation_id);
  return relation_d;
}

// Remove this relation from the mesh.
inline void remove_relation( STKGenericMesh::relation_descriptor relation_d, STKGenericMesh & mesh )
{
  mesh.destroy_relation(*relation_d.m_entity_from, *relation_d.m_entity_to, relation_d.m_relation_id );
}


// Get a const reference to the Entity from the entity_local_id.
inline const STKGenericMesh::entity_value & get_entity(
    generic_mesh_traits<Mesh>::entity_local_id entity_lid,
    const Mesh & Mesh
    )
{
  return *entity_lid;
}


// Generic API:  Get a const reference to the Entity from the entity_descriptor.
template <typename Mesh>
inline const generic_mesh_traits<Mesh>::entity_value & get_entity(
    generic_mesh_traits<Mesh>::entity_descriptor entity_d,
    const Mesh & Mesh
    )
{
  return *entity_d;
}


// Generic API:  Get an entity_descriptor from an entity_local_id
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_local_id entity_descriptor_to_local_id(
    typename generic_mesh_traits<Mesh>::entity_descriptor entity_d,
    const Mesh & mesh
    )
{
  return entity_d;
}


// Generic API:  Get an entity_local_id from an entity_descriptor.
template <typename Mesh>
inline typename generic_mesh_traits<Mesh>::entity_descriptor entity_local_id_to_descriptor(
    typename generic_mesh_traits<Mesh>::entity_local_id entity_lid,
    const Mesh & mesh
    )
{
  return entity_lid;
}



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


// Not supported by stk_mesh:
inline
std::pair<
          typename STKGenericMesh::bucket_descriptor_iterator,
          typename STKGenericMesh::bucket_descriptor_iterator
         >
get_buckets( const STKGenericMesh & mesh )
{
  ThrowRequireMsg(false, "Error!  STK_Mesh only supports get_buckets with an EntityRank specification!");
  return std::pair<STKGenericMesh::bucket_descriptor_iterator,STKGenericMesh::bucket_descriptor_iterator>(NULL,NULL);
}

// An API that is almost the same as the Generic Mesh API:
inline
std::pair<
          typename STKGenericMesh::bucket_descriptor_iterator,
          typename STKGenericMesh::bucket_descriptor_iterator
         >
get_buckets( EntityRank entity_rank, const STKGenericMesh & mesh )
{
  const BucketVector & buckets = mesh.buckets(entity_rank);
  return std::make_pair(buckets.begin(),buckets.end());
}

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
  STKGenericMesh::entity_value & entity = get_entity(entity_descriptor, mesh);
  mesh.change_entity_parts( entity, add_parts, remove_parts );
  STKGenericMesh::bucket_value & bucket = entity.bucket();
  STKGenericMesh::bucket_descriptor bucket_d = &bucket;
  return bucket_d;
}



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


} // namespace mesh
} // namespace stk

#endif // STK_MESH_GENERIC_MESH_FUNCTIONS_HPP

