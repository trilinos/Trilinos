#ifndef CSRMESH_API_HPP
#define CSRMESH_API_HPP

#include <sierra/mesh/csr_mesh_traits.hpp>
#include <sierra/mesh/details/constant_size_field.hpp>
#include <sierra/mesh/csr/csr_mesh.hpp>

namespace sierra {
namespace mesh {

// Given a relation, get the targeted entity
inline
mesh_traits<csr_mesh>::entity_descriptor
target_entity( const mesh_traits<csr_mesh>::relation_descriptor& relation,
        const csr_mesh & mesh )
{
  return mesh.target_entity(relation);
}

// Generic API:  Get an entity_descriptor from an entity_key
inline
mesh_traits<csr_mesh>::entity_descriptor
entity_key_to_entity_descriptor(const mesh_traits<csr_mesh>::entity_key& entity_key,
                                const csr_mesh & mesh)
{
  return mesh.get_entity_descriptor(entity_key);
}

// Generic API:  Get an entity_key from an entity_descriptor
inline
mesh_traits<csr_mesh>::entity_key
entity_descriptor_to_entity_key(const mesh_traits<csr_mesh>::entity_descriptor& entity_d,
                                const csr_mesh & mesh)
{
  return mesh.get_entity_key(entity_d);
}

// Generic API: Get all entities in a bucket
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<csr_mesh>::bucket_entity_range
get_entities(const mesh_traits<csr_mesh>::bucket_key& bucket_key,
             const csr_mesh & mesh)
{
  return mesh.get_entities(bucket_key);
}

// Generic API:  Get a range to selected relations for this entity_local_id
// selector is a unary-predicate that takes a relation_descriptor and returns true/false
// Concept: ??? TODO
template <class RelationSelector>
inline
mesh_traits<csr_mesh>::relation_range
get_relations(const mesh_traits<csr_mesh>::entity_descriptor& entity,
              const RelationSelector & select,
              const csr_mesh & mesh)
{
  return mesh.get_relations(entity, select);
}

// Generic API:  Get a bucket for an entity_descriptor
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<csr_mesh>::bucket_key
get_bucket(const mesh_traits<csr_mesh>::entity_descriptor& entity,
           const csr_mesh & mesh)
{
  return mesh.get_bucket_location(entity).bucket();
}

// Generic API:  Get a bucket for an entity_descriptor
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<csr_mesh>::bucket_location
get_bucket_location(const mesh_traits<csr_mesh>::entity_descriptor& entity,
                    const csr_mesh & mesh)
{
  return mesh.get_bucket_location(entity);
}

// Generic API:  Get buckets associated with a selector.
// selector is a unary-predicate that takes a bucket_descriptor and returns true/false
// Concept: HomogeneousSubsetMesh + selector
inline
mesh_traits<csr_mesh>::selected_bucket_range
get_buckets( const mesh_traits<csr_mesh>::selector & select,
             const csr_mesh & mesh )
{
  return get_selected_buckets(select, mesh);
}

// Generic API:  Get all parts on the mesh.
// Concept: Subsets
inline
mesh_traits<csr_mesh>::part_range
get_parts( const csr_mesh & mesh )
{
  return mesh.get_parts();
}

// Generic API:  Get all Parts associated with a bucket.
// Concept: HomogeneousSubsetMesh
inline
mesh_traits<csr_mesh>::bucket_part_range
get_parts(const mesh_traits<csr_mesh>::bucket_key& bucket,
          const csr_mesh & mesh)
{
  const csr_mesh::part_set& parts = mesh.get_parts(bucket);
  return std::make_pair(parts.begin(), parts.end());
}

// Generic API:  Query if we are in a modification cycle
inline
bool
is_modifiable( const csr_mesh & mesh )
{
  return false;
}

// Parallel support

inline
int
get_owner_proc(const mesh_traits<csr_mesh>::entity_descriptor& entity,
               const csr_mesh & mesh)
{
  return mesh.get_owner_proc(entity);
}

inline
mesh_traits<csr_mesh>::sharing_proc_range
get_sharing_procs(const mesh_traits<csr_mesh>::entity_descriptor& entity,
               const csr_mesh & mesh)
{
  return mesh.get_sharing_procs(entity);
}

inline
int
get_owner_proc(const mesh_traits<csr_mesh>::entity_key& entity,
               const csr_mesh & mesh)
{
  return mesh.get_owner_proc(entity);
}

inline
mesh_traits<csr_mesh>::sharing_proc_range
get_sharing_procs(const mesh_traits<csr_mesh>::entity_key& entity,
               const csr_mesh & mesh)
{
  return mesh.get_sharing_procs(entity);
}

// Access field data
template<class Field>
inline
typename field_traits<Field>::data_type*
get_field_data( Field& field, const mesh_traits<csr_mesh>::entity_descriptor& entity, const csr_mesh& mesh )
{
  return field[mesh.get_bucket_location(entity)];
}

//inline
//bool
//is_selected(const mesh_traits<csr_mesh>::bucket_key& bucket,
//            const mesh_traits<csr_mesh>::selector& select,
//            const csr_mesh& mesh)
//{
//!!!!!!!!
//This function has been replaced by the is_selected function in
//generic_mesh_utilities.hpp
//!!!!!!!!
//  return select(bucket, mesh);
//}

// Access field data
template<class Field>
inline
typename field_traits<Field>::data_type*
get_field_data( Field& field, const mesh_traits<csr_mesh>::bucket_key& bucket, const csr_mesh& mesh )
{
  return field[bucket];
}

} // namespace mesh
} // namespace sierra

#include <sierra/mesh/generic_mesh_utilities.hpp>

#endif // CSRMESH_API_HPP
