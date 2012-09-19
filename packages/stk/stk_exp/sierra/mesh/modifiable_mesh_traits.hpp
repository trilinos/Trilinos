#ifndef SIERRA_MESH_MODIFIABLE_MESH_TRAITS_HPP
#define SIERRA_MESH_MODIFIABLE_MESH_TRAITS_HPP

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>
#include <sierra/mesh/mesh_traits.hpp>

namespace sierra {
namespace mesh {

// mesh_traits specialization for modifiable_mesh
template <>
struct mesh_traits<modifiable_mesh>
{
  typedef modifiable_mesh::part_identifier         part_identifier;

  typedef modifiable_mesh::entity_key              entity_key;
  typedef modifiable_mesh::entity_descriptor       entity_descriptor;
  typedef modifiable_mesh::entity_property         entity_property;
  typedef modifiable_mesh::entity_rank             entity_rank;
  typedef modifiable_mesh::entity_descriptor_range entity_descriptor_range;
  typedef modifiable_mesh::bucket_entity_range     bucket_entity_range;

  typedef modifiable_mesh::relation_descriptor  relation_descriptor;
  typedef modifiable_mesh::out_relation_range       out_relation_range;
  typedef modifiable_mesh::in_relation_range       in_relation_range;
  typedef modifiable_mesh::relation_position    relation_position;
  typedef modifiable_mesh::out_relation_iterator    out_relation_iterator;
  typedef modifiable_mesh::in_relation_iterator    in_relation_iterator;

  typedef modifiable_mesh::bucket_key           bucket_key;
  typedef modifiable_mesh::bucket_iterator      bucket_iterator;
  typedef modifiable_mesh::bucket_range         bucket_range;
  typedef modifiable_mesh::bucket_location      bucket_location;

  typedef modifiable_mesh::selected_bucket_range selected_bucket_range;
  typedef modifiable_mesh::selected_bucket_iterator selected_bucket_iterator;

  typedef modifiable_mesh::part_bucket_range   part_bucket_range;

  typedef modifiable_mesh::bucket_entity_iterator bucket_entity_iterator;

  typedef modifiable_mesh::part_key             part_key;
  typedef modifiable_mesh::part_range           part_range;
  typedef modifiable_mesh::part_property        part_property;

  typedef modifiable_mesh::bucket_part_iterator bucket_part_iterator;
  typedef modifiable_mesh::bucket_part_range   bucket_part_range;

  typedef modifiable_mesh::selector            selector;
};


} // namespace mesh
} // namespace sierra

#endif // SIERRA_MESH_GENERIC_MESH_TRAITS_HPP
