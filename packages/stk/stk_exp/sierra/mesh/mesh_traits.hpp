#ifndef SIERRA_MESH_GENERIC_MESH_TRAITS_HPP
#define SIERRA_MESH_GENERIC_MESH_TRAITS_HPP

namespace sierra {
namespace mesh {

// Completely general version, all associated types are expected to be defined
// in the mesh class.
template <class Mesh>
struct mesh_traits
{
  typedef typename Mesh::part_identifier         part_identifier;

  typedef typename Mesh::entity_key              entity_key;
  typedef typename Mesh::entity_descriptor       entity_descriptor;
  typedef typename Mesh::entity_property         entity_property;
  typedef typename Mesh::entity_rank             entity_rank;
  typedef typename Mesh::entity_descriptor_range entity_descriptor_range;

  typedef typename Mesh::bucket_entity_iterator  bucket_entity_iterator;
  typedef typename Mesh::bucket_entity_range     bucket_entity_range;

  typedef typename Mesh::relation_descriptor  relation_descriptor;
//  typedef typename Mesh::relation_range       relation_range;
  typedef typename Mesh::relation_position    relation_position;
//  typedef typename Mesh::relation_iterator    relation_iterator;

  typedef typename Mesh::bucket_key           bucket_key;
  typedef typename Mesh::bucket_range         bucket_range;
  typedef typename Mesh::bucket_location      bucket_location;
  typedef typename Mesh::bucket_iterator bucket_iterator;

  typedef typename Mesh::selected_bucket_iterator selected_bucket_iterator;
  typedef typename Mesh::selected_bucket_range selected_bucket_range;

//  typedef typename Mesh::selected_entity_iterator selected_entity_iterator;
//  typedef typename Mesh::selected_entity_range selected_entity_range;

  typedef typename Mesh::part_bucket_range   part_bucket_range;

  typedef typename Mesh::part_key             part_key;
  typedef typename Mesh::part_range           part_range;
  typedef typename Mesh::part_property        part_property;
  typedef typename Mesh::bucket_part_range   bucket_part_range;

  typedef typename Mesh::sharing_proc_iterator sharing_proc_iterator;
  typedef typename Mesh::sharing_proc_range    sharing_proc_range;

  typedef typename Mesh::selector            selector;
};


} // namespace mesh
} // namespace sierra

#endif // SIERRA_MESH_GENERIC_MESH_TRAITS_HPP
