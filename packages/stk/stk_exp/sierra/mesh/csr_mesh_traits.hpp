#ifndef SIERRA_MESH_CSR_MESH_TRAITS_HPP
#define SIERRA_MESH_CSR_MESH_TRAITS_HPP

#include <sierra/mesh/csr/csr_mesh.hpp>

namespace sierra {
namespace mesh {

template <>
struct mesh_traits<csr_mesh>
{
  typedef csr_mesh::part_identifier         part_identifier;

  typedef csr_mesh::entity_key              entity_key;
  typedef csr_mesh::entity_descriptor       entity_descriptor;
  typedef csr_mesh::entity_property         entity_property;
  typedef csr_mesh::entity_rank             entity_rank;
  typedef csr_mesh::entity_descriptor_range entity_descriptor_range;

  typedef csr_mesh::bucket_entity_iterator  bucket_entity_iterator;
  typedef csr_mesh::bucket_entity_range     bucket_entity_range;

  typedef csr_mesh::relation_descriptor  relation_descriptor;
  typedef csr_mesh::relation_range       relation_range;
  typedef csr_mesh::relation_position    relation_position;
  typedef csr_mesh::relation_iterator    relation_iterator;

  typedef csr_mesh::bucket_key           bucket_key;
  typedef csr_mesh::bucket_range         bucket_range;
  typedef csr_mesh::bucket_location      bucket_location;
  typedef csr_mesh::bucket_iterator bucket_iterator;

  typedef csr_mesh::selected_bucket_iterator selected_bucket_iterator;
  typedef csr_mesh::selected_bucket_range selected_bucket_range;

  typedef csr_mesh::selected_entity_iterator selected_entity_iterator;
  typedef csr_mesh::selected_entity_range selected_entity_range;

  typedef csr_mesh::part_bucket_range   part_bucket_range;

  typedef csr_mesh::part_key             part_key;
  typedef csr_mesh::part_range           part_range;
  typedef csr_mesh::part_property        part_property;

  typedef csr_mesh::bucket_part_iterator   bucket_part_iterator;
  typedef csr_mesh::bucket_part_range   bucket_part_range;

  typedef csr_mesh::sharing_proc_iterator sharing_proc_iterator;
  typedef csr_mesh::sharing_proc_range    sharing_proc_range;

  typedef csr_mesh::selector            selector;
};

//Provide the mesh_traits types outside the mesh_traits struct, to clean up
//syntax for non-templated client code:

  typedef mesh_traits< csr_mesh >::part_identifier         part_identifier;

  typedef mesh_traits< csr_mesh >::entity_key              entity_key;
  typedef mesh_traits< csr_mesh >::entity_descriptor       entity_descriptor;
  typedef mesh_traits< csr_mesh >::entity_property         entity_property;
  typedef mesh_traits< csr_mesh >::entity_rank             entity_rank;
  typedef mesh_traits< csr_mesh >::entity_descriptor_range entity_descriptor_range;

  typedef mesh_traits< csr_mesh >::bucket_entity_iterator  bucket_entity_iterator;
  typedef mesh_traits< csr_mesh >::bucket_entity_range     bucket_entity_range;

  typedef mesh_traits< csr_mesh >::relation_descriptor  relation_descriptor;
  typedef mesh_traits< csr_mesh >::relation_range       relation_range;
  typedef mesh_traits< csr_mesh >::relation_position    relation_position;
  typedef mesh_traits< csr_mesh >::relation_iterator    relation_iterator;

  typedef mesh_traits< csr_mesh >::bucket_key           bucket_key;
  typedef mesh_traits< csr_mesh >::bucket_range         bucket_range;
  typedef mesh_traits< csr_mesh >::bucket_location      bucket_location;
  typedef mesh_traits< csr_mesh >::bucket_iterator bucket_iterator;

  typedef mesh_traits< csr_mesh >::selected_bucket_iterator selected_bucket_iterator;
  typedef mesh_traits< csr_mesh >::selected_bucket_range selected_bucket_range;

  typedef mesh_traits< csr_mesh >::selected_entity_iterator selected_entity_iterator;
  typedef mesh_traits< csr_mesh >::selected_entity_range selected_entity_range;

  typedef mesh_traits< csr_mesh >::part_bucket_range   part_bucket_range;

  typedef mesh_traits< csr_mesh >::part_key             part_key;
  typedef mesh_traits< csr_mesh >::part_range           part_range;
  typedef mesh_traits< csr_mesh >::part_property        part_property;
  typedef mesh_traits< csr_mesh >::bucket_part_range   bucket_part_range;

  typedef mesh_traits< csr_mesh >::sharing_proc_iterator sharing_proc_iterator;
  typedef mesh_traits< csr_mesh >::sharing_proc_range    sharing_proc_range;

  typedef mesh_traits< csr_mesh >::selector            selector;

} // namespace mesh
} // namespace sierra

#endif // SIERRA_MESH_GENERIC_MESH_TRAITS_HPP
