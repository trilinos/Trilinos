#ifndef GenericGather_hpp
#define GenericGather_hpp

#include <sierra/mesh/mesh_traits.hpp>

#include <stk_performance_test_includes/calculate_centroid.hpp>

#include <boost/foreach.hpp>
#include <boost/range.hpp>

namespace stk {
namespace performance_tests {

// TODO: Should we further generalize this function so that it's not centroid specific?
template <class Mesh, class Field, class selector>
void gather_centroid_algorithm(const selector & element_select,
                               const Field & coordinates,
                               const Mesh & mesh,
                               std::vector<double> & avg_centroid)
{
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_key bucket_key;
  typedef typename sierra::mesh::mesh_traits<Mesh>::selected_bucket_range selected_bucket_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_entity_range bucket_entity_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_entity_iterator bucket_entity_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::selected_bucket_iterator selected_bucket_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_descriptor entity_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_range relation_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_descriptor relation_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_iterator relation_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_rank entity_rank;

  // hardcoded to 3d
  const unsigned spatial_dim = 3;

  const entity_rank node_rank(0);

  std::vector<double> elem_centroid(spatial_dim, 0.0);
  std::vector<double> elem_node_coords;

  selected_bucket_range buckets = sierra::mesh::get_buckets(element_select, mesh);

  for(selected_bucket_iterator b_iter = boost::begin(buckets), b_end = boost::end(buckets);  b_iter != b_end; ++b_iter) {
    bucket_entity_range entities = sierra::mesh::get_entities(*b_iter, mesh);

    for(bucket_entity_iterator ent_iter = boost::begin(entities), ent_end = boost::end(entities); ent_iter != ent_end; ++ent_iter) {
      entity_descriptor elem = *ent_iter;
      relation_range node_relations = sierra::mesh::get_relations(elem, node_rank, mesh);
      const size_t num_nodes = boost::distance(node_relations);
      elem_node_coords.resize(num_nodes * spatial_dim, 0.0);
      size_t offset = 0;

      for(relation_iterator rel_iter = boost::begin(node_relations), rel_end = boost::end(node_relations); rel_iter != rel_end; ++rel_iter) {
        entity_descriptor node = sierra::mesh::target_entity(*rel_iter, mesh);
        double * node_coords = sierra::mesh::get_field_data(coordinates, node, mesh);
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }
      stk::performance_tests::calculate_centroid_3d(num_nodes,&elem_node_coords[0],&elem_centroid[0]);

      //add this element-centroid to the avg_centroid vector, and
      //re-zero the element-centroid vector:
      avg_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      avg_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      avg_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
}

}
}

#endif
