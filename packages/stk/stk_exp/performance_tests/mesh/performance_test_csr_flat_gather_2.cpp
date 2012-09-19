#ifndef __IBMCPP__
#include <gtest/gtest.h>

#include <sierra/mesh/fixture/hex_fixture.hpp>
#include <sierra/mesh/fixture/flat_array_hex_fixture.hpp>
#include <sierra/mesh/details/selected_buckets.hpp>
#include <sierra/mesh/csr/csr_mesh.hpp>
#include <sierra/mesh/fixture/csr_mesh_factory.hpp>

#include <performance_tests/mesh/calculate_centroid.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

using namespace sierra::mesh;
using namespace sierra::mesh::details;

namespace {

void setup_flat_array_mesh(fixture::flat_array_mesh& fmesh,
                           const selector& element_select,
                           const fixture::hex_fixture::CoordinateField& coordinates,
                           const csr_mesh& mesh)
{
  std::vector<std::pair<int,int> > elem_blocks(1, std::make_pair(1000000, 8));
  fmesh.allocate(elem_blocks);
  int counter = 0;
  BOOST_FOREACH(bucket_key bucket, get_selected_buckets(element_select, mesh))
  {
    BOOST_FOREACH(entity_descriptor elem, mesh.get_entities(bucket))
    {
      BOOST_FOREACH(csr_mesh::relation_descriptor node_rel, mesh.get_relations(elem,entity_rank(0)))
      {
        entity_descriptor node = mesh.target_entity(node_rel);
        bucket_location bucket_loc = mesh.get_bucket_location(node);
        fmesh.connectivity_table[counter++] = coordinates.offset(bucket_loc);
      }
    }
  }
}

void flat_gather(const fixture::flat_array_mesh& mesh,
                 const std::vector<double>& coord_field,
                 std::vector<double>& avg_centroid)
{
  std::vector<double> elem_centroid(3,0);
  std::vector<double> elem_node_coords;

  const std::vector<int>& conn = mesh.connectivity_table;

  int conn_offset = 0;
  for(size_t b=0; b<mesh.elem_blocks.size(); ++b) {
    int num_elems = mesh.elem_blocks[b].first;
    int num_nodes = mesh.elem_blocks[b].second;

    elem_node_coords.resize(num_elems*num_nodes*3);
    elem_centroid.resize(num_elems*3);
//    std::fill(elem_centroid.begin(),elem_centroid.end(), 0.0);

    int offset = 0;
    for( int i=0; i<num_elems; ++i) {
      for( int n=0; n<num_nodes; ++n) {
        const double * node_coords = &coord_field[conn[conn_offset++]];
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }
    }

    performance_tests::calculate_centroid_3d(num_elems, num_nodes, &elem_node_coords[0], &elem_centroid[0]);

    for( int i=0; i<num_elems; ++i) {
      //add this element-centroid to the avg_centroid vector, and
      //re-zero the element-centroid vector:
      avg_centroid[0] += elem_centroid[i*3 + 0]; elem_centroid[i*3 + 0] = 0.0;
      avg_centroid[1] += elem_centroid[i*3 + 1]; elem_centroid[i*3 + 1] = 0.0;
      avg_centroid[2] += elem_centroid[i*3 + 2]; elem_centroid[i*3 + 2] = 0.0;
    }
  }
}

}

TEST( Hex, CSR_flat_gather_2)
{
  double start_time = stk::cpu_time();
  unsigned ex=100, ey=100, ez=100;
  unsigned num_elems = ex*ey*ez;
  fixture::hex_fixture fixture(ex,ey,ez);
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "hex_fixture: "<<std::endl;
  std::cout << "\tNum Nodes: " << fixture.m_num_nodes<<std::endl;
  std::cout << "\tNum Elements: " << fixture.m_num_elements << std::endl;
  std::cout << "Time to create hex(" << ex << "," << ey << "," << ez << ") mesh: " << end_time << std::endl;

  start_time = stk::cpu_time();
  boost::shared_ptr<csr_mesh> csr = csr_mesh_factory::create_from_modifiable(fixture.m_mesh);
  end_time = stk::cpu_time() - start_time;

  std::cout << "Time to convert mesh: " << end_time << std::endl;
  std::cout << "sizeof(int): "<<sizeof(int)<<std::endl;
  std::cout << "sizeof(double*): "<<sizeof(double*)<<std::endl;

  std::vector<double> avg_centroid(3,0);
  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2.0;

  start_time = stk::cpu_time();

  fixture::flat_array_mesh fmesh;
  setup_flat_array_mesh(fmesh, fixture.m_hex_part&!fixture.m_node_part, fixture.m_coordinates, *csr);

  double end_setup_time = stk::cpu_time()-start_time;

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    flat_gather(fmesh, fixture.m_coordinates.field_data_flat_array(), avg_centroid);

    for(size_t d=0; d<3u; ++d) {
      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - expected), tolerance);
    }

    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num centroid iterations: " << num_iters << std::endl;
  std::cout << "Time to compute centroids: " << test_time << std::endl;
  std::cout << "Setup-gather time: " << end_setup_time << std::endl;
}
#endif
