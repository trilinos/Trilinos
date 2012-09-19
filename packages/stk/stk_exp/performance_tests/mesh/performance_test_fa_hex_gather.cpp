#ifndef __IBMCPP__
#include <gtest/gtest.h>

#include <sierra/mesh/fixture/flat_array_hex_fixture.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <performance_tests/mesh/calculate_centroid.hpp>

#include <cmath>
#include <iostream>

using namespace sierra::mesh;

namespace {

void flat_array_gather( const fixture::flat_array_mesh & mesh,
                        const std::vector<double>& coord_field,
                        std::vector<double> & avg_centroid )
{
  std::vector<double> elem_centroid(3,0);
  std::vector<double> elem_node_coords(8*3,0);

  const std::vector<int>& conn = mesh.connectivity_table;

  int conn_offset = 0;
  for(size_t b=0; b<mesh.elem_blocks.size(); ++b) {
    int num_elems = mesh.elem_blocks[b].first;
    const int num_nodes = mesh.elem_blocks[b].second;

    elem_node_coords.resize(num_nodes*3);

    for( int i=0; i<num_elems; ++i) {
      int offset = 0;
      for( int n=0; n<num_nodes; ++n) {
        const double * node_coords = &coord_field[conn[conn_offset++]*3];
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }
      performance_tests::calculate_centroid_3d(num_nodes,&elem_node_coords[0], &elem_centroid[0]);

      //add this element-centroid to the avg_centroid vector, and
      //re-zero the element-centroid vector:
      avg_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      avg_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      avg_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
}

}

TEST( FlatArray, gather_centroid_hex_elem)
{
  double start_time = stk::cpu_time();
#ifndef NDEBUG
  unsigned ex=50,  ey=50,  ez=50; // make things smaller in debug so tests don't timeout
#else
  unsigned ex=100, ey=100, ez=100;
#endif
  unsigned num_elems = ex*ey*ez;
  fixture::flat_array_hex_fixture fixture(ex,ey,ez);
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "flat_array_hex_fixture: "<<std::endl;
  std::cout << "\tNum Nodes: " << fixture.m_num_nodes<<std::endl;
  std::cout << "\tNum Elements: " << fixture.m_num_elements << std::endl;
  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3,0);
  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    flat_array_gather(fixture.m_mesh, fixture.m_coord_field, avg_centroid);

    for(size_t d=0; d<3u; ++d) {
      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - expected), tolerance);
    }

    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Time to compute centroids ("<<num_iters<<" iters): " << test_time << std::endl;
}

// TEST( ArrayMesh, gather_centroid_hex_elem)
//{
//  double start_time = stk::cpu_time();
//#ifndef NDEBUG
//  unsigned ex=2,  ey=2,  ez=2; // make things smaller in debug so tests don't timeout
//#else
//  unsigned ex=100, ey=100, ez=100;
//#endif
//  unsigned num_elems = ex*ey*ez;
//  fixture::array_mesh_hex_fixture fixture(ex,ey,ez, true/*create upward node-to-elem co    nnectivity*/);
//  fixture.generate_mesh();
//  double end_time = stk::cpu_time() - start_time;
//
//  std::cout << "array_mesh_hex_fixture: "<<std::endl;
//  std::cout << "\tNum Nodes: " << fixture.m_num_nodes<<std::endl;
//  std::cout << "\tNum Elements: " << fixture.m_num_elements << std::endl;
//  std::cout << "Time to create hex mesh: " << end_time << std::endl;
//
//  std::vector<double> avg_centroid(3,0);
//  const double tolerance = 1.e-6;
//  const double expected = ((double)ex)/2;
//
//  start_time = stk::cpu_time();
//
//  const int num_iters = 100;
//  for(int t=0; t<num_iters; ++t) {
//    array_mesh_gather(fixture.m_mesh, fixture.m_coord_field, avg_centroid);
//
//    for(size_t d=0; d<3u; ++d) {
//      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - expected), tolerance);
//    }
//
//    avg_centroid[0] = 0;
//    avg_centroid[1] = 0;
//    avg_centroid[2] = 0;
//  }
//
//  end_time = stk::cpu_time();
//
//  double test_time = end_time - start_time;
//  std::cout << "Time to compute centroids ("<<num_iters<<" iters): " << test_time << st    d::endl;
//}

#endif
