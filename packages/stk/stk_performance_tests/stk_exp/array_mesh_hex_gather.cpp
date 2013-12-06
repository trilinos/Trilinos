#include <sierra/mesh/array_mesh/array_mesh.hpp>
#include <sierra/io/array_mesh_reader.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/memory_util.hpp>

#include <stk_performance_test_includes/calculate_centroid.hpp>

#include <cmath>
#include <iostream>

using namespace sierra::mesh;

namespace {

const unsigned spatial_dim = 3;


void block_gather(int num_elems,
                  const std::vector<int>& connectivity,
                  const std::vector<double>& coord_field,
                  std::vector<double>& avg_centroid, 
                  stk::topology topo)
{
  double elem_centroid[spatial_dim] = {0, 0, 0};
  const int num_nodes_per_elem = topo.num_nodes(); 
  double elem_node_coords[num_nodes_per_elem*spatial_dim];

  //gather nodal coordinates for each element, call calculate_centroid.
  int elem_offset = 0;
  for( int i=0; i<num_elems; ++i) {
    int offset = 0;
    for( int n=0; n<num_nodes_per_elem; ++n) {
      const double * node_coords = &coord_field[connectivity[elem_offset]*spatial_dim];
      elem_offset++;
      elem_node_coords[offset++] = node_coords[0];
      elem_node_coords[offset++] = node_coords[1];
      elem_node_coords[offset++] = node_coords[2];
    }

    stk::performance_tests::calculate_centroid_3d(num_nodes_per_elem,&elem_node_coords[0], elem_centroid);

    //add this element-centroid to the avg_centroid vector, and
    //re-zero the element-centroid vector:
    avg_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
    avg_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
    avg_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
  }
}

void array_mesh_gather( const array_mesh & mesh,
                  const std::vector<double>& coord_field,
                  std::vector<double> & avg_centroid )
{
  array_mesh::BlockRange blocks = mesh.get_blocks();

  for(array_mesh::BlockIterator b_it=blocks.first, b_end=blocks.second; b_it!=b_end; ++b_it) {
    int num_elems = mesh.get_num_elems(*b_it);
    stk::topology topo = mesh.get_topology(*b_it);

    switch(topo) {
      case stk::topology::HEX_8:
       block_gather(num_elems, mesh.get_block_connectivity(*b_it),
                            coord_field, avg_centroid, topo);
       break;
      case stk::topology::NODE: continue; break;
      case stk::topology::TET_4:
     default:
       std::cout<<"Unsupported topology!"<<std::endl;
    }
  }
}

}

STKUNIT_UNIT_TEST(array_mesh, gather_centroid_hex_elem_genmesh)
{
  double start_time = stk::cpu_time();
#ifndef NDEBUG
  unsigned ex=2,  ey=2,  ez=2; // make things smaller in debug
#else
  unsigned ex=100, ey=100, ez=100;
#endif
  unsigned num_elems = ex*ey*ez;

  std::ostringstream oss;
  oss << ex<<"x"<<ey<<"x"<<ez;
  std::string file_name = oss.str();
  std::string mesh_type("generated");

  bool create_upward_connectivity = true;
  sierra::mesh::array_mesh mesh(create_upward_connectivity);
  sierra::mesh::io::array_mesh_reader reader(MPI_COMM_WORLD, mesh_type, file_name, mesh);

  std::vector<double> mesh_coords;
  reader.read_nodal_field(mesh_coords, "mesh_model_coordinates");

  double end_time = stk::cpu_time() - start_time;

  std::cout << "\tNum Nodes: " << mesh.get_num_nodes()<<std::endl;
  std::cout << "\tNum Elements: " << mesh.get_num_elements() << std::endl;
  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3,0);
  const double tolerance = 1.e-6;
  const double expected = ex/2.0;

  start_time = stk::cpu_time();

  const int spatial_dim = 3;
  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    array_mesh_gather(mesh, mesh_coords, avg_centroid);

    for(int d=0; d<spatial_dim; ++d) {
      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - expected), tolerance);
    }

    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Time to compute centroids ("<<num_iters<<" iters): " << test_time << std::endl;

  size_t now = 0;
  size_t hwm = 0;
  stk::get_memory_usage(now, hwm);
  std::cout<<"Memory high-water-mark: "<<stk::human_bytes(hwm)<<std::endl;
}
