/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <sstream>
#include <gtest/gtest.h>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/memory_util.hpp>
#include <stk_util/util/perf_util.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
#include <stk_performance_test_includes/calculate_centroid.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MemoryUsage.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/Comm.hpp>

namespace stk {
namespace performance_tests {

namespace {

//This very simple test will visit all local elements, gather coordinates,
//compute element-centroid (simple average of nodal coords) for each, and
//store the sum of the centroids in sum_centroid.
void do_stk_gather_test(stk::mesh::BulkData& bulk, std::vector<double>& sum_centroid)
{
  using namespace stk::mesh;
  typedef Field<double,Cartesian> VectorField;

  MetaData& meta = MetaData::get(bulk);
  const unsigned spatial_dim = meta.spatial_dimension();
  for(unsigned d=0; d<spatial_dim; ++d) sum_centroid[d] = 0;

  std::vector<double> elem_centroid(spatial_dim, 0);

  const VectorField& coord_field = *meta.get_field<VectorField>("coordinates");

  Selector local = meta.locally_owned_part();

  BucketVector buckets;
  get_buckets(local, bulk.buckets(stk::mesh::MetaData::ELEMENT_RANK), buckets);

  std::vector<double> elem_node_coords;

  size_t num_elems = 0;
  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];
    num_elems += b.size();
    const size_t num_nodes = b.topology().num_nodes();
    const size_t len = num_nodes*spatial_dim;
    if (elem_node_coords.size() != len) elem_node_coords.resize(len);

    for(size_t i=0; i<b.size(); ++i) {
      Entity elem = b[i];
      Entity const* node_rels = bulk.begin_nodes(elem);
      const int num_elem_nodes = bulk.num_nodes(elem);

      //here's the gather:

      unsigned offset = 0;
      for(int n = 0; n < num_elem_nodes; ++n) {
        Entity node = node_rels[n];
        double* node_coords = bulk.field_data(coord_field, node);
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }

      stk::performance_tests::calculate_centroid_3d(num_elem_nodes, &elem_node_coords[0], &elem_centroid[0]);

      //add this element-centroid to the sum_centroid vector, and
      //re-zero the element-centroid vector:
      sum_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      sum_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      sum_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
}

} // empty namespace

TEST(hex_edges, hex_edges)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
  int proc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  int numprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=100; //num_elems_x
  mesh_dims[1]=100; //num_elems_y
  mesh_dims[2]=100*numprocs; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2];

  double start_time = stk::cpu_time();

  size_t num_nodes = mesh_dims[0]+1;
  num_nodes *= mesh_dims[1]+1;
  num_nodes *= mesh_dims[2]+1;

  //if num_nodes is greater than the int limit (2.1B) then we need to
  //use the 64-bit api in the IO (which lies beneath the gmesh fixture).
  bool use_64bit_IO_api = num_nodes > 2150000000;

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str(), use_64bit_IO_api);
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  start_time = stk::cpu_time();

  stk::mesh::create_edges(fixture.getBulkData());

  double create_edges_time = stk::cpu_time() - start_time;
  double total_time = mesh_create_time + create_edges_time;

  std::vector<size_t> mesh_counts;
  stk::mesh::comm_mesh_counts(fixture.getBulkData(), mesh_counts);

  if (proc == 0) {
    std::cout<< "\nnum nodes: "<<mesh_counts[stk::topology::NODE_RANK]<<std::endl;
    std::cout<< "num edges: "<<mesh_counts[stk::topology::EDGE_RANK]<<std::endl;
    std::cout<< "num faces: "<<mesh_counts[stk::topology::FACE_RANK]<<std::endl;
    std::cout<< "num elems: "<<mesh_counts[stk::topology::ELEMENT_RANK]<<std::endl;
  }

  size_t now=0, hwm=0;
  stk::get_memory_usage(now,hwm);
  size_t global_hwm=0;
  MPI_Allreduce(&hwm, &global_hwm, 1, MPI_LONG_LONG, MPI_MAX, fixture.getBulkData().parallel());

  if (proc == 0) {
    static const int NUM_TIMERS = 3;
    const double timers[NUM_TIMERS] = {mesh_create_time, create_edges_time, total_time};
    const char* timer_names[NUM_TIMERS] = {"Create mesh", "Create edges", "Total time"};

    stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);
    std::cout<<"Global HWM: "<<stk::human_bytes(global_hwm)<<std::endl;
  }
}

TEST(hex_edges, minimal_hex_edges)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
  int proc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  int numprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=100; //num_elems_x
  mesh_dims[1]=100; //num_elems_y
  mesh_dims[2]=100*numprocs; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2];

  double start_time = stk::cpu_time();

  size_t num_nodes = mesh_dims[0]+1;
  num_nodes *= mesh_dims[1]+1;
  num_nodes *= mesh_dims[2]+1;

  //if num_nodes is greater than the int limit (2.1B) then we need to
  //use the 64-bit api in the IO (which lies beneath the gmesh fixture).
  bool use_64bit_IO_api = num_nodes > 2150000000;

  stk::mesh::ConnectivityMap connectivity_map = stk::mesh::ConnectivityMap::minimal_back_relations_map();
  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str(), use_64bit_IO_api, 1000, &connectivity_map );
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  start_time = stk::cpu_time();

  stk::mesh::create_edges(fixture.getBulkData());

  double create_edges_time = stk::cpu_time() - start_time;
  double total_time = mesh_create_time + create_edges_time;

  std::vector<size_t> mesh_counts;
  stk::mesh::comm_mesh_counts(fixture.getBulkData(), mesh_counts);

  if (proc == 0) {
    std::cout<< "\nnum nodes: "<<mesh_counts[stk::topology::NODE_RANK]<<std::endl;
    std::cout<< "num edges: "<<mesh_counts[stk::topology::EDGE_RANK]<<std::endl;
    std::cout<< "num faces: "<<mesh_counts[stk::topology::FACE_RANK]<<std::endl;
    std::cout<< "num elems: "<<mesh_counts[stk::topology::ELEMENT_RANK]<<std::endl;
  }

  size_t now=0, hwm=0;
  stk::get_memory_usage(now,hwm);
  size_t global_hwm=0;
  MPI_Allreduce(&hwm, &global_hwm, 1, MPI_LONG_LONG, MPI_MAX, fixture.getBulkData().parallel());

  if (proc == 0) {
    static const int NUM_TIMERS = 3;
    const double timers[NUM_TIMERS] = {mesh_create_time, create_edges_time, total_time};
    const char* timer_names[NUM_TIMERS] = {"Create mesh", "Create edges", "Total time"};

    stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);
    std::cout<<"Global HWM: "<<stk::human_bytes(global_hwm)<<std::endl;
  }
}

} //namespace performance_tests
} //namespace stk
