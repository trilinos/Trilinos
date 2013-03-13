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

namespace stk {
namespace performance_tests {

namespace {

//This very simple test will visit all local elements, gather coordinates,
//compute element-centroid (simple average of nodal coords) for each, and
//store the sum of the centroids in sum_centroid.
void copy_out_flat_connectivity_and_coords(stk::mesh::BulkData& bulk, std::vector<int>& connectivity, std::vector<double>& coords)
{
  using namespace stk::mesh;
  typedef Field<double,Cartesian> VectorField;

  MetaData& meta = MetaData::get(bulk);
  const unsigned spatial_dim = meta.spatial_dimension();

  const VectorField& coord_field = *meta.get_field<VectorField>("coordinates");

  Selector local = meta.locally_owned_part();

  BucketVector buckets;
  get_buckets(local, bulk.buckets(stk::mesh::MetaData::ELEMENT_RANK), buckets);

  size_t num_elems = 0;
  int nodes_per_elem = 0;
  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];
    if (nodes_per_elem == 0) nodes_per_elem = b.topology().num_nodes();
    num_elems += b.size();
  }

  connectivity.clear();
  connectivity.reserve(num_elems*nodes_per_elem);

  unsigned min_lid = 99;
  unsigned max_lid = 99;
  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];

    for(size_t i=0; i<b.size(); ++i) {
      Entity elem = b[i];
      PairIterRelation node_rels = elem.node_relations();

      //here's the gather:

      for(int n=0; n<8; ++n) {
        Entity node = node_rels[n].entity();
        if (min_lid > node.local_id()) min_lid = node.local_id();
        if (max_lid < node.local_id()) max_lid = node.local_id();
        connectivity.push_back(node.local_id());
      }
    }
  }

  get_buckets(local, bulk.buckets(stk::mesh::MetaData::NODE_RANK), buckets);

  size_t num_nodes = 0;
  for(size_t ib=0; ib<buckets.size(); ++ib) {
    num_nodes += buckets[ib]->size();
  }

  coords.resize(num_nodes*spatial_dim);

  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];
    const double* data = b.field_data(coord_field, b[0]);
    for(size_t i=0; i<b.size(); ++i) {
      int index = b[i].local_id()*spatial_dim;
      coords[index] = *data++;
      coords[index+1] = *data++;
      coords[index+2] = *data++;
    }
  }
}

void compute_centroid(int num_elems, int nodes_per_elem, int spatial_dim, const std::vector<int>& connectivity, const std::vector<double>& coords, std::vector<double>& sum_centroid)
{
  std::vector<double> elem_centroid(spatial_dim, 0);

  std::vector<double> elem_node_coords;

  for(int d=0; d<spatial_dim; ++d) sum_centroid[d] = 0;

  const size_t len = nodes_per_elem*spatial_dim;
  elem_node_coords.resize(len);

  const int* conn_ptr = &connectivity[0];
  for(int ie=0; ie<num_elems; ++ie) {

    unsigned offset = 0;
    for(int i=0; i<8; ++i) {
      int index = *conn_ptr++ * spatial_dim;
      elem_node_coords[offset++] = coords[index];
      elem_node_coords[offset++] = coords[index+1];
      elem_node_coords[offset++] = coords[index+2];
    }

    stk::performance_tests::calculate_centroid_3d(nodes_per_elem, &elem_node_coords[0], &elem_centroid[0]);

    //add this element-centroid to the sum_centroid vector, and
    //re-zero the element-centroid vector:
    sum_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
    sum_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
    sum_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
  }
}

} // empty namespace

TEST(hex_gather, hex_gather_copy)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
#ifndef NDEBUG
  mesh_dims[0]=30; //num_elems_x
  mesh_dims[1]=30; //num_elems_y
  mesh_dims[2]=30; //num_elems_z
#else
  mesh_dims[0]=100; //num_elems_x
  mesh_dims[1]=100; //num_elems_y
  mesh_dims[2]=100; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2];

  double start_time = stk::cpu_time();

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str());
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  const size_t spatial_dim = fixture.getMetaData().spatial_dimension();

  //compute total number of elements in the mesh:
  size_t num_elems = mesh_dims[0];
  for(size_t d=1; d<spatial_dim; ++d) {
    num_elems *= mesh_dims[d];
  }

  std::cout << "num_elems: " << num_elems << std::endl;

  std::vector<double> avg_centroid(spatial_dim, 0);

  start_time = stk::cpu_time();

  std::vector<int> connectivity;
  std::vector<double> coords;

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {

    copy_out_flat_connectivity_and_coords(fixture.getBulkData(), connectivity, coords);
    compute_centroid(num_elems, 8, spatial_dim, connectivity, coords, avg_centroid);

    //Next, divide by num_elems to turn the sum into an average,
    //then insist that the average matches our expected value.
    //NOTE: This error-check won't work in parallel. It is expecting
    //the local-average centroid to be a global-average.

    double tolerance = 1.e-6;

    for(size_t d=0; d<spatial_dim; ++d) {
      avg_centroid[d] /= num_elems;

      double expected = mesh_dims[d]/2.0;
      EXPECT_LT(std::abs(avg_centroid[d] - expected), tolerance);
    }
  }

  double gather_time = stk::cpu_time() - start_time;
  double total_time = mesh_create_time + gather_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {mesh_create_time, gather_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Create mesh", "Gather", "Total time"};

  stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);
}

} //namespace performance_tests
} //namespace stk
