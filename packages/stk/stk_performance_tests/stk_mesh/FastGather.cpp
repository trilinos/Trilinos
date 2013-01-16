
/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/perf_util.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MemoryUsage.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_performance_test_includes/calculate_centroid.hpp>

#include <sstream>

namespace stk {
namespace performance_tests {

namespace {

//This very simple test will visit all local elements, gather coordinates,
//compute element-centroid (simple average of nodal coords) for each, and
//store the sum of the centroids in sum_centroid.
void do_stk_fast_gather_test(stk::mesh::BulkData& bulk, std::vector<double>& sum_centroid)
{
  using namespace stk::mesh;
  typedef Field<double,Cartesian> VectorField;

  MetaData& meta = MetaData::get(bulk);
  const unsigned spatial_dim = meta.spatial_dimension();
  for(unsigned d=0; d<spatial_dim; ++d) sum_centroid[d] = 0;

  std::vector<double> elem_centroid(spatial_dim, 0);

  ElementNodePointerField& coord_gather = *meta.get_field<ElementNodePointerField>("coord_gather");

  Selector local = meta.locally_owned_part();

  BucketVector buckets;
  get_buckets(local, bulk.buckets(stk::mesh::MetaData::ELEMENT_RANK), buckets);

  std::vector<double> elem_node_coords;

  size_t num_elems = 0;
  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];
    num_elems += b.size();

    for(size_t i=0; i<b.size(); ++i) {
      Entity elem = b[i];
      double** elem_node_coord_ptrs = field_data(coord_gather, elem);

      size_t num_nodes = 8;//!!!HARD-WIRED for HEX-8!!!!!
      size_t len = num_nodes*spatial_dim;
      if (elem_node_coords.size() != len) elem_node_coords.resize(len);

      //here's the gather:

      unsigned offset = 0;
      for(size_t j=0; j<num_nodes; ++j) {
        double* node_coords = *elem_node_coord_ptrs++;
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }

      stk::performance_tests::calculate_centroid_3d(num_nodes, &elem_node_coords[0], &elem_centroid[0]);

      //add this element-centroid to the sum_centroid vector, and
      //re-zero the element-centroid vector:
      sum_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      sum_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      sum_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
}

} // empty namespace

STKUNIT_UNIT_TEST(fast_gather, fast_gather)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=150; //num_elems_x
  mesh_dims[1]=150; //num_elems_y
  mesh_dims[2]=150; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2] ;

  double start_time = stk::cpu_time();

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str());

  stk::mesh::MetaData& meta = fixture.getMetaData();
  stk::mesh::Field<double,stk::mesh::Cartesian>* coord_field = meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");
  stk::mesh::ElementNodePointerField& coord_gather = stk::mesh::declare_element_node_pointer_field(meta, "coord_gather",*coord_field);
  stk::mesh::put_field(coord_gather, stk::mesh::MetaData::ELEMENT_RANK, meta.locally_owned_part(), 8/*hard-wired for hex-8*/);

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

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {

    do_stk_fast_gather_test(fixture.getBulkData(), avg_centroid);

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

} // namespace performance_tests
} // namespace stk
