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

#include <stk_mesh/fixtures/Gear.hpp>
#include <stk_mesh/fixtures/GearsFixture.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_performance_test_includes/calculate_centroid.hpp>

#include <sstream>

namespace stk {
namespace performance_tests {

namespace {

//This very simple test will visit all local elements, gather coordinates,
//compute element-centroid (simple average of nodal coords) for each, and
//store the sum of the centroids in sum_centroid.
void do_stk_gather_gears_test(stk::mesh::BulkData& bulk, std::vector<double>& sum_centroid)
{
  using namespace stk::mesh;
  typedef Field<double,Cartesian> VectorField;

  MetaData& meta = MetaData::get(bulk);
  const unsigned spatial_dim = meta.spatial_dimension();
  for(unsigned d=0; d<spatial_dim; ++d) sum_centroid[d] = 0;

  std::vector<double> elem_centroid(spatial_dim, 0);

  const VectorField& coord_field = *meta.get_field<VectorField>("coordinates");
  ThrowAssert(&coord_field != NULL);

  Selector local = meta.locally_owned_part();

  BucketVector const& buckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, local);

  std::vector<double> elem_node_coords;

  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];

    for(size_t i=0; i<b.size(); ++i) {

      Entity const *node_rels = b.begin_nodes(i);
      size_t num_nodes = b.num_nodes(i);
      size_t len = num_nodes*spatial_dim;
      if (elem_node_coords.size() != len) elem_node_coords.resize(len);

      //here's the gather:

      unsigned offset = 0;
      for(size_t n=0; n<num_nodes; ++n) {
        Entity node = node_rels[n];
        double* node_coords = stk::mesh::field_data(coord_field, node);
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

STKUNIT_UNIT_TEST(gather_gears, gather_gears)
{
  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, 1,
      stk::mesh::fixtures::GearParams(0.01, 0.4, 1.5, -0.4, 0.4));
  fixture.meta_data.commit();

  double start_time = stk::cpu_time();

  fixture.generate_mesh();

  double mesh_create_time = stk::cpu_time() - start_time;

  const size_t spatial_dim = fixture.meta_data.spatial_dimension();

  std::vector<double> avg_centroid(spatial_dim, 0);

  start_time = stk::cpu_time();

  const int num_iters = 150;
  for(int t=0; t<num_iters; ++t) {

    do_stk_gather_gears_test(fixture.bulk_data, avg_centroid);

    double tolerance = 1.e-6;

    for(size_t d=0; d<spatial_dim; ++d) {
      double expected = 0;
      EXPECT_LT(std::abs(avg_centroid[d] - expected), tolerance);
    }
  }

  double gather_time = stk::cpu_time() - start_time;

  std::cout << "Gear: ";
  std::cout << "\tNum Nodes: " << fixture.get_gear(0).num_nodes;
  std::cout << "\tNum Elements: " << fixture.get_gear(0).num_elements << std::endl;

  double total_time = mesh_create_time + gather_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {mesh_create_time, gather_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Create mesh", "Gather", "Total time"};

  stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}

} // namespace performance_tests
} // namespace stk
