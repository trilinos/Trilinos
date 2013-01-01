#ifndef __IBMCPP__
#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_performance_test_includes/calculate_centroid.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

#include <boost/range.hpp>

using namespace stk::mesh;

namespace stk {
namespace performance_tests {

namespace {

void gather_hex( const Selector & element_select,
                 const fixtures::HexFixture::CoordFieldType & coordinates,
                 const BulkData & mesh,
                 std::vector<double> & avg_centroid )
{
  std::vector<double> elem_centroid(3,0);
  std::vector<double> elem_node_coords;

  const EntityRank node_rank(0);

  AllSelectedBucketsRange buckets = get_buckets( element_select, mesh );
  for (AllSelectedBucketsIterator b_iter = boost::begin(buckets), b_end = boost::end(buckets); b_iter != b_end; ++b_iter) {
    BucketIterator ent_begin = (*b_iter)->begin();
    BucketIterator ent_end   = (*b_iter)->end();

    for (BucketIterator ent_iter = ent_begin; ent_iter != ent_end; ++ent_iter) {
      Entity elem = *ent_iter;
      PairIterRelation node_relations = elem.relations(node_rank);
      const size_t num_nodes = node_relations.size();
      elem_node_coords.resize(num_nodes*3, 0.0);
      size_t offset = 0;

      for ( ; !node_relations.empty() ; ++node_relations ) {
        Entity node = node_relations->entity();
        double * node_coords = field_data(coordinates, node);
        elem_node_coords[offset++] = node_coords[0];
        elem_node_coords[offset++] = node_coords[1];
        elem_node_coords[offset++] = node_coords[2];
      }
      calculate_centroid_3d(num_nodes,&elem_node_coords[0],&elem_centroid[0]);

      //add this element-centroid to the avg_centroid vector, and
      //re-zero the element-centroid vector:
      avg_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
      avg_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
      avg_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
    }
  }
}

}

TEST( hex_gather, old_hex_gather )
{
  double start_time = stk::cpu_time();
  unsigned ex=100, ey=100, ez=100;
  unsigned num_elems = ex*ey*ez;
  fixtures::HexFixture fixture(MPI_COMM_SELF, ex, ey, ez);
  fixture.m_fem_meta.commit();
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3, 0.0);
  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  const int num_iters = 100;
  for (int t=0; t<num_iters; ++t) {
    Selector hex_elem_selector(fixture.m_hex_part & !fixture.m_node_part);
    gather_hex(hex_elem_selector, fixture.m_coord_field, fixture.m_bulk_data, avg_centroid);

    for (size_t d=0; d<3u; ++d) {
      const bool invert = d == 2;
      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - (invert ? -expected : expected)), tolerance);
    }

    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num centroid iterations: " << num_iters << std::endl;
  std::cout << "Time to compute centroids: " << test_time << std::endl;
}

} // namespace performance_tests
} // namespace stk

#endif
