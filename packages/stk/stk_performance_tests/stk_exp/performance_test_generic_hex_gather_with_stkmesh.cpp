#ifndef __IBMCPP__
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <sierra/mesh/stkmesh_mesh_traits.hpp>
#include <sierra/mesh/stkmesh_mesh.hpp>

#include <stk_performance_test_includes/generic_gather.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

#include <boost/range.hpp>

namespace stk {
namespace performance_tests {

STKUNIT_UNIT_TEST( stk_mesh, generic_hex_gather )
{
  double start_time = stk::cpu_time();
  unsigned ex=100, ey=100, ez=100;
  unsigned num_elems = ex*ey*ez;
  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_SELF, ex, ey, ez);
  fixture.m_fem_meta.commit();
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "Time to create hex(" << ex << "," << ey << "," << ez << ") mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3, 0.0);
  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  const int num_iters = 100;
  for (int t=0; t<num_iters; ++t) {
    stk::mesh::Selector hex_elem_selector(fixture.m_hex_part & !fixture.m_node_part);
    gather_centroid_algorithm(hex_elem_selector, fixture.m_coord_field, fixture.m_bulk_data, avg_centroid);

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

}
}
#endif
