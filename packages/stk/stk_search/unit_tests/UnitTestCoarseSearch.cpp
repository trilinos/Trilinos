
#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

#include <stk_util/diag/Timer.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <unit_tests/UnitTestStkSearchFixture.hpp>

namespace stk_search_unit {

STKUNIT_UNIT_TEST(UnitTestingOfSearchCoarse, octtree)
{
  MPI_Barrier( MPI_COMM_WORLD );

  stk::search::FactoryOrder order;
  order.m_algorithm = stk::search::FactoryOrder::OCTREE;
  order.m_communicator = MPI_COMM_WORLD;

  OverlappingBoxes box_test;
  OverlappingBoxes::IdentProcRelation relation;

  stk::search::coarse_search(relation, box_test.range(), box_test.domain(), order);
  STKUNIT_ASSERT_TRUE(box_test.check_results(relation));
}

STKUNIT_UNIT_TEST(UnitTestingOfSearchCoarse, bihtree)
{
  MPI_Barrier( MPI_COMM_WORLD );

  stk::search::FactoryOrder order;
  order.m_algorithm = stk::search::FactoryOrder::BIHTREE;
  order.m_communicator = MPI_COMM_WORLD;

  OverlappingBoxes box_test;
  OverlappingBoxes::IdentProcRelation relation;

  stk::search::coarse_search(relation, box_test.range(), box_test.domain(), order);
  STKUNIT_ASSERT_TRUE(box_test.check_results(relation));
}

} // namespace stk_search_unit

