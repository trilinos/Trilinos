/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

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

  stk_classic::search::FactoryOrder order;
  order.m_algorithm = stk_classic::search::FactoryOrder::OCTREE;
  order.m_communicator = MPI_COMM_WORLD;

  OverlappingBoxes box_test;
  OverlappingBoxes::IdentProcRelation relation;

  stk_classic::search::coarse_search(relation, box_test.range(), box_test.domain(), order);
  STKUNIT_ASSERT_TRUE(box_test.check_results(relation));
}

STKUNIT_UNIT_TEST(UnitTestingOfSearchCoarse, bihtree)
{
  MPI_Barrier( MPI_COMM_WORLD );

  stk_classic::search::FactoryOrder order;
  order.m_algorithm = stk_classic::search::FactoryOrder::BIHTREE;
  order.m_communicator = MPI_COMM_WORLD;

  OverlappingBoxes box_test;
  OverlappingBoxes::IdentProcRelation relation;

  stk_classic::search::coarse_search(relation, box_test.range(), box_test.domain(), order);
  STKUNIT_ASSERT_TRUE(box_test.check_results(relation));
}

} // namespace stk_search_unit

