
#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

namespace stk_transfer_unit_tests {

//The following ::testUnit() function is where the actual unit-test is:
//(This one doesn't do much, it's mainly a template to demonstrate how
// to create a unit-test.)
//
//To create another unit-test, copy this file, change all occurrences of
//'UnitTest_1' to something else, and write your test code in its body.

void UnitTest_1( MPI_Comm comm )
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  STKUNIT_ASSERT(mpi_rank < mpi_size);
}

} // namespace stk_transfer_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfTransfer, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_transfer_unit_tests::UnitTest_1 ( MPI_COMM_WORLD );
}

