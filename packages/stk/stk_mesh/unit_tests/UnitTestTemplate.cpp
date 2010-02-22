
#include <stk_util/parallel/Parallel.hpp>
#include <unit_tests/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(UnitTestTemplate, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
#ifdef STK_HAS_MPI
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
#endif
  
  STKUNIT_ASSERT(mpi_rank < mpi_size);
}
