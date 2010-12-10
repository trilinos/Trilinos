/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <cstdlib>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(UnitTestParallel, testUnit)
{
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  std::string s;
  std::ostringstream strout;

  std::cout << "all_write_string " << std::flush;
  
//  for (size_t i = 0; i < 250000; ++i) {
  for (size_t i = 0; i < 100; ++i) {
    if (mpi_rank == 0 && i%1000 == 0)
      std::cout << "." << std::flush;
    
    stk::all_write_string(MPI_COMM_WORLD, strout, s);
  }
  
  STKUNIT_ASSERT_LT(mpi_rank, mpi_size);
}
