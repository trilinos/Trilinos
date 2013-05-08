/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_transfer/TransferP2P.hpp>

namespace stk_transfer_unit_tests {

void UnitTest_Linterp( MPI_Comm comm )
{
  int mpi_rank = 0;
  int mpi_size = 1;
  
#ifdef STK_HAS_MPI
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
#endif
  
  STKUNIT_ASSERT(mpi_rank < mpi_size);

  const double   EPSILON  = 0.000001;
  const unsigned    RUNS  =      100;
  const double   rand_max = RAND_MAX;

  STK_TransferP2P::MDArray A(3,3);
  std::vector<double> B(3);

  for (unsigned n=0; n<RUNS; ++n) {
    std::vector<double> Y(3,0);
    for (unsigned i=0; i<3; ++i) {
      B[i] = rand()/rand_max;
      for (unsigned j=0; j<3; ++j) A(i,j) = rand()/rand_max;
    }
    const std::vector<double> X = STK_TransferP2P::solve_3_by_3_with_LU(A, B);
    for (unsigned i=0; i<3; ++i) {
      for (unsigned j=0; j<3; ++j) Y[i] += A(i,j)*X[j] ;
    }
    for (unsigned i=0; i<3; ++i) STKUNIT_ASSERT(abs(Y[i]-B[i])<EPSILON);
  }
}

} // namespace stk_transfer_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfTransferLinterp, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_transfer_unit_tests::UnitTest_Linterp ( MPI_COMM_WORLD );
}

