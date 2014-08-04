/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <ostream>                      // for basic_ostream::operator<<
#include <gtest/gtest.h>
#include "gtest/gtest.h"                // for AssertHelper
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine

TEST(UnitTestTemplate, testUnit)
{
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  
  ASSERT_TRUE(mpi_rank < mpi_size);
}
