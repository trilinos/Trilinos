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
#include "stk_util/stk_config.h"        // for STK_HAS_MPI
#if defined ( STK_HAS_MPI )
#  include <mpi.h>                        // for MPI_Comm
#endif

TEST(UnitTestTemplate, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 1;
  
#ifdef STK_HAS_MPI
  ASSERT_EQ(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  ASSERT_EQ(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
#endif
  
  ASSERT_TRUE(mpi_rank < mpi_size);
}
