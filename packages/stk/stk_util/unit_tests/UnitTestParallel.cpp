/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll, CommBuffer
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_write_string
#include <gtest/gtest.h>
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD


TEST(UnitTestParallel, testUnit)
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
  
  ASSERT_LT(mpi_rank, mpi_size);
}

TEST(UnitTestParallel, testCommAll)
{
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  ASSERT_LT(mpi_rank, mpi_size);

  stk::CommAll comm_all(MPI_COMM_WORLD);

  ASSERT_EQ(comm_all.parallel_size(), mpi_size);
  ASSERT_EQ(comm_all.parallel_rank(), mpi_rank);

  for (int p = 0; p<mpi_size; ++p) {
    if (p - mpi_rank <= 10 || mpi_rank - p <= 10) {
      stk::CommBuffer& buf = comm_all.send_buffer(p);
      buf.skip<unsigned long>(1);
    }
  }

  bool symmetric = false;
  bool local_error = false;
  bool global_bad_input = comm_all.allocate_buffers(mpi_size / 4, symmetric, local_error);
  ASSERT_EQ(global_bad_input, false);
}
