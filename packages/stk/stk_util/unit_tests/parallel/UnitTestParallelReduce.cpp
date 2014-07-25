/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/MPI.hpp>
#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_util/parallel/ParallelReduce.hpp>


#if defined ( STK_HAS_MPI )
TEST(ParallelComm, AllReduceLoc)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    stk::CommAll comm_all(comm);

    int myProcId = comm_all.parallel_rank();
    int numProcs = comm_all.parallel_size();

    int nvalues=5;
    uint64_t limit_32bit_integer = std::numeric_limits<uint32_t>::max();

    std::vector<double> values(nvalues);
    std::vector<uint64_t> locations(nvalues);
    std::vector<double> globalValues(nvalues);
    std::vector<uint64_t> globalLocations(nvalues);

    for(int n = 0; n < nvalues; n++) {
      values[n] = myProcId*numProcs+n;
      locations[n] = limit_32bit_integer + myProcId*numProcs*numProcs+n; //Want to test when outside 32-bit range for index.
    }

    stk::all_reduce_maxloc(comm, &values[0], &locations[0], &globalValues[0], &globalLocations[0], nvalues);

    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ((numProcs-1)*numProcs+n, globalValues[n]);
      EXPECT_EQ(limit_32bit_integer + (numProcs-1)*numProcs*numProcs+n, globalLocations[n]);
    }

    stk::all_reduce_minloc(comm, &values[0], &locations[0], &globalValues[0], &globalLocations[0], nvalues);

    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n, globalValues[n]);
      EXPECT_EQ(limit_32bit_integer+n, globalLocations[n]);
    }
}

TEST(ParallelComm, AllReduce)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    stk::CommAll comm_all(comm);

    int myProcId = comm_all.parallel_rank();
    int numProcs = comm_all.parallel_size();

    int nvalues=5;

    std::vector<double> values(nvalues);
    std::vector<double> globalValues(nvalues);

    for(int n = 0; n < nvalues; n++) {
      values[n] = myProcId*numProcs+n;
    }

    stk::all_reduce_max(comm, &values[0], &globalValues[0], nvalues);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ((numProcs-1)*numProcs+n, globalValues[n]);
    }

    stk::all_reduce_min(comm, &values[0], &globalValues[0], nvalues);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n, globalValues[n]);
    }

    int alpha = 0;
    for(int n = 0; n < numProcs; ++n) {
      alpha += n*numProcs;
    }

    stk::all_reduce_sum(comm, &values[0], &globalValues[0], nvalues);
    for(int n = 0; n < nvalues; n++) {
      EXPECT_EQ(n*numProcs + alpha, globalValues[n]);
    }

}
#endif

