// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AllReduce.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <array>

namespace krino
{

//--------------------------------------------------------------------------
// Akri_AllReduce.hpp provides thin variadic/array wrappers over
// stk::all_reduce_{sum,min,max}. These exercise the reductions across 1-4
// ranks, where each rank contributes (rank+1) so the expected reduced values
// depend on the communicator size.
//--------------------------------------------------------------------------

TEST(AllReduce, make_array_deducesCommonTypeAndSize)
{
  const auto a = make_array(1, 2, 3);
  static_assert(std::is_same<decltype(a), const std::array<int, 3>>::value,
      "make_array should produce a std::array of the common type and arity");
  EXPECT_EQ(3u, a.size());
  EXPECT_EQ(1, a[0]);
  EXPECT_EQ(2, a[1]);
  EXPECT_EQ(3, a[2]);
}

TEST(AllReduce, variadicSumMinMax_reduceScalarsAcrossRanks)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int numRanks = stk::parallel_machine_size(comm);
  const int myRank = stk::parallel_machine_rank(comm);

  // Each rank contributes (rank+1): 1, 2, ..., numRanks.
  const int expectedSum = numRanks * (numRanks + 1) / 2;
  const int expectedMin = 1;
  const int expectedMax = numRanks;

  double sum = myRank + 1;
  all_reduce_sum(comm, sum);
  EXPECT_DOUBLE_EQ(static_cast<double>(expectedSum), sum);

  double minVal = myRank + 1;
  all_reduce_min(comm, minVal);
  EXPECT_DOUBLE_EQ(static_cast<double>(expectedMin), minVal);

  double maxVal = myRank + 1;
  all_reduce_max(comm, maxVal);
  EXPECT_DOUBLE_EQ(static_cast<double>(expectedMax), maxVal);
}

TEST(AllReduce, variadicSum_reducesMultipleScalarsIndependently)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int numRanks = stk::parallel_machine_size(comm);
  const int myRank = stk::parallel_machine_rank(comm);

  // Reduce two values in one call: (rank+1) and 2*(rank+1).
  const int expectedSum = numRanks * (numRanks + 1) / 2;

  double a = myRank + 1;
  double b = 2 * (myRank + 1);
  all_reduce_sum(comm, a, b);
  EXPECT_DOUBLE_EQ(static_cast<double>(expectedSum), a);
  EXPECT_DOUBLE_EQ(static_cast<double>(2 * expectedSum), b);
}

TEST(AllReduce, sumArray_reducesEachComponent)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int numRanks = stk::parallel_machine_size(comm);
  const int myRank = stk::parallel_machine_rank(comm);

  const int rankContribution = myRank + 1;
  const std::array<double, 3> local{ 1.0 * rankContribution, 10.0 * rankContribution, 100.0 * rankContribution };

  const std::array<double, 3> global = all_reduce_sum_array(comm, local);

  const int expectedSum = numRanks * (numRanks + 1) / 2;
  EXPECT_DOUBLE_EQ(1.0 * expectedSum, global[0]);
  EXPECT_DOUBLE_EQ(10.0 * expectedSum, global[1]);
  EXPECT_DOUBLE_EQ(100.0 * expectedSum, global[2]);
}

TEST(AllReduce, minMaxArray_reduceEachComponent)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int numRanks = stk::parallel_machine_size(comm);
  const int myRank = stk::parallel_machine_rank(comm);

  const std::array<double, 2> local{ static_cast<double>(myRank + 1), static_cast<double>(-(myRank + 1)) };

  const std::array<double, 2> globalMin = all_reduce_min_array(comm, local);
  EXPECT_DOUBLE_EQ(1.0, globalMin[0]);
  EXPECT_DOUBLE_EQ(static_cast<double>(-numRanks), globalMin[1]);

  const std::array<double, 2> globalMax = all_reduce_max_array(comm, local);
  EXPECT_DOUBLE_EQ(static_cast<double>(numRanks), globalMax[0]);
  EXPECT_DOUBLE_EQ(-1.0, globalMax[1]);
}

}
