// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/OldCommSplitting.hpp>
#include <stdexcept>
#include <algorithm>
#include <vector>

#ifndef STK_HIDE_DEPRECATED_CODE  // delete October 2022
namespace {

TEST(UnitTestSplitComm, has_split_comm_false_when_same)
{
  EXPECT_FALSE(stk::coupling::has_split_comm(MPI_COMM_WORLD, MPI_COMM_WORLD));
}

TEST(UnitTestSplitComm, has_split_comm_true_when_different)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  EXPECT_TRUE(stk::coupling::has_split_comm(MPI_COMM_WORLD, MPI_COMM_SELF));
}

TEST(UnitTestSplitComm, split_comm_np1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);

  MPI_Comm splitComm = stk::coupling::split_comm(MPI_COMM_WORLD, myColor);
  EXPECT_FALSE(stk::coupling::has_split_comm(MPI_COMM_WORLD, splitComm));
  EXPECT_EQ(1, stk::parallel_machine_size(splitComm));
}

TEST(UnitTestSplitComm, split_comm_same_color)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  
  int myColor = 0;

  MPI_Comm splitComm = stk::coupling::split_comm(MPI_COMM_WORLD, myColor);
  EXPECT_FALSE(stk::coupling::has_split_comm(MPI_COMM_WORLD, splitComm));
  EXPECT_EQ(2, stk::parallel_machine_size(splitComm));
}

TEST(UnitTestSplitComm, split_comm)
{
  int numWorldProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numWorldProcs <= 1) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = myRank == 0 ? 0 : 1;

  MPI_Comm splitComm = stk::coupling::split_comm(MPI_COMM_WORLD, myColor);
  EXPECT_TRUE(stk::coupling::has_split_comm(MPI_COMM_WORLD, splitComm));
  int expectedSize = myRank == 0 ? 1 : (numWorldProcs - 1);
  EXPECT_EQ(expectedSize, stk::parallel_machine_size(splitComm));
}

TEST(UnitTestSplitComm, calc_my_root_and_other_root_ranks_both_comm_world)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  std::pair<int,int> rootRanks = 
      stk::coupling::calc_my_root_and_other_root_ranks(MPI_COMM_WORLD, MPI_COMM_WORLD);
  EXPECT_EQ(rootRanks.first, rootRanks.second);
}
  
TEST(UnitTestSplitComm, calc_my_root_and_other_root_ranks_both_split_comms)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);
  MPI_Comm splitComm = stk::coupling::split_comm(MPI_COMM_WORLD, myColor);

  std::pair<int,int> rootRanks = 
      stk::coupling::calc_my_root_and_other_root_ranks(splitComm, splitComm);
  EXPECT_EQ(rootRanks.first, rootRanks.second);
}

TEST(UnitTestSplitComm, calc_my_root_and_other_root_ranks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);
  MPI_Comm splitComm = stk::coupling::split_comm(MPI_COMM_WORLD, myColor);

  std::pair<int,int> rootRanks = 
      stk::coupling::calc_my_root_and_other_root_ranks(MPI_COMM_WORLD, splitComm);

  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int otherRank = 1 - myRank;
  EXPECT_EQ(myRank, rootRanks.first);
  EXPECT_EQ(otherRank, rootRanks.second);
}

TEST(UnitTestSplitComm, calc_my_root_and_other_root_ranks_non_contig_comm)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) < 3) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = myRank == 1 ? 1 : 0;
  MPI_Comm splitComm = stk::coupling::split_comm(MPI_COMM_WORLD, myColor);

  int expectedOtherRootRank = myRank == 1 ? 0 : 1;
  std::pair<int,int> rootRanks = 
      stk::coupling::calc_my_root_and_other_root_ranks(MPI_COMM_WORLD, splitComm);
  EXPECT_EQ(expectedOtherRootRank, rootRanks.second);
}

}

#endif