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
#include <stk_coupling/SplitComms.hpp>
#include "TestCompatibilityMode.hpp"
#include <stdexcept>
#include <algorithm>
#include <vector>

namespace {

TEST(UnitTestSplitComm, string_to_color_empty_string_throw)
{
  std::string empty;
  EXPECT_THROW(stk::coupling::string_to_color(empty), std::logic_error);
}

TEST(UnitTestSplitComm, string_to_color_is_not_random)
{
  std::string str("arbitraryString");
  int firstColor = stk::coupling::string_to_color(str);
  EXPECT_EQ(firstColor, stk::coupling::string_to_color(str));
}

TEST(UnitTestSplitComm, string_to_color_nominal)
{
  std::string appString1("myApp"), appString2("myOtherApp");
  EXPECT_NE(stk::coupling::string_to_color(appString1), stk::coupling::string_to_color(appString2));
  EXPECT_TRUE(0 <= stk::coupling::string_to_color(appString1));
  EXPECT_TRUE(0 <= stk::coupling::string_to_color(appString2));
}

TEST(UnitTestSplitComm, split_comms_my_comm_np1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  MPI_Comm splitComm = splitComms.get_split_comm();
  EXPECT_EQ(1, stk::parallel_machine_size(splitComm));
}

TEST(UnitTestSplitComm, split_comms_my_comm_same_color)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  int myColor = 0;

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  MPI_Comm splitComm = splitComms.get_split_comm();
  EXPECT_EQ(2, stk::parallel_machine_size(splitComm));
}

TEST(UnitTestSplitComm, split_comms_my_comm)
{
  int numWorldProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numWorldProcs <= 1) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = myRank == 0 ? 0 : 1;

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  MPI_Comm splitComm = splitComms.get_split_comm();
  int expectedSize = myRank == 0 ? 1 : (numWorldProcs - 1);
  EXPECT_EQ(expectedSize, stk::parallel_machine_size(splitComm));
}

TEST(UnitTestSplitComm, get_other_colors_2_colors)
{
  int numWorldProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numWorldProcs <= 1) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = myRank == 0 ? 0 : 1;

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  std::vector<int> otherColors = splitComms.get_other_colors();
  EXPECT_EQ(1u, otherColors.size());
  EXPECT_EQ((myRank == 0) ? 1 : 0, otherColors[0]);
}

TEST(UnitTestSplitComm, get_other_colors_3_colors)
{
  int numWorldProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (numWorldProcs != 3) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = myRank;

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  std::vector<int> otherColors = splitComms.get_other_colors();
  EXPECT_EQ(2u, otherColors.size());
  std::sort(otherColors.begin(), otherColors.end());
  if(myRank == 0)
  {
    EXPECT_EQ(1, otherColors[0]);
    EXPECT_EQ(2, otherColors[1]);
  }
  else if(myRank == 1)
  {
    EXPECT_EQ(0, otherColors[0]);
    EXPECT_EQ(2, otherColors[1]);
  }
  else
  {
    EXPECT_EQ(0, otherColors[0]);
    EXPECT_EQ(1, otherColors[1]);
  }
}

TEST(UnitTestSplitComm, split_comms_comm_world)
{
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  MPI_Comm commWorld = splitComms.get_parent_comm();

  int commCompareResult;
  MPI_Comm_compare(MPI_COMM_WORLD, commWorld, &commCompareResult);
  EXPECT_EQ(MPI_IDENT, commCompareResult);
}

int find_pairwise_comm_size(const std::vector<int>& colors, int color1, int color2)
{
  int count = 0;
  for(unsigned i = 0; i < colors.size(); i++) {
    if(colors[i] == color1 || colors[i] == color2) {
      count++;
    }
  }
  return count;
}

void check_round_robin_communicate(MPI_Comm comm)
{
  MPI_Request request;
  MPI_Status status;

  int recvValue = -1;
  int commSize = stk::parallel_machine_size(comm);
  int myLocalRank = stk::parallel_machine_rank(comm);
  int recvRank = (myLocalRank - 1 + commSize) % commSize;
  int sendRank = (myLocalRank + 1) % commSize;

  MPI_Irecv(&recvValue, 1, MPI_INT, recvRank, MPI_ANY_TAG, comm, &request);
  MPI_Send(&myLocalRank, 1, MPI_INT, sendRank, 0, comm);
  MPI_Wait(&request, &status);

  EXPECT_EQ(recvValue, recvRank);
}

void check_pairwise_comms(const stk::coupling::SplitComms& splitComms, const std::vector<int>& colors, int myColor)
{
  std::vector<int> uniqueColors(colors.size());
  std::copy(colors.begin(), colors.end(), uniqueColors.begin());
  stk::util::sort_and_unique(uniqueColors);

  for(int color : uniqueColors) {
    if(color == myColor) {
      EXPECT_ANY_THROW(splitComms.get_pairwise_comm(color));
    } else {
      MPI_Comm pairwiseComm = splitComms.get_pairwise_comm(color);
      int pairwiseCommSize = find_pairwise_comm_size(colors, myColor, color);
      EXPECT_EQ(stk::parallel_machine_size(pairwiseComm), pairwiseCommSize);

      check_round_robin_communicate(pairwiseComm);
    }
  }
}

TEST(UnitTestSplitComm, split_comms_2colors_pairwise_comm_equals_comm_world)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  int otherColor = 1 - myColor;
  MPI_Comm pairwiseComm = splitComms.get_pairwise_comm(otherColor);
  int result;
  MPI_Comm_compare(MPI_COMM_WORLD, pairwiseComm, &result);
  EXPECT_EQ(MPI_IDENT, result);
  MPI_Comm_compare(splitComms.get_parent_comm(), pairwiseComm, &result);
  EXPECT_EQ(MPI_IDENT, result);
}

TEST(UnitTestSplitComm, split_comms_pairwise_comm)
{
  int commSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<int> colors;

  for(int i = 0; i < commSize; i++) {
    colors.push_back(i);
  }

  int myColor = colors[myRank];
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  check_pairwise_comms(splitComms, colors, myColor);
}

TEST(UnitTestSplitComm, split_comms_pairwise_comm_two_ranks_per_color)
{
  int commSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (commSize < 4 || commSize % 2 != 0) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<int> colors;

  for(int i = 0; i < commSize; i++) {
    colors.push_back(i / 2);
  }

  int myColor = colors[myRank];
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  check_pairwise_comms(splitComms, colors, myColor);
}

TEST(UnitTestSplitComm, split_comms_pairwise_comm_random_colors_per_rank)
{
  int commSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (commSize != 3) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<int> colors = {5,8,3};

  int myColor = colors[myRank];
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  check_pairwise_comms(splitComms, colors, myColor);
}

TEST(UnitTestSplitComm, split_comms_pairwise_comm_non_consecutive_ranks_per_color)
{
  int commSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (commSize != 8) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<int> colors = {0,1,1,2,3,3,2,0};

  int myColor = colors[myRank];
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  check_pairwise_comms(splitComms, colors, myColor);
}

TEST(UnitTestSplitComm, split_comms_pairwise_comm_different_rank_counts_per_color)
{
  int commSize = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (commSize != 8) { GTEST_SKIP(); }
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<int> colors = {0,1,1,1,1,1,2,2};

  int myColor = colors[myRank];
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  check_pairwise_comms(splitComms, colors, myColor);
}

TEST(UnitTestSplitComm, get_pairwise_root_ranks_2procs)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int otherColor = 1 - myColor;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  stk::coupling::PairwiseRanks ranks = splitComms.get_pairwise_root_ranks(otherColor);
  EXPECT_EQ(myColor, ranks.localColorRoot);
  EXPECT_EQ(otherColor, ranks.otherColorRoot);
}

TEST(UnitTestSplitComm, get_pairwise_root_ranks_3proc)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  std::vector<int> otherColors = splitComms.get_other_colors();

  for(int otherColor : otherColors)
  {
    stk::coupling::PairwiseRanks pairwiseRanks = splitComms.get_pairwise_root_ranks(otherColor);

    MPI_Comm pairedComm = splitComms.get_pairwise_comm(otherColor);
    int myRank = stk::parallel_machine_rank(pairedComm);
    int otherRank = (myRank == 0) ? 1 : 0;

    EXPECT_EQ(myRank, pairwiseRanks.localColorRoot);
    EXPECT_EQ(otherRank, pairwiseRanks.otherColorRoot);
  }
}

TEST(UnitTestSplitComm, get_pairwise_root_ranks_3colors_6ranks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 6) { GTEST_SKIP(); }
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD)/2;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  std::vector<int> otherColors = splitComms.get_other_colors();

  for(int otherColor : otherColors)
  {
    stk::coupling::PairwiseRanks pairwiseRanks = splitComms.get_pairwise_root_ranks(otherColor);
    int myRoot = (myColor < otherColor) ? 0 : 2;
    int otherRoot = 2 - myRoot;

    EXPECT_EQ(myRoot, pairwiseRanks.localColorRoot);
    EXPECT_EQ(otherRoot, pairwiseRanks.otherColorRoot);
  }
}

TEST(UnitTestSplitComm, get_pairwise_root_ranks_3colors_unequal_rank_distribution)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 6) { GTEST_SKIP(); }
  std::vector<int> numRanksPerColor = {2, 1, 3};
  std::vector<int> colorOfRank = {0, 0, 1, 2, 2, 2};
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = colorOfRank[myRank];

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  std::vector<int> otherColors = splitComms.get_other_colors();

  for(int otherColor : otherColors)
  {
    stk::coupling::PairwiseRanks pairwiseRanks = splitComms.get_pairwise_root_ranks(otherColor);
    int myRoot, otherRoot;
    if (myColor < otherColor) {
      myRoot = 0;
      otherRoot = numRanksPerColor[myColor];
    } else {
      myRoot = numRanksPerColor[otherColor];
      otherRoot = 0;
    }

    EXPECT_EQ(myRoot, pairwiseRanks.localColorRoot);
    EXPECT_EQ(otherRoot, pairwiseRanks.otherColorRoot);
  }
}

TEST(UnitTestSplitComm, get_pairwise_root_ranks_3colors_unequal_noncontiguous_rank_distribution)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 6) { GTEST_SKIP(); }
  std::vector<int> colorOfRank = {2, 1, 0, 2, 0, 2};
  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int myColor = colorOfRank[myRank];

  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);
  std::vector<int> otherColors = splitComms.get_other_colors();

  for(int otherColor : otherColors)
  {
    stk::coupling::PairwiseRanks pairwiseRanks = splitComms.get_pairwise_root_ranks(otherColor);
    int myRoot = (myColor > otherColor) ? 0 : 1;
    int otherRoot = 1 - myRoot;

    EXPECT_EQ(myRoot, pairwiseRanks.localColorRoot);
    EXPECT_EQ(otherRoot, pairwiseRanks.otherColorRoot);
  }
}

TEST(UnitTestSplitComm, get_pairwise_root_ranks_2procs_invalid_color)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  int myColor = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int otherColor = 100;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, myColor);

  EXPECT_ANY_THROW(splitComms.get_pairwise_root_ranks(otherColor));
}

TEST(UnitTestSplitComm, is_coupling_version_incompatible)
{
  TestCompatibilityMode testMode(stk::coupling::impl::Incompatible);
  int color = 0;
  EXPECT_ANY_THROW(stk::coupling::SplitComms(MPI_COMM_WORLD, color));
}

TEST(UnitTestSplitComm, is_coupling_version_deprecated_current)
{
  TestCompatibilityMode testMode(stk::coupling::impl::Current);
  int color = 0;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, color);
  EXPECT_FALSE(splitComms.is_coupling_version_deprecated());
}

TEST(UnitTestSplitComm, is_coupling_version_deprecated_BackwardsCompatible)
{
  TestCompatibilityMode testMode(stk::coupling::impl::BackwardsCompatible);
  int color = 0;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, color);
  EXPECT_FALSE(splitComms.is_coupling_version_deprecated());
}

TEST(UnitTestSplitComm, is_coupling_version_deprecated_deprecated)
{
  TestCompatibilityMode testMode(stk::coupling::impl::Deprecated);
  int color = 0;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, color);
  EXPECT_TRUE(splitComms.is_coupling_version_deprecated());
}

}
