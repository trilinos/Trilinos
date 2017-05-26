// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <unit_tests/UnitTestUtils.hpp>
#include <unit_tests/MeshUtilsForBoundingVolumes.hpp>
#include <stk_search/CommonSearchUtilsInstrumented.hpp>
#include <stk_search/PrototypeSearchUtilsInstrumented.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <fstream>

namespace std {
template <typename Ident, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::IdentProc<Ident,Proc>,stk::search::IdentProc<Ident,Proc> > const& ip)
{
  return out << "[" << ip.first << ":" << ip.second << "]";
}
} // namespace std


namespace {


void testGhostingSearchLinearAdjacentCase(stk::search::SearchMethod searchMethod,
                                                             int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(comm);
    const int p_size = stk::parallel_machine_size(comm);

    typedef std::vector< std::pair<Sphere,Ident> > SphereVector;
    SphereVector source_bbox_vector;

    Point coords(1.0 + p_rank, 1.0, 1.0);
    double radius = 0.6;
    Sphere node(coords, radius);
    uint64_t global_id = 1000 + p_rank;
    Ident id = Ident(global_id, p_rank);

    source_bbox_vector.push_back(std::make_pair(node, id));

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<Ident> rangeGhostIdentifiers;

    stk::search::SplitTimer timer;
    stk::search::instrumented::GhostingSearchTimeBreakdown timeBreakdown;
    timeBreakdown.reset();

    double totalTime = 0;

    for (int count = 0; count < numLoops; ++count) {

      const bool lastLoop = (count + 1 == numLoops);

      rangeGhostIdentifiers.clear();
      timer.split();

      stk::search::instrumented::ComputeRangeWithGhostsForCoarseSearch(
                                    source_bbox_vector, source_bbox_vector, p_size,
                                    rangeBoxes, rangeGhostIdentifiers, comm, timeBreakdown);
      totalTime += timer.split();

      if (lastLoop) {
        if (p_size == 1) {
          EXPECT_EQ(0, static_cast<int>(rangeGhostIdentifiers.size()));
        }
        else {
          const int expectedSize = ((p_rank == 0) || (p_rank == p_size - 1) ? 1 : 2);
          EXPECT_EQ(expectedSize, static_cast<int>(rangeGhostIdentifiers.size()));
        }
      }

    }

    if (p_rank == 0) {
      double avgTime =  totalTime / numLoops;
      std::cout << "Average ghosting search time measured on rank 0 of " << p_size << " ranks is "
                << avgTime << " through " << numLoops << " loops" << std::endl;
      std::cout << "Accumulated splits:" << std::endl;
      timeBreakdown.streamit(std::cout);
    }
}

TEST(InstrumentedGhostingSearch, LinearAdjacentCase_KDTREE)
{
  testGhostingSearchLinearAdjacentCase(stk::search::KDTREE, 1000);
}


struct Tiling2D {
  int tileXSize;
  int tileYSize;
  int xCount;
  int yCount;
};

Tiling2D computeTileSizeForPow2Ranks(int xDim, int yDim, int numRanks)
{
  int whatLg = 0;
  int currProduct = 1;
  while (currProduct >= 0 && currProduct < numRanks) {
    currProduct *= 2;
    ++ whatLg;
  }

  if (currProduct != numRanks) {
    throw("Problem size is not a power of 2.");
  }

  Tiling2D tiling = {0,0,0,0};
  if ((whatLg % 2) == 0) {
    tiling.xCount = tiling.yCount = 1 << (whatLg / 2);
  }
  else
  {
    tiling.yCount = 1 << (whatLg / 2);
    tiling.xCount = 2 * tiling.yCount;
  }

  double xDim_d = xDim;
  double yDim_d = yDim;
  tiling.tileXSize = ceil(xDim_d / tiling.xCount);
  tiling.tileYSize = ceil(yDim_d / tiling.yCount);

  return tiling;

}

enum ScalingArray2DRankType {
  INTERIOR_RANK,
  HORIZONTAL_EDGE_RANK,
  SIDE_EDGE_RANK,
  CORNER_RANK,
  ONLY_RANK,
  HORIZONTAL_ARRAY_END_RANK,
  HORIZONTAL_ARRAY_INTERIOR_RANK,
  VERTICAL_ARRAY_END_RANK,
  VERTICAL_ARRAY_INTERIOR_RANK
};

ScalingArray2DRankType computeScalingArray2DRankType(const Tiling2D &tiling, const int procRank)
{
  if (tiling.xCount == 1) {
    if (tiling.yCount == 1) {
      return ONLY_RANK;
    }
    else if ((procRank == 0) || (procRank == tiling.yCount - 1)){
      return VERTICAL_ARRAY_END_RANK;
    }
    else {
      return VERTICAL_ARRAY_INTERIOR_RANK;
    }
  }
  else if (tiling.yCount == 1) {
    if ((procRank == 0) || (procRank == tiling.xCount - 1)){
      return HORIZONTAL_ARRAY_END_RANK;
    }
    else {
      return HORIZONTAL_ARRAY_INTERIOR_RANK;
    }
  }
  else {
    bool leftEdgeOrRightEdge =
        ((procRank < tiling.yCount) || ((procRank / tiling.yCount) == (tiling.xCount - 1)) );

    bool topEdgeOrBottomEdge =
        (((procRank % tiling.yCount ) == 0)
             || ((procRank % tiling.yCount) == (tiling.yCount - 1)));

    if (leftEdgeOrRightEdge && topEdgeOrBottomEdge) {
      return CORNER_RANK;
    }
    else if (topEdgeOrBottomEdge) {
      return HORIZONTAL_EDGE_RANK;
    }
    else if (leftEdgeOrRightEdge) {
      return SIDE_EDGE_RANK;
    }
    else {
      return INTERIOR_RANK;
    }
  }
}

int computeExpectedRangeGhostIdentifiersSize(const Tiling2D &tiling, const int procRank)
{
  ScalingArray2DRankType rankType = computeScalingArray2DRankType(tiling, procRank);

  switch (rankType) {
    case INTERIOR_RANK:
      return 2 *( tiling.tileXSize + tiling.tileYSize + 2);
    case HORIZONTAL_EDGE_RANK:
      return 2 * tiling.tileYSize + tiling.tileXSize + 2;
    case SIDE_EDGE_RANK:
      return 2 * tiling.tileXSize + tiling.tileYSize + 2;
    case CORNER_RANK:
      return tiling.tileXSize + tiling.tileYSize + 1;
    case ONLY_RANK:
      return 0;
    case HORIZONTAL_ARRAY_END_RANK:
      return tiling.tileYSize;
    case VERTICAL_ARRAY_END_RANK:
      return tiling.tileXSize;
    case HORIZONTAL_ARRAY_INTERIOR_RANK:
      return 2 * tiling.tileYSize;
    case VERTICAL_ARRAY_INTERIOR_RANK:
      return 2 * tiling.tileXSize;
    default:
      return -1;
  }
}

typedef std::vector< std::pair<Sphere,Ident> > SphereVector;

enum GhostingSearchScalableTestType {
  SINGLE_PLATE_TEST,
  TWO_PLATE_TEST
};

void setupProblemForSinglePlateGhostingSearchScalingTest(int xDim, int yDim,
                                                      SphereVector &source_bbox_vector,
                                                      Tiling2D &tiling)
{
  const stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank(comm);
  const int p_size = stk::parallel_machine_size(comm);

  const double spacing = 1.1;
  const double sphereRadius = 1.0;

  tiling = computeTileSizeForPow2Ranks(xDim, yDim, p_size);
  const int myFirstGid  = p_size * tiling.tileXSize * tiling.tileYSize;
  double myFirstX       = p_rank / tiling.yCount * tiling.tileXSize * spacing;
  double myFirstY       = p_rank % tiling.yCount * tiling.tileYSize * spacing;

  int global_id = myFirstGid;
  for (int i = 0; i < tiling.tileXSize; ++i) {
    double x = i * spacing + myFirstX;
    for (int j = 0; j < tiling.tileYSize; ++j)
    {
      double y = j * spacing + myFirstY;
      Point coords(x, y, 10.0);
      Sphere node(coords, sphereRadius);
      Ident id = Ident(global_id, p_rank);

      source_bbox_vector.push_back(std::make_pair(node, id));
      ++global_id;
    }
  }
}

void printGhostingSearchTestTiming(
    double totalTime, int numLoops, const int p_size, int xDim, int yDim,
    const Tiling2D& tiling,
    stk::search::instrumented::GhostingSearchTimeBreakdown timeBreakdown) {
  double avgTime = totalTime / numLoops;
  std::cout << "Average ghosting search time measured on rank 0 of " << p_size
            << " ranks is " << avgTime << " through " << numLoops
            << " loops on problem size " << xDim << "x" << yDim << " with "
            << tiling.tileXSize << "x" << tiling.tileYSize << " domains"
            << std::endl;
  std::cout << "Accumulated splits:" << std::endl;
  timeBreakdown.streamit(std::cout);
}

int computeExpectedTwoPlateRangeGhostIdentifiersSize(const Tiling2D &tiling, const int procRankDiv2)
{
  const ScalingArray2DRankType rankType = computeScalingArray2DRankType(tiling, procRankDiv2);
  const int numSpheresFacing = tiling.tileXSize *tiling.tileYSize;

  switch (rankType) {
    case INTERIOR_RANK:
      return 4 *( tiling.tileXSize + tiling.tileYSize + 2) + numSpheresFacing;
    case HORIZONTAL_EDGE_RANK:
      return 4 * tiling.tileYSize + 2 * tiling.tileXSize + 4 + numSpheresFacing;
    case SIDE_EDGE_RANK:
      return 4 * tiling.tileXSize + 2 * tiling.tileYSize + 4 + numSpheresFacing;
    case CORNER_RANK:
      return 2 * tiling.tileXSize + 2 * tiling.tileYSize + 2 + numSpheresFacing;
    case ONLY_RANK:
      return numSpheresFacing;
    case HORIZONTAL_ARRAY_END_RANK:
      return 2 * tiling.tileYSize + numSpheresFacing;
    case VERTICAL_ARRAY_END_RANK:
      return 2 * tiling.tileXSize + numSpheresFacing;
    case HORIZONTAL_ARRAY_INTERIOR_RANK:
      return 4 * tiling.tileYSize + numSpheresFacing;
    case VERTICAL_ARRAY_INTERIOR_RANK:
      return 4 * tiling.tileXSize + numSpheresFacing;
    default:
      return -1;
  }
}

void setupProblemForTwoPlateGhostingSearchScalingTest(int xDim, int yDim,
                                                      SphereVector &source_bbox_vector,
                                                      Tiling2D &tiling)
{
  const stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank(comm);
  const int p_size = stk::parallel_machine_size(comm);

  if (p_size < 2) {
    std::cout << "Trivial pass!  testGhostingTwoPlateSearchStrongScalable(.) requires # of MPI ranks >= 2"
              << " and to be an odd power of 2, i.e., 2, 8, 32,..." << std::endl;
    return;
  }

  const int pRankDiv2 = p_rank / 2;
  const int pSizeDiv2 = p_size / 2;

  const double spacing = 1.1;
  const double sphereRadius  = 1.0;

  tiling = computeTileSizeForPow2Ranks(xDim, yDim, pSizeDiv2);
  const int myFirstGid  = p_size * tiling.tileXSize * tiling.tileYSize;
  double myFirstX       = pRankDiv2 / tiling.yCount * tiling.tileXSize * spacing;
  double myFirstY       = pRankDiv2 % tiling.yCount * tiling.tileYSize * spacing;
  double myZ            = spacing * (p_rank % 2);

  int global_id = myFirstGid;
  for (int i = 0; i < tiling.tileXSize; ++i) {
    double x = i * spacing + myFirstX;
    for (int j = 0; j < tiling.tileYSize; ++j)
    {
      double y = j * spacing + myFirstY;
      Point coords(x, y, myZ);
      Sphere node(coords, sphereRadius);
      Ident id = Ident(global_id, p_rank);

      source_bbox_vector.push_back(std::make_pair(node, id));
      ++global_id;
    }
  }
}

void testGhostingSearchStrongScalable(GhostingSearchScalableTestType testType,
                                      stk::search::SearchMethod searchMethod,
                                      int xDim, int yDim,
                                      int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int pRank = stk::parallel_machine_rank(comm);
    const int pSize = stk::parallel_machine_size(comm);
    const int pRankDiv2 = pRank / 2;

    SphereVector source_bbox_vector;
    Tiling2D tiling = {0, 0, 0, 0};
    if (testType == SINGLE_PLATE_TEST) {
      setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }
    else if (testType == TWO_PLATE_TEST) {
      setupProblemForTwoPlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<Ident> rangeGhostIdentifiers;

    stk::search::SplitTimer timer;
    stk::search::instrumented::GhostingSearchTimeBreakdown timeBreakdown;
    timeBreakdown.reset();

    double totalTime = 0;

    for (int count = 0; count < numLoops; ++count) {
      const bool lastLoop = (count + 1 == numLoops);
      rangeGhostIdentifiers.clear();
      timer.split();
      stk::search::instrumented::ComputeRangeWithGhostsForCoarseSearch(
                                    source_bbox_vector, source_bbox_vector, pSize,
                                    rangeBoxes, rangeGhostIdentifiers, comm, timeBreakdown);
      totalTime += timer.split();
      if (lastLoop) {
        if (pSize == 1) {
          EXPECT_EQ(0, static_cast<int>(rangeGhostIdentifiers.size()));
        }
        else {
          int expectedSize = 0;
          if (testType == SINGLE_PLATE_TEST) {
            expectedSize = computeExpectedRangeGhostIdentifiersSize(tiling, pRank);
          }
          else if (testType == TWO_PLATE_TEST) {
            expectedSize = computeExpectedTwoPlateRangeGhostIdentifiersSize(tiling, pRankDiv2);
          }
          EXPECT_EQ(expectedSize, static_cast<int>(rangeGhostIdentifiers.size()));
        }
      }
    }
    if (pRank == 0) {
      printGhostingSearchTestTiming(totalTime, numLoops, pSize, xDim, yDim, tiling, timeBreakdown);
    }
}

TEST(InstrumentedGhostingSearch, StrongScalable_KDTREE)
{
  testGhostingSearchStrongScalable(SINGLE_PLATE_TEST, stk::search::KDTREE, 1024, 1024);
}


TEST(InstrumentedGhostingSearch, TwoPlateStrongScalable_KDTREE)
{
  testGhostingSearchStrongScalable(TWO_PLATE_TEST, stk::search::KDTREE, 1024, 1024);
}

void testPrototypeGhostingSearchStrongScalable(GhostingSearchScalableTestType testType,
                                               stk::search::SearchMethod searchMethod,
                                               int xDim, int yDim,
                                               int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int pRank = stk::parallel_machine_rank(comm);
    const int pSize = stk::parallel_machine_size(comm);
    const int pRankDiv2 = pRank / 2;

    SphereVector source_bbox_vector;
    Tiling2D tiling = {0, 0, 0, 0};
    if (testType == SINGLE_PLATE_TEST) {
      setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }
    else if (testType == TWO_PLATE_TEST) {
      setupProblemForTwoPlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<Ident> rangeGhostIdentifiers;

    stk::search::SplitTimer timer;
    stk::search::instrumented::GhostingSearchTimeBreakdown timeBreakdown;
    timeBreakdown.reset();

    stk::search::experimental::GhostingSearcher<Ident, Ident, Sphere, Sphere>
      ghostSearcher(source_bbox_vector, source_bbox_vector,
                    rangeBoxes, rangeGhostIdentifiers, comm, timeBreakdown);

    double totalTime = 0;

    for (int count = 0; count < numLoops; ++count) {
      const bool lastLoop = (count + 1 == numLoops);
      rangeGhostIdentifiers.clear();
      timer.split();
      ghostSearcher.searchFromScratch(source_bbox_vector, source_bbox_vector,
                                      rangeBoxes, rangeGhostIdentifiers,timeBreakdown);
      totalTime += timer.split();
      if (lastLoop) {
        if (pSize == 1) {
          EXPECT_EQ(0, static_cast<int>(rangeGhostIdentifiers.size()));
        }
        else {
          int expectedSize = 0;
          if (testType == SINGLE_PLATE_TEST) {
            expectedSize = computeExpectedRangeGhostIdentifiersSize(tiling, pRank);
          }
          else if (testType == TWO_PLATE_TEST) {
            expectedSize = computeExpectedTwoPlateRangeGhostIdentifiersSize(tiling, pRankDiv2);
          }
          EXPECT_EQ(expectedSize, static_cast<int>(rangeGhostIdentifiers.size()));
        }
      }
    }
    if (pRank == 0) {
      printGhostingSearchTestTiming(totalTime, numLoops, pSize, xDim, yDim, tiling, timeBreakdown);
    }
}

TEST(PrototypeGhostingSearch, StrongScalable_KDTREE)
{
  testPrototypeGhostingSearchStrongScalable(SINGLE_PLATE_TEST, stk::search::KDTREE, 1024, 1024);
}


TEST(PrototypeGhostingSearch, TwoPlateStrongScalable_KDTREE)
{
  testPrototypeGhostingSearchStrongScalable(TWO_PLATE_TEST, stk::search::KDTREE, 1024, 1024);
}

}
