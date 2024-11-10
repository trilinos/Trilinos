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

#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
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

#define TEST_GRID_DIAMETER 1024
// #define TEST_GRID_DIAMETER 8

namespace std {
template <typename Ident, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::IdentProc<Ident,Proc>,stk::search::IdentProc<Ident,Proc> > const& ip)
{
  return out << "[" << ip.first << ":" << ip.second << "]";
}
} // namespace std


namespace {

const double TestSphereSpacing = 1.0;
const double TestSphereRadius  = 0.85;

typedef std::vector< std::pair<Sphere,IdentProc> > SphereVector;


void checkTile3DExpected(const stk::search::experimental::SuperTile3D &tile, double xLen, double yLen, double zLen,
                         int numXTiles, int numYTiles, int numZTiles) {
  EXPECT_FLOAT_EQ(xLen, tile.xLen);
  EXPECT_FLOAT_EQ(yLen, tile.yLen);
  EXPECT_FLOAT_EQ(zLen, tile.zLen);
  EXPECT_EQ(numXTiles, tile.numXTiles);
  EXPECT_EQ(numYTiles, tile.numYTiles);
  EXPECT_EQ(numZTiles, tile.numZTiles);
}

TEST(PrototypeGhostingSearch, OptimizeTile3D) {

  stk::search::experimental::SuperTile3D tile_4_1_1 = {4, 1, 1, 1, 1, 1};
  optimizeTiling3D(tile_4_1_1, 4);
  checkTile3DExpected(tile_4_1_1, 1, 1, 1, 4, 1, 1);

  stk::search::experimental::SuperTile3D tile_1_4_1 = {1, 4, 1, 1, 1, 1};
  optimizeTiling3D(tile_1_4_1, 4);
  checkTile3DExpected(tile_1_4_1, 1, 1, 1, 1, 4, 1);

  stk::search::experimental::SuperTile3D tile_1_1_4 = {1, 1, 4, 1, 1, 1};
  optimizeTiling3D(tile_1_1_4, 4);
  checkTile3DExpected(tile_1_1_4, 1, 1, 1, 1, 1, 4);

  stk::search::experimental::SuperTile3D tile_3_4_5 = {3, 4, 5, 1, 1, 1};
  optimizeTiling3D(tile_3_4_5, 60);
  checkTile3DExpected(tile_3_4_5, 1, 1, 1, 3, 4, 5);

  stk::search::experimental::SuperTile3D tile_6_4_5 = {6, 4, 5, 1, 1, 1};
  optimizeTiling3D(tile_6_4_5, 120);
  checkTile3DExpected(tile_6_4_5, 1, 1, 1, 6, 4, 5);

  stk::search::experimental::SuperTile3D tile_3_4_2 = {3, 4, 2, 1, 1, 1};
  optimizeTiling3D(tile_3_4_2, 24);
  checkTile3DExpected(tile_3_4_2, 1, 1, 1, 3, 4, 2);
}


void testGhostingSearchLinearAdjacentCase(stk::search::SearchMethod searchMethod,
                                                             int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(comm);
    const int p_size = stk::parallel_machine_size(comm);

    SphereVector source_bbox_vector;

    Point coords(1.0 + p_rank, 1.0, 1.0);
    double radius = 0.6;
    Sphere node(coords, radius);
    uint64_t global_id = 1000 + p_rank;
    IdentProc id = IdentProc(global_id, p_rank);

    source_bbox_vector.push_back(std::make_pair(node, id));

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<IdentProc> rangeGhostIdentifiers;

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


bool myIsPow2(int num) {
  if (num < 1) {
    return false;
  }
  int product = 1;
  while (product >= 0 && product < num) {
    product *= 2;
  }
  return (num == product);
}

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


enum GhostingSearchScalableTestType {
  SINGLE_PLATE_TEST,
  TWO_PLATE_TEST
};

void setupProblemForSinglePlateGhostingSearchScalingTest(int xDim, int yDim,
                                                      SphereVector &source_bbox_vector,
                                                      Tiling2D &tiling,
                                                      Point offset = Point(),
                                                      int idOffset = 0)
{
  const stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank(comm);
  const int p_size = stk::parallel_machine_size(comm);

  const double spacing      = TestSphereSpacing;
  const double sphereRadius = TestSphereRadius;

  tiling = computeTileSizeForPow2Ranks(xDim, yDim, p_size);
  const int myFirstGid  = p_size * tiling.tileXSize * tiling.tileYSize + idOffset;
  double myFirstX       = p_rank / tiling.yCount * tiling.tileXSize * spacing + offset[0];
  double myFirstY       = p_rank % tiling.yCount * tiling.tileYSize * spacing + offset[1];

  // if ((p_rank == 0) || (p_rank == 63)) {
  //   std::cout << "p_" << p_rank << ": tiling xCount=" << tiling.xCount << " yCount=" << tiling.xCount
  //             << "   LL sph ctr_xy is (" << myFirstX << " " << myFirstY
  //             << ")  UR sph ctr_xy is (" << (myFirstX + spacing * (tiling.tileXSize - 1))
  //             << " " << (myFirstY + spacing * (tiling.tileYSize - 1)) << ")" << std::endl;
  // }

  int global_id = myFirstGid;
  for (int i = 0; i < tiling.tileXSize; ++i) {
    double x = i * spacing + myFirstX;
    for (int j = 0; j < tiling.tileYSize; ++j)
    {
      double y = j * spacing + myFirstY;
      Point coords(x, y, 10.0);
      Sphere node(coords, sphereRadius);
      IdentProc id = IdentProc(global_id, p_rank);

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

void printGhostingSearchTestTiming(
    double totalTime, int numLoops, const int p_size, int xDim, int yDim,
    const Tiling2D& tiling,
    stk::search::experimental::GhostingSearchTimeBreakdown timeBreakdown) {
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
                                                      Tiling2D &tiling,
                                                      Point offset = Point(),
                                                      int idOffset = 0)
{
  const stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank(comm);
  const int p_size = stk::parallel_machine_size(comm);

  if (p_size < 2) {
    std::cout << "Trivial pass!  testGhostingTwoPlateSearchStrongScalable(.) requires # of MPI ranks >= 2"
              << std::endl;
    return;
  }

  const int pRankDiv2 = p_rank / 2;
  const int pSizeDiv2 = p_size / 2;

  const double spacing       = TestSphereSpacing;
  const double sphereRadius  = TestSphereRadius;

  tiling = computeTileSizeForPow2Ranks(xDim, yDim, pSizeDiv2);
  const int myFirstGid  = p_size * tiling.tileXSize * tiling.tileYSize + idOffset;
  double myFirstX       = pRankDiv2 / tiling.yCount * tiling.tileXSize * spacing + offset[0];
  double myFirstY       = pRankDiv2 % tiling.yCount * tiling.tileYSize * spacing + offset[1];
  double myZ            = spacing * (p_rank % 2) + offset[2];

  int global_id = myFirstGid;
  for (int i = 0; i < tiling.tileXSize; ++i) {
    double x = i * spacing + myFirstX;
    for (int j = 0; j < tiling.tileYSize; ++j)
    {
      double y = j * spacing + myFirstY;
      Point coords(x, y, myZ);
      Sphere node(coords, sphereRadius);
      IdentProc id = IdentProc(global_id, p_rank);

      source_bbox_vector.push_back(std::make_pair(node, id));
      ++global_id;
    }
  }
}

void testGhostingSearchFindRangeStrongScalable(GhostingSearchScalableTestType testType,
                                      stk::search::SearchMethod searchMethod,
                                      int xDim, int yDim,
                                      int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int pRank = stk::parallel_machine_rank(comm);
    const int pSize = stk::parallel_machine_size(comm);
    const int pRankDiv2 = pRank / 2;

    bool parmsArePowersOfTwo = (myIsPow2(pSize) && myIsPow2(xDim) && myIsPow2(yDim));
    EXPECT_TRUE(parmsArePowersOfTwo);
    if (!parmsArePowersOfTwo) {
      // Bail because other things will throw.
      return;
    }

    SphereVector source_bbox_vector;
    Tiling2D tiling = {0, 0, 0, 0};
    if (testType == SINGLE_PLATE_TEST) {
      setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }
    else if (testType == TWO_PLATE_TEST) {
      setupProblemForTwoPlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<IdentProc> rangeGhostIdentifiers;

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

TEST(InstrumentedGhostingSearch, OnePlateFindRangeStrongScalable_KDTREE)
{
  testGhostingSearchFindRangeStrongScalable(SINGLE_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER);
}


TEST(InstrumentedGhostingSearch, TwoPlateFindRangeStrongScalable_KDTREE)
{
  testGhostingSearchFindRangeStrongScalable(TWO_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER);
}

void testPrototypeGhostingSearchFindRangeStrongScalable(GhostingSearchScalableTestType testType,
                                               stk::search::SearchMethod searchMethod,
                                               int xDim, int yDim,
                                               stk::search::experimental::FindNeighborsAlgorithmChoice findNeighborsAlg,
                                               int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int pRank = stk::parallel_machine_rank(comm);
    const int pSize = stk::parallel_machine_size(comm);
    const int pRankDiv2 = pRank / 2;

    bool parmsArePowersOfTwo = (myIsPow2(pSize) && myIsPow2(xDim) && myIsPow2(yDim));
    EXPECT_TRUE(parmsArePowersOfTwo);
    if (!parmsArePowersOfTwo) {
      // Bail because other things will throw.
      return;
    }

    SphereVector source_bbox_vector;
    Tiling2D tiling = {0, 0, 0, 0};
    if (testType == SINGLE_PLATE_TEST) {
      setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }
    else if (testType == TWO_PLATE_TEST) {
      setupProblemForTwoPlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
    }

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<IdentProc> rangeGhostIdentifiers;

    stk::search::SplitTimer timer;
    stk::search::experimental::GhostingSearchTimeBreakdown timeBreakdown;
    timeBreakdown.reset();

    stk::search::experimental::GhostingSearcher<IdentProc, IdentProc, Sphere, Sphere>
      ghostSearcher(source_bbox_vector, source_bbox_vector,
                    rangeBoxes, rangeGhostIdentifiers, comm, timeBreakdown);

    double totalTime = 0;

    for (int count = 0; count < numLoops; ++count) {
      const bool lastLoop = (count + 1 == numLoops);
      rangeGhostIdentifiers.clear();
      timer.split();
      ghostSearcher.searchFromScratch(source_bbox_vector, source_bbox_vector,
                                      rangeBoxes, rangeGhostIdentifiers, timeBreakdown, findNeighborsAlg);
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

TEST(PrototypeGhostingSearch, OnePlateFindRangeStrongScalable_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeStrongScalable(SINGLE_PLATE_TEST, stk::search::KDTREE,
                                                     TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                     SCALABLE_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, TwoPlateFindFangeStrongScalable_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeStrongScalable(TWO_PLATE_TEST, stk::search::KDTREE,
                                                     TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                     SCALABLE_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, OnePlateFindRangeStrongScalableUseDDRobust_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeStrongScalable(SINGLE_PLATE_TEST, stk::search::KDTREE,
                                                     TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                     SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, TwoPlateFindFangeStrongScalableUseDDRobust_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeStrongScalable(TWO_PLATE_TEST, stk::search::KDTREE,
                                                     TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                     SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS);
}


void testPrototypeGhostingSearchFindRangeDecompRobust(GhostingSearchScalableTestType testType,
                                                          stk::search::SearchMethod searchMethod,
                                                          int xDim, int yDim,
                                                          stk::search::experimental::FindNeighborsAlgorithmChoice findNeighborsAlg,
                                                          int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int pRank = stk::parallel_machine_rank(comm);
    const int pSize = stk::parallel_machine_size(comm);
    const int pRankDiv2 = pRank / 2;

    bool parmsArePowersOfTwo = (myIsPow2(pSize) && myIsPow2(xDim) && myIsPow2(yDim));
    EXPECT_TRUE(parmsArePowersOfTwo);
    if (!parmsArePowersOfTwo) {
      // Bail because other things will throw.
      return;
    }

    if (pSize <= 8) {
      if (pRank == 0) {
        std::cout << "p_0:  Test body skipped (all ranks): > 8 mpi ranks required for real test." << std::endl;
      }
      return;
    }

    const Point cloneOffset(xDim * TestSphereSpacing + 10, yDim + 10, 0.0);
    const int   cloneIdOffset = 2 * xDim * yDim + 1000000;

    SphereVector source_bbox_vector;
    Tiling2D tiling = {0, 0, 0, 0};
    if (testType == SINGLE_PLATE_TEST) {
      setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
      setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling,
                                                          cloneOffset, cloneIdOffset);
    }
    else if (testType == TWO_PLATE_TEST) {
      setupProblemForTwoPlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling);
      setupProblemForTwoPlateGhostingSearchScalingTest(xDim, yDim, source_bbox_vector, tiling,
                                                       cloneOffset, cloneIdOffset);
    }

    std::vector<Sphere> rangeBoxes( source_bbox_vector.size() );
    std::vector<IdentProc> rangeGhostIdentifiers;

    stk::search::SplitTimer timer;
    stk::search::experimental::GhostingSearchTimeBreakdown timeBreakdown;
    timeBreakdown.reset();

    stk::search::experimental::GhostingSearcher<IdentProc, IdentProc, Sphere, Sphere>
      ghostSearcher(source_bbox_vector, source_bbox_vector,
                    rangeBoxes, rangeGhostIdentifiers, comm, timeBreakdown);

    double totalTime = 0;

    for (int count = 0; count < numLoops; ++count) {
      const bool lastLoop = (count + 1 == numLoops);
      rangeGhostIdentifiers.clear();
      timer.split();
      ghostSearcher.searchFromScratch(source_bbox_vector, source_bbox_vector,
                                      rangeBoxes, rangeGhostIdentifiers,timeBreakdown, findNeighborsAlg);
      totalTime += timer.split();
      if (lastLoop) {
        if (pSize == 1) {
          EXPECT_EQ(0, static_cast<int>(rangeGhostIdentifiers.size()));
        }
        else {
          int expectedSize = 0;
          if (testType == SINGLE_PLATE_TEST) {
            expectedSize = 2 * computeExpectedRangeGhostIdentifiersSize(tiling, pRank);
          }
          else if (testType == TWO_PLATE_TEST) {
            expectedSize = 2 * computeExpectedTwoPlateRangeGhostIdentifiersSize(tiling, pRankDiv2);
          }
          if (findNeighborsAlg == stk::search::experimental::SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS) {
            EXPECT_EQ(expectedSize, static_cast<int>(rangeGhostIdentifiers.size()));
          }
          else { // Verify that non-robust algorithm finds more neighbors than necessary.
            EXPECT_LT(expectedSize, static_cast<int>(rangeGhostIdentifiers.size()));
          }
        }
      }
    }
    if (pRank == 0) {
      printGhostingSearchTestTiming(totalTime, numLoops, pSize, xDim, yDim, tiling, timeBreakdown);
    }
}

TEST(PrototypeGhostingSearch, OnePlateFindRangeDecompRobust_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeDecompRobust(SINGLE_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                   SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, TwoPlateFindRangeDecompRobust_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeDecompRobust(TWO_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                   SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, OnePlateFindRangeNonDecompRobustFindsMoreNeighbors_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeDecompRobust(SINGLE_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                   SCALABLE_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, TwoPlateFindRangeNonDecompRobustFindsMoreNeighbors_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeDecompRobust(TWO_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                   SCALABLE_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, OnePlateFindRangeNonDecompRobustFindsMoreNeighborsOrigAlg_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeDecompRobust(SINGLE_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                   ORIGINAL_FIND_NEIGHBORS);
}

TEST(PrototypeGhostingSearch, TwoPlateFindRangeNonDecompRobustFindsMoreNeighborsOrigAlg_KDTREE)
{
  using namespace stk::search::experimental;
  testPrototypeGhostingSearchFindRangeDecompRobust(TWO_PLATE_TEST, stk::search::KDTREE, TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                   ORIGINAL_FIND_NEIGHBORS);
}

int computeExpectedOddEvenNeighbors(const Tiling2D &tiling, const int procRank)
{
  ScalingArray2DRankType rankType = computeScalingArray2DRankType(tiling, procRank);

  switch (rankType) {
    case INTERIOR_RANK:
      return 6;
    case HORIZONTAL_EDGE_RANK:
      return 3;
    case SIDE_EDGE_RANK:
      return 4;
    case CORNER_RANK:
      return 2;
    case ONLY_RANK:
      return 0;
    case HORIZONTAL_ARRAY_END_RANK:
      return 1;
    case VERTICAL_ARRAY_END_RANK:
      return 1;
    case HORIZONTAL_ARRAY_INTERIOR_RANK:
      return 2;
    case VERTICAL_ARRAY_INTERIOR_RANK:
      return 2;
    default:
      return -1;
  }
}


void testPrototypeFindGhostingNeighborsFromsScratchSourceAndRange(
    int xDim, int yDim,
    stk::search::experimental::FindNeighborsAlgorithmChoice findNeighbors)
{
  const stk::ParallelMachine comm = MPI_COMM_WORLD;
  const int pRank = stk::parallel_machine_rank(comm);
  const int pSize = stk::parallel_machine_size(comm);

  bool parmsArePowersOfTwo = (myIsPow2(pSize) && myIsPow2(xDim) && myIsPow2(yDim));
  EXPECT_TRUE(parmsArePowersOfTwo);
  if (!parmsArePowersOfTwo) {
    // Bail because other things will throw.
    return;
  }

  SphereVector source_bbox_vector;

  Tiling2D tiling = {0, 0, 0, 0};
  SphereVector localDomain, localRange;
  if ((pRank % 2) == 0) {
    setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, localDomain, tiling);
  }
  else {
    setupProblemForSinglePlateGhostingSearchScalingTest(xDim, yDim, localRange, tiling);
  }

  typedef typename Sphere::value_type domainValueType;
  typedef stk::search::Box<domainValueType>  DomainBox;
  typedef typename Sphere::value_type  rangeValueType;
  typedef stk::search::Box<rangeValueType>   RangeBox;

  std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > boxA_proc_box_array(pSize);
  std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > boxB_proc_box_array(pSize);

  stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc;
  boxA_proc.set_object_number(pRank);
  const unsigned numBoxDomain = localDomain.size();
  for(unsigned iboxA = 0; iboxA < numBoxDomain; ++iboxA) {
    stk::search::add_to_box(boxA_proc.GetBox(), localDomain[iboxA].first);
  }
  boxA_proc_box_array[pRank] = boxA_proc;

  stk::search::ObjectBoundingBox_T<RangeBox> boxB_proc;
  boxB_proc.set_object_number(pRank);
  const unsigned numBoxRange = localRange.size();
  for(unsigned iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    stk::search::add_to_box(boxB_proc.GetBox(), localRange[iboxB].first);
  }
  boxB_proc_box_array[pRank] = boxB_proc;

  stk::search::SplitTimer timer;
  stk::search::experimental::GhostingSearchTimeBreakdown timeBreakdown;
  timeBreakdown.reset();

  if ((findNeighbors == stk::search::experimental::SCALABLE_FIND_NEIGHBORS)
      || (findNeighbors == stk::search::experimental::SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS)) {
    if (findNeighbors == stk::search::experimental::SCALABLE_FIND_NEIGHBORS) {
      stk::search::experimental::findGhostingNeighborsFromScratch(comm, pRank, pSize, boxA_proc, boxB_proc,
                                                                  boxA_proc_box_array, boxB_proc_box_array,
                                                                  timeBreakdown, timer);
      int expectedOppPolarityNeighbors = computeExpectedOddEvenNeighbors(tiling, pRank);
      int expectedNumSrcNeighbors = ((pRank % 2) ? expectedOppPolarityNeighbors : 0);
      int expectedNumRngNeighbors = ((pRank % 2) ? 0 : expectedOppPolarityNeighbors);
      EXPECT_EQ(expectedNumSrcNeighbors, static_cast<int>(boxA_proc_box_array.size()));
      EXPECT_EQ(expectedNumRngNeighbors, static_cast<int>(boxB_proc_box_array.size()));
    }
    else {
      std::vector<DomainBox> localDomainBoxes;
      std::vector<RangeBox>  localRangeBoxes;
      stk::search::experimental::ComputeBoxVector(localDomain, localDomainBoxes);
      stk::search::experimental::ComputeBoxVector(localRange, localRangeBoxes);
      findGhostingNeighborsFromScratchDDEfficient(comm, pRank, pSize, boxA_proc, boxB_proc,
                                                  localDomainBoxes, localRangeBoxes,
                                                  boxA_proc_box_array, boxB_proc_box_array,
                                                  timeBreakdown, timer);
      if ((pRank %2) == 0) {
       EXPECT_EQ(0u, boxA_proc_box_array.size());
        bool rangeNeighborsAllOdd = true;
        for (auto &obb : boxB_proc_box_array) {
          int nbr = obb.get_object_number();
          if ((nbr % 2) == 0) {
            rangeNeighborsAllOdd = false;
          }
          EXPECT_TRUE(rangeNeighborsAllOdd);
        }
      }
      else {
        EXPECT_EQ(0u, boxB_proc_box_array.size());
         bool sourceNeighborsAllEven = true;
         for (auto &obb : boxA_proc_box_array) {
           int nbr = obb.get_object_number();
           if ((nbr % 2) != 0) {
             sourceNeighborsAllEven= false;
           }
           EXPECT_TRUE(sourceNeighborsAllEven);
         }
      }
    }
  }
  else {
    std::cerr << "Unexpected findNeighbors choice: " << findNeighbors << std::endl;
  }

}

TEST(PrototypeGhostingSearch, FindGhostingNeighborsFromScratchSrcRngScalable)
{
  using namespace stk::search::experimental;
  testPrototypeFindGhostingNeighborsFromsScratchSourceAndRange(TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                               SCALABLE_FIND_NEIGHBORS);
}

// Coming soon!
TEST(PrototypeGhostingSearch, FindGhostingNeighborsFromsScratchSrcRngScalableDDEfficient)
{
  using namespace stk::search::experimental;
  testPrototypeFindGhostingNeighborsFromsScratchSourceAndRange(TEST_GRID_DIAMETER, TEST_GRID_DIAMETER,
                                                               SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS);
}

}
