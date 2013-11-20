#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_util/util/memory_util.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <mpi.h>
#include <gtest/gtest.h>
#include <vector>

#include <unit_tests/UnitTestUtils.hpp>

namespace
{

void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm);

TEST(Performance, ofAxisAlignedBoundingBoxesUsingOctTree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::OCTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

TEST(Performance, ofAxisAlignedBoundingBoxesUsingBoostRtree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::BOOST_RTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm)
{
  int proc = 0;
  MPI_Comm_rank(comm, &proc);
  int numProcs = 1;
  MPI_Comm_size(comm, &numProcs);

  size_t numColumnsPerProcessor = 1000;
  double boxSize = 1.0;

  double startTime = stk::wall_time();

  std::vector< std::pair<Box,Ident> > smallBoxVector;
  std::vector< std::pair<Box,Ident> > bigBoxVector;

  std::vector<std::pair<Ident, Ident> > boxIdPairResults;

  if(proc % 2 == 0)
  {
    double startX = numColumnsPerProcessor * boxSize * proc / 2;

    for(size_t x = 0; x < numColumnsPerProcessor; x++)
    {
      for(size_t y = 0; y < numColumnsPerProcessor; y++)
      {
        double radius = boxSize / 2;
        double centerX = startX + x * boxSize + radius;
        double centerY = y * boxSize + radius;
        double centerZ = radius;

        int id = x * numColumnsPerProcessor + y;

        smallBoxVector.push_back(generateBoundingVolume<Box>(centerX,
              centerY,
              centerZ,
              radius,
              id,
              proc));
      }
    }
  }
  else
  {
    double radius = numColumnsPerProcessor * boxSize / 2;
    double startX = numColumnsPerProcessor * boxSize * (proc - 1) / 2;
    double centerX = startX + radius;
    double centerY = radius;
    double centerZ = boxSize / 2;
    int id = 1;
    bigBoxVector.push_back(generateBoundingVolume<Box>(centerX,
          centerY,
          centerZ,
          radius - boxSize / 2,
          id,
          proc));
  }
  stk::search::coarse_search(smallBoxVector, bigBoxVector, searchMethod, comm, boxIdPairResults);

  size_t numExpectedResults = numColumnsPerProcessor * numColumnsPerProcessor;
  bool lastProcessorWithOddNumberOfProcs = numProcs % 2 != 0 && proc == numProcs - 1;
  if(lastProcessorWithOddNumberOfProcs)
  {
    numExpectedResults = 0;
  }

  EXPECT_EQ(numExpectedResults, boxIdPairResults.size());

  long int maxHwm = 0, minHwm = 0;
  double avgHwm = 0;
  stk::get_memory_high_water_mark_across_processors(comm, maxHwm, minHwm, avgHwm);

  double elapsedTime = stk::wall_time() - startTime;

  double minTime = 0, maxTime = 0, avgTime = 0;
  MPI_Allreduce(&elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&elapsedTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, comm);
  double elapsedTimeDivided = elapsedTime/numProcs;
  MPI_Allreduce(&elapsedTimeDivided, &avgTime, 1, MPI_DOUBLE, MPI_SUM, comm);

  if (proc == 0)
  {
    double bytesInMegabyte = 1024*1024;
    std::cout << "Max time: "  << maxTime << ", Min time: " << minTime << ", Avg time: " << avgTime << std::endl;
    std::cout << std::setw(6) << std::fixed << std::setprecision(1) << "Max HWM: "<<double(maxHwm)/double(bytesInMegabyte)
      <<", Min HWM: "<<double(minHwm)/double(bytesInMegabyte)<<", Avg HWM: "<<avgHwm/bytesInMegabyte<<std::endl;
    std::cout<<"### Total Number of Steps Taken ###: 1"<<std::endl;
    std::cout<<"### Total Wall Clock Run Time Used ###: "<< maxTime <<std::endl;
  }
}


} //namespace
