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

#include <gtest/gtest.h>
#include <vector>
#include <fstream>

#include <unit_tests/UnitTestUtils.hpp>
#include <unit_tests/MeshUtilsForBoundingVolumes.hpp>

#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_unit_test_utils/getOption.h>

namespace
{

void runStkSearchTestUsingStkAABoxes(stk::search::SearchMethod search);
void runStkSearchTestUsingFloatAABoxes(stk::search::SearchMethod search);
void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm);
void testStkSearchUsingStkAABoxes(MPI_Comm comm, std::vector<FloatBox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults);
void testStkSearchUsingFloatAABoxes(MPI_Comm comm, std::vector<FloatBox> &domainBoxes,
                stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults);

TEST(Performance, ofAxisAlignedBoundingBoxesUsingOctTree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::KDTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

TEST(Performance, ofAxisAlignedBoundingBoxesUsingBoostRtree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::BOOST_RTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

TEST(Performance, ofAxisAlignedBoundingBoxesUsingKdtree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::KDTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm)
{
    int proc = stk::parallel_machine_rank(comm);
    int numProcs = stk::parallel_machine_size(comm);

    size_t numColumnsPerProcessor = 1000;
    double boxSize = 1.0;

    StkBoxVector smallBoxVector;
    StkBoxVector bigBoxVector;
    SearchResults boxIdPairResults;

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

                smallBoxVector.push_back(generateBoundingVolume<StkBox>(centerX,
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
        bigBoxVector.push_back(generateBoundingVolume<StkBox>(centerX,
                centerY,
                centerZ,
                radius - boxSize / 2,
                id,
                proc));
    }

    double startTime = stk::wall_time();

    stk::search::coarse_search(smallBoxVector, bigBoxVector, searchMethod, comm, boxIdPairResults);

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, comm);

    size_t numExpectedResults = numColumnsPerProcessor * numColumnsPerProcessor;
    bool lastProcessorWithOddNumberOfProcs = numProcs % 2 != 0 && proc == numProcs - 1;
    if(lastProcessorWithOddNumberOfProcs)
    {
        numExpectedResults = 0;
    }

    EXPECT_EQ(numExpectedResults, boxIdPairResults.size());
}

////////////////////////////////////////////////////////////

TEST(Performance, stkSearchUsingBoostUsingStkAABoxes)
{
    runStkSearchTestUsingStkAABoxes(stk::search::BOOST_RTREE);
}

TEST(Performance, stkSearchUsingKdtreeUsingStkAABoxes)
{
    runStkSearchTestUsingStkAABoxes(stk::search::KDTREE);
}

TEST(Performance, stkSearchUsingBoostUsingFloatAABoxes)
{
    runStkSearchTestUsingFloatAABoxes(stk::search::BOOST_RTREE);
}

TEST(Performance, stkSearchUsingKdtreeUsingFloatAABoxes)
{
    runStkSearchTestUsingFloatAABoxes(stk::search::KDTREE);
}

void runStkSearchTestUsingFloatAABoxes(stk::search::SearchMethod searchMethod)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<FloatBox> domainBoxes( fillDomainBoxes(comm) );

    SearchResults boxIdPairResults;
    testStkSearchUsingFloatAABoxes(comm, domainBoxes, searchMethod, boxIdPairResults);
}

void runStkSearchTestUsingStkAABoxes(stk::search::SearchMethod searchMethod)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<FloatBox> domainBoxes( fillDomainBoxes(comm) );

    SearchResults boxIdPairResults;
    testStkSearchUsingStkAABoxes(comm, domainBoxes, searchMethod, boxIdPairResults);
}

void testStkSearchUsingStkAABoxes(MPI_Comm comm, std::vector<FloatBox> &domainBoxes, stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults)
{
    int procId = stk::parallel_machine_rank(comm);

    StkBoxVector stkBoxes(domainBoxes.size());
    fillStkBoxesUsingFloatBoxes(domainBoxes, procId, stkBoxes);

    std::string rangeBoxComm = stk::unit_test_util::get_option("-rangeBoxComm", "yes");
    bool rangeResultsCommunicated = ( rangeBoxComm == "yes" );

    double startTime = stk::wall_time();
    stk::search::coarse_search(stkBoxes, stkBoxes, searchMethod, comm, boxIdPairResults, rangeResultsCommunicated);
    double elapsedTime = stk::wall_time() - startTime;

    printPeformanceStats(elapsedTime, comm);

    gatherResultstoProcZero(comm, boxIdPairResults);
    size_t goldValueNumber = getGoldValueForTest();
    if ( procId == 0 )
    {
        if ( goldValueNumber != 0u)
        {
            EXPECT_EQ(goldValueNumber, boxIdPairResults.size());
        }
        std::cerr << "Number of interactions: " << boxIdPairResults.size() << std::endl;
    }
}

void testStkSearchUsingFloatAABoxes(MPI_Comm comm, std::vector<FloatBox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults)
{
    int procId = stk::parallel_machine_rank(comm);

    FloatBoxVector searchBoxPairs(domainBoxes.size());
    for(size_t i = 0; i < domainBoxes.size(); i++)
    {
        Ident domainBoxId(i, procId);
        searchBoxPairs[i] = std::make_pair(domainBoxes[i], domainBoxId);
    }

    std::string rangeBoxComm = stk::unit_test_util::get_option("-rangeBoxComm", "yes");
    bool rangeResultsCommunicated = ( rangeBoxComm == "yes" );

    double startTime = stk::wall_time();
    stk::search::coarse_search(searchBoxPairs, searchBoxPairs, searchMethod, comm, boxIdPairResults, rangeResultsCommunicated);
    double elapsedTime = stk::wall_time() - startTime;

    printPeformanceStats(elapsedTime, comm);

    gatherResultstoProcZero(comm, boxIdPairResults);
    size_t goldValueNumber = getGoldValueForTest();
    if ( procId == 0 )
    {
        if ( goldValueNumber != 0u)
        {
            EXPECT_EQ(goldValueNumber, boxIdPairResults.size());
        }
        std::cerr << "Number of interactions: " << boxIdPairResults.size() << std::endl;
    }
}

}
