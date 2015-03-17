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

#include <exodusMeshInterface.h>

#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/unit_test_support/perf_unit_util.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace
{

void runStkSearchTestUsingStkAABoxes(stk::search::SearchMethod search);
void runStkSearchTestUsingGtkAABoxes(stk::search::SearchMethod search);
void testGtkSearch(MPI_Comm comm, std::vector<GtkBox>&domainBoxes, SearchResults& boxIdPairResults);
void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm);
void testStkSearchUsingStkAABoxes(MPI_Comm comm, std::vector<GtkBox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults);
void testStkSearchUsingGtkAABoxes(MPI_Comm comm, std::vector<GtkBox> &domainBoxes,
                stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults);

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

TEST(Performance, stkSearchUsingOcttreeUsingStkAABoxes)
{
    runStkSearchTestUsingStkAABoxes(stk::search::OCTREE);
}

TEST(Performance, stkSearchUsingBoostUsingGtkAABoxes)
{
    runStkSearchTestUsingGtkAABoxes(stk::search::BOOST_RTREE);
}

TEST(Performance, stkSearchUsingOcttreeUsingGtkAABoxes)
{
    runStkSearchTestUsingGtkAABoxes(stk::search::OCTREE);
}

TEST(Performance, gtkSearchUsingOcttreeUsingGtkAABoxes)
{
    runStkSearchTestUsingGtkAABoxes(stk::search::OCTREE);
}

TEST(Performance, gtkSearch)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    std::vector<GtkBox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;
    testGtkSearch(comm, domainBoxes, boxIdPairResults);
}

void runStkSearchTestUsingGtkAABoxes(stk::search::SearchMethod searchMethod)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<GtkBox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;
    testStkSearchUsingGtkAABoxes(comm, domainBoxes, searchMethod, boxIdPairResults);
}

void runStkSearchTestUsingStkAABoxes(stk::search::SearchMethod searchMethod)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<GtkBox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;
    testStkSearchUsingStkAABoxes(comm, domainBoxes, searchMethod, boxIdPairResults);
}

void testGtkSearch(MPI_Comm comm, std::vector<GtkBox>&domainBoxes, SearchResults& searchResults)
{
  // This is an unusual situation where these tests are both performance unit tests
  // and normal tests; therefore, they need to work regardless of whether we have
  // callgrind available or not.
#if __VALGRIND_MAJOR__
    stk::check_valgrind_version();
    CALLGRIND_START_INSTRUMENTATION;
#endif

    int proc_id = stk::parallel_machine_rank(comm);
    std::vector<int> procThatOwnsBox;

    for(size_t i=0;i<domainBoxes.size();i++)
    {
        procThatOwnsBox.push_back(proc_id);
    }

    std::vector<GtkBox> rangeBoxes(domainBoxes);

#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
#endif

    double startTime = stk::wall_time();


    std::vector<std::pair<int, int> > interaction_list;
    std::vector<int> first_interaction;
    std::vector<int> last_interaction;

    ACME::BoxA_BoxB_Par_Search(domainBoxes, rangeBoxes, comm, interaction_list, first_interaction, last_interaction);

    /*
    std::vector<int> ghost_indices;
    std::vector<int> ghost_procs;
    ACME::BoxA_BoxB_Ghost(domainBoxes, rangeBoxes, comm, ghost_indices, ghost_procs);

    int num_procs = stk::parallel_machine_size(comm);
    std::vector< std::vector<GtkBox> > send_list(num_procs);
    std::vector< std::vector<GtkBox> > recv_list(num_procs);

    for (size_t i=0;i<ghost_indices.size();i++)
    {
        send_list[ghost_procs[i]].push_back(rangeBoxes[ghost_indices[i]]);
    }

    ACME::Parallel_Data_Exchange(send_list, recv_list, comm);

    ASSERT_EQ(static_cast<size_t>(num_procs), recv_list.size());
    for (size_t i=0;i<recv_list.size();i++)
    {
        for (size_t j=0;j<recv_list[i].size();j++)
        {
            rangeBoxes.push_back(recv_list[i][j]);
            procThatOwnsBox.push_back(i);
        }
    }
    std::vector<int> interaction_list;
    std::vector<int> first_interaction;
    std::vector<int> last_interaction;

    gtk::BoxA_BoxB_Search(domainBoxes, rangeBoxes, interaction_list, first_interaction, last_interaction);
    */



    double elapsedTime = stk::wall_time() - startTime;



#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif

    printPeformanceStats(elapsedTime, comm);
    stk::print_debug_skip(comm);

    EXPECT_EQ(domainBoxes.size(), first_interaction.size());
    EXPECT_EQ(domainBoxes.size(), last_interaction.size());
    typedef std::set <std::pair<Ident,Ident> > localJunk;
    localJunk localResults;

    // Ident box1, Ident box2
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Ident box1(i, procThatOwnsBox[i]);
        for (int j=first_interaction[i];j<last_interaction[i];j++)
        {
	  Ident box2(interaction_list[j].first, interaction_list[j].second);
            localResults.insert(std::make_pair(box1, box2));
        }
    }
    std::string rangeBoxComm = unitTestUtils::getOption("-rangeBoxComm", "yes");
    bool rangeResultsCommunicated = ( rangeBoxComm == "yes" );

    if ( rangeResultsCommunicated )
    {
        localJunk tmp;
        stk::search::communicate< std::pair<Ident,Ident>, std::pair<Ident,Ident> >(comm, localResults, tmp);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(searchResults));
    }
    else
    {
        std::copy(localResults.begin(), localResults.end(), std::back_inserter(searchResults));
    }

    gatherResultstoProcZero(comm, searchResults);
    size_t goldValueNumber=getGoldValueForTest();
    if ( proc_id == 0 )
    {
        if ( goldValueNumber != 0u)
        {
            EXPECT_EQ(goldValueNumber, searchResults.size());
        }
        else
        {
            std::cerr << "Number of interactions: " << searchResults.size() << std::endl;
        }
    }
}

void testStkSearchUsingStkAABoxes(MPI_Comm comm, std::vector<GtkBox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults)
{
#if __VALGRIND_MAJOR__
    stk::check_valgrind_version();
    CALLGRIND_START_INSTRUMENTATION;
#endif

    int procId = stk::parallel_machine_rank(comm);

    StkBoxVector stkBoxes(domainBoxes.size());
    fillStkBoxesUsingGtkBoxes(domainBoxes, procId, stkBoxes);

    std::string rangeBoxComm = unitTestUtils::getOption("-rangeBoxComm", "yes");
    bool rangeResultsCommunicated = ( rangeBoxComm == "yes" );

#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
#endif

    double startTime = stk::wall_time();
    stk::search::coarse_search(stkBoxes, stkBoxes, searchMethod, comm, boxIdPairResults, rangeResultsCommunicated);
    double elapsedTime = stk::wall_time() - startTime;

#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif

    printPeformanceStats(elapsedTime, comm);
    stk::print_debug_skip(comm);

    gatherResultstoProcZero(comm, boxIdPairResults);
    size_t goldValueNumber=getGoldValueForTest();
    if ( procId == 0 )
    {
        if ( goldValueNumber != 0u)
        {
            EXPECT_EQ(goldValueNumber, boxIdPairResults.size());
        }
        else
        {
            std::cerr << "Number of interactions: " << boxIdPairResults.size() << std::endl;
        }
    }
}

void testStkSearchUsingGtkAABoxes(MPI_Comm comm, std::vector<GtkBox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults)
{
#if __VALGRIND_MAJOR__
    stk::check_valgrind_version();
    CALLGRIND_START_INSTRUMENTATION;
#endif

    int procId = stk::parallel_machine_rank(comm);

    GtkBoxVector searchBoxPairs(domainBoxes.size());
    for(size_t i=0;i<domainBoxes.size();i++)
    {
        Ident domainBoxId(i, procId);
        searchBoxPairs[i] = std::make_pair(domainBoxes[i], domainBoxId);
    }

    std::string rangeBoxComm = unitTestUtils::getOption("-rangeBoxComm", "yes");
    bool rangeResultsCommunicated = ( rangeBoxComm == "yes" );

#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
#endif

    double startTime = stk::wall_time();
    stk::search::coarse_search(searchBoxPairs, searchBoxPairs, searchMethod, comm, boxIdPairResults, rangeResultsCommunicated);
    double elapsedTime = stk::wall_time() - startTime;

#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif

    printPeformanceStats(elapsedTime, comm);
    stk::print_debug_skip(comm);

    gatherResultstoProcZero(comm, boxIdPairResults);
    size_t goldValueNumber=getGoldValueForTest();
    if ( procId == 0 )
    {
        if ( goldValueNumber != 0u)
        {
            EXPECT_EQ(goldValueNumber, boxIdPairResults.size());
        }
        else
        {
            std::cerr << "Number of interactions: " << boxIdPairResults.size() << std::endl;
        }
    }
}

TEST(Performance, getGoldResults)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = stk::parallel_machine_rank(comm);

    std::vector<GtkBox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;

    StkBoxVector stkBoxes(domainBoxes.size());
    fillStkBoxesUsingGtkBoxes(domainBoxes, procId, stkBoxes);

    double startTime = stk::wall_time();
    for (size_t i=0;i<stkBoxes.size();++i)
    {
        for (size_t j=0;j<stkBoxes.size();++j)
        {
            if ( stk::search::intersects(stkBoxes[i].first, stkBoxes[j].first) )
            {
                boxIdPairResults.push_back(std::make_pair(stkBoxes[i].second, stkBoxes[j].second));
            }
        }
    }

    std::sort(boxIdPairResults.begin(), boxIdPairResults.end());
    SearchResults::iterator iter_end = std::unique(boxIdPairResults.begin(), boxIdPairResults.end());
    boxIdPairResults.erase(iter_end, boxIdPairResults.end());

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, comm);

    std::cerr << "Number of boxes: " << boxIdPairResults.size() << std::endl;
}

}
