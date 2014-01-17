#include <valgrind/callgrind.h>

#include <mpi.h>
#include <gtest/gtest.h>
#include <vector>
#include <fstream>

#include <unit_tests/UnitTestUtils.hpp>
#include <unit_tests/MeshUtilsForBoundingVolumes.hpp>

#include <exodusMeshInterface.h>

namespace
{

void check_valgrind_version()
{
    ASSERT_EQ(__VALGRIND_MAJOR__, 3);
    ASSERT_EQ(__VALGRIND_MINOR__, 8);
}

void print_debug_skip(stk::ParallelMachine pm)
{
#ifndef NDEBUG
  // We're in debug; need to tell test script not to validate cycle count
  const size_t p_rank = stk::parallel_machine_rank(pm);
  if (p_rank == 0) {
    std::cout << "\nSTKPERF SKIP VALIDATION" << std::endl;
  }
#endif
}

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
    int proc = 0;
    MPI_Comm_rank(comm, &proc);
    int numProcs = 1;
    MPI_Comm_size(comm, &numProcs);

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
    check_valgrind_version();
    CALLGRIND_START_INSTRUMENTATION;
#endif

    int num_procs = -1;
    int proc_id   = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);
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
    std::vector<int> ghost_indices;
    std::vector<int> ghost_procs;
    ACME::BoxA_BoxB_Ghost(domainBoxes, rangeBoxes, comm, ghost_indices, ghost_procs);

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

    geometry::BoxA_BoxB_Search(domainBoxes, rangeBoxes, interaction_list, first_interaction, last_interaction);

    double elapsedTime = stk::wall_time() - startTime;

#if __VALGRIND_MAJOR__
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif

    printPeformanceStats(elapsedTime, comm);
    print_debug_skip(comm);

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
            Ident box2(interaction_list[j], procThatOwnsBox[interaction_list[j]]);
            localResults.insert(std::make_pair(box1, box2));
        }
    }
    std::string rangeBoxComm = getOption("-rangeBoxComm", "yes");
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
    check_valgrind_version();
    CALLGRIND_START_INSTRUMENTATION;
#endif

    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    int numProc=-1;
    MPI_Comm_size(comm, &numProc);

    StkBoxVector stkBoxes(domainBoxes.size());
    fillStkBoxesUsingGtkBoxes(domainBoxes, procId, stkBoxes);

    std::string rangeBoxComm = getOption("-rangeBoxComm", "yes");
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
    print_debug_skip(comm);

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
    check_valgrind_version();
    CALLGRIND_START_INSTRUMENTATION;
#endif

    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    int numProc=-1;
    MPI_Comm_size(comm, &numProc);

    GtkBoxVector searchBoxPairs(domainBoxes.size());
    for(size_t i=0;i<domainBoxes.size();i++)
    {
        Ident domainBoxId(i, procId);
        searchBoxPairs[i] = std::make_pair(domainBoxes[i], domainBoxId);
    }

    std::string rangeBoxComm = getOption("-rangeBoxComm", "yes");
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
    print_debug_skip(comm);

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
    int procId=-1;
    MPI_Comm_rank(comm, &procId);

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
