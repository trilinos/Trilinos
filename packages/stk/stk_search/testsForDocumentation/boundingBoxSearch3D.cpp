#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <gtest/gtest.h>

namespace
{
typedef stk::search::Box<double> Box;
typedef stk::search::IdentProc<int, int> BoxId;
void assertPairInResults(BoxId a, BoxId b, const std::vector<std::pair<BoxId, BoxId> > &searchResults);
TEST(StkSearchHowTo, useBoostRtreeSearch)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int myProcId = stk::parallel_machine_rank(comm);
    std::vector<std::pair<Box, BoxId> > firstList, secondList;
    Box unitBox(Box::point_type(0, 0, 0), Box::point_type(1, 1, 1));
    BoxId firstBoxId(0, myProcId);
    firstList.push_back(std::make_pair(unitBox, firstBoxId));
    BoxId secondBoxId(1, myProcId);
    secondList.push_back(std::make_pair(unitBox, secondBoxId));

    std::vector<std::pair<BoxId, BoxId> > searchResults;
    stk::search::coarse_search(firstList, secondList, stk::search::BOOST_RTREE, comm, searchResults);

    int numProc = stk::parallel_machine_size(comm);
    const size_t numMatchesFromOtherProcs = 2 * (numProc - 1);
    const size_t numMatchesThisProc = 1;
    ASSERT_EQ(numMatchesFromOtherProcs + numMatchesThisProc, searchResults.size());

    for(int procId = 0; procId < numProc; procId++)
    {
        assertPairInResults(BoxId(0, myProcId), BoxId(1, procId), searchResults);
        assertPairInResults(BoxId(0, procId), BoxId(1, myProcId), searchResults);
    }
}
void assertPairInResults(BoxId a, BoxId b, const std::vector<std::pair<BoxId, BoxId> > &searchResults)
{
    std::pair<BoxId, BoxId> expectedIntersectionPair(a, b);
    std::vector<std::pair<BoxId, BoxId> >::const_iterator resultsIter =
            std::find(searchResults.begin(), searchResults.end(), expectedIntersectionPair);
    bool foundExpectedPairInResults = resultsIter != searchResults.end();
    ASSERT_TRUE(foundExpectedPairInResults);
}
}
