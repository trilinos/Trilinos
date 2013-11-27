#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <gtest/gtest.h>

namespace
{
typedef stk::search::Box<double> Box;
typedef stk::search::IdentProc<int, int> Id;
void assertPairInResults(Id a, Id b, const std::vector<std::pair<Id, Id> > &searchResults);
TEST(StkSearchHowTo, useBoostRtreeSearch)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int myProcId = stk::parallel_machine_rank(comm);
    std::vector<std::pair<Box, Id> > firstList, secondList;
    Box unitBox(Box::point_type(0, 0, 0), Box::point_type(1, 1, 1));
    Id firstId(0, myProcId);
    firstList.push_back(std::make_pair(unitBox, firstId));
    Id secondId(1, myProcId);
    secondList.push_back(std::make_pair(unitBox, secondId));

    std::vector<std::pair<Id, Id> > searchResults;
    stk::search::coarse_search(firstList, secondList, stk::search::BOOST_RTREE, comm, searchResults);

    int numProc = stk::parallel_machine_size(comm);
    for(int procId = 0; procId < numProc; procId++)
    {
        assertPairInResults(Id(0, myProcId), Id(1, procId), searchResults);
        assertPairInResults(Id(0, procId), Id(1, myProcId), searchResults);
    }
}
TEST(StkSearchHowTo, useSphereAndPointBoundingVolumes)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int myProcId = stk::parallel_machine_rank(comm);
    std::vector<std::pair<stk::search::Sphere<double>, Id> > firstList;
    stk::search::Point<double> center(0, 0, 0);
    const double radius = 0.5;
    stk::search::Sphere<double> unitSphere(center, radius);
    Id firstId(0, myProcId);
    firstList.push_back(std::make_pair(unitSphere, firstId));
    std::vector<std::pair<stk::search::Point<double>, Id> > secondList;
    stk::search::Point<double> point(0.1, 0.2, 0.3);
    Id secondId(1, myProcId);
    secondList.push_back(std::make_pair(point, secondId));

    std::vector<std::pair<Id, Id> > searchResults;
    stk::search::coarse_search(firstList, secondList, stk::search::BOOST_RTREE, comm, searchResults);

    int numProc = stk::parallel_machine_size(comm);
    for(int procId = 0; procId < numProc; procId++)
    {
        assertPairInResults(Id(0, myProcId), Id(1, procId), searchResults);
        assertPairInResults(Id(0, procId), Id(1, myProcId), searchResults);
    }
}
void assertPairInResults(Id a, Id b, const std::vector<std::pair<Id, Id> > &searchResults)
{
    std::pair<Id, Id> expectedIntersectionPair(a, b);
    std::vector<std::pair<Id, Id> >::const_iterator resultsIter =
            std::find(searchResults.begin(), searchResults.end(), expectedIntersectionPair);
    bool foundExpectedPairInResults = resultsIter != searchResults.end();
    ASSERT_TRUE(foundExpectedPairInResults);
}
}
