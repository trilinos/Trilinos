#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <gtest/gtest.h>

namespace
{
typedef stk::search::Box<double> Box;
typedef stk::search::IdentProc<int, int> BoxId;
TEST(StkSearchHowTo, useBoostRtreeSearch)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = stk::parallel_machine_rank(comm);
    std::vector< std::pair<Box, BoxId> > localDomain, localRange;
    Box unitBox(Box::point_type(0, 0, 0), Box::point_type(1, 1, 1));
    localDomain.push_back(std::make_pair(unitBox, BoxId(0, procId)));
    localRange.push_back(std::make_pair(unitBox, BoxId(1, procId)));

    std::vector< std::pair<BoxId, BoxId> > searchResults;
    stk::search::coarse_search(localDomain, localRange,
                               stk::search::BOOST_RTREE,
                               comm, searchResults);

    int numProc = stk::parallel_machine_size(comm);
    const size_t numMatchesFromOtherProcs = 2 * (numProc - 1);
    const size_t numMatchesThisProc = 1;
    ASSERT_EQ(numMatchesFromOtherProcs + numMatchesThisProc,
              searchResults.size());
}
}
