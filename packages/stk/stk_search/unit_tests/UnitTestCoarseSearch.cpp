#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <algorithm>
#include <vector>
#include <iterator>
#include <cstdlib>

namespace std {
template <typename Key, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::ident::IdentProc<Key,Proc>,stk::search::ident::IdentProc<Key,Proc> > const& ident)
{
  return out << "[" << ident.first << ":" << ident.second << "]";
}
} // namespace std


namespace {

typedef stk::search::ident::IdentProc<uint64_t, unsigned> Ident;
typedef std::vector<std::pair<Ident,Ident> > SearchResults;

void checkSearchResults(const int proc_id, SearchResults &goldResults, SearchResults& searchResults, const size_t numGoldResults)
{
    Ident unUsedBox1(987654321, 1000000);
    Ident unUsedBox2(987654322, 1000000);
    size_t numResultsMatchingGoldResults = 0;

    for (size_t i = 0; i < searchResults.size(); i++ )
    {
//        if ( proc_id == 0 )
//        {
//            std::cerr << searchResults[i].first << "\t" << searchResults[i].second << std::endl;
//        }
        for (size_t j = 0; j < goldResults.size(); j++)
        {
            if ( searchResults[i] == goldResults[j] )
            {
                goldResults[j] = std::make_pair(unUsedBox1, unUsedBox2);
                numResultsMatchingGoldResults++;
            }
        }
    }
    EXPECT_EQ(numGoldResults, numResultsMatchingGoldResults) << "proc id = " << proc_id << std::endl;
}

void testCoarseSearchAABBForAlgorithm(stk::search::FactoryOrder &searchParameters, MPI_Comm comm)
{
    typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
    typedef std::vector<Box> BoxVector;

    int num_procs = stk::parallel_machine_size(comm);
    int proc_id   = stk::parallel_machine_rank(comm);

    double data[6];

    BoxVector local_domain, local_range;
    // what if identifier is NOT unique
    // x_min <= x_max
    // y_min <= y_max
    // z_min <= z_max

    data[0] = proc_id + 0.1; data[1] = 0.0; data[2] = 0.0;
    data[3] = proc_id + 0.9; data[4] = 1.0; data[5] = 1.0;

    Ident domainBox1(proc_id*4, proc_id);
    local_domain.push_back(Box(data, domainBox1));

    data[0] = proc_id + 0.1; data[1] = 2.0; data[2] = 0.0;
    data[3] = proc_id + 0.9; data[4] = 3.0; data[5] = 1.0;

    Ident domainBox2(proc_id*4+1, proc_id);
    local_domain.push_back(Box(data, domainBox2));

    data[0] = proc_id + 0.6; data[1] = 0.5; data[2] = 0.0;
    data[3] = proc_id + 1.4; data[4] = 1.5; data[5] = 1.0;

    Ident rangeBox1(proc_id*4+2, proc_id);
    local_range.push_back(Box(data, rangeBox1));

    data[0] = proc_id + 0.6; data[1] = 2.5; data[2] = 0.0;
    data[3] = proc_id + 1.4; data[4] = 3.5; data[5] = 1.0;

    Ident rangeBox2(proc_id*4+3, proc_id);
    local_range.push_back(Box(data, rangeBox2));

    SearchResults searchResults;

    stk::search::coarse_search(searchResults, local_range, local_domain, searchParameters);
    SearchResults goldResults;

    if (num_procs == 1) {
        size_t numGoldResults = 2u;
        EXPECT_TRUE(searchResults.size() >= numGoldResults);
        goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
        goldResults.push_back(std::make_pair(domainBox2, rangeBox2));

        checkSearchResults(proc_id, goldResults, searchResults, numGoldResults);
    }
    else {
        if (proc_id == 0) {
            size_t numGoldResults = 4u;
            EXPECT_TRUE(searchResults.size() >= numGoldResults);
            Ident domainBox1OnProcessor1(4,1);
            Ident domainBox2OnProcessor1(5,1);
            goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox1OnProcessor1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2));
            goldResults.push_back(std::make_pair(domainBox2OnProcessor1, rangeBox2));

            checkSearchResults(proc_id, goldResults, searchResults, numGoldResults);
        }
        else if (proc_id == num_procs - 1) {
            size_t numGoldResults = 4u;
            EXPECT_TRUE(searchResults.size() >= numGoldResults);

            Ident rangeBox1OnPreviousProcessor1((proc_id-1)*4 + 2, proc_id - 1);
            Ident rangeBox2OnPreviousProcessor1((proc_id-1)*4 + 3, proc_id - 1);

            goldResults.push_back(std::make_pair(domainBox1, rangeBox1OnPreviousProcessor1));
            goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2OnPreviousProcessor1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2));

            checkSearchResults(proc_id, goldResults, searchResults, numGoldResults);
        }
        else {
            size_t numGoldResults = 6u;
            EXPECT_TRUE(searchResults.size() >= numGoldResults);

            Ident rangeBox1OnPreviousProcessor((proc_id-1)*4 + 2, proc_id - 1);
            Ident rangeBox2OnPreviousProcessor((proc_id-1)*4 + 3, proc_id - 1);
            Ident domainBox1OnNextProcessor((proc_id+1)*4,     proc_id + 1);
            Ident domainBox2OnNextProcessor((proc_id+1)*4 + 1, proc_id + 1);

            goldResults.push_back(std::make_pair(domainBox1, rangeBox1OnPreviousProcessor));
            goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2OnPreviousProcessor));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2));
            goldResults.push_back(std::make_pair(domainBox1OnNextProcessor, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2OnNextProcessor, rangeBox2));

            checkSearchResults(proc_id, goldResults, searchResults, numGoldResults);
        }
    }
}

//  axis aligned bounding box search

STKUNIT_UNIT_TEST(stk_search_not_boost, coarse_search_3D_oct_tree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::FactoryOrder searchParameters;
  searchParameters.m_algorithm = stk::search::FactoryOrder::OCTREE;
  searchParameters.m_communicator = comm;
  testCoarseSearchAABBForAlgorithm(searchParameters, comm);
}

STKUNIT_UNIT_TEST(stk_search_not_boost, coarse_search_3D_bih_tree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::FactoryOrder searchParameters;
  searchParameters.m_algorithm = stk::search::FactoryOrder::BIHTREE;
  searchParameters.m_communicator = comm;
  testCoarseSearchAABBForAlgorithm(searchParameters, comm);
}

} //namespace
