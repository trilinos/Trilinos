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

//  axis aligned bounding box search

STKUNIT_UNIT_TEST(stk_search, coarse_search_3D)
{
  typedef stk::search::ident::IdentProc<uint64_t, unsigned> Ident;
  typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
  typedef std::vector<Box> BoxVector;
  typedef std::vector<std::pair<Ident,Ident> > SearchResults;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
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

  stk::search::FactoryOrder searchParameters;
  searchParameters.m_communicator = comm;
  stk::search::coarse_search(searchResults, local_range, local_domain, searchParameters);

  if (num_procs == 1) {
    STKUNIT_ASSERT_EQ(searchResults.size(), 2u);

    STKUNIT_EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
    STKUNIT_EXPECT_EQ(searchResults[1], std::make_pair(domainBox2, rangeBox2));
  }
  else {
    if (proc_id == 0) {
      STKUNIT_ASSERT_EQ(searchResults.size(), 4u);

      Ident domainBox1OnProcessor1(4,1);
      Ident domainBox2OnProcessor1(5,1);
      STKUNIT_EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
      STKUNIT_EXPECT_EQ(searchResults[1], std::make_pair(domainBox2, rangeBox2));
      STKUNIT_EXPECT_EQ(searchResults[2], std::make_pair(domainBox1OnProcessor1, rangeBox1));
      STKUNIT_EXPECT_EQ(searchResults[3], std::make_pair(domainBox2OnProcessor1, rangeBox2));
    }
    else if (proc_id == num_procs - 1) {
      STKUNIT_ASSERT_EQ(searchResults.size(), 4u);

      Ident rangeBox1OnPreviousProcessor1((proc_id-1)*4 + 2, proc_id - 1);
      Ident rangeBox2OnPreviousProcessor1((proc_id-1)*4 + 3, proc_id - 1);
      STKUNIT_EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1OnPreviousProcessor1));
      STKUNIT_EXPECT_EQ(searchResults[1], std::make_pair(domainBox1, rangeBox1));
      STKUNIT_EXPECT_EQ(searchResults[2], std::make_pair(domainBox2, rangeBox2OnPreviousProcessor1));
      STKUNIT_EXPECT_EQ(searchResults[3], std::make_pair(domainBox2, rangeBox2));
    }
    else {
      STKUNIT_ASSERT_EQ(searchResults.size(), 6u);

      Ident rangeBox1OnPreviousProcessor((proc_id-1)*4 + 2, proc_id - 1);
      Ident rangeBox2OnPreviousProcessor((proc_id-1)*4 + 3, proc_id - 1);
      Ident domainBox1OnNextProcessor((proc_id+1)*4,     proc_id + 1);
      Ident domainBox2OnNextProcessor((proc_id+1)*4 + 1, proc_id + 1);

      STKUNIT_EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1OnPreviousProcessor));
      STKUNIT_EXPECT_EQ(searchResults[1], std::make_pair(domainBox1, rangeBox1));
      STKUNIT_EXPECT_EQ(searchResults[2], std::make_pair(domainBox2, rangeBox2OnPreviousProcessor));
      STKUNIT_EXPECT_EQ(searchResults[3], std::make_pair(domainBox2, rangeBox2));
      STKUNIT_EXPECT_EQ(searchResults[4], std::make_pair(domainBox1OnNextProcessor, rangeBox1));
      STKUNIT_EXPECT_EQ(searchResults[5], std::make_pair(domainBox2OnNextProcessor, rangeBox2));
    }
  }
}

} //namespace
