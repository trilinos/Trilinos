#include <stk_search/Box.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <unit_tests/UnitTestUtils.hpp>

#include <Geom_AxisAlignedBB.h>
#include <stk_search/CoarseSearch.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <algorithm>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <fstream>

namespace std {
template <typename Ident, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::IdentProc<Ident,Proc>,stk::search::IdentProc<Ident,Proc> > const& ip)
{
  return out << "[" << ip.first << ":" << ip.second << "]";
}
} // namespace std


namespace {

void testCoarseSearchForAlgorithm(stk::search::SearchMethod algorithm, MPI_Comm comm)
{
  typedef stk::search::Point<double> Point;
  typedef stk::search::Box<double> Box;
  typedef stk::search::IdentProc<int,int> Ident;
  typedef std::vector<std::pair<Ident,Ident> > SearchResults;
  typedef std::vector< std::pair<Box,Ident> > BoxVector;

  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  BoxVector local_domain, local_range;
  // what if identifier is NOT unique

  Box box;
  Ident id;

  box = Box( Point(proc_id + 0.1, 0.0, 0.0), Point(proc_id + 0.9, 1.0, 1.0));
  id = Ident(proc_id * 4, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box = Box( Point(proc_id + 0.1, 2.0, 0.0), Point(proc_id + 0.9, 3.0, 1.0));
  id = Ident(proc_id * 4+1, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box = Box( Point(proc_id + 0.6, 0.5, 0.0), Point(proc_id + 1.4, 1.5, 1.0));
  id = Ident(proc_id * 4+2, proc_id);
  local_range.push_back(std::make_pair(box,id));

  box = Box( Point(proc_id + 0.6, 2.5, 0.0), Point(proc_id + 1.4, 3.5, 1.0));
  id = Ident(proc_id * 4+3, proc_id);
  local_range.push_back(std::make_pair(box,id));

  SearchResults searchResults;

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResults);

  if (num_procs == 1) {
    STKUNIT_ASSERT_EQ( searchResults.size(), 2u);
    STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
    STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
  }
  else {
    if (proc_id == 0) {
      STKUNIT_ASSERT_EQ( searchResults.size(), 4u);
      STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
      STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
      STKUNIT_EXPECT_EQ( searchResults[2], std::make_pair( Ident(4,1), Ident(2,0)) );
      STKUNIT_EXPECT_EQ( searchResults[3], std::make_pair( Ident(5,1), Ident(3,0)) );
    }
    else if (proc_id == num_procs - 1) {
      STKUNIT_ASSERT_EQ( searchResults.size(), 4u);
      int prev = proc_id -1;
      STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
    }
    else {
      STKUNIT_ASSERT_EQ( searchResults.size(), 6u);
      int prev = proc_id -1;
      int next = proc_id + 1;
      STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[4], std::make_pair( Ident(next*4,next), Ident(proc_id*4+2,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[5], std::make_pair( Ident(next*4+1,next), Ident(proc_id*4+3,proc_id)) );
    }
  }
}

void testCoarseSearchForAlgorithmUsingGtkAABoxes(stk::search::SearchMethod algorithm, MPI_Comm comm)
{
  typedef geometry::AxisAlignedBB Box;
  typedef stk::search::IdentProc<int,int> Ident;
  typedef std::vector<std::pair<Ident,Ident> > SearchResults;
  typedef std::vector< std::pair<Box,Ident> > BoxVector;

  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  BoxVector local_domain, local_range;
  // what if identifier is NOT unique

  Box box;
  Ident id;

  box.set_box(proc_id + 0.1, 0.0, 0.0, proc_id + 0.9, 1.0, 1.0);
  id = Ident(proc_id * 4, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box.set_box(proc_id + 0.1, 2.0, 0.0, proc_id + 0.9, 3.0, 1.0);
  id = Ident(proc_id * 4+1, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box.set_box(proc_id + 0.6, 0.5, 0.0, proc_id + 1.4, 1.5, 1.0);
  id = Ident(proc_id * 4+2, proc_id);
  local_range.push_back(std::make_pair(box,id));

  box.set_box(proc_id + 0.6, 2.5, 0.0, proc_id + 1.4, 3.5, 1.0);
  id = Ident(proc_id * 4+3, proc_id);
  local_range.push_back(std::make_pair(box,id));

  SearchResults searchResults;

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResults);

  if (num_procs == 1) {
    STKUNIT_ASSERT_EQ( searchResults.size(), 2u);
    STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
    STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
  }
  else {
    if (proc_id == 0) {
      STKUNIT_ASSERT_EQ( searchResults.size(), 4u);
      STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
      STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
      STKUNIT_EXPECT_EQ( searchResults[2], std::make_pair( Ident(4,1), Ident(2,0)) );
      STKUNIT_EXPECT_EQ( searchResults[3], std::make_pair( Ident(5,1), Ident(3,0)) );
    }
    else if (proc_id == num_procs - 1) {
      STKUNIT_ASSERT_EQ( searchResults.size(), 4u);
      int prev = proc_id -1;
      STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
    }
    else {
      STKUNIT_ASSERT_EQ( searchResults.size(), 6u);
      int prev = proc_id -1;
      int next = proc_id + 1;
      STKUNIT_EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      STKUNIT_EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[4], std::make_pair( Ident(next*4,next), Ident(proc_id*4+2,proc_id)) );
      STKUNIT_EXPECT_EQ( searchResults[5], std::make_pair( Ident(next*4+1,next), Ident(proc_id*4+3,proc_id)) );
    }
  }
}

STKUNIT_UNIT_TEST(stk_search, coarse_search_boost_rtree)
{
  testCoarseSearchForAlgorithm(stk::search::BOOST_RTREE, MPI_COMM_WORLD);
}

STKUNIT_UNIT_TEST(stk_search, coarse_search_octree)
{
  testCoarseSearchForAlgorithm(stk::search::OCTREE, MPI_COMM_WORLD);
}

STKUNIT_UNIT_TEST(stk_search, coarse_search_boost_rtree_using_gtk_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingGtkAABoxes(stk::search::BOOST_RTREE, MPI_COMM_WORLD);
}

STKUNIT_UNIT_TEST(stk_search, coarse_search_octree_using_gtk_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingGtkAABoxes(stk::search::OCTREE, MPI_COMM_WORLD);
}


STKUNIT_UNIT_TEST(stk_search, coarse_search_one_point)
{
  typedef stk::search::IdentProc<uint64_t, unsigned> Ident;
  typedef stk::search::Point<double> Point;
  typedef stk::search::Box<double> Box;
  typedef std::vector<std::pair<Box,Ident> > BoxVector;
  typedef std::vector<std::pair<Ident,Ident> > SearchResults;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  //int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  Point min_corner, max_corner;

  BoxVector local_domain, local_range;
  // what if identifier is NOT unique
  // x_min <= x_max
  // y_min <= y_max
  // z_min <= z_max

  min_corner[0] = 0.0; min_corner[1] = 0.0; min_corner[2] = 0.0;
  max_corner[0] = 1.0; max_corner[1] = 1.0; max_corner[2] = 1.0;

  // One bounding box on processor 0 with the label:  0
  // All other processors have empty domain.
  Ident domainBox1(0, 0);
  if (proc_id == 0) {
    local_domain.push_back(std::make_pair(Box(min_corner, max_corner), domainBox1));
  }

  min_corner[0] = 0.5; min_corner[1] = 0.5; min_corner[2] = 0.5;
  max_corner[0] = 0.5; max_corner[1] = 0.5; max_corner[2] = 0.5;

  // One range target on processor 0 with the label:  1
  // All other processors have empty range.
  Ident rangeBox1(1, 0);
  if (proc_id == 0) {
    local_range.push_back(std::make_pair(Box(min_corner, max_corner), rangeBox1));
  }

  SearchResults searchResults;

  stk::search::coarse_search(local_domain, local_range, stk::search::OCTREE, comm, searchResults);

  if (proc_id == 0) {
    STKUNIT_ASSERT_EQ(searchResults.size(), 1u);
    STKUNIT_EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
  } else {
    STKUNIT_ASSERT_EQ(searchResults.size(), 0u);
  }
}

} //namespace
