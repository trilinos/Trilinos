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

#include <unit_tests/UnitTestUtils.hpp>
#include <unit_tests/MeshUtilsForBoundingVolumes.hpp>

#include <gtest/gtest.h>

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
  typedef stk::search::Box<double> StkBox;
  typedef std::vector< std::pair<StkBox,Ident> > BoxVector;

  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  BoxVector local_domain, local_range;
  // what if identifier is NOT unique

  StkBox box;
  Ident id;

  box = StkBox( Point(proc_id + 0.1, 0.0, 0.0), Point(proc_id + 0.9, 1.0, 1.0));
  id = Ident(proc_id * 4, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box = StkBox( Point(proc_id + 0.1, 2.0, 0.0), Point(proc_id + 0.9, 3.0, 1.0));
  id = Ident(proc_id * 4+1, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box = StkBox( Point(proc_id + 0.6, 0.5, 0.0), Point(proc_id + 1.4, 1.5, 1.0));
  id = Ident(proc_id * 4+2, proc_id);
  local_range.push_back(std::make_pair(box,id));

  box = StkBox( Point(proc_id + 0.6, 2.5, 0.0), Point(proc_id + 1.4, 3.5, 1.0));
  id = Ident(proc_id * 4+3, proc_id);
  local_range.push_back(std::make_pair(box,id));

  SearchResults searchResults;

  stk::search::coarse_search(local_domain, local_range, algorithm, comm, searchResults);

  if (num_procs == 1) {
    ASSERT_EQ( searchResults.size(), 2u);
    EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
    EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
  }
  else {
    if (proc_id == 0) {
      ASSERT_EQ( searchResults.size(), 4u);
      EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
      EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
      EXPECT_EQ( searchResults[2], std::make_pair( Ident(4,1), Ident(2,0)) );
      EXPECT_EQ( searchResults[3], std::make_pair( Ident(5,1), Ident(3,0)) );
    }
    else if (proc_id == num_procs - 1) {
      ASSERT_EQ( searchResults.size(), 4u);
      int prev = proc_id -1;
      EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
    }
    else {
      ASSERT_EQ( searchResults.size(), 6u);
      int prev = proc_id -1;
      int next = proc_id + 1;
      EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
      EXPECT_EQ( searchResults[4], std::make_pair( Ident(next*4,next), Ident(proc_id*4+2,proc_id)) );
      EXPECT_EQ( searchResults[5], std::make_pair( Ident(next*4+1,next), Ident(proc_id*4+3,proc_id)) );
    }
  }
}


int IdentForTest(int id, int proc_id)
{
  return id + 1000 * proc_id;
}

void testCoarseSearchForAlgorithm_IntsForIdents(stk::search::SearchMethod algorithm, MPI_Comm comm)
{
  typedef stk::search::Point<double> Point;
  typedef stk::search::Box<double> StkBox;
  typedef std::vector< std::pair<StkBox, int > > BoxVector;

  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  BoxVector local_domain, local_range;
  // what if identifier is NOT unique

  StkBox box;
  int id;

  box = StkBox( Point(proc_id + 0.1, 0.0, 0.0), Point(proc_id + 0.9, 1.0, 1.0));
  id = IdentForTest(proc_id * 4, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box = StkBox( Point(proc_id + 0.1, 2.0, 0.0), Point(proc_id + 0.9, 3.0, 1.0));
  id = IdentForTest(proc_id * 4+1, proc_id);
  local_domain.push_back(std::make_pair(box,id));

  box = StkBox( Point(proc_id + 0.6, 0.5, 0.0), Point(proc_id + 1.4, 1.5, 1.0));
  id = IdentForTest(proc_id * 4+2, proc_id);
  local_range.push_back(std::make_pair(box,id));

  box = StkBox( Point(proc_id + 0.6, 2.5, 0.0), Point(proc_id + 1.4, 3.5, 1.0));
  id = IdentForTest(proc_id * 4+3, proc_id);
  local_range.push_back(std::make_pair(box,id));

  std::vector<std::pair<int, int> > searchResults;

  stk::search::coarse_search_nonIdentProc<StkBox, int, StkBox, int>(local_domain, local_range, algorithm, comm, searchResults);

  if (num_procs == 1) {
    ASSERT_EQ( searchResults.size(), 2u);
    EXPECT_EQ( searchResults[0].first, IdentForTest(0,0) );
    EXPECT_EQ( searchResults[1].first, IdentForTest(1,0) );
    EXPECT_EQ( searchResults[0].second,IdentForTest(2,0) );
    EXPECT_EQ( searchResults[1].second,IdentForTest(3,0) );
  }
  else {
    if (proc_id == 0) {
      ASSERT_EQ( searchResults.size(), 2u);
      EXPECT_EQ( searchResults[0].first, IdentForTest(0,0) );
      EXPECT_EQ( searchResults[1].first, IdentForTest(1,0) );
      EXPECT_EQ( searchResults[0].second, IdentForTest(2,0) );
      EXPECT_EQ( searchResults[1].second, IdentForTest(3,0) );
    }
    else if (proc_id == num_procs - 1) {
      ASSERT_EQ( searchResults.size(), 4u);
      int prev = proc_id -1;
      EXPECT_EQ( searchResults[0].first, IdentForTest(proc_id*4,proc_id) );
      EXPECT_EQ( searchResults[1].first, IdentForTest(proc_id*4,proc_id) );
      EXPECT_EQ( searchResults[2].first, IdentForTest(proc_id*4+1,proc_id) );
      EXPECT_EQ( searchResults[3].first, IdentForTest(proc_id*4+1,proc_id) );
      EXPECT_EQ( searchResults[0].second, IdentForTest(prev*4+2,prev) );
      EXPECT_EQ( searchResults[1].second, IdentForTest(proc_id*4+2,proc_id) );
      EXPECT_EQ( searchResults[2].second, IdentForTest(prev*4+3,prev) );
      EXPECT_EQ( searchResults[3].second, IdentForTest(proc_id*4+3,proc_id) );

    }
    else {
      ASSERT_EQ( searchResults.size(), 4u);
      int prev = proc_id -1;
      EXPECT_EQ( searchResults[0].first, IdentForTest(proc_id*4,proc_id) );
      EXPECT_EQ( searchResults[1].first, IdentForTest(proc_id*4,proc_id) );
      EXPECT_EQ( searchResults[2].first, IdentForTest(proc_id*4+1,proc_id) );
      EXPECT_EQ( searchResults[3].first, IdentForTest(proc_id*4+1,proc_id) );
      EXPECT_EQ( searchResults[0].second, IdentForTest(prev*4+2,prev) );
      EXPECT_EQ( searchResults[1].second, IdentForTest(proc_id*4+2,proc_id) );
      EXPECT_EQ( searchResults[2].second, IdentForTest(prev*4+3,prev) );
      EXPECT_EQ( searchResults[3].second, IdentForTest(proc_id*4+3,proc_id) );
    }
  }
}

void testCoarseSearchForAlgorithmUsingGtkAABoxes(NewSearchMethod algorithm, MPI_Comm comm)
{
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  GtkBoxVector local_domain, local_range;
  // what if identifier is NOT unique

  GtkBox box;
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

  coarse_search_new(local_domain, local_range, algorithm, comm, searchResults);

  if (num_procs == 1) {
    ASSERT_EQ( searchResults.size(), 2u);
    EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
    EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
  }
  else {
    if (proc_id == 0) {
      ASSERT_EQ( searchResults.size(), 4u);
      EXPECT_EQ( searchResults[0], std::make_pair( Ident(0,0), Ident(2,0)) );
      EXPECT_EQ( searchResults[1], std::make_pair( Ident(1,0), Ident(3,0)) );
      EXPECT_EQ( searchResults[2], std::make_pair( Ident(4,1), Ident(2,0)) );
      EXPECT_EQ( searchResults[3], std::make_pair( Ident(5,1), Ident(3,0)) );
    }
    else if (proc_id == num_procs - 1) {
      ASSERT_EQ( searchResults.size(), 4u);
      int prev = proc_id -1;
      EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
    }
    else {
      ASSERT_EQ( searchResults.size(), 6u);
      int prev = proc_id -1;
      int next = proc_id + 1;
      EXPECT_EQ( searchResults[0], std::make_pair( Ident(proc_id*4,proc_id), Ident(prev*4+2,prev)) );
      EXPECT_EQ( searchResults[1], std::make_pair( Ident(proc_id*4,proc_id), Ident(proc_id*4+2,proc_id)) );
      EXPECT_EQ( searchResults[2], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(prev*4+3,prev)) );
      EXPECT_EQ( searchResults[3], std::make_pair( Ident(proc_id*4+1,proc_id), Ident(proc_id*4+3,proc_id)) );
      EXPECT_EQ( searchResults[4], std::make_pair( Ident(next*4,next), Ident(proc_id*4+2,proc_id)) );
      EXPECT_EQ( searchResults[5], std::make_pair( Ident(next*4+1,next), Ident(proc_id*4+3,proc_id)) );
    }
  }
}

TEST(stk_search, coarse_search_noIdentProc_boost_rtree)
{
  testCoarseSearchForAlgorithm_IntsForIdents(stk::search::BOOST_RTREE, MPI_COMM_WORLD);
}

#if 0
TEST(stk_search, basic_coarse_search_octree)
{
  testCoarseSearchForAlgorithm_IntsForIdents(stk::search::OCTREE, MPI_COMM_WORLD);
}
#endif

TEST(stk_search, coarse_search_boost_rtree)
{
  testCoarseSearchForAlgorithm(stk::search::BOOST_RTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_octree)
{
  testCoarseSearchForAlgorithm(stk::search::OCTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_boost_rtree_using_gtk_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingGtkAABoxes(BOOST_RTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_octree_using_gtk_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingGtkAABoxes(OCTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_gtk_using_gtk_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingGtkAABoxes(GTK, MPI_COMM_WORLD);
}

void testIdentProcWithSearch(stk::search::SearchMethod searchMethod)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId=-1;
    MPI_Comm_rank(comm, &procId);
    int numProcs = -1;
    MPI_Comm_size(comm, &numProcs);

    if ( numProcs != 1 )
    {
        GtkBox box1(0,0,0,1,1,1);
        GtkBox box2(0.5, 0.5, 0.5, 1.5, 1.5, 1.5);
        Ident id1(1, -1);
        Ident id2(1, -1);

        GtkBoxVector boxes;
        if ( procId == 0 )
        {
          boxes.push_back(std::make_pair(box1, id1));
        }
        else if ( procId == 1 )
        {
          boxes.push_back(std::make_pair(box2, id2));
        }

        SearchResults searchResults;

        coarse_search(boxes, boxes, searchMethod, comm, searchResults);

        SearchResults goldResults;

        Ident goldId1(1, 0);
        Ident goldId2(1, 1);

        if (procId == 0 )
        {
            goldResults.push_back(std::make_pair(goldId1, goldId1));
            goldResults.push_back(std::make_pair(goldId1, goldId2));
            goldResults.push_back(std::make_pair(goldId2, goldId1));
            ASSERT_EQ(goldResults.size(), searchResults.size());
        }
        else if ( procId == 1 )
        {
            goldResults.push_back(std::make_pair(goldId1, goldId2));
            goldResults.push_back(std::make_pair(goldId2, goldId1));
            goldResults.push_back(std::make_pair(goldId2, goldId2));
            ASSERT_EQ(3u, searchResults.size());
        }

        for (size_t i=0;i<goldResults.size();i++)
        {
            EXPECT_EQ(goldResults[i], searchResults[i]) << "Test comparison for proc " << procId << " failed for comparsion #" << i << std::endl;
        }
    }
}

//TEST(stk_search, coarse_search_boost_ident_proc_switch)
//{
//    stk::search::SearchMethod sm = stk::search::BOOST_RTREE;
//    testIdentProcWithSearch(sm);
//}

TEST(stk_search, coarse_search_octree_ident_proc_switch)
{
    stk::search::SearchMethod sm = stk::search::OCTREE;
    testIdentProcWithSearch(sm);
}

TEST(stk_search, coarse_search_one_point)
{
  typedef stk::search::IdentProc<uint64_t, unsigned> Ident;
  typedef stk::search::Point<double> Point;
  typedef stk::search::Box<double> StkBox;
  typedef std::vector<std::pair<StkBox,Ident> > BoxVector;
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
    local_domain.push_back(std::make_pair(StkBox(min_corner, max_corner), domainBox1));
  }

  min_corner[0] = 0.5; min_corner[1] = 0.5; min_corner[2] = 0.5;
  max_corner[0] = 0.5; max_corner[1] = 0.5; max_corner[2] = 0.5;

  // One range target on processor 0 with the label:  1
  // All other processors have empty range.
  Ident rangeBox1(1, 0);
  if (proc_id == 0) {
    local_range.push_back(std::make_pair(StkBox(min_corner, max_corner), rangeBox1));
  }

  SearchResults searchResults;

  stk::search::coarse_search(local_domain, local_range, stk::search::OCTREE, comm, searchResults);

  if (proc_id == 0) {
    ASSERT_EQ(searchResults.size(), 1u);
    EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
  } else {
    ASSERT_EQ(searchResults.size(), 0u);
  }
}

} //namespace
