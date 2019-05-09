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
#include <stk_util/parallel/Parallel.hpp>

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

void expect_search_results(int num_procs, int proc_id, const SearchResults&  searchResults)
{
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

  expect_search_results(num_procs, proc_id, searchResults);
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
  int id = 0;

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

void testCoarseSearchForAlgorithmUsingFloatAABoxes(NewSearchMethod algorithm, MPI_Comm comm)
{
  int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  FloatBoxVector local_domain, local_range;
  // what if identifier is NOT unique

  FloatBox box;
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

  expect_search_results(num_procs, proc_id, searchResults);
}

TEST(stk_search, coarse_search_noIdentProc_boost_rtree)
{
  testCoarseSearchForAlgorithm_IntsForIdents(stk::search::BOOST_RTREE, MPI_COMM_WORLD);
}

//  
//  rdjamis, Dec 12, 2016: Added this test w/ KDTREE and the unit test failed,
//  needs to be addressed. Leaving it off for the time being
//
// TEST(stk_search, coarse_search_noIdentProc_kdtree)
// {
//   testCoarseSearchForAlgorithm_IntsForIdents(stk::search::KDTREE, MPI_COMM_WORLD);
// }

TEST(stk_search, coarse_search_boost_rtree)
{
  testCoarseSearchForAlgorithm(stk::search::BOOST_RTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_kdtree)
{
  testCoarseSearchForAlgorithm(stk::search::KDTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_boost_rtree_using_float_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingFloatAABoxes(BOOST_RTREE, MPI_COMM_WORLD);
}

TEST(stk_search, coarse_search_kdtree_using_float_aa_boxes)
{
    testCoarseSearchForAlgorithmUsingFloatAABoxes(KDTREE, MPI_COMM_WORLD);
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
        FloatBox box1(0,0,0,1,1,1);
        FloatBox box2(0.5, 0.5, 0.5, 1.5, 1.5, 1.5);
        Ident id1(1, 0);
        Ident id2(1, 1);

        FloatBoxVector boxes;
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

TEST(stk_search, coarse_search_boost_ident_proc_switch)
{
    testIdentProcWithSearch(stk::search::BOOST_RTREE);
}

TEST(stk_search, coarse_search_kdtree_ident_proc_switch)
{
    testIdentProcWithSearch(stk::search::KDTREE);
}

void testCoarseSearchOnePoint(stk::search::SearchMethod searchMethod)
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

    stk::search::coarse_search(local_domain, local_range, searchMethod, comm, searchResults);

    if (proc_id == 0) {
      ASSERT_EQ(searchResults.size(), 1u);
      EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
    } else {
      ASSERT_EQ(searchResults.size(), 0u);
    }
}

TEST(stk_search, coarse_search_one_point_BOOST_RTREE)
{
    testCoarseSearchOnePoint(stk::search::BOOST_RTREE);
}

TEST(stk_search, coarse_search_one_point_KDTREE)
{
    testCoarseSearchOnePoint(stk::search::KDTREE);
}

void testCoarseSearchForDeterminingSharingAllAllCase(stk::search::SearchMethod searchMethod)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(comm);
    const int p_size = stk::parallel_machine_size(comm);

    typedef std::vector< std::pair<Sphere,Ident> > SphereVector;

    SphereVector source_bbox_vector;

    Point coords(1.0, 1.0, 1.0);
    double radius = 1.0e-6;
    Sphere node(coords, radius);
    uint64_t global_id = 1000 + p_rank;
    Ident id = Ident(global_id, p_rank);

    source_bbox_vector.push_back(std::make_pair(node, id));

    const bool communicateRangeBoxInfo = true;
    const bool dontCommunicateRangeBoxInfo = false;

    SearchResults searchResults;
    stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                               communicateRangeBoxInfo);
    EXPECT_EQ(2*p_size - 1, static_cast<int>(searchResults.size()));

    searchResults.clear();
    stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                               dontCommunicateRangeBoxInfo);
    EXPECT_EQ(p_size, static_cast<int>(searchResults.size()));

    std::set<int> procs;

    for(size_t i=0;i<searchResults.size();++i)
    {
        procs.insert(searchResults[i].second.proc());
        procs.insert(searchResults[i].first.proc());
    }

    std::set<int>::iterator iter = procs.begin();

    int procCounter = 0;
    for (;iter!=procs.end();++iter)
    {
        EXPECT_EQ(procCounter, *iter);
        procCounter++;
    }
}

TEST(CoarseSearch, forDeterminingSharingAllAllCase_BOOST_RTREE)
{
  testCoarseSearchForDeterminingSharingAllAllCase(stk::search::BOOST_RTREE);
}

TEST(CoarseSearch, forDeterminingSharingAllAllCase_KDTREE)
{
  testCoarseSearchForDeterminingSharingAllAllCase(stk::search::KDTREE);
}


void testCoarseSearchForDeterminingSharingLinearAdjacentCase(stk::search::SearchMethod searchMethod,
                                                             int numLoops = 1)
{
    const stk::ParallelMachine comm = MPI_COMM_WORLD;
    const int p_rank = stk::parallel_machine_rank(comm);
    const int p_size = stk::parallel_machine_size(comm);

    typedef std::vector< std::pair<Sphere,Ident> > SphereVector;

    SphereVector source_bbox_vector;

    Point coords(1.0 + p_rank, 1.0, 1.0);
    double radius = 0.6;
    Sphere node(coords, radius);
    uint64_t global_id = 1000 + p_rank;
    Ident id = Ident(global_id, p_rank);

    source_bbox_vector.push_back(std::make_pair(node, id));

    const bool communicateRangeBoxInfo = true;
    const bool dontCommunicateRangeBoxInfo = false;

    SearchResults searchResults;

    double markTime = 0.0;
    double totalTime = 0.0;

    for (int count = 0; count < numLoops; ++count) {

      const bool lastLoop = (count + 1 == numLoops);

      searchResults.clear();
      markTime = stk::wall_time();
      stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                                 communicateRangeBoxInfo);
      totalTime += stk::wall_time() - markTime;

      if (lastLoop) {
        if (p_size == 1) {
          EXPECT_EQ(1, static_cast<int>(searchResults.size()));
        }
        else {
          const int expectedSize = ((p_rank == 0) || (p_rank == p_size - 1) ? 3 : 5);
          EXPECT_EQ(expectedSize, static_cast<int>(searchResults.size()));
        }
      }

      searchResults.clear();
      markTime = stk::wall_time();
      stk::search::coarse_search(source_bbox_vector, source_bbox_vector, searchMethod, comm, searchResults,
                                 dontCommunicateRangeBoxInfo);
      totalTime += stk::wall_time() - markTime;

      if (lastLoop) {
        if (p_size == 1) {
          EXPECT_EQ(1, static_cast<int>(searchResults.size()));
        }
        else {
          const int expectedSize = ((p_rank == 0) || (p_rank == p_size - 1) ? 2 : 3);
          EXPECT_EQ(expectedSize, static_cast<int>(searchResults.size()));
        }
      }
    }

    if (p_rank == 0) {
      double avgTime =  totalTime / numLoops;
      std::cout << "Average search time measured on rank 0 of " << p_size << " ranks is " << avgTime << " through " << numLoops << " loops" << std::endl;
    }

    // Do more correctness checking after the last loop.
    std::set<int> procs;
    for(size_t i=0;i<searchResults.size();++i)
    {
      procs.insert(searchResults[i].second.proc());
      procs.insert(searchResults[i].first.proc());
    }
    std::set<int>::iterator iter = procs.begin();

    int minNeighbor = p_rank - 1;
    int maxNeighbor = p_rank + 1;

    int procCounter = 0;
    for (;iter!=procs.end();++iter)
    {
      EXPECT_LE(minNeighbor, *iter);
      EXPECT_GE(maxNeighbor, *iter);
      ++procCounter;
    }
    if (p_size == 1) {
      EXPECT_EQ(1, procCounter);
    }
    else {
      const int expectedCount = ((p_rank == 0) || (p_rank == p_size - 1) ? 2 : 3);
      EXPECT_EQ(expectedCount, static_cast<int>(procCounter));
    }
}

TEST(CoarseSearch, forDeterminingSharingLinearAdjacentCase_BOOST_RTREE)
{
  testCoarseSearchForDeterminingSharingLinearAdjacentCase(stk::search::BOOST_RTREE);
}

TEST(CoarseSearch, forDeterminingSharingLinearAdjacentCase_KDTREE)
{
  testCoarseSearchForDeterminingSharingLinearAdjacentCase(stk::search::KDTREE);
}

TEST(CoarseSearchScaling, forDeterminingSharingLinearAdjacentCase_KDTREE)
{
  testCoarseSearchForDeterminingSharingLinearAdjacentCase(stk::search::KDTREE, 1000);
}

} //namespace
