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

#include <stk_search/CoarseSearchBoostRTree.hpp>
#include <gtest/gtest.h>

namespace {

struct CompareSecond
{
  template <typename First, typename Second>
  bool operator()( std::pair<First,Second> const& a, std::pair<First,Second> const& b) const
  { return a.second < b.second; }
};

TEST(stk_search, boost_global_spatial_index)
{
  int parallel_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (parallel_size == 1) return;

  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;
  typedef std::pair<Box,int> BoxProc;
  const unsigned MaxVolumesPerNode = 16;
  typedef bgi::rtree< BoxProc, bgi::quadratic<MaxVolumesPerNode> > Rtree;

  int size = parallel_size;
  int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

  Point min_corner(rank - 0.2, 0, 0);
  Point max_corner(rank + 1.2, 1, 1);
  Box box(min_corner, max_corner);

  Rtree tree;
  stk::search::create_global_spatial_index(tree, box, MPI_COMM_WORLD);

  EXPECT_EQ(tree.size(), size_t(size));

  Box global_bounds = tree.bounds();
  EXPECT_EQ(bg::distance(global_bounds.min_corner(), Point(-0.2,0,0)), 0.0);
  EXPECT_EQ(bg::distance(global_bounds.max_corner(), Point(size + 0.2, 1, 1)), 0.0);

  std::vector<BoxProc> intersections;
  bgi::query(tree, bgi::intersects(box), std::back_inserter(intersections));

  std::sort(intersections.begin(), intersections.end(), CompareSecond());

  if (rank > 0 && rank < size-1) {
    EXPECT_EQ(intersections.size(), 3u);
    EXPECT_EQ(intersections[0].second, rank-1);
    EXPECT_EQ(intersections[1].second, rank);
    EXPECT_EQ(intersections[2].second, rank+1);
  }
  else if (size > 1) {
    EXPECT_EQ(intersections.size(), 2u);
    if (rank == 0) {
      EXPECT_EQ(intersections[0].second, rank);
      EXPECT_EQ(intersections[1].second, rank+1);
    }
    else {
      EXPECT_EQ(intersections[0].second, rank-1);
      EXPECT_EQ(intersections[1].second, rank);
    }
  }
}

TEST(stk_search, boost_bounding_box)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;
  typedef std::pair<Box,int> BoxProc;

  Point min_corner(-1,-2,-3);
  Point max_corner(1,2,3);
  Box box(min_corner, max_corner);

  double coords[6] = {};
  stk::search::impl::fill_array(box, coords);

  EXPECT_EQ(coords[0], -1.0);
  EXPECT_EQ(coords[1], -2.0);
  EXPECT_EQ(coords[2], -3.0);
  EXPECT_EQ(coords[3], 1.0);
  EXPECT_EQ(coords[4], 2.0);
  EXPECT_EQ(coords[5], 3.0);

  Box temp;
  stk::search::impl::set_box(temp, coords);

  EXPECT_EQ(bg::distance(temp.min_corner(), min_corner), 0.0);
  EXPECT_EQ(bg::distance(temp.max_corner(), max_corner), 0.0);
}

} // namespace
