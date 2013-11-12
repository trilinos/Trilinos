#include <stk_search/CoarseSearch2.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/TrackingAllocator.hpp>
#include <stk_util/util/memory_util.hpp>

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

struct CompareSecond
{
  template <typename First, typename Second>
  bool operator()( std::pair<First,Second> const& a, std::pair<First,Second> const& b) const
  { return a.second < b.second; }
};

STKUNIT_UNIT_TEST(stk_search, bounding_box)
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

  STKUNIT_EXPECT_EQ(coords[0], -1.0);
  STKUNIT_EXPECT_EQ(coords[1], -2.0);
  STKUNIT_EXPECT_EQ(coords[2], -3.0);
  STKUNIT_EXPECT_EQ(coords[3], 1.0);
  STKUNIT_EXPECT_EQ(coords[4], 2.0);
  STKUNIT_EXPECT_EQ(coords[5], 3.0);

  Box temp;
  stk::search::impl::set_box(temp, coords);

  STKUNIT_EXPECT_EQ(bg::distance(temp.min_corner(), min_corner), 0.0);
  STKUNIT_EXPECT_EQ(bg::distance(temp.max_corner(), max_corner), 0.0);
}

STKUNIT_UNIT_TEST(stk_search, global_spatial_index)
{
  int parallel_size = stk::parallel_machine_size(MPI_COMM_WORLD);
  if (parallel_size == 1) return;

  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;
  typedef std::pair<Box,int> BoxProc;
  const unsigned MaxVolumesPerNode = 16;
  typedef bgi::rtree< BoxProc, bgi::quadratic<MaxVolumesPerNode>, bgi::indexable<BoxProc>, bgi::equal_to<BoxProc>, stk::tracking_allocator<BoxProc, BoxProc> > Rtree;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Point min_corner(rank - 0.2, 0, 0);
  Point max_corner(rank + 1.2, 1, 1);
  Box box(min_corner, max_corner);

  Rtree tree;
  stk::search::create_global_spatial_index(tree, box, MPI_COMM_WORLD);

  STKUNIT_EXPECT_EQ(tree.size(), size_t(size));

  Box global_bounds = tree.bounds();
  STKUNIT_EXPECT_EQ(bg::distance(global_bounds.min_corner(), Point(-0.2,0,0)), 0.0);
  STKUNIT_EXPECT_EQ(bg::distance(global_bounds.max_corner(), Point(size + 0.2, 1, 1)), 0.0);

  std::vector<BoxProc> intersections;
  bgi::query(tree, bgi::intersects(box), std::back_inserter(intersections));

  std::sort(intersections.begin(), intersections.end(), CompareSecond());

  if (rank > 0 && rank < size-1) {
    STKUNIT_EXPECT_EQ(intersections.size(), 3u);
    STKUNIT_EXPECT_EQ(intersections[0].second, rank-1);
    STKUNIT_EXPECT_EQ(intersections[1].second, rank);
    STKUNIT_EXPECT_EQ(intersections[2].second, rank+1);
  }
  else if (size > 1) {
    STKUNIT_EXPECT_EQ(intersections.size(), 2u);
    if (rank == 0) {
      STKUNIT_EXPECT_EQ(intersections[0].second, rank);
      STKUNIT_EXPECT_EQ(intersections[1].second, rank+1);
    }
    else {
      STKUNIT_EXPECT_EQ(intersections[0].second, rank-1);
      STKUNIT_EXPECT_EQ(intersections[1].second, rank);
    }
  }
}

#if 0
STKUNIT_UNIT_TEST(stk_search, coarse_search2_2D)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef std::pair<int, int> Ident; // id, proc
  typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int size = stk::parallel_machine_size(comm);
  int rank = stk::parallel_machine_rank(comm);

  Box domain_box1(Point(double(rank) +0.1, 0.0), Point(double(rank) + 0.9, 1.0));
  Box domain_box2(Point(double(rank) +0.1, 2.0), Point(double(rank) + 0.9, 3.0));

  Box range_box1(Point(double(rank) +0.1 + 0.5, 0.5), Point(double(rank) + 0.9 + 0.5, 1.5));
  Box range_box2(Point(double(rank) +0.1 + 0.5, 2.5), Point(double(rank) + 0.9 + 0.5, 3.5));

  std::vector<std::pair<Box, Ident> > local_domain, local_range;

  local_domain.push_back(std::make_pair(domain_box1, std::make_pair(rank * 4, rank)));
  local_domain.push_back(std::make_pair(domain_box2, std::make_pair(rank * 4 + 1, rank)));

  local_range.push_back(std::make_pair(range_box1, std::make_pair(rank * 4 + 2, rank)));
  local_range.push_back(std::make_pair(range_box2, std::make_pair(rank * 4 + 3, rank)));

  std::vector<std::pair<Ident, Ident> > output;
  stk::search::coarse_search2(local_domain, local_range, comm, output);

  if (size == 1) {
    STKUNIT_EXPECT_EQ(output.size(), 2u);

    STKUNIT_EXPECT_EQ(output[0], std::make_pair(std::make_pair(0, 0), std::make_pair(2, 0)));
    STKUNIT_EXPECT_EQ(output[1], std::make_pair(std::make_pair(1, 0), std::make_pair(3, 0)));
  }
  else {
    if (rank == 0) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(std::make_pair(0, 0), std::make_pair(2, 0)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(std::make_pair(1, 0), std::make_pair(3, 0)));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(std::make_pair(4, 1), std::make_pair(2, 0)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(std::make_pair(5, 1), std::make_pair(3, 0)));
    }
    else if (rank == size - 1) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(std::make_pair(rank*4,     rank), std::make_pair(rank*4 - 2, rank - 1)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(std::make_pair(rank*4,     rank), std::make_pair(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(std::make_pair(rank*4 + 1, rank), std::make_pair(rank*4 - 1, rank - 1)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(std::make_pair(rank*4 + 1, rank), std::make_pair(rank*4 + 3, rank   )));
    }
    else {
      STKUNIT_EXPECT_EQ(output.size(), 6u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(std::make_pair(rank*4,         rank    ), std::make_pair(rank*4 - 2, rank - 1)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(std::make_pair(rank*4,         rank    ), std::make_pair(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(std::make_pair(rank*4 + 1,     rank    ), std::make_pair(rank*4 - 1, rank - 1)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(std::make_pair(rank*4 + 1,     rank    ), std::make_pair(rank*4 + 3, rank   )));
      STKUNIT_EXPECT_EQ(output[4], std::make_pair(std::make_pair((rank+1)*4,     rank + 1), std::make_pair(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[5], std::make_pair(std::make_pair((rank+1)*4 + 1, rank + 1), std::make_pair(rank*4 + 3, rank   )));
    }
  }
}
#endif

STKUNIT_UNIT_TEST(stk_search, coarse_search2_3D)
{
  typedef stk::search::ident::IdentProc<uint64_t, unsigned> Ident;
  typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
  typedef std::vector<Box> BoxVector;
  typedef std::vector<std::pair<Ident,Ident> > Output;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int size = stk::parallel_machine_size(comm);
  int rank = stk::parallel_machine_rank(comm);

  double data[6];

  BoxVector local_domain, local_range;
  data[0] = rank + 0.1; data[1] = 0.0; data[2] = 0.0;
  data[3] = rank + 0.9; data[4] = 1.0; data[5] = 1.0;
  local_domain.push_back(Box(data, Ident(rank*4, rank)));

  data[0] = rank + 0.1; data[1] = 2.0; data[2] = 0.0;
  data[3] = rank + 0.9; data[4] = 3.0; data[5] = 1.0;
  local_domain.push_back(Box(data, Ident(rank*4+1, rank)));

  data[0] = rank + 0.6; data[1] = 0.5; data[2] = 0.0;
  data[3] = rank + 1.4; data[4] = 1.5; data[5] = 1.0;
  local_range.push_back(Box(data, Ident(rank*4+2, rank)));

  data[0] = rank + 0.6; data[1] = 2.5; data[2] = 0.0;
  data[3] = rank + 1.4; data[4] = 3.5; data[5] = 1.0;
  local_range.push_back(Box(data, Ident(rank*4+3, rank)));

  Output output;
  stk::search::coarse_search2(local_domain, local_range, comm, output);

  if (size == 1) {
    STKUNIT_EXPECT_EQ(output.size(), 2u);

    STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(0, 0), Ident(2, 0)));
    STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(1, 0), Ident(3, 0)));
  }
  else {
    if (rank == 0) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(0, 0), Ident(2, 0)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(1, 0), Ident(3, 0)));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(Ident(4, 1), Ident(2, 0)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(Ident(5, 1), Ident(3, 0)));
    }
    else if (rank == size - 1) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(rank*4,     rank), Ident(rank*4 - 2, rank - 1)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(rank*4,     rank), Ident(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(Ident(rank*4 + 1, rank), Ident(rank*4 - 1, rank - 1)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(Ident(rank*4 + 1, rank), Ident(rank*4 + 3, rank   )));
    }
    else {
      STKUNIT_EXPECT_EQ(output.size(), 6u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(rank*4,         rank    ), Ident(rank*4 - 2, rank - 1)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(rank*4,         rank    ), Ident(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(Ident(rank*4 + 1,     rank    ), Ident(rank*4 - 1, rank - 1)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(Ident(rank*4 + 1,     rank    ), Ident(rank*4 + 3, rank   )));
      STKUNIT_EXPECT_EQ(output[4], std::make_pair(Ident((rank+1)*4,     rank + 1), Ident(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[5], std::make_pair(Ident((rank+1)*4 + 1, rank + 1), Ident(rank*4 + 3, rank   )));
    }
  }
}

} //namespace
