#include <stk_search/PeriodicBCSearch.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/TrackingAllocator.hpp>
#include <stk_util/util/memory_util.hpp>

#include <vector>
#include <iterator>
#include <cstdlib>

namespace {

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

  // Check intersections. TODO
  // std::vector<BoxProc> intersections;
  // bgi::query(tree, bgi::intersects(box), std::back_inserter(intersections));

  // if (rank != 0 && rank != size-1) {

  //   STKUNIT_EXPECT_EQ(intersections.size(), 2u);

  //   if (intersections[0].second < intersections[1].second) {
  //     STKUNIT_EXPECT_EQ(intersections[0].second, rank-1);
  //     STKUNIT_EXPECT_EQ(intersections[1].second, rank+1);
  //   } else {
  //     STKUNIT_EXPECT_EQ(intersections[0].second, rank+1);
  //     STKUNIT_EXPECT_EQ(intersections[1].second, rank-1);
  //   }
  // }
  // else {
  //   STKUNIT_EXPECT_EQ(intersections.size(), 2u);

  //   if (rank == 0) {
  //     STKUNIT_EXPECT_EQ(intersections[0].second, 0);
  //   }
  //   else {
  //     STKUNIT_EXPECT_EQ(intersections[0].second, rank - 1);
  //   }
  // }
}

STKUNIT_UNIT_TEST(stk_search, create_parallel_domain)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef std::pair<int, int> IdType; // <proc, id>
  typedef std::pair<Point, IdType > SpatialIndexValue;
  typedef bgi::rtree< SpatialIndexValue, bgi::quadratic<16> > Rtree;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int rank = stk::parallel_machine_rank(comm);

  const int num_points = 10;
  Rtree local_domain;
  for (int i = 0; i < num_points; ++i) {
    IdType id = std::make_pair(rank, i);
    Point p(rank + i, rank + i, rank + i);
    local_domain.insert(std::make_pair(p, id));
  }

  stk::search::create_parallel_domain(local_domain, comm);

  // TODO - complete
}

std::ostream& operator<<(std::ostream& out, std::pair<std::pair<int, int>, std::pair<int, int> > const& data)
{
  out << "( (" << data.first.first << ", " << data.first.second << ") -> (" << data.second.first << ", " << data.second.second << ") )";
  return out;
}

// For testing purposes
void check_equal(std::pair<std::pair<int, int>, std::pair<int, int> > const& lhs,
                 std::pair<std::pair<int, int>, std::pair<int, int> > const& rhs)
{
  STKUNIT_EXPECT_EQ(lhs.first.first,   rhs.first.first);
  STKUNIT_EXPECT_EQ(lhs.first.second,  rhs.first.second);
  STKUNIT_EXPECT_EQ(lhs.second.first,  rhs.second.first);
  STKUNIT_EXPECT_EQ(lhs.second.second, rhs.second.second);
}

STKUNIT_UNIT_TEST(stk_search, periodic_search)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef std::pair<int, int> Ident; // id, proc
  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int size = stk::parallel_machine_size(comm);
  int rank = stk::parallel_machine_rank(comm);

  Box domain_box1(Point(double(rank) +0.1, 0.0, 0.0), Point(double(rank) + 0.9, 1.0, 0.0));
  Box domain_box2(Point(double(rank) +0.1, 2.0, 0.0), Point(double(rank) + 0.9, 3.0, 0.0));

  Box range_box1(Point(double(rank) +0.1 + 0.5, 0.5, 0.0), Point(double(rank) + 0.9 + 0.5, 1.5, 0.0));
  Box range_box2(Point(double(rank) +0.1 + 0.5, 2.5, 0.0), Point(double(rank) + 0.9 + 0.5, 3.5, 0.0));

  std::vector<std::pair<Box, Ident> > local_domain, local_range;

  local_domain.push_back(std::make_pair(domain_box1, std::make_pair(rank * 4, rank)));
  local_domain.push_back(std::make_pair(domain_box2, std::make_pair(rank * 4 + 1, rank)));

  local_range.push_back(std::make_pair(range_box1, std::make_pair(rank * 4 + 2, rank)));
  local_range.push_back(std::make_pair(range_box2, std::make_pair(rank * 4 + 3, rank)));

  std::vector<std::pair<Ident, Ident> > output;
  stk::search::periodic_search(local_domain, local_range, comm, output);

  if (size == 1) {
    STKUNIT_EXPECT_EQ(output.size(), 2u);

    check_equal(output[0], std::make_pair(std::make_pair(0, 0), std::make_pair(2, 0)));
    check_equal(output[1], std::make_pair(std::make_pair(1, 0), std::make_pair(3, 0)));
  }
  else {
    if (rank == 0) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      check_equal(output[0], std::make_pair(std::make_pair(0, 0), std::make_pair(2, 0)));
      check_equal(output[1], std::make_pair(std::make_pair(1, 0), std::make_pair(3, 0)));
      check_equal(output[2], std::make_pair(std::make_pair(4, 1), std::make_pair(2, 0)));
      check_equal(output[3], std::make_pair(std::make_pair(5, 1), std::make_pair(3, 0)));
    }
    else if (rank == size - 1) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      check_equal(output[0], std::make_pair(std::make_pair(rank*4,     rank), std::make_pair(rank*4 - 2, rank - 1)));
      check_equal(output[1], std::make_pair(std::make_pair(rank*4 + 1, rank), std::make_pair(rank*4 - 1, rank - 1)));
      check_equal(output[2], std::make_pair(std::make_pair(rank*4,     rank), std::make_pair(rank*4 + 2, rank   )));
      check_equal(output[3], std::make_pair(std::make_pair(rank*4 + 1, rank), std::make_pair(rank*4 + 3, rank   )));
    }
    else {
      STKUNIT_EXPECT_EQ(output.size(), 6u);

      check_equal(output[0], std::make_pair(std::make_pair(rank*4,         rank    ), std::make_pair(rank*4 - 2, rank - 1)));
      check_equal(output[1], std::make_pair(std::make_pair(rank*4 + 1,     rank    ), std::make_pair(rank*4 - 1, rank - 1)));
      check_equal(output[2], std::make_pair(std::make_pair(rank*4,         rank    ), std::make_pair(rank*4 + 2, rank   )));
      check_equal(output[3], std::make_pair(std::make_pair(rank*4 + 1,     rank    ), std::make_pair(rank*4 + 3, rank   )));
      check_equal(output[4], std::make_pair(std::make_pair((rank+1)*4,     rank + 1), std::make_pair(rank*4 + 2, rank   )));
      check_equal(output[5], std::make_pair(std::make_pair((rank+1)*4 + 1, rank + 1), std::make_pair(rank*4 + 3, rank   )));
    }
  }
}

STKUNIT_UNIT_TEST(stk_search, boost_geometry)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;
  typedef std::pair<Point,unsigned> PointAndId;

  {
    const unsigned MaxVolumesPerNode = 16;
    typedef bgi::rtree< PointAndId, bgi::quadratic<MaxVolumesPerNode>, bgi::indexable<PointAndId>, bgi::equal_to<PointAndId>, stk::tracking_allocator<PointAndId, PointAndId> > Rtree;

    Rtree domain;

    for (unsigned i=0; i<10; ++i) {
      domain.insert( std::make_pair( Point(i,i,i), i) );
    }

    std::vector<PointAndId> matches;

    bgi::query( domain, bgi::intersects(Box(Point(0.5,0.5,0.5),Point(9.5,9.5,9.5))), std::back_inserter(matches));

    STKUNIT_EXPECT_EQ(matches.size(), 9u);
    for (size_t i=0, e=matches.size(); i<e; ++i) {
      STKUNIT_EXPECT_EQ(matches[i].second, i+1);
    }

    std::cout << "peak: " << stk::human_bytes(stk::allocator_memory_usage<PointAndId>::peak_memory()) << std::endl;
    std::cout << "current: " << stk::human_bytes(stk::allocator_memory_usage<PointAndId>::current_memory()) << std::endl;
  }

  std::cout << "num alloc: " << stk::allocator_memory_usage<PointAndId>::num_allocations() << std::endl;
  std::cout << "num dealloc: " << stk::allocator_memory_usage<PointAndId>::num_deallocations() << std::endl;
}

STKUNIT_UNIT_TEST(stk_search, index_building)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;
  typedef std::pair<Point,unsigned> PointAndId;

  {
    const unsigned MaxVolumesPerNode = 16;
    typedef bgi::rtree< PointAndId, bgi::quadratic<MaxVolumesPerNode>, bgi::indexable<PointAndId>, bgi::equal_to<PointAndId>, stk::tracking_allocator<PointAndId, PointAndId> > Rtree;

    Rtree domain;
    const unsigned num_points = 1000000;

    boost::timer t;
    std::srand(0);
    for (unsigned i=0; i<num_points; ++i) {
      double x = std::rand() / double(RAND_MAX);
      double y = std::rand() / double(RAND_MAX);
      double z = std::rand() / double(RAND_MAX);
      domain.insert( std::make_pair( Point(x,y,z), i) );
    }
    std::cout << "time for: " << num_points << " insertions is " << t.elapsed() << std::endl;

    std::cout << "peak: " << stk::human_bytes(stk::allocator_memory_usage<PointAndId>::peak_memory()) << std::endl;
    std::cout << "current: " << stk::human_bytes(stk::allocator_memory_usage<PointAndId>::current_memory()) << std::endl;
  }

  std::cout << "num alloc: " << stk::allocator_memory_usage<PointAndId>::num_allocations() << std::endl;
  std::cout << "num dealloc: " << stk::allocator_memory_usage<PointAndId>::num_deallocations() << std::endl;
}

}
