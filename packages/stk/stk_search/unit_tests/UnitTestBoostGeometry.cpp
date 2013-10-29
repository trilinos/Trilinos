#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>
#include <iterator>


STKUNIT_UNIT_TEST(stk_search, boost_geometry)
{
  namespace bg = boost::geometry;
  namespace bgi = boost::geometry::index;

  typedef bg::model::point<double, 3, bg::cs::cartesian> Point;
  typedef bg::model::box<Point> Box;
  typedef std::pair<Point,unsigned> PointAndId;

  const unsigned MaxVolumesPerNode = 16;
  typedef bgi::rtree< PointAndId, bgi::quadratic<MaxVolumesPerNode> > Rtree;

  Rtree domain;

  for (unsigned i=0; i<10; ++i) {
    domain.insert( std::make_pair( Point(i,i,i), i) );
  }


  std::vector<PointAndId> matches;

  bgi::query( domain, bgi::intersects(Box(Point(0.5,0.5,0.5),Point(9.5,9.5,9.5))), std::back_inserter(matches));


  for (int i=0, e=matches.size(); i<e; ++i) {
    std::cout << bg::wkt(matches[i].first) << ", " << matches[i].second << std::endl;
  }
}
