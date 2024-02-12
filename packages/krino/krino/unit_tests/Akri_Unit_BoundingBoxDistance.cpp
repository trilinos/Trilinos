#include <Akri_BoundingBox.hpp>
#include <Akri_BoundingBoxDistance.hpp>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <random>

namespace krino {

template<class REAL, unsigned DIM>
REAL compute_closest_distance_squared_to_points(const stk::math::Vec<REAL,DIM> & queryPt, const std::vector<stk::math::Vec<REAL,DIM>> & searchPts)
{
  REAL minSqrDist = std::numeric_limits<REAL>::max();
  for (auto & searchPt : searchPts)
    minSqrDist = std::min(minSqrDist, (queryPt-searchPt).length_squared());
  return minSqrDist;
}

template<class REAL, unsigned DIM>
REAL compute_max_closest_distance_squared_between_points(const std::vector<stk::math::Vec<REAL,DIM>> & queryPts, const std::vector<stk::math::Vec<REAL,DIM>> & searchPts)
{
  REAL maxMinSqrDist = 0;
  for (auto & queryPt : queryPts)
    maxMinSqrDist = std::max(maxMinSqrDist, compute_closest_distance_squared_to_points(queryPt, searchPts));
  return maxMinSqrDist;
}

template<class REAL, unsigned DIM>
BoundingBox_T<REAL,DIM> compute_point_bounding_box(const std::vector<stk::math::Vec<REAL,DIM>> & pts)
{
  BoundingBox_T<REAL,DIM> bbox;
  for(auto & pt : pts)
    bbox.accommodate(pt);
  return bbox;
}

template<class REAL, unsigned DIM>
void test_min_and_max_possible_closest_distance_between_each_query_pt_and_search_pts_via_bounding_box(const std::vector<stk::math::Vec<REAL,DIM>> & queryPts, const std::vector<stk::math::Vec<REAL,DIM>> & searchPts)
{
  const BoundingBox_T<REAL,DIM> searchBbox = compute_point_bounding_box(searchPts);

  for (auto & queryPt : queryPts)
  {
    const auto lowerBoundMinSqrDist = min_possible_closest_squared_distance(searchBbox, queryPt);
    const auto upperBoundMinSqrDist = max_possible_closest_squared_distance(searchBbox, queryPt);
    const auto actualMinSqrDist = compute_closest_distance_squared_to_points(queryPt, searchPts);

    EXPECT_LE(lowerBoundMinSqrDist, actualMinSqrDist);
    EXPECT_GE(upperBoundMinSqrDist, actualMinSqrDist);
  }
}

template<class REAL, unsigned DIM>
stk::math::Vec<REAL,DIM> create_point_in_box(const stk::math::Vec<REAL,DIM> & param, const BoundingBox_T<REAL,DIM> & ptBbox)
{
  const auto & boxMin = ptBbox.get_min();
  const auto & boxMax = ptBbox.get_max();
  stk::math::Vec<REAL,DIM> pt = stk::math::Vec<REAL,DIM>::ZERO;
  for (unsigned i=0; i<DIM; ++i)
    pt[i] = (1.-param[i])*boxMin[i]+param[i]*boxMax[i];
  return pt;
}

template<class REAL, unsigned DIM>
std::vector<stk::math::Vec<REAL,DIM>> create_random_points_in_bounding_box(const size_t numPts, const BoundingBox_T<REAL,DIM> & ptBbox)
{
  std::mt19937 mt(std::mt19937::default_seed);
  std::uniform_real_distribution<REAL> dist(0.,1.);

  std::vector<stk::math::Vec<REAL,DIM>> pts;
  pts.reserve(numPts);
  for (size_t i=0; i<numPts; ++i)
  {
    const stk::math::Vec<REAL,DIM> param(dist(mt),dist(mt),dist(mt));
    pts.push_back(create_point_in_box(param, ptBbox));
  }

  return pts;
}

TEST(PointToBoxDistanceBounds, randomPointsDistributedInBoxToRandomPointsInSmallerBox_eachQueryWithinBoundsOfQueryToSearchBox)
{
  const BoundingBox queryPtBbox{stk::math::Float3d{-3,-3,-3}, stk::math::Float3d{3,3,3}};
  const BoundingBox searchPtBBox{stk::math::Float3d{-1,-1,-1}, stk::math::Float3d{1,1,1}};
  const std::vector<stk::math::Float3d> queryPts = create_random_points_in_bounding_box(3*3*3*20, queryPtBbox);
  const std::vector<stk::math::Float3d> searchPts  = create_random_points_in_bounding_box(20, searchPtBBox);
  test_min_and_max_possible_closest_distance_between_each_query_pt_and_search_pts_via_bounding_box(queryPts, searchPts);
}

template<class REAL>
std::array<stk::math::Vec<REAL,3>,2> get_box_corners(const unsigned config, const BoundingBox_T<REAL,3> & box)
{
  const auto & boxMin = box.get_min();
  const auto & boxMax = box.get_max();
  if (0 == config)
    return {stk::math::Vec<REAL,3>{boxMin[0],boxMin[1],boxMin[2]}, stk::math::Vec<REAL,3>{boxMax[0],boxMax[1],boxMax[2]}};
  else if (1 == config)
    return {stk::math::Vec<REAL,3>{boxMin[0],boxMin[1],boxMax[2]}, stk::math::Vec<REAL,3>{boxMax[0],boxMax[1],boxMin[2]}};
  else if (2 == config)
    return {stk::math::Vec<REAL,3>{boxMin[0],boxMax[1],boxMin[2]}, stk::math::Vec<REAL,3>{boxMax[0],boxMin[1],boxMax[2]}};
  return {stk::math::Vec<REAL,3>{boxMax[0],boxMin[1],boxMin[2]}, stk::math::Vec<REAL,3>{boxMin[0],boxMax[1],boxMax[2]}};
}

template<class REAL, unsigned DIM>
void test_max_possible_closest_distance_between_each_minimal_configuration_via_bounding_boxes(const BoundingBox_T<REAL,DIM> & box1, const BoundingBox_T<REAL,DIM> & box2)
{
  const auto upperBoundMinSqrDist = max_possible_closest_squared_distance_between_contained_points(box1, box2);

  REAL overallMaxMinSqrDist = 0.;
  for (int box1Config=0; box1Config<4; ++box1Config)
  {
    const auto box1Corners = get_box_corners(box1Config, box1);
    for (int box2Config=0; box2Config<4; ++box2Config)
    {
      const auto box2Corners = get_box_corners(box2Config, box2);
      const REAL actualMaxMinSqrDist = compute_max_closest_distance_squared_between_points(std::vector<stk::math::Vec<REAL,DIM>>{box1Corners[0], box1Corners[1]}, std::vector<stk::math::Vec<REAL,DIM>>{box2Corners[0], box2Corners[1]});
      overallMaxMinSqrDist = std::max(overallMaxMinSqrDist, actualMaxMinSqrDist);
      EXPECT_LE(actualMaxMinSqrDist, upperBoundMinSqrDist);
    }
  }

  EXPECT_NEAR(upperBoundMinSqrDist, overallMaxMinSqrDist, 1.e-4) << "Upper bound is not tight.";
}

TEST(BoxToBoxMaxDistance,eachPossibleConfigurationOfMinimalPointsInBoundingBoxes_minDistanceQuerySatisfiesUpperBoundOfQueryBoxToSearchBox)
{
  const BoundingBox box1{stk::math::Float3d{-1,-1,-1}, stk::math::Float3d{1,1,1}};

  {
    // box1 corners: Vec3d: -1  1 -1  Vec3d:  1 -1  1
    // box2 corners: Vec3d: -2  1 -2  Vec3d: -1  2 -1
    // maxSqrDist = 2^2 + 3^2 + 2^2 = 17
    const BoundingBox box2_XmYpZm{stk::math::Float3d{-2,+1,-2}, stk::math::Float3d{-1,+2,-1}};
    test_max_possible_closest_distance_between_each_minimal_configuration_via_bounding_boxes(box1, box2_XmYpZm);
    const double goldUpperBoundSqrDist = 17;
    EXPECT_NEAR(goldUpperBoundSqrDist, max_possible_closest_squared_distance_between_contained_points(box1, box2_XmYpZm), 1.e-6);
  }

  {
    // box1 corners: Vec3d: -1  1 -1  Vec3d:  1 -1  1
    // box2 corners: Vec3d: -2  1 -1  Vec3d: -1  2  1
    // maxSqrDist = 2^2 + 3^2 + 0^2 = 13
    const BoundingBox box2_XmYpZ0{stk::math::Float3d{-2,+1,-1}, stk::math::Float3d{-1,+2,+1}};
    test_max_possible_closest_distance_between_each_minimal_configuration_via_bounding_boxes(box1, box2_XmYpZ0);
    const double goldUpperBoundSqrDist = 13;
    EXPECT_NEAR(goldUpperBoundSqrDist, max_possible_closest_squared_distance_between_contained_points(box1, box2_XmYpZ0), 1.e-6);
  }

  {
    // box1 corners: Vec3d: -1 -1  1  Vec3d: 1  1 -1
    // box2 corners: Vec3d: -1 -1 -1  Vec3d: 1  1  1
    // maxSqrDist = 0^2 + 0^2 + 2^2 = 4
    const BoundingBox box2_X0Y0Z0{stk::math::Float3d{-1,-1,-1}, stk::math::Float3d{+1,+1,+1}};
    test_max_possible_closest_distance_between_each_minimal_configuration_via_bounding_boxes(box1, box2_X0Y0Z0);
    const double goldUpperBoundSqrDist = 4;
    EXPECT_NEAR(goldUpperBoundSqrDist, max_possible_closest_squared_distance_between_contained_points(box1, box2_X0Y0Z0), 1.e-6);
  }

  {
    // box1 corners: Vec3d:  1 -1 -1  Vec3d: -1 1 1
    // box2 corners: Vec3d: -2 -1 -1  Vec3d: -1 1 1
    // maxSqrDist = 3^2 + 0^2 + 0^2 = 9
    const BoundingBox box2_XmY0Z0{stk::math::Float3d{-2,-1,-1}, stk::math::Float3d{-1,+1,+1}};
    test_max_possible_closest_distance_between_each_minimal_configuration_via_bounding_boxes(box1, box2_XmY0Z0);
    const double goldUpperBoundSqrDist = 9;
    EXPECT_NEAR(goldUpperBoundSqrDist, max_possible_closest_squared_distance_between_contained_points(box1, box2_XmY0Z0), 1.e-6);
  }
}

BoundingBox compute_search_box(const BoundingBox & fullSearchBox, const unsigned numDx, const unsigned i, const unsigned j, const unsigned k)
{
  const stk::math::Float3d dx = (fullSearchBox.get_max()-fullSearchBox.get_min())/numDx;
  const stk::math::Float3d minCorner = fullSearchBox.get_min() + stk::math::Float3d{i*dx[0], j*dx[1], k*dx[2]};
  const stk::math::Float3d maxCorner = fullSearchBox.get_min() + stk::math::Float3d{(i+1)*dx[0], (j+1)*dx[1], (k+1)*dx[2]};
  return BoundingBox{minCorner, maxCorner};
}

TEST(BoxToBoxMaxDistance, gridOfSearchBoxes_eachPossibleConfigurationOfMinimalPointsInBoundingBoxes_eachQueryBelowUpperBound)
{
  const BoundingBox box1{stk::math::Float3d{-0.5,-1,-2}, stk::math::Float3d{1,1,1}};

  const BoundingBox search{stk::math::Float3d{-1.5,-2,-3}, stk::math::Float3d{2,2,2}};
  const int numDx = 4;
  for (int i=0; i<numDx; ++i)
    for (int j=0; j<numDx; ++j)
      for (int k=0; k<numDx; ++k)
        test_max_possible_closest_distance_between_each_minimal_configuration_via_bounding_boxes(box1, compute_search_box(search, numDx, i,j,k));
}

template<class REAL, unsigned DIM>
void test_max_possible_closest_distance_between_points_via_bounding_boxes(const std::vector<stk::math::Vec<REAL,DIM>> & box1Pts, const std::vector<stk::math::Vec<REAL,DIM>> & box2Pts)
{
  const REAL actualMaxMinSqrDist = compute_max_closest_distance_squared_between_points(box1Pts, box2Pts);

  const auto box1 = compute_point_bounding_box(box1Pts);
  const auto box2 = compute_point_bounding_box(box2Pts);
  const auto upperBoundMinSqrDist = max_possible_closest_squared_distance_between_contained_points(box1, box2);

  EXPECT_LE(actualMaxMinSqrDist, upperBoundMinSqrDist);
}

TEST(BoxToBoxMaxDistance, gridOfSearchBoxes_randomPointsDistributedInBoundingBoxes_eachQueryBelowUpperBound)
{
  const size_t numPts = 1000;
  const BoundingBox box1{stk::math::Float3d{-0.5,-1,-2}, stk::math::Float3d{1,1,1}};
  const std::vector<stk::math::Float3d> box1Pts  = create_random_points_in_bounding_box(numPts, box1);

  const BoundingBox search{stk::math::Float3d{-1.5,-2,-3}, stk::math::Float3d{2,2,2}};
  const int numDx = 4;
  for (int i=0; i<numDx; ++i)
    for (int j=0; j<numDx; ++j)
      for (int k=0; k<numDx; ++k)
        test_max_possible_closest_distance_between_points_via_bounding_boxes(box1Pts, create_random_points_in_bounding_box(numPts, compute_search_box(search, numDx, i,j,k)));
}


}
