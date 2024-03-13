#include <Akri_BoundingBoxDistance.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_BoundingBox.hpp>

namespace krino {

template<class REAL, unsigned DIM>
REAL min_possible_closest_squared_distance(const BoundingBox_T<REAL,DIM> & bbox, const stk::math::Vec<REAL,3> & queryPt)
{
  STK_ThrowAssert( bbox.valid() );

  const stk::math::Vec<REAL,3> & bboxMin = bbox.get_min();
  const stk::math::Vec<REAL,3> & bboxMax = bbox.get_max();

  REAL minSqrDist = 0.;

  for ( unsigned i = 0; i < DIM; i++ )
  {
    if ( queryPt[i] < bboxMin[i] )
    {
      const REAL delta = bboxMin[i] - queryPt[i];
      minSqrDist += delta * delta;
    }
    else if ( queryPt[i] > bboxMax[i] )
    {
      const REAL delta = queryPt[i] - bboxMax[i];
      minSqrDist += delta * delta;
    }
  }
  return ( minSqrDist );
}

template<class REAL, unsigned DIM>
std::pair<stk::math::Vec<REAL,3>, stk::math::Vec<REAL,3>> get_close_and_far_corners_of_bounding_box(const BoundingBox_T<REAL,DIM> & bbox, const stk::math::Vec<REAL,3> & queryPt)
{
  STK_ThrowAssert( bbox.valid() );
  const stk::math::Vec<REAL,3> & bboxMin = bbox.get_min();
  const stk::math::Vec<REAL,3> & bboxMax = bbox.get_max();

  stk::math::Vec<REAL,3> closePt;
  stk::math::Vec<REAL,3> farPt;

  for ( unsigned i = 0; i < DIM; i++ )
  {
    if ( queryPt[i] < bboxMin[i] )
    {
      closePt[i] = bboxMin[i];
      farPt[i] = bboxMax[i];
    }
    else if ( queryPt[i] > bboxMax[i] )
    {
      closePt[i] = bboxMax[i];
      farPt[i] = bboxMin[i];
    }
    else
    {
      if (queryPt[i]-bboxMin[i] < bboxMax[i]-queryPt[i])
      {
        closePt[i] = bboxMin[i];
        farPt[i] = bboxMax[i];
      }
      else
      {
        closePt[i] = bboxMax[i];
        farPt[i] = bboxMin[i];
      }
    }
  }

  return {closePt, farPt};
}

template<class REAL, unsigned DIM>
REAL max_possible_closest_squared_distance(const BoundingBox_T<REAL,DIM> & bbox, const stk::math::Vec<REAL,3> & queryPt)
{
  typedef stk::math::Vec<REAL,3> VecType;

  // We are guaranteed that there is a point on the surface on each face of the
  // bounding box.  So we know that the upper bound for the distance to a face is
  // the distance to the farthest point on that face.  So the upper bound for this
  // bounding box is the minimum of the upper bounds for each face.  In other words,
  // the upper bound is the minimum distance to the farthest point on each face.

  const auto [closePt, farPt] = get_close_and_far_corners_of_bounding_box(bbox, queryPt);

  REAL minSqrDist;
  if (3 == DIM)
  {
    minSqrDist = (queryPt-VecType(closePt[0],farPt[1],farPt[2])).length_squared();
    minSqrDist = std::min(minSqrDist,(queryPt-VecType(farPt[0],closePt[1],farPt[2])).length_squared());
    minSqrDist = std::min(minSqrDist,(queryPt-VecType(farPt[0],farPt[1],closePt[2])).length_squared());
  }
  else
  {
    STK_ThrowAssert(2 == DIM);
    minSqrDist = (queryPt-VecType(closePt[0],farPt[1],0.0)).length_squared();
    minSqrDist = std::min(minSqrDist,(queryPt-VecType(farPt[0],closePt[1],0.0)).length_squared());
  }

  return minSqrDist;
}

template<class REAL>
std::array<stk::math::Vec<REAL,2>,2> get_box_corners(const unsigned config, const BoundingBox_T<REAL,2> & box)
{
  const auto & boxMin = box.get_min();
  const auto & boxMax = box.get_max();
  if (0 == config)
    return {stk::math::Vec<REAL,2>{boxMin[0],boxMin[1]}, stk::math::Vec<REAL,2>{boxMax[0],boxMax[1]}};
  return {stk::math::Vec<REAL,2>{boxMin[0],boxMax[1]}, stk::math::Vec<REAL,2>{boxMax[0],boxMin[1]}};
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
REAL max_possible_closest_squared_distance_between_contained_points(const BoundingBox_T<REAL,DIM> & box1, const BoundingBox_T<REAL,DIM> & box2)
{
  if( !box1.valid() || !box2.valid() )
    return std::numeric_limits<REAL>::max();

  constexpr unsigned numConfig = 1 << (DIM-1);
  REAL maxMinSqrDist = 0.;
  for (unsigned box1Config=0; box1Config<numConfig; ++box1Config)
  {
    const auto box1Corners = get_box_corners(box1Config, box1);
    for (unsigned box2Config=0; box2Config<numConfig; ++box2Config)
    {
      const auto box2Corners = get_box_corners(box2Config, box2);
      maxMinSqrDist = std::max(maxMinSqrDist, std::max(
          std::min((box1Corners[0]-box2Corners[0]).length_squared(),(box1Corners[0]-box2Corners[1]).length_squared()),
          std::min((box1Corners[1]-box2Corners[0]).length_squared(),(box1Corners[1]-box2Corners[1]).length_squared())));
    }
  }
  return maxMinSqrDist;
}

// Explicit template instantiation
template float min_possible_closest_squared_distance(const BoundingBox_T<float,3> & bbox, const stk::math::Vec<float,3> & queryPt);
template float min_possible_closest_squared_distance(const BoundingBox_T<float,2> & bbox, const stk::math::Vec<float,3> & queryPt);
template float max_possible_closest_squared_distance(const BoundingBox_T<float,3> & bbox, const stk::math::Vec<float,3> & queryPt);
template float max_possible_closest_squared_distance(const BoundingBox_T<float,2> & bbox, const stk::math::Vec<float,3> & queryPt);
template float max_possible_closest_squared_distance_between_contained_points(const BoundingBox_T<float,3> & bbox1, const BoundingBox_T<float,3> & bbox2);
template float max_possible_closest_squared_distance_between_contained_points(const BoundingBox_T<float,2> & bbox1, const BoundingBox_T<float,2> & bbox2);

}

