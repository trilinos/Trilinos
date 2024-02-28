#include <Akri_SurfaceIntersectionFromSignedDistance.hpp>

#include <Akri_MathUtil.hpp>
#include <Akri_Sign.hpp>
#include <Akri_Surface.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <functional>

namespace krino {

// Note:  Seems pretty strange, but the boost::toms748_solve in find_root appears to do better on the interval -1->1 than from 0->1

static std::function<double(const double)> build_edge_distance_function(const Surface & surface, const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1)
{
  std::function<double(const double)> distanceFunction =
    [&surface, &pt0, &pt1](const double x)
    {
      return surface.point_signed_distance(0.5*(1.-x)*pt0 + 0.5*(1.+x)*pt1);
    };
  return distanceFunction;
}

static double find_edge_crossing_position(const std::function<double(const double)> & edgeDistanceFunc, const double dist0, const double dist1, const double edgeTol)
{
  const int maxIters = 100;
  const auto result = find_root(edgeDistanceFunc, -1., 1., dist0, dist1, maxIters, edgeTol);
  STK_ThrowRequire(result.first);

  return result.second;
}

std::pair<int, double> compute_surface_intersection_with_crossed_segment_from_signed_distance(const Surface & surface, const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double dist0, const double dist1, const double edgeCrossingTol)
{
  STK_ThrowAssert(sign_change(dist0, dist1));
  const double location = find_edge_crossing_position(build_edge_distance_function(surface, pt0, pt1), dist0, dist1, 2*edgeCrossingTol); // 2x edgeTol because of 2x range (-1->1)
  return {sign(dist1), 0.5*(1.+location)};
}

std::pair<int, double> compute_surface_intersection_with_segment_from_signed_distance(const Surface & surface, const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol)
{
  const double phi0 = surface.point_signed_distance(pt0);
  const double phi1 = surface.point_signed_distance(pt1);
  if (sign_change(phi0, phi1))
    return compute_surface_intersection_with_crossed_segment_from_signed_distance(surface, pt0, pt1, phi0, phi1, edgeCrossingTol);
  return {0, -1.};
}

}

