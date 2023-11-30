#include <Akri_SurfaceIntersectionFromSignedDistance.hpp>

#include <Akri_MathUtil.hpp>
#include <Akri_Sign.hpp>
#include <Akri_Surface.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <functional>

namespace krino {

// Note:  Seems pretty strange, but the boost::toms748_solve in find_root appears to do better on the interval -1->1 than from 0->1

static std::function<double(const double)> build_edge_distance_function(const Surface & surface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords)
{
  std::function<double(const double)> distanceFunction =
    [&surface, &edgeNodeCoords](const double x)
    {
      return surface.point_signed_distance(0.5*(1.-x)*edgeNodeCoords[0] + 0.5*(1.+x)*edgeNodeCoords[1]);
    };
  return distanceFunction;
}

static double find_crossing_position(const Surface & surface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords, const double edgeTol)
{
  const double phi0 = surface.point_signed_distance(edgeNodeCoords[0]);
  const double phi1 = surface.point_signed_distance(edgeNodeCoords[1]);
  const int maxIters = 100;
  const auto result = find_root(build_edge_distance_function(surface, edgeNodeCoords), -1., 1., phi0, phi1, maxIters, 2.*edgeTol);
  STK_ThrowRequire(result.first);

  return 0.5*(1.+result.second);
}

std::pair<int, double> compute_surface_intersection_with_segment_from_signed_distance(const Surface & surface, const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol)
{
  const double phi0 = surface.point_signed_distance(pt0);
  const double phi1 = surface.point_signed_distance(pt1);
  if (sign_change(phi0, phi1))
  {
    const std::array<stk::math::Vector3d,2> edgeNodeCoords{pt0, pt1};
    const double location = find_crossing_position(surface, edgeNodeCoords, edgeCrossingTol);
    return {sign(phi1), location};
  }
  return {0, -1.};
}

}

