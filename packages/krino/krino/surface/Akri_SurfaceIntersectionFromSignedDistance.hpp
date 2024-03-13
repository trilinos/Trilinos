#ifndef KRINO_KRINO_SURFACE_AKRI_SURFACEINTERSECTIONFROMSIGNEDDISTANCE_HPP_
#define KRINO_KRINO_SURFACE_AKRI_SURFACEINTERSECTIONFROMSIGNEDDISTANCE_HPP_

#include <stk_math/StkVector.hpp>

namespace krino {
  class Surface;

  std::pair<int, double> compute_surface_intersection_with_crossed_segment_from_signed_distance(const Surface & surface, const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double dist0, const double dist1, const double edgeCrossingTol);
  std::pair<int, double> compute_surface_intersection_with_segment_from_signed_distance(const Surface & surface, const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol);
}

#endif /* KRINO_KRINO_SURFACE_AKRI_SURFACEINTERSECTIONFROMSIGNEDDISTANCE_HPP_ */
