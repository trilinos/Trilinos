#ifndef KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVEELEMENTCONTOUR_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVEELEMENTCONTOUR_HPP_

#include <stk_math/StkVector.hpp>

#include <array>
#include <functional>

namespace krino {

class FacetedSurfaceBase;

void adaptively_append_facets_for_tri_using_interpolated_distance(const std::array<stk::math::Vector3d,6> & tri6Coords,
  const std::array<double,6> & tri6Dist,
  const double lengthScale,
  const int currentDepth,
  const int maxDepth,
  FacetedSurfaceBase & facets);

void adaptively_append_facets_for_tri_using_semilagrangian_distance(const std::array<stk::math::Vector3d,3> & coords,
  const std::array<double,3> & distance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double distTol,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVEELEMENTCONTOUR_HPP_ */
