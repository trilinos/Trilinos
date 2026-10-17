#ifndef KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVECONTOURTET_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVECONTOURTET_HPP_

#include <stk_math/StkVector.hpp>

#include <array>
#include <functional>

namespace krino {

class FacetedSurfaceBase;

void adaptively_append_facets_for_tet_using_semilagrangian_distance(const std::array<stk::math::Vector3d,4> & coords,
  const std::array<stk::math::Vector3d,4> & departureCoords,
  const std::array<double,4> & distance,
  const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
  const double distTol,
  FacetedSurfaceBase & facets,
  const int currentDepth,
  const int minDepth,
  const int maxDepth);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVECONTOURTET_HPP_ */
