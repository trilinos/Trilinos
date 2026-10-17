#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBER_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBER_HPP_

#include <vector>
#include <stk_math/StkVector.hpp>

namespace krino {

double compute_facet_winding_number(const stk::math::Vector3d & facetX0, const stk::math::Vector3d & facetX1, const stk::math::Vector3d & facetX2, const stk::math::Vector3d & queryLoc);
double compute_faceted_surface_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc);

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBER_HPP_ */
