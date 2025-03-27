/*
 * Akri_IntersectionUtils.hpp
 *
 *  Created on: Dec 10, 2024
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_SURFACE_AKRI_INTERSECTIONUTILS_HPP_
#define KRINO_KRINO_SURFACE_AKRI_INTERSECTIONUTILS_HPP_
#include <array>
#include <Akri_BoundingBox.hpp>

namespace krino {

bool does_bounding_box_intersect_tetrahedron(const BoundingBox & bbox, const std::array<stk::math::Vector3d,4> & tetNodes);

bool does_bounding_box_intersect_triangle_3d(const BoundingBox & bbox, const std::array<stk::math::Vector3d,3> & triNodes);

bool does_bounding_box_intersect_edge_2d(const BoundingBox & bbox, const std::array<stk::math::Vector3d,2> & edgeNodes);

bool does_bounding_box_intersect_tetrahedron(const BoundingBox & bbox, const std::array<stk::math::Vector3d,4> & tetNodes);

void append_tetrahedron_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,4> & tetNodes, const stk::math::Vector3d & candidate, std::vector<stk::math::Vector3d> & intParamCoords);

void append_triangle_edge_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,3> & triNodes,
    const double triArea, // passed in, rather than recomputed, for performance
    const stk::math::Vector3d & triNormal, // passed in, rather than recomputed, for performance
    const stk::math::Vector3d & edgePt0,
    const stk::math::Vector3d & edgePt1,
    std::vector<stk::math::Vector3d> & intParamCoords);

}



#endif /* KRINO_KRINO_SURFACE_AKRI_INTERSECTIONUTILS_HPP_ */
