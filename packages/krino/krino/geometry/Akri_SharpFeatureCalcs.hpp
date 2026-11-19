// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_SHARPFEATURECALCS_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_SHARPFEATURECALCS_HPP_

#include <stk_math/StkVector.hpp>

namespace krino
{

// The geometry of sharp features, for surfaces given by any kind of mesh.  The angle at a feature
// is the interior angle between the two faces (3D) or segments (2D) that meet there: 180 degrees
// where they are flat (coplanar or collinear), 0 where they fold back onto each other.  A feature
// is sharp when its angle is less than the feature angle, that is when its cosine is greater than
// the cosine of the feature angle.

// The cosine of the angle between two triangles sharing the edge from edgePt0 to edgePt1, given by
// their vertices off the edge.
double cosine_of_angle_between_triangles_at_edge(const stk::math::Vector3d & edgePt0,
    const stk::math::Vector3d & edgePt1,
    const stk::math::Vector3d & triVertexOffEdge0,
    const stk::math::Vector3d & triVertexOffEdge1);

// The cosine of the angle between two segments sharing the vertex, given by their other vertices.
double cosine_of_angle_between_segments_at_vertex(const stk::math::Vector3d & vertex,
    const stk::math::Vector3d & otherSegmentVertex0,
    const stk::math::Vector3d & otherSegmentVertex1);

inline bool is_angle_sharper_than_feature_angle(const double cosAngle, const double cosFeatureAngle)
{
  return cosAngle > cosFeatureAngle;
}

} // namespace krino

#endif /* KRINO_KRINO_GEOMETRY_AKRI_SHARPFEATURECALCS_HPP_ */
