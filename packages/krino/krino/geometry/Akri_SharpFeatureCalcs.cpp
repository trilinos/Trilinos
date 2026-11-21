// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SharpFeatureCalcs.hpp>

namespace krino
{

double cosine_of_angle_between_triangles_at_edge(const stk::math::Vector3d & edgePt0,
    const stk::math::Vector3d & edgePt1,
    const stk::math::Vector3d & triVertexOffEdge0,
    const stk::math::Vector3d & triVertexOffEdge1)
{
  // The dihedral angle (https://en.wikipedia.org/wiki/Dihedral_angle), from the tangents of the
  // triangles normal to the edge.
  const stk::math::Vector3d edgeVec = edgePt1 - edgePt0;
  const stk::math::Vector3d crossEdgeTri0 = Cross(edgeVec, triVertexOffEdge0 - edgePt0);
  const stk::math::Vector3d crossEdgeTri1 = Cross(edgeVec, triVertexOffEdge1 - edgePt0);
  return Dot(crossEdgeTri0, crossEdgeTri1) / (crossEdgeTri0.length() * crossEdgeTri1.length());
}

double cosine_of_angle_between_segments_at_vertex(const stk::math::Vector3d & vertex,
    const stk::math::Vector3d & otherSegmentVertex0,
    const stk::math::Vector3d & otherSegmentVertex1)
{
  return Dot(
      (otherSegmentVertex0 - vertex).unit_vector(), (otherSegmentVertex1 - vertex).unit_vector());
}

} // namespace krino
