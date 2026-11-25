// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_FIELD_CALCS_AKRI_MESHTOPOLOGYSHARPFEATURES_HPP_
#define KRINO_KRINO_FIELD_CALCS_AKRI_MESHTOPOLOGYSHARPFEATURES_HPP_

#include <optional>
#include <vector>

#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

// How a node on a topological feature of the mesh may move when the mesh is smoothed: within the
// plane of its surface (direction is the unit normal), along the curve where surfaces meet
// (direction is the unit tangent), or not at all (direction is zero).
struct SmoothingDirection
{
  enum Type
  {
    SMOOTH_IN_PLANE,
    SMOOTH_ALONG_TANGENT,
    FIXED
  };
  stk::mesh::Entity node;
  Type type;
  stk::math::Vector3d direction;
};

// The smoothing directions of the locally owned and shared nodes of the selected elements that are on
// a topological feature of the mesh, sorted by node id.  The features come from the parts, and, with a
// feature angle, also from where the sides of a surface meet at an angle sharper than it (see
// Akri_SharpFeatureCalcs.hpp: the interior angle, 180 degrees where flat):
// - A surface is made of the sides that are in the same set of side-rank parts (other than stk's
//   topology root parts) and have a selected element.  A node is on a feature if it has such a side,
//   or if it is in a nodeset (an IO part of node rank; other node-rank parts, such as krino's own, do
//   not count).
// - In 3D, a feature edge from the node is one that is not shared by exactly two of its sides of the
//   same surface (where surfaces meet, or at the open edge of a surface), or that is an edge in an
//   edge-rank part, or, with a feature angle, whose two sides meet at a sharp angle.  With no feature
//   edge the node moves within the plane of its surface (the unit area-weighted normal of its sides);
//   with two, along the curve (the unit tangent from its two neighbors on it); otherwise (corners,
//   ends and junctions of curves) it is fixed.
// - In 2D, the sides are segments: a node with exactly two segments of the same surface, meeting at an
//   angle that is not sharp if there is a feature angle, moves along the curve (the unit tangent from
//   its two neighbors on it); otherwise it is fixed.
// - A node in a nodeset is fixed.
// Linear sides (their vertices).  On more than one processor this needs the aura, for the sides of
// the shared nodes.
std::vector<SmoothingDirection> compute_smoothing_directions_from_mesh_topology(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const std::optional<double> cosFeatureAngle = std::nullopt);

}

#endif /* KRINO_KRINO_FIELD_CALCS_AKRI_MESHTOPOLOGYSHARPFEATURES_HPP_ */
