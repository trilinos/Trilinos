// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_FIELD_CALCS_AKRI_NODALMOTIONCONSTRAINT_HPP_
#define KRINO_KRINO_FIELD_CALCS_AKRI_NODALMOTIONCONSTRAINT_HPP_
#include <functional>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino
{

class FieldRef;

// How a node may move: within a surface (direction is its unit normal), along a line (direction is
// the line's unit tangent), or freely.
struct NodalMotionConstraint
{
  enum Type
  {
    FREE,
    SURFACE,
    LINE
  };
  Type type{FREE};
  stk::math::Vector3d direction{stk::math::Vector3d::ZERO};
};

// Where a node with zero level set on the boundary moves as the level set increases, so that the
// interface stays on the boundary.  The node's boundary sides (in boundarySideSelector, with an
// element in elementSelector) are grouped into walls: two sides sharing an edge from the node (in
// 2D, any two sides) are in the same wall if areSidesInSameWall says so.  The interface touches the
// boundary at the node along contact segments: toward the cut on a side's opposite edge (3D), or
// along a side edge whose other node also has zero level set, which is in each wall containing that
// edge, or in those of them with a negative neighbor of the node if any.
//   SURFACE  one wall contains every contact: the node moves within it (its unit area-weighted
//            normal, out of the elements)
//   LINE     otherwise: the node moves along the line where two of the contacts' walls meet
//   FREE     the node has no boundary side, or the interface does not touch the boundary there
// In 2D the contact is the node itself: SURFACE in the one wall with a negative neighbor, or the
// only wall, else FREE.  For a planar interface on flat walls this is the motion of the interface's
// contact with the boundary.  Linear line or tri sides of tri or tet elements.  Local: needs the
// node's sides and their nodes (on a shared node, the aura).
NodalMotionConstraint
compute_zero_levelset_boundary_node_constraint(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef levelSetField,
    const stk::mesh::Entity node,
    const stk::mesh::Selector & boundarySideSelector,
    const stk::mesh::Selector & elementSelector,
    const std::function<bool(stk::mesh::Entity side0, stk::mesh::Entity side1)> &
        areSidesInSameWall);

// The part of v that satisfies the constraint (removing the SURFACE normal component, or keeping
// the LINE component), scaled so that its dot product with v is v.v.  For v along an interface
// normal this keeps the interface's normal speed.  For FREE, or if the constrained part is nearly
// orthogonal to v, returns v.
stk::math::Vector3d constrain_keeping_normal_speed(
    const NodalMotionConstraint & constraint, const stk::math::Vector3d & v);

} // namespace krino

#endif /* KRINO_KRINO_FIELD_CALCS_AKRI_NODALMOTIONCONSTRAINT_HPP_ */
