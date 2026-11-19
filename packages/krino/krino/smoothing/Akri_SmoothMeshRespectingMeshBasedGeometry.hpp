// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_SMOOTHMESHRESPECTINGMESHBASEDGEOMETRY_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_SMOOTHMESHRESPECTINGMESHBASEDGEOMETRY_HPP_

#include <optional>

#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

// The nodes on any boundary that the mesh describes: those of every part whose primary rank is face,
// edge or node (sidesets, nodesets, boundary and interface parts), other than stk's topology root parts,
// whatever their IO role.
stk::mesh::Selector select_nodes_on_any_boundary(const stk::mesh::MetaData & meta);

// Smooths the selected elements by optimizing the mean-ratio norm of all of their node locations at
// once (krino's L-BFGS, with the default tolerances), except the nodes on any boundary
// (select_nodes_on_any_boundary), which do not move.  The objective is a norm over the elements, so the
// worst element can get worse: if the worst mean ratio drops, the node locations are restored.
// Returns whether the mesh changed.
bool smooth_mesh_keeping_boundary_nodes_fixed(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector);

// The same, except that the nodes on the topological features of the mesh, and, with a feature angle,
// also on its sharp features (compute_smoothing_directions_from_mesh_topology), move within the plane
// of their surface, or along the curve where surfaces meet, or not at all.  The search directions are
// projected onto those planes and lines, but the moved nodes are not projected back onto the
// surfaces: on curved surfaces or curves they move off them by up to the curvature times the move
// squared.  On more than one processor this needs the aura.
bool smooth_mesh_respecting_mesh_based_geometry(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const std::optional<double> cosFeatureAngle = std::nullopt);

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_SMOOTHMESHRESPECTINGMESHBASEDGEOMETRY_HPP_ */
