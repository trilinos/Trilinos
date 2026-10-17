// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_METRICEDGEMARKING_HPP_
#define KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_METRICEDGEMARKING_HPP_

#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/Types.hpp>

#include <cstddef>

namespace krino {

class EdgeBasedEdgeMarker;
class MeshSizeMetric;
class Refinement;

// Refinement of edges that are too long for a target mesh size (isotropic or anisotropic, see MeshSizeMetric).

// Marks in edgeMarker every edge of the locally owned non-parent elements whose relative length in sizeMetric
// exceeds maxRelativeEdgeLength (ngs refines edges longer than sqrt(2)).  Returns the number of edges marked on
// this processor.
size_t mark_edges_too_long(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const FieldRef coordsField,
    const MeshSizeMetric & sizeMetric,
    const double maxRelativeEdgeLength,
    EdgeBasedEdgeMarker & edgeMarker);

// Sets the size of sizeMetric at each node created by the last refinement from its parent nodes, and makes it
// consistent on ghost nodes.
void update_size_for_new_nodes(const stk::mesh::BulkData & mesh,
    const Refinement & refinement,
    const MeshSizeMetric & sizeMetric);

// Refines the edges that are too long in sizeMetric (see mark_edges_too_long) and updates the size at the new nodes,
// repeating until no edge is too long or maxNumPasses passes have been done.  Returns the global number of edges
// still too long (0 if converged).
size_t refine_edges_too_long(const stk::mesh::BulkData & mesh,
    Refinement & refinement,
    const FieldRef coordsField,
    const MeshSizeMetric & sizeMetric,
    const double maxRelativeEdgeLength,
    const unsigned maxNumPasses);

} // namespace krino

#endif /* KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_METRICEDGEMARKING_HPP_ */
