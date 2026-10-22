// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_KINEMATICMESHUTILS_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_KINEMATICMESHUTILS_HPP_

#include <Akri_FieldRef.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <stk_mesh/base/Types.hpp>

#include <vector>

namespace krino {

enum class SmoothReference;

// Per-element data for the kinematic mesh objectives, computed from the mesh for the given elements,
// in order (e.g. the owned elements of a mesh objective, from get_owned_elements()).

// Ideal size of each of the given tet4 elements from its node coordinates, in the form used by
// KinematicIsotropicMeshObjective::set_element_ref_sizes.
std::vector<double> compute_element_ref_sizes_from_reference(
    const stk::mesh::BulkData & mesh,
    const FieldRef referenceCoordsField,
    const std::vector<stk::mesh::Entity> & elems,
    const SmoothReference & reference);

// Metric mapping F_M of each of the given tri3 (2D) or tet4 (3D) elements from a nodal metric field,
// in the form used by KinematicAnisotropicMeshObjective::set_element_metric_map_data.  The field holds
// the metric M = F_M^T F_M at each node as a full spatialDim*spatialDim row-major matrix, as used for
// the metric tensor size fields of ngs (e.g. improve_mesh).  The metric of each element is the
// log-Euclidean mean of the metrics at its nodes.
// Each nodal metric must be symmetric and positive definite (an error names the node otherwise).  A
// field read as a full tensor by IOSS (components xx, yy, zz, xy, ...) differs from row-major only by
// swapping each off-diagonal pair, so it gives the same symmetric metric.
std::vector<double> compute_element_metric_mappings_from_nodal_metric_field(
    const stk::mesh::BulkData & mesh,
    const FieldRef nodalMetricField,
    const std::vector<stk::mesh::Entity> & elems);

// As above, but with the nodal metrics given directly: nodalMetrics holds the metric of each of the
// nodes (spatialDim*spatialDim row-major values per node, in the order of nodes), and elemNodeIndices
// gives the indices of each element's nodes in nodes (as for a mesh objective, from get_nodes() and its
// element node indices).
std::vector<double> compute_element_metric_mappings_from_nodal_metrics(
    const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    const MeshElemNodeIndices & elemNodeIndices,
    const std::vector<double> & nodalMetrics);

} // namespace krino

#endif /* KRINO_KRINO_SMOOTHING_AKRI_KINEMATICMESHUTILS_HPP_ */
