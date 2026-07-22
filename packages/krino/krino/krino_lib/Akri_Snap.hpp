// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_SNAP_H_
#define KRINO_INCLUDE_AKRI_SNAP_H_
#include <Akri_FieldRef.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_math/StkVector.hpp>

namespace krino
{
class InterfaceGeometry;
class SharpFeatureInfo;
class IntersectionPoint;
class QualityMetric;

typedef std::map<stk::mesh::Entity, std::vector<std::pair<size_t,bool>>> mapFromEntityToIntPtIndexAndSnapAllowed;

NodeToCapturedDomainsMap snap_as_much_as_possible_while_maintaining_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & potentialParentElementSelector,
    const stk::mesh::Selector & decomposedParentElementSelector,
    const FieldRef coordsField,
    const FieldSet & interpolationFields,
    const InterfaceGeometry & geometry,
    const bool globalIDsAreParallelConsistent,
    const double snappingSharpFeatureAngleInDegrees,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxSnapForEdges);

void undo_previous_snaps_using_interpolation(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField, const FieldSet & snapFields);
void snap_fields_using_interpolation(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField, const FieldSet & interpFields);

stk::math::Vector3d compute_intersection_point_location(const int dim, const FieldRef coordsField, const IntersectionPoint & intersectionPoint);

mapFromEntityToIntPtIndexAndSnapAllowed get_node_to_intersection_point_indices_and_which_snaps_allowed(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const double maxSnapForEdges,
    const std::vector<IntersectionPoint> & intersectionPoints);

std::map<std::vector<int>, std::map<stk::mesh::EntityId,double>> determine_quality_per_node_per_domain(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed,
    const QualityMetric &qualityMetric,
    const double minIntPtWeightForEstimatingCutQuality,
    const bool globalIDsAreParallelConsistent);

}



#endif /* KRINO_INCLUDE_AKRI_SNAP_H_ */
