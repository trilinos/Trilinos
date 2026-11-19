// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_SNAP_H_
#define KRINO_INCLUDE_AKRI_SNAP_H_
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CutQualityEstimator.hpp>
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

NodeToCapturedDomainsMap
snap_as_much_as_possible_while_maintaining_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & potentialParentElementSelector,
    const stk::mesh::Selector & decomposedParentElementSelector,
    const FieldRef coordsField,
    const FieldSet & interpolationFields,
    const InterfaceGeometry & geometry,
    const bool globalIDsAreParallelConsistent,
    const double snappingSharpFeatureAngleInDegrees,
    const double minIntPtWeightForEstimatingCutQuality,
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const Snap_Selection snapSelection);

void undo_previous_snaps_using_interpolation(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField, const FieldSet & snapFields);
void snap_fields_using_interpolation(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const FieldRef coordsField, FieldRef cdfemSnapField, const FieldSet & interpFields);

mapFromEntityToIntPtIndexAndSnapAllowed get_node_to_intersection_point_indices_and_which_snaps_allowed(const stk::mesh::BulkData & mesh,
    const SharpFeatureInfo * sharpFeatureInfo,
    const double maxEdgeSnap,
    const double minWeightAtNode,
    const std::vector<IntersectionPoint> & intersectionPoints);

}



#endif /* KRINO_INCLUDE_AKRI_SNAP_H_ */
