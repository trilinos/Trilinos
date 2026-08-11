// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SMOOTHMESH_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SMOOTHMESH_HPP_
#include <stk_mesh/base/Types.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_OptimizeTypes.hpp>
#include <Akri_Objective_Mesh.hpp>

namespace krino {

void smooth_mesh_by_optimization(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const OptimizerFunction & optimizer,
    const MeshObjectiveBase & smoothingObj,
    const bool doPrintDiagnostics = false);

void smooth_mesh_by_optimization(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizerFunction & optimizer,
    const ElementObjective & elemSmoothingObj,
    const GradientFilter & gradientFilter,
    const SolutionProjector & solutionProjector,
    const bool doPrintDiagnostics = false);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SMOOTHMESH_HPP_ */
