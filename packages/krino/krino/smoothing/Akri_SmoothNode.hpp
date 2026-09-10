// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SMOOTHNODE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SMOOTHNODE_HPP_
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_OptimizeTypes.hpp>

namespace krino {

void smooth_mesh_by_ODT_of_interior_node_locations(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Selector & interiorNodeSelector,
    const size_t maxNumSmoothIterations = 50);

void smooth_mesh_by_optimizing_interior_node_locations(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const stk::mesh::Selector & interiorNodeSelector,
    const Optimizer3DFunction & optimizer,
    const NodeCoordinateObjectiveInterface & smoothingObj,
    const size_t maxNumSmoothIterations = 50);

void smooth_mesh_by_optimizing_node_locations(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const Optimizer3DFunction & optimizer,
    const NodeCoordinateObjectiveInterface & smoothingObj,
    const GradientFilter & node_search_direction_filter,
    const size_t maxNumSmoothIterations = 50);
}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SMOOTHNODE_HPP_ */
