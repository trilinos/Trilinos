// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_ELEMENTCUTTERUTILS_H_
#define KRINO_INCLUDE_AKRI_ELEMENTCUTTERUTILS_H_
#include <Akri_Intersection_Points.hpp>
#include <Akri_MasterElement.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <vector>

#include "../interface_geometry_interface/Akri_InterfaceGeometry.hpp"

namespace krino {

void append_intersection_points_from_element_interior(const MasterElement & masterElement,
    const std::vector<stk::mesh::Entity> & nodes,
    const ElementCutter & elementCutter,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints);

void append_intersection_points_from_within_elements_and_owned_faces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    const std::vector<stk::mesh::Entity> & elements,
    const InterfaceGeometry & geometry,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints);

}

#endif /* KRINO_INCLUDE_AKRI_ELEMENTCUTTERUTILS_H_ */
