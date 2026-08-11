// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_NODEMARKER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_NODEMARKER_HPP_
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

class IntersectionPoint;

void mark_nodes_near_edge_intersections(const stk::mesh::BulkData& mesh,
    const std::vector<IntersectionPoint> & edgeIntersections,
    FieldRef node_marker_field);

bool element_has_marked_node(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem, const FieldRef nodeMarkerField);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_NODEMARKER_HPP_ */
