// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_CDMESH_DEBUG_H_
#define KRINO_INCLUDE_AKRI_CDMESH_DEBUG_H_
#include <Akri_Element.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace krino {
  void debug_elem_parts_and_relations(const stk::mesh::BulkData & mesh, const Mesh_Element & elem);
  void debug_nodal_parts_and_fields(const stk::mesh::BulkData & mesh, const SubElementNode * node);
  void debug_sides(const stk::mesh::BulkData & mesh, stk::mesh::Part & active_part);
}

#endif /* KRINO_INCLUDE_AKRI_CDMESH_DEBUG_H_ */
