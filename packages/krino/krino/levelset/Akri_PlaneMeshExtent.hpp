// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_LEVELSET_AKRI_PLANEMESHEXTENT_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_PLANEMESHEXTENT_HPP_

#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>

#include <array>
#include <vector>

namespace stk { namespace mesh { class BulkData; } }

namespace krino {

class Surface;

// Minimum and maximum of direction.x over the nodes of the mesh, on all processors.
std::array<double,2> compute_min_max_of_nodes_along_direction(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::math::Vector3d & direction);

// Sets the offset of each plane in the surfaces (or among their subsurfaces) that is placed at a
// fraction of the mesh extent along its normal.  The extent is computed once for each distinct normal.
void set_plane_offsets_from_mesh_extent(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const std::vector<Surface*> & surfaces);

// Does the above for the bounding surfaces and level set initialization surfaces of the mesh, using
// its current coordinates.
void set_plane_offsets_from_mesh_extent(const stk::mesh::BulkData & mesh);

}

#endif /* KRINO_KRINO_LEVELSET_AKRI_PLANEMESHEXTENT_HPP_ */
