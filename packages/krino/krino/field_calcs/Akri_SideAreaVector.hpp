// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_FIELD_CALCS_AKRI_SIDEAREAVECTOR_HPP_
#define KRINO_KRINO_FIELD_CALCS_AKRI_SIDEAREAVECTOR_HPP_
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino
{

class FieldRef;

// The length (2D) or area (3D) times the unit normal of a linear line or tri side, pointing out of
// the given tri or tet element of the side.
stk::math::Vector3d
compute_line_or_tri_side_area_vector_out_of_element(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Entity side,
    const stk::mesh::Entity element);

} // namespace krino

#endif /* KRINO_KRINO_FIELD_CALCS_AKRI_SIDEAREAVECTOR_HPP_ */
