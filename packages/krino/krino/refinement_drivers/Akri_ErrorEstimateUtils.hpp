// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_ERRORESTIMATEUTILS_HPP_
#define KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_ERRORESTIMATEUTILS_HPP_

#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>

#include <vector>

namespace krino {

// Element error estimates from the normals of an interface in elements, indexed by element offset (zero where an
// element has no interface normal).  Each adds to errorByElemOffset, so calling it for each of several surfaces sums
// their contributions.

// For each owned selected element with an interface normal, adds the mean angle between its normal and the normals of
// its selected side neighbors that have one.
void accumulate_facet_angle_errors_from_interface_normals(const stk::mesh::BulkData& mesh,
    const stk::mesh::Selector & elemSelector,
    const std::vector<stk::math::Vector3d> & interfaceNormals,
    std::vector<double> & errorByElemOffset);

// For each owned selected element with an interface normal, adds the chordal error estimated from the jump between its
// normal and the normals of its selected side neighbors that have one, scaled by the element size.
void accumulate_facet_chordal_errors_from_interface_normals(const stk::mesh::BulkData& mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const std::vector<stk::math::Vector3d> & interfaceNormals,
    std::vector<double> & errorByElemOffset);

} // namespace krino

#endif /* KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_ERRORESTIMATEUTILS_HPP_ */
