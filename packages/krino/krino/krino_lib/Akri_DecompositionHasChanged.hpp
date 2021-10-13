// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_DECOMPOSITIONHASCHANGED_H_
#define KRINO_INCLUDE_AKRI_DECOMPOSITIONHASCHANGED_H_
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "../interface_geometry_interface/Akri_InterfaceGeometry.hpp"
namespace krino {

bool decomposition_has_changed(const stk::mesh::BulkData & mesh,
    const InterfaceGeometry & interfaceGeometry,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport);

}

#endif /* KRINO_INCLUDE_AKRI_DECOMPOSITIONHASCHANGED_H_ */
