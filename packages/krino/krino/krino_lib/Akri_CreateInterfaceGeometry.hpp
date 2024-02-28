// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_CREATEINTERFACEGEOMETRY_H_
#define KRINO_INCLUDE_AKRI_CREATEINTERFACEGEOMETRY_H_

#include <stk_mesh/base/Part.hpp>
#include <Akri_InterfaceGeometry.hpp>

namespace krino {

class CDFEM_Support;
class Phase_Support;
class LS_Field;
class Surface_Manager;

std::unique_ptr<InterfaceGeometry> create_interface_geometry(const stk::mesh::MetaData & meta);

std::unique_ptr<InterfaceGeometry> create_bounding_surface_geometry(Surface_Manager & manager,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport);

std::unique_ptr<InterfaceGeometry> create_levelset_geometry(const int dim,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields);

}

#endif /* KRINO_INCLUDE_AKRI_CREATEINTERFACEGEOMETRY_H_ */
