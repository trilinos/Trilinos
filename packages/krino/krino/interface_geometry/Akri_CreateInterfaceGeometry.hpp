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
// TEMPORARY COMPATIBILITY INCLUDE -- REMOVE ME.
// This header only needs a forward declaration of LS_Field for its own
// signatures, but external callers (platoengine third_party_integration/krino)
// include it specifically "//LS_Field" and rely on the COMPLETE ::krino::LS_Field
// type -- some of them construct/copy std::vector<LS_Field> (e.g. the
// KrinoTestFixture / Test_Utilities helpers), for which a forward declaration is
// not enough.  Both this header and Akri_LS_Field.hpp live in
// krino_interface_geometry, so pulling it in here is layering-safe.  The
// permanent fix is for those callers to include <Akri_LS_Field.hpp> directly,
// then restore the forward declaration below.
#include <Akri_LS_Field.hpp>

namespace krino {

class CDFEM_Support;
class Phase_Support;
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
