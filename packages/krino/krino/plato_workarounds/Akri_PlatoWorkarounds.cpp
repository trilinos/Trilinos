// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_PlatoWorkarounds.hpp>

#include <Akri_CSG_Tree.hpp>
#include <Akri_LevelSet_Phase_Utils.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Surface.hpp>

namespace krino {

// ============================================================================
// TEMPORARY COMPATIBILITY SHIM -- DO NOT BUILD ON THIS.  REMOVE ME.
// ----------------------------------------------------------------------------
// get_levelset_fields() used to be a static member of Phase_Support.  During
// the krino library modularization it was extracted into a FREE function,
//     krino::get_levelset_fields(const stk::mesh::MetaData &)
// (declared in Akri_LevelSet_Phase_Utils.hpp, in krino_interface_geometry).
//
// External callers (platoengine's third_party_integration/krino, e.g.
// Utilities.cpp:117,126) still spell it the old way,
// "::krino::Phase_Support::get_levelset_fields(meta)", so their build breaks:
//     no member named 'get_levelset_fields' in 'krino::Phase_Support'
//
// The member is re-declared in Akri_Phase_Support.hpp (a static member must be
// declared inside its class, which lives in krino_phase, so the declaration
// cannot move out of that header).  Its DEFINITION lives here in
// krino_plato_workarounds -- the home for platoengine compatibility shims --
// so the whole workaround can be deleted by removing this one library's
// contribution once platoengine migrates.  This requires krino_plato_workarounds
// to depend on krino_interface_geometry (for the complete LS_Field type and the
// free function); there is no cycle since interface_geometry does not depend on
// plato_workarounds.
//
// It simply forwards to the free function so callers get identical behavior to
// the real (pre-modularization) API.  The correct permanent fix is for callers
// to use ::krino::get_levelset_fields(meta) directly, then delete this shim.
// ============================================================================
std::vector<LS_Field> Phase_Support::get_levelset_fields(const stk::mesh::MetaData & meta)
{
  return krino::get_levelset_fields(meta);
}

void
Composite_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  // Unite the accumulated primitives into the CSG_Surface base, then defer to
  // its prepare_to_compute for the narrow-band filtering and wrong-sign
  // reduction.  CSGTree::unite requires a non-empty operand list, so leave the
  // tree empty (CSG_Surface then returns the far-field value) when nothing was
  // added, mirroring the old empty Composite_Surface.
  if (!mySubSurfaces.empty())
    set_tree(CSGTree::unite(mySubSurfaces));

  CSG_Surface::prepare_to_compute(time, point_bbox, truncation_length);
}

} // namespace krino
