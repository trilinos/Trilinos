// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_LevelSetInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_BoundingSurface.hpp>
#include <Akri_Surface_Identifier.hpp>
#include <Akri_Surface_Manager.hpp>

namespace krino {

std::unique_ptr<InterfaceGeometry> create_interface_geometry(const stk::mesh::MetaData & meta)
{
  Surface_Manager & manager = Surface_Manager::get(meta);
  const auto & boundingSurfaces = manager.get_bounding_surfaces();
  if (boundingSurfaces.size() > 0)
  {
    ThrowRequire(1 == boundingSurfaces.size());
    const Surface_Identifier surfaceIdentifier = Surface_Manager::get(meta).get_identifier(boundingSurfaces[0]->name());
    return std::make_unique<AnalyticSurfaceInterfaceGeometry>(surfaceIdentifier, boundingSurfaces[0]->surface(), AuxMetaData::get(meta).active_part(), CDFEM_Support::get(meta), Phase_Support::get(meta));
  }
  return create_levelset_geometry(AuxMetaData::get(meta).active_part(), CDFEM_Support::get(meta), Phase_Support::get(meta), Phase_Support::get_levelset_fields(meta));
}

std::unique_ptr<InterfaceGeometry> create_levelset_geometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
{
  return std::make_unique<LevelSetInterfaceGeometry>(activePart, cdfemSupport, phaseSupport, LSFields);
}

}
