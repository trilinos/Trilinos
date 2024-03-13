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
#include "Akri_LevelSetSurfaceInterfaceGeometry.hpp"

namespace krino {

std::unique_ptr<InterfaceGeometry> create_interface_geometry(const stk::mesh::MetaData & meta)
{
  const stk::mesh::Part & activePart = AuxMetaData::get(meta).active_part();
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
  const Phase_Support & phaseSupport = Phase_Support::get(meta);

  Surface_Manager & manager = Surface_Manager::get(meta);
  if (manager.get_bounding_surfaces().size() > 0)
    return create_bounding_surface_geometry(manager, activePart, cdfemSupport, phaseSupport);

  return create_levelset_geometry(meta.spatial_dimension(), activePart, cdfemSupport, phaseSupport, Phase_Support::get_levelset_fields(meta));
}

std::unique_ptr<InterfaceGeometry> create_bounding_surface_geometry(Surface_Manager & manager,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  std::unique_ptr<AnalyticSurfaceInterfaceGeometry> geom = std::make_unique<AnalyticSurfaceInterfaceGeometry>(activePart, cdfemSupport, phaseSupport);
  for (auto && boundingSurface : manager.get_bounding_surfaces())
    geom->add_surface(manager.get_identifier(boundingSurface->name()), boundingSurface->surface());
  return geom;
}

std::unique_ptr<InterfaceGeometry> create_levelset_geometry(const int dim,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
{
  const double facetsRequested = cdfemSupport.use_facets_instead_of_levelset_fields() ||
      (SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE == cdfemSupport.get_cdfem_edge_degeneracy_handling() &&
       (RESNAP_USING_INTERFACE_ON_PREVIOUS_SNAPPED_MESH == cdfemSupport.get_resnap_method() ||
        RESNAP_USING_INTERPOLATION == cdfemSupport.get_resnap_method()));
  if (facetsRequested && !phaseSupport.has_one_levelset_per_phase())
    return std::make_unique<LevelSetSurfaceInterfaceGeometry>(dim, activePart, cdfemSupport, phaseSupport, LSFields);
  return std::make_unique<LevelSetInterfaceGeometry>(activePart, cdfemSupport, phaseSupport, LSFields);
}

}
