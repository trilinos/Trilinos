// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Surface_Manager.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_BoundingSurface.hpp>
#include <Akri_LevelSet.hpp>

namespace krino {

Surface_Manager &
Surface_Manager::get(const stk::mesh::MetaData & meta)
{
  Surface_Manager * mgr = const_cast<Surface_Manager *>(meta.get_attribute<Surface_Manager>());
  STK_ThrowRequireMsg(nullptr != mgr, "No Surface_Manager found for MetaData.");
  return *mgr;
}

Surface_Manager &
Surface_Manager::get(stk::mesh::MetaData & meta)
{
  Surface_Manager * mgr = const_cast<Surface_Manager *>(meta.get_attribute<Surface_Manager>());
  if (nullptr == mgr)
  {
    mgr = new Surface_Manager();
    meta.declare_attribute_with_delete<Surface_Manager>(mgr);
  }
  return *mgr;
}

void Surface_Manager::set_surface_names(const SurfaceNames & surfaceNames)
{
  STK_ThrowRequireMsg(mySurfaceNames.get_names().empty(), "Surface_Manager::set_surface_names called after surfaces were already named.");
  mySurfaceNames = surfaceNames;
}

Surface_Identifier Surface_Manager::get_identifier(const std::string & name)
{
  return mySurfaceNames.get_identifier(name);
}

const std::string & Surface_Manager::get_name(const Surface_Identifier identifier) const
{
  return mySurfaceNames.get_names().at(identifier.get());
}

bool Surface_Manager::has_bounding_surface(const std::string & surfName) const
{
  for (auto&& boundingSurface : myBoundingSurfaces)
    if (boundingSurface->name() == surfName)
      return true;

  return false;
}

void Surface_Manager::add_bounding_surface(BoundingSurface * surf)
{
  myBoundingSurfaces.emplace_back(surf);
}

bool Surface_Manager::has_levelset(const std::string & surfName) const
{
  for (auto&& ls : myLevelSetSurfaces)
    if (ls->name() == surfName || ls->get_composite_name() == surfName)
      return true;

  return false;
}

void Surface_Manager::clear()
{
  myLevelSetSurfaces.clear();
  myBoundingSurfaces.clear();
}

BoundingSurface & Surface_Manager::get_bounding_surface(const std::string & surfName) const
{
  for (auto&& boundingSurface : myBoundingSurfaces)
    if (boundingSurface->name() == surfName)
      return *boundingSurface;

  STK_ThrowErrorMsg("No bounding surface with name " << surfName);
  return *myBoundingSurfaces.front(); // unreachable; STK_ThrowErrorMsg does not return
}

LevelSet & Surface_Manager::get_levelset(const std::string & surfName) const
{
  for (auto&& ls : myLevelSetSurfaces)
    if (ls->name() == surfName || ls->get_composite_name() == surfName)
      return *ls;

  STK_ThrowErrorMsg("No levelset with name " << surfName);
  return *myLevelSetSurfaces.front(); // unreachable; STK_ThrowErrorMsg does not return
}

void Surface_Manager::add_levelset(LevelSet * ls)
{
  myLevelSetSurfaces.emplace_back(ls);
}

} // namespace krino
