// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Surface_Manager.hpp>
#include <algorithm>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_BoundingSurface.hpp>
#include <Akri_LevelSet.hpp>

namespace krino {

std::map<std::string,std::unique_ptr<Surface_Manager>> Surface_Manager::theModeltoManagerMap;

Surface_Manager & Surface_Manager::get_or_create(const std::string & FEModelName)
{
  std::string upperName = FEModelName;
  std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
  auto iter = theModeltoManagerMap.find(upperName);

  if (iter != theModeltoManagerMap.end())
    return *(iter->second);

  Surface_Manager * mgr = new Surface_Manager();
  theModeltoManagerMap[upperName] = std::unique_ptr<Surface_Manager>(mgr);
  return *mgr;
}

void Surface_Manager::associate_FEModel_and_metadata(const std::string & FEModelName, stk::mesh::MetaData & meta)
{
  Surface_Manager * existingManagerOnMeta = const_cast<Surface_Manager *>(meta.get_attribute<Surface_Manager>());
  Surface_Manager * mgr = &get_or_create(FEModelName);
  STK_ThrowRequireMsg(nullptr == existingManagerOnMeta || existingManagerOnMeta == mgr,
    "krino::Surface_Manager already set on stk::mesh::MetaData and it doesn't match the one associated with the FEModel " << FEModelName);
  meta.declare_attribute_no_delete<Surface_Manager>(mgr);
}

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

Surface_Identifier Surface_Manager::get_identifier(const std::string & name)
{
  std::string upper_name = name;
  std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(), ::toupper);
  for (unsigned i=0; i<mySurfaceNames.size(); ++i)
    if (upper_name == mySurfaceNames[i])
      return Surface_Identifier(i);

  const unsigned id = mySurfaceNames.size();
  mySurfaceNames.push_back(upper_name);
  return Surface_Identifier(id);
}


const std::string & Surface_Manager::get_name(const Surface_Identifier identifier) const
{
  return mySurfaceNames.at(identifier.get());
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

void Surface_Manager::add_levelset(LevelSet * ls)
{
  myLevelSetSurfaces.emplace_back(ls);
}

} // namespace krino
