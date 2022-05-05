// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_LIB_AKRI_SURFACE_MANAGER_HPP_
#define KRINO_LIB_AKRI_SURFACE_MANAGER_HPP_
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <Akri_Surface_Identifier.hpp>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class BoundingSurface;
class LevelSet;

class Surface_Manager
{
public:
  static Surface_Manager & get_or_create(const std::string & FEModelName);
  static void associate_FEModel_and_metadata(const std::string & FEModelName, stk::mesh::MetaData & meta); // for typical parsing scenario where Surface_Manager is originally created with the FEModel and now is being associated with MetaData
  static Surface_Manager & get(const stk::mesh::MetaData & meta);
  static Surface_Manager & get(stk::mesh::MetaData & meta); // for unit test usage where the Surface_Manager is created with the MetaData

  Surface_Identifier get_identifier(const std::string & name);
  const std::string & get_name(const Surface_Identifier identifier) const;

  const std::vector< std::unique_ptr<BoundingSurface> > & get_bounding_surfaces() const { return myBoundingSurfaces; }
  bool has_bounding_surface(const std::string & name) const;
  void add_bounding_surface(BoundingSurface * surf);

  const std::vector< std::unique_ptr<LevelSet> > & get_levelsets() const { return myLevelSetSurfaces; }
  bool has_levelset(const std::string & name) const;
  void add_levelset(LevelSet * ls);

private:
  static std::map<std::string,std::unique_ptr<Surface_Manager>> theModeltoManagerMap;
  std::vector<std::string> mySurfaceNames;
  std::vector< std::unique_ptr<BoundingSurface> > myBoundingSurfaces;
  std::vector< std::unique_ptr<LevelSet> > myLevelSetSurfaces;
};

} // namespace krino



#endif /* KRINO_LIB_AKRI_SURFACE_MANAGER_HPP_ */
