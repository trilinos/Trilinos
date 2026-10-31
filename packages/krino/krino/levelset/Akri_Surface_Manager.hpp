// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_LEVELSET_AKRI_SURFACE_MANAGER_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_SURFACE_MANAGER_HPP_
#include <memory>
#include <string>
#include <vector>
#include <Akri_PhaseInput.hpp>
#include <Akri_Surface_Identifier.hpp>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class BoundingSurface;
class LevelSet;

class Surface_Manager
{
public:
  static Surface_Manager & get(const stk::mesh::MetaData & meta);
  static Surface_Manager & get(stk::mesh::MetaData & meta); // creates the Surface_Manager for the MetaData if needed

  // Sets the surface names parsed for this mesh's finite element model, so the surface identifiers
  // used in its parsed phases are valid here.  Must be called before any other surface is named.
  void set_surface_names(const SurfaceNames & surfaceNames);
  Surface_Identifier get_identifier(const std::string & name);
  const std::string & get_name(const Surface_Identifier identifier) const;

  const std::vector< std::unique_ptr<BoundingSurface> > & get_bounding_surfaces() const { return myBoundingSurfaces; }
  bool has_bounding_surface(const std::string & name) const;
  BoundingSurface & get_bounding_surface(const std::string & surfName) const;
  void add_bounding_surface(BoundingSurface * surf);

  const std::vector< std::unique_ptr<LevelSet> > & get_levelsets() const { return myLevelSetSurfaces; }
  bool has_levelset(const std::string & name) const;
  LevelSet & get_levelset(const std::string & surfName) const;
  void add_levelset(LevelSet * ls);
  void clear();

private:
  SurfaceNames mySurfaceNames;
  std::vector< std::unique_ptr<BoundingSurface> > myBoundingSurfaces;
  std::vector< std::unique_ptr<LevelSet> > myLevelSetSurfaces;
};

} // namespace krino



#endif /* KRINO_KRINO_LEVELSET_AKRI_SURFACE_MANAGER_HPP_ */
