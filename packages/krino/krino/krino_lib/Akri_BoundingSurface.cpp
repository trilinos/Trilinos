/*
 * Akri_BoundingSurface.cpp
 *
 *  Created on: Jan 24, 2022
 *      Author: drnoble
 */
#include <Akri_BoundingSurface.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_Surface_Manager.hpp>

namespace krino {

BoundingSurface::BoundingSurface(const std::string & surfName, Surface * surf)
: myName(surfName),
  mySurface(surf)
{
  std::transform(myName.begin(), myName.end(), myName.begin(), ::toupper);
}

BoundingSurface &
BoundingSurface::build(
    stk::mesh::MetaData & meta,
    const std::string & surfName,
    Surface * surf)
{
  Surface_Manager & manager = Surface_Manager::get(meta);

  STK_ThrowRequire(!manager.has_bounding_surface(surfName));
  BoundingSurface * boundingSurf = new BoundingSurface(surfName, surf);
  manager.add_bounding_surface(boundingSurf);
  return *boundingSurf;
}

}

