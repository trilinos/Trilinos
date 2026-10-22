// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingSurface.hpp>
#include <Akri_Surface_Manager.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <algorithm>
#include <cctype>
#include <memory>

namespace krino {

static std::string get_uppercase_name(const std::string & name)
{
  std::string upperName = name;
  std::transform(upperName.begin(), upperName.end(), upperName.begin(),
      [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return upperName;
}

BoundingSurface::BoundingSurface(const std::string & surfName, std::shared_ptr<Surface> surf)
: myName(get_uppercase_name(surfName)),
  mySurface(std::move(surf))
{
}

BoundingSurface &
BoundingSurface::build(
    stk::mesh::MetaData & meta,
    const std::string & surfName,
    std::shared_ptr<Surface> surf)
{
  STK_ThrowRequire(surf != nullptr);

  Surface_Manager & manager = Surface_Manager::get(meta);

  const std::string normalizedName = get_uppercase_name(surfName);

  STK_ThrowRequire(!manager.has_bounding_surface(normalizedName));

  std::unique_ptr<BoundingSurface> boundingSurf(new BoundingSurface(normalizedName, std::move(surf)));
  BoundingSurface & result = *boundingSurf;

  manager.add_bounding_surface(boundingSurf.release());

  return result;
}

} // namespace krino
