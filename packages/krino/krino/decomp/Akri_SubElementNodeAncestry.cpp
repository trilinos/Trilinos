// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>

namespace krino {

void SubElementNodeAncestry::print(std::ostream & os) const
{
  if (my_node->is_mesh_node())
  {
    os << my_node->entityId();
  }
  else
  {
    os << "{ ";
    for (auto&& parent : get_parents())
    {
      parent.print(os);
      os << " ";
    }
    os << "}";
  }
}

}
