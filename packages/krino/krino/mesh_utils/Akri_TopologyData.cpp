// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_TopologyData.hpp>

#include <stdexcept>
#include <string>

namespace krino {

stk::topology interface_topology(const stk::topology elementTopology)
{
  if (elementTopology.is_shell())
  {
    switch(elementTopology())
    {
    case stk::topology::SHELL_QUAD_4:
    case stk::topology::SHELL_TRI_3:
    case stk::topology::SHELL_TRI_3_ALL_FACE_SIDES:
    case stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES:
        return stk::topology::SHELL_SIDE_BEAM_2;
    case stk::topology::SHELL_QUAD_9:
    case stk::topology::SHELL_TRI_6:
    case stk::topology::SHELL_TRI_6_ALL_FACE_SIDES:
    case stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES:
      return stk::topology::SHELL_SIDE_BEAM_3;
    default:
        throw std::runtime_error("Shell topology not found in interface_topology(): " + elementTopology.name());
    }
  }
  return elementTopology.side_topology();
}

} // namespace krino
