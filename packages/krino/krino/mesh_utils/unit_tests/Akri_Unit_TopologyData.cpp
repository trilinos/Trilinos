// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_TopologyData.hpp>

#include <stk_topology/topology.hpp>

#include <stdexcept>

namespace krino {

//--------------------------------------------------------------------
// interface_topology
//--------------------------------------------------------------------

TEST(Interface_Topology, NonShell_ReturnsSideTopology)
{
  // For non-shell topologies the interface is just the element's side topology.
  EXPECT_EQ(stk::topology(stk::topology::TRI_3).side_topology(),
      interface_topology(stk::topology::TRI_3));
  EXPECT_EQ(stk::topology(stk::topology::TET_4).side_topology(),
      interface_topology(stk::topology::TET_4));
  EXPECT_EQ(stk::topology(stk::topology::HEX_8).side_topology(),
      interface_topology(stk::topology::HEX_8));
  EXPECT_EQ(stk::topology(stk::topology::QUAD_4_2D).side_topology(),
      interface_topology(stk::topology::QUAD_4_2D));
}

TEST(Interface_Topology, LinearShells_ReturnBeam2)
{
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_2, interface_topology(stk::topology::SHELL_QUAD_4));
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_2, interface_topology(stk::topology::SHELL_TRI_3));
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_2, interface_topology(stk::topology::SHELL_TRI_3_ALL_FACE_SIDES));
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_2, interface_topology(stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES));
}

TEST(Interface_Topology, QuadraticShells_ReturnBeam3)
{
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_3, interface_topology(stk::topology::SHELL_QUAD_9));
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_3, interface_topology(stk::topology::SHELL_TRI_6));
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_3, interface_topology(stk::topology::SHELL_TRI_6_ALL_FACE_SIDES));
  EXPECT_EQ(stk::topology::SHELL_SIDE_BEAM_3, interface_topology(stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES));
}

TEST(Interface_Topology, UnrecognizedShell_Throws)
{
  EXPECT_THROW(interface_topology(stk::topology::SHELL_LINE_2), std::runtime_error);
}

} // namespace krino
