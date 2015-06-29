// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>

namespace {
static const unsigned num_test_topologies = 30;
static const stk::topology test_topologies[num_test_topologies] = {
                                                       stk::topology::NODE
                                                     //EDGE_RANK
                                                     , stk::topology::LINE_2
                                                     , stk::topology::LINE_3
                                                     //FACE_RANK
                                                     , stk::topology::TRI_3
                                                     , stk::topology::TRI_4
                                                     , stk::topology::TRI_6
                                                     , stk::topology::QUAD_4
                                                     , stk::topology::QUAD_8
                                                     , stk::topology::QUAD_9
                                                     //ELEMENT_RANK
                                                     , stk::topology::PARTICLE
                                                     , stk::topology::BEAM_2
                                                     , stk::topology::BEAM_3
                                                     , stk::topology::SHELL_TRI_3
                                                     // Cannot create SHELL_TRI_4 entities because Shards does not support!
                                                     // , stk::topology::SHELL_TRI_4
                                                     , stk::topology::SHELL_TRI_6
                                                     , stk::topology::SHELL_QUAD_4
                                                     , stk::topology::SHELL_QUAD_8
                                                     , stk::topology::SHELL_QUAD_9
                                                     , stk::topology::TET_4
                                                     , stk::topology::TET_8
                                                     , stk::topology::TET_10
                                                     , stk::topology::TET_11
                                                     , stk::topology::PYRAMID_5
                                                     , stk::topology::PYRAMID_13
                                                     , stk::topology::PYRAMID_14
                                                     , stk::topology::WEDGE_6
                                                     , stk::topology::WEDGE_15
                                                     , stk::topology::WEDGE_18
                                                     , stk::topology::HEX_8
                                                     , stk::topology::HEX_20
                                                     , stk::topology::HEX_27};

}

TEST(UnitTestRootTopology, TestRootTopologyPartGetters)
{
  const int spatial_dim = 3;
  stk::mesh::MetaData meta(spatial_dim);
  meta.commit();

  for (unsigned i = 0; i < num_test_topologies; ++i)
  {
    stk::mesh::Part &root_part1 = meta.get_topology_root_part(test_topologies[i]);

    stk::mesh::CellTopology cell_topo = stk::mesh::get_cell_topology(test_topologies[i]);
    stk::mesh::Part &root_part2 = meta.get_cell_topology_root_part(cell_topo);

    EXPECT_TRUE(root_part1 == root_part2);

    stk::topology root_topo = root_part1.topology();

    // The root_topology_part has the same topology information as the original topology.
    EXPECT_TRUE(test_topologies[i] == root_topo);

    const stk::mesh::PartVector &subsets = root_part1.subsets();
    const stk::mesh::PartVector &supersets = root_part1.supersets();

    EXPECT_TRUE(subsets.empty());
    EXPECT_EQ(1u, supersets.size());
    EXPECT_EQ(&meta.universal_part(), supersets[0]);

    stk::mesh::Part &root_part3 = meta.get_part(root_part1.mesh_meta_data_ordinal());
    EXPECT_TRUE(root_part1 == root_part3);

    stk::mesh::Part *root_part4 = meta.get_part(root_part1.name());
    EXPECT_TRUE(root_part1 == *root_part4);

    const stk::mesh::PartVector &all_parts = meta.get_parts();
    const stk::mesh::PartVector non_internal_parts = meta.get_mesh_parts();

    EXPECT_TRUE(std::find(all_parts.begin(), all_parts.end(), &root_part1) != all_parts.end());
    EXPECT_TRUE(std::find(non_internal_parts.begin(), non_internal_parts.end(), &root_part1)
                == non_internal_parts.end());
  }
}

TEST(UnitTestRootTopology, TestRootTopologySubsets)
{
  const int spatial_dim = 3;
  stk::mesh::MetaData meta(spatial_dim);

  for (unsigned i = 0; i < num_test_topologies; ++i)
  {
    stk::mesh::Part &root_part1 = meta.get_topology_root_part(test_topologies[i]);
    stk::topology root_topo = root_part1.topology();

    EXPECT_FALSE(root_part1.force_no_induce());

    // The root_topology_part has the same topology information as the original topology.
    EXPECT_TRUE(test_topologies[i] == root_topo);

    const stk::mesh::PartVector &rootTopoPartSubsets = root_part1.subsets();
    const stk::mesh::PartVector &rootTopoPartSupersets = root_part1.supersets();

    EXPECT_TRUE(rootTopoPartSubsets.empty());
    EXPECT_EQ(1u, rootTopoPartSupersets.size());
    EXPECT_EQ(&meta.universal_part(), rootTopoPartSupersets[0]);

    stk::mesh::Part & newPart = meta.declare_part_with_topology( std::string("topo_part") + root_part1.name(), root_topo );

    EXPECT_EQ(1u, rootTopoPartSubsets.size());
    EXPECT_EQ(&newPart, rootTopoPartSubsets[0]);

    const stk::mesh::PartVector & newPartSupersets = newPart.supersets();

    EXPECT_EQ(2u, newPartSupersets.size());
    EXPECT_EQ(&meta.universal_part(), newPartSupersets[0]);
    EXPECT_EQ(&root_part1, newPartSupersets[1]);
  }


  meta.commit();

}
