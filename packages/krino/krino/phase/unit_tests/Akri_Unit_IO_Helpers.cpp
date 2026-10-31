// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_IO_Helpers.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <memory>
#include <set>

namespace krino {

namespace {
using OrdSet = std::set<stk::mesh::PartOrdinal>;

// Named ordinals; blocks and surfaces share the ordinal space but are kept
// in disjoint ranges here for readability.
constexpr stk::mesh::PartOrdinal BLOCK_1 = 10;
constexpr stk::mesh::PartOrdinal BLOCK_2 = 11;
constexpr stk::mesh::PartOrdinal SURF_1 = 20;
constexpr stk::mesh::PartOrdinal SURF_2 = 21;
}

//--------------------------------------------------------------------
// Block_Surface_Connectivity
//--------------------------------------------------------------------

TEST(Block_Surface_Connectivity, DefaultConstructed_ReturnsEmptyForAnyOrdinal)
{
  const Block_Surface_Connectivity connectivity;
  EXPECT_TRUE(connectivity.get_surfaces_touching_block(BLOCK_1).empty());
  EXPECT_TRUE(connectivity.get_blocks_touching_surface(SURF_1).empty());
}

TEST(Block_Surface_Connectivity, AddSurface_IsQueryableFromBothDirections)
{
  Block_Surface_Connectivity connectivity;
  connectivity.add_surface(SURF_1, {BLOCK_1, BLOCK_2});

  EXPECT_EQ(OrdSet({BLOCK_1, BLOCK_2}), connectivity.get_blocks_touching_surface(SURF_1));
  EXPECT_EQ(OrdSet({SURF_1}), connectivity.get_surfaces_touching_block(BLOCK_1));
  EXPECT_EQ(OrdSet({SURF_1}), connectivity.get_surfaces_touching_block(BLOCK_2));
}

TEST(Block_Surface_Connectivity, UnknownOrdinal_ReturnsEmpty)
{
  Block_Surface_Connectivity connectivity;
  connectivity.add_surface(SURF_1, {BLOCK_1});

  EXPECT_TRUE(connectivity.get_surfaces_touching_block(BLOCK_2).empty());
  EXPECT_TRUE(connectivity.get_blocks_touching_surface(SURF_2).empty());
}

TEST(Block_Surface_Connectivity, BlockTouchingMultipleSurfaces_AccumulatesSurfaces)
{
  Block_Surface_Connectivity connectivity;
  connectivity.add_surface(SURF_1, {BLOCK_1});
  connectivity.add_surface(SURF_2, {BLOCK_1});

  EXPECT_EQ(OrdSet({SURF_1, SURF_2}), connectivity.get_surfaces_touching_block(BLOCK_1));
}

TEST(Block_Surface_Connectivity, RepeatedAddOfSameSurface_MergesTouchingBlocks)
{
  Block_Surface_Connectivity connectivity;
  connectivity.add_surface(SURF_1, {BLOCK_1});
  connectivity.add_surface(SURF_1, {BLOCK_2});

  // The two adds accumulate rather than replace.
  EXPECT_EQ(OrdSet({BLOCK_1, BLOCK_2}), connectivity.get_blocks_touching_surface(SURF_1));
  EXPECT_EQ(OrdSet({SURF_1}), connectivity.get_surfaces_touching_block(BLOCK_1));
  EXPECT_EQ(OrdSet({SURF_1}), connectivity.get_surfaces_touching_block(BLOCK_2));
}

TEST(Block_Surface_Connectivity, AddSurfaceWithNoTouchingBlocks_HasEmptyBlockSet)
{
  Block_Surface_Connectivity connectivity;
  connectivity.add_surface(SURF_1, {});

  EXPECT_TRUE(connectivity.get_blocks_touching_surface(SURF_1).empty());
}

//--------------------------------------------------------------------
// Ioss sideblock naming
//--------------------------------------------------------------------

TEST(Ioss_Sideblock_Naming, GenerateName_CanonicalSidesetName_InsertsBlockAndTopologyBeforeId)
{
  EXPECT_EQ("surface_block_1_tri3_1", generate_ioss_sideblock_name("surface_1", "block_1", "tri3"));
  EXPECT_EQ("surface_tetra4_tri3_12", generate_ioss_sideblock_name("surface_12", "tetra4", "tri3"));
}

TEST(Ioss_Sideblock_Naming, GenerateName_NonCanonicalSidesetName_AppendsBlockAndTopology)
{
  EXPECT_EQ("my_surf_block_1_tri3", generate_ioss_sideblock_name("my_surf", "block_1", "tri3"));
  EXPECT_EQ("surface_1_A_block_1_A_tri3", generate_ioss_sideblock_name("surface_1_A", "block_1_A", "tri3"));
  EXPECT_EQ("surface_x_block_1_tri3", generate_ioss_sideblock_name("surface_x", "block_1", "tri3"));
}

TEST(Ioss_Sideblock_Naming, GetIossTopologyName_MapsStkTopologies)
{
  EXPECT_EQ("tri3", get_ioss_topology_name(stk::topology::TRIANGLE_3));
  EXPECT_EQ("tetra4", get_ioss_topology_name(stk::topology::TETRAHEDRON_4));
}

TEST(Ioss_Sideblock_Naming, IsSplitByTopology_DetectsBlockTopologyInName)
{
  EXPECT_TRUE(is_sideblock_split_by_topology("surface_tetra4_tri3_1", stk::topology::TETRAHEDRON_4));
  EXPECT_FALSE(is_sideblock_split_by_topology("surface_block_1_tri3_1", stk::topology::TETRAHEDRON_4));
}

//--------------------------------------------------------------------
// get_parent_sideset
//--------------------------------------------------------------------

class GetParentSidesetFixture : public ::testing::Test
{
protected:
  GetParentSidesetFixture()
    : bulk(stk::mesh::MeshBuilder(MPI_COMM_SELF).set_spatial_dimension(3).create()),
      meta(bulk->mesh_meta_data())
  {
  }

  stk::mesh::Part & declare_io_side_part(const std::string & name)
  {
    stk::mesh::Part & part = meta.declare_part_with_topology(name, stk::topology::TRIANGLE_3);
    stk::io::put_io_part_attribute(part);
    return part;
  }

  std::unique_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData & meta;
};

TEST_F(GetParentSidesetFixture, SideblockOfIoSideset_ReturnsSideset)
{
  stk::mesh::Part & sideset = meta.declare_part("surface_1", meta.side_rank());
  stk::io::put_io_part_attribute(sideset);
  stk::mesh::Part & sideblock = declare_io_side_part("surface_block_1_tri3_1");
  meta.declare_part_subset(sideset, sideblock);

  EXPECT_EQ(&sideset, get_parent_sideset(sideblock));
  EXPECT_EQ(nullptr, get_parent_sideset(sideset));
}

TEST_F(GetParentSidesetFixture, AssemblySuperset_IsNotAParentSideset)
{
  stk::mesh::Part & assembly = meta.declare_part("assembly_1", meta.side_rank());
  stk::io::put_assembly_io_part_attribute(assembly);
  stk::mesh::Part & sideset = declare_io_side_part("surface_1");
  meta.declare_part_subset(assembly, sideset);

  EXPECT_EQ(nullptr, get_parent_sideset(sideset));
}

TEST_F(GetParentSidesetFixture, NonSideRankPart_HasNoParentSideset)
{
  stk::mesh::Part & block = meta.declare_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
  stk::io::put_io_part_attribute(block);
  EXPECT_EQ(nullptr, get_parent_sideset(block));
}

} // namespace krino
