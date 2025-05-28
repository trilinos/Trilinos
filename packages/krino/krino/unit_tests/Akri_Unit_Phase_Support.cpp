// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Phase_Support.hpp>

#include <Akri_Unit_Part_Decomposition_Fixture.hpp>
#include <Ioss_Initializer.h>
#include <Ioss_ConcreteVariableType.h>

Ioss::StorageInitializer initialize_storage;
Ioss::Initializer        initialize_topologies;

namespace krino
{

TEST_F(Part_Decomposition_Fixture, One_Block_LS)
{
  Block_Surface_Connectivity block_surface_info;
  performDecomposition({findPart("block_1")}, block_surface_info, false);

  assert_conformal_part_exists("block_1_A", "block_1_nonconformal");
  assert_conformal_part_exists("block_1_B", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_A_B", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_B_A", "block_1_nonconformal");
}

TEST_F(Part_Decomposition_Fixture, One_Block_Death)
{
  Block_Surface_Connectivity block_surface_info;
  performDecomposition({findPart("block_1")}, block_surface_info, true);

  assert_conformal_part_exists("block_1", "block_1_nonconformal");
  assert_conformal_part_exists("block_1_dead", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_test", "block_1_nonconformal");
}

TEST_F(Part_Decomposition_Fixture, OneSidedSideset_LS)
{
  Block_Surface_Connectivity block_surface_info = addOneSidedSideset();
  performDecomposition({findPart("block_1")}, block_surface_info, false);

  assert_conformal_part_exists("surface_1_A", "surface_1_nonconformal");
  assert_conformal_part_exists("surface_1_B", "surface_1_nonconformal");
}

TEST_F(Part_Decomposition_Fixture, OneSidedSideset_Death)
{
  Block_Surface_Connectivity block_surface_info = addOneSidedSideset();
  performDecomposition({findPart("block_1")}, block_surface_info, true);

  assert_conformal_part_exists("surface_1", "surface_1_nonconformal");
  assert_conformal_part_exists("surface_1_dead", "surface_1_nonconformal");
}

TEST_F(Part_Decomposition_Fixture, TwoSidedSideset_LS)
{
  Block_Surface_Connectivity block_surface_info = addTwoSidedSideset();
  performDecomposition({findPart("block_1")}, block_surface_info, false);

  //const stk::mesh::Part * surf_1_A_block_1 = findPart("surface_1_A_block_1"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_A_block_1 = findPart("surface_1_A_block_1_A_tri3");
  ASSERT_TRUE( surf_1_A_block_1 != NULL );
  EXPECT_TRUE( findSuperset("surface_1_A", surf_1_A_block_1) != NULL );

  //const stk::mesh::Part * surf_1_B_block_1 = findPart("surface_1_B_block_1"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_B_block_1 = findPart("surface_1_B_block_1_B_tri3");
  ASSERT_TRUE( surf_1_B_block_1 != NULL );
  EXPECT_TRUE( findSuperset("surface_1_B", surf_1_B_block_1) != NULL );

  //const stk::mesh::Part * surf_1_A_block_2 = findPart("surface_1_A_block_2"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_A_block_2 = findPart("surface_1_A_block_2_tri3");
  ASSERT_TRUE( surf_1_A_block_2 != NULL );
  EXPECT_TRUE( findSuperset("surface_1_A", surf_1_A_block_2) != NULL );

  //const stk::mesh::Part * surf_1_B_block_2 = findPart("surface_1_B_block_2"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_B_block_2 = findPart("surface_1_B_block_2_tri3");
  ASSERT_TRUE( surf_1_B_block_2 != NULL );
  EXPECT_TRUE( findSuperset("surface_1_B", surf_1_B_block_2) != NULL );
}

TEST_F(Part_Decomposition_Fixture, TwoSidedSideset_Death)
{
  Block_Surface_Connectivity block_surface_info = addTwoSidedSideset();
  performDecomposition({findPart("block_1")}, block_surface_info, true);

  //const stk::mesh::Part * surf_1_block_1 = findPart("surface_1_block_1"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_block_1 = findPart("surface_block_1_tri3_1");
  ASSERT_TRUE( surf_1_block_1 != NULL );
  EXPECT_TRUE( findSuperset("surface_1", surf_1_block_1) != NULL );

  //const stk::mesh::Part * surf_1_dead_block_1 = findPart("surface_1_dead_block_1"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_dead_block_1 = findPart("surface_1_dead_block_1_dead_tri3");
  ASSERT_TRUE( surf_1_dead_block_1 != NULL );
  EXPECT_TRUE( findSuperset("surface_1_dead", surf_1_dead_block_1) != NULL );

  //const stk::mesh::Part * surf_1_block_2 = findPart("surface_1_block_2"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_block_2 = findPart("surface_block_2_tri3_1");
  ASSERT_TRUE( surf_1_block_2 != NULL );
  EXPECT_TRUE( findSuperset("surface_1", surf_1_block_2) != NULL );

  //const stk::mesh::Part * surf_1_dead_block_2 = findPart("surface_1_dead_block_2"); // No support for aliases in stk yet
  const stk::mesh::Part * surf_1_dead_block_2 = findPart("surface_1_dead_block_2_tri3");
  ASSERT_TRUE( surf_1_dead_block_2 != NULL );
  EXPECT_TRUE( findSuperset("surface_1_dead", surf_1_dead_block_2) != NULL );
}

TEST_F(Part_Decomposition_Fixture, Multiple_LS_Decomposition)
{
  Block_Surface_Connectivity block_surface_info;
  performDecomposition({findPart("block_1")}, block_surface_info, false, 2);

  assert_conformal_part_exists("block_1_A", "block_1_nonconformal");
  assert_conformal_part_exists("block_1_B", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_A_B", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_B_A", "block_1_nonconformal");

  assert_conformal_part_exists("block_1_C", "block_1_nonconformal");
  assert_conformal_part_exists("block_1_D", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_C_D", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_D_C", "block_1_nonconformal");

  assert_conformal_part_exists("surface_block_1_A_C", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_C_A", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_A_D", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_D_A", "block_1_nonconformal");

  assert_conformal_part_exists("surface_block_1_B_C", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_C_B", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_B_D", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_D_B", "block_1_nonconformal");
}

TEST_F(Part_Decomposition_Fixture, find_conformal_io_part)
{
  std::vector<std::string> decomposed_blocks;
  decomposed_blocks.push_back("block_1");
  Block_Surface_Connectivity block_surface_info;

  stk::mesh::Part * block_1 = findPart("block_1");
  performDecomposition({block_1}, block_surface_info, false, 2);

  // Test volume conformal io part lookup
  const stk::mesh::Part & block_1_A = phase_support().find_conformal_io_part(*block_1, get_phase_A());
  EXPECT_EQ( "block_1_A", block_1_A.name() );

  const stk::mesh::Part & block_1_B = phase_support().find_conformal_io_part(*block_1, get_phase_B());
  EXPECT_EQ( "block_1_B", block_1_B.name() );

  const stk::mesh::Part & block_1_C = phase_support().find_conformal_io_part(*block_1, get_phase_C());
  EXPECT_EQ( "block_1_C", block_1_C.name() );

  const stk::mesh::Part & block_1_D = phase_support().find_conformal_io_part(*block_1, get_phase_D());
  EXPECT_EQ( "block_1_D", block_1_D.name() );

  const stk::mesh::Part * surface_block_1_A_B = phase_support().find_interface_part(block_1_A, block_1_B);
  ASSERT_TRUE( surface_block_1_A_B != NULL );
  EXPECT_EQ( "surface_block_1_A_B", surface_block_1_A_B->name() );

  const stk::mesh::Part * surface_block_1_A_C = phase_support().find_interface_part(block_1_A, block_1_C);
  ASSERT_TRUE( surface_block_1_A_C != NULL );
  EXPECT_EQ( "surface_block_1_A_C", surface_block_1_A_C->name() );

  const stk::mesh::Part * surface_block_1_A_D = phase_support().find_interface_part(block_1_A, block_1_D);
  ASSERT_TRUE( surface_block_1_A_D != NULL );
  EXPECT_EQ( "surface_block_1_A_D", surface_block_1_A_D->name() );

  const stk::mesh::Part * surface_block_1_C_D = phase_support().find_interface_part(block_1_C, block_1_D);
  ASSERT_TRUE( surface_block_1_C_D != NULL );
  EXPECT_EQ( "surface_block_1_C_D", surface_block_1_C_D->name() );
}

TEST_F(Part_Decomposition_Fixture, get_blocks_touching_surface)
{
  Block_Surface_Connectivity block_surface_info = addTwoSidedSideset();
  performDecomposition({findPart("block_1")}, block_surface_info, false);

  std::vector<const stk::mesh::Part*> blocks;

  blocks = meta_data().get_blocks_touching_surface(meta_data().get_part("surface_1_A"));
  EXPECT_EQ( 2u, blocks.size() );
  EXPECT_TRUE( blocks.end() != std::find(blocks.begin(), blocks.end(), meta_data().get_part("block_1_A")));
  EXPECT_TRUE( blocks.end() != std::find(blocks.begin(), blocks.end(), meta_data().get_part("block_2")));

  blocks = meta_data().get_blocks_touching_surface(meta_data().get_part("surface_1_A_block_1_A_tri3"));
  EXPECT_EQ( 1u, blocks.size() );
  EXPECT_TRUE( blocks.end() != std::find(blocks.begin(), blocks.end(), meta_data().get_part("block_1_A")));

  blocks = meta_data().get_blocks_touching_surface(meta_data().get_part("surface_1_A_block_2_tri3"));
  EXPECT_EQ( 1u, blocks.size() );
  EXPECT_TRUE( blocks.end() != std::find(blocks.begin(), blocks.end(), meta_data().get_part("block_2")));

  blocks = meta_data().get_blocks_touching_surface(meta_data().get_part("surface_block_1_A_B"));
  EXPECT_EQ( 1u, blocks.size() );
  EXPECT_TRUE( blocks.end() != std::find(blocks.begin(), blocks.end(), meta_data().get_part("block_1_A")));

  blocks = meta_data().get_blocks_touching_surface(meta_data().get_part("surface_block_1_B_A"));
  EXPECT_EQ( 1u, blocks.size() );
  EXPECT_TRUE( blocks.end() != std::find(blocks.begin(), blocks.end(), meta_data().get_part("block_1_B")));

}

TEST_F(Part_Decomposition_Fixture, One_Block_Two_LS_One_LS_Per_Phase)
{
  Block_Surface_Connectivity block_surface_info;
  performDecomposition({findPart("block_1")}, block_surface_info, false, 2, true);

  assert_conformal_part_exists("block_1_LS1", "block_1_nonconformal");
  assert_conformal_part_exists("block_1_LS2", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_LS1_LS2", "block_1_nonconformal");
  assert_conformal_part_exists("surface_block_1_LS2_LS1", "block_1_nonconformal");
}

} // namespace krino
