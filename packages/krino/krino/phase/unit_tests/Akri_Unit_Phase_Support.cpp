// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Phase_Support.hpp>
#include <Akri_PhaseInput.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <Akri_Unit_Part_Decomposition_Fixture.hpp>
#include <Ioss_Initializer.h>
#include <Ioss_ConcreteVariableType.h>

#include <set>
#include <string>

Ioss::StorageInitializer initialize_storage;
Ioss::Initializer        initialize_topologies;

namespace krino
{

TEST_F(Part_Decomposition_Fixture, SetPhaseInput_CopiesParsedPhases)
{
  PhaseTag negative;
  negative.add(Surface_Identifier(0), -1);
  PhaseInput phaseInput;
  phaseInput.meshPhases = {NamedPhase("A", negative), NamedPhase("B")};
  phaseInput.blockPhasesByName["block_1"] = {"A", "B"};
  phaseInput.oneLevelSetPerPhase = true;

  phase_support().set_phase_input(phaseInput);

  ASSERT_EQ(2u, phase_support().get_mesh_phases().size());
  EXPECT_EQ("A", phase_support().get_mesh_phases()[0].name());
  EXPECT_TRUE(phase_support().get_mesh_phases()[0].tag() == negative);
  EXPECT_EQ("B", phase_support().get_mesh_phases()[1].name());
  EXPECT_EQ(phaseInput.blockPhasesByName, phase_support().get_block_phases_by_name());
  EXPECT_TRUE(phase_support().has_one_levelset_per_phase());
}

TEST_F(Part_Decomposition_Fixture, SetPhaseInput_AfterPhasesAreSet_Throws)
{
  PhaseInput phaseInput;
  phaseInput.meshPhases = {NamedPhase("A")};
  phase_support().set_phase_input(phaseInput);
  EXPECT_ANY_THROW(phase_support().set_phase_input(phaseInput));
}

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

static bool contains(const std::string & text, const std::string & substring)
{
  return text.find(substring) != std::string::npos;
}

TEST_F(Part_Decomposition_Fixture, TwoSidedSideset_Two_Death_Criteria_Validation_Reports_Doomed)
{
  // surface_1 touches block_1 and block_2, which are decomposed by distinct death criteria.  Both
  // criteria name their phases "" (alive) and "dead", so surface_1 touches blocks with different
  // criteria for phases with the same name.  This is reported at parse time.
  Block_Surface_Connectivity block_surface_info = addTwoSidedSideset();
  EXPECT_FALSE( validate(make_two_death_decompositions(false), block_surface_info) );
  EXPECT_EQ(1u, stk::get_doomed_count());

  const std::string & msg = get_reported_doomed_messages();
  EXPECT_TRUE( contains(msg, "Sideset surface_1 touches blocks that are decomposed with different criteria") ) << msg;
  EXPECT_TRUE( contains(msg, "Phase \"\" (the unnamed alive phase of CDFEM death)") ) << msg;
  EXPECT_TRUE( contains(msg, "Phase \"dead\"") ) << msg;
  EXPECT_TRUE( contains(msg, "(CDFEM death death_by_temp_1): blocks block_1\n") ) << msg;
  EXPECT_TRUE( contains(msg, "(CDFEM death death_by_temp_2): blocks block_2\n") ) << msg;
  EXPECT_FALSE( contains(msg, "not decomposed") ) << msg;
  EXPECT_TRUE( contains(msg, "replaced by separate sidesets") ) << msg;
}

TEST_F(Part_Decomposition_Fixture, ThreeSidedSideset_Two_Death_Criteria_Validation_Reports_Undecomposed_Block)
{
  // As above, with surface_1 also touching the undecomposed block_3.
  Block_Surface_Connectivity block_surface_info = addThreeSidedSideset();
  EXPECT_FALSE( validate(make_two_death_decompositions(false), block_surface_info) );
  EXPECT_EQ(1u, stk::get_doomed_count());

  const std::string & msg = get_reported_doomed_messages();
  EXPECT_TRUE( contains(msg, "Sideset surface_1 touches blocks") ) << msg;
  EXPECT_TRUE( contains(msg, "surface_1 also touches blocks that are not decomposed: block_3\n") ) << msg;
}

TEST_F(Part_Decomposition_Fixture, Subset_Touching_Two_Death_Criteria_Validation_Reports_Doomed)
{
  // A single topology-split subset of surface_1 touches both blocks.
  Block_Surface_Connectivity block_surface_info = addSidesetWithSubsetTouchingTwoBlocks();
  EXPECT_FALSE( validate(make_two_death_decompositions(false), block_surface_info) );
  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_TRUE( contains(get_reported_doomed_messages(), "Sideset surface_1 touches blocks") ) << get_reported_doomed_messages();
}

TEST_F(Part_Decomposition_Fixture, OneSidedSideset_Two_Death_Criteria_Validation_Passes)
{
  // surface_1 touches only block_1, so the two criteria never meet on a sideset.
  Block_Surface_Connectivity block_surface_info = addOneSidedSideset();
  EXPECT_TRUE( validate(make_two_death_decompositions(false), block_surface_info) );
  EXPECT_EQ(0u, stk::get_doomed_count());
  EXPECT_TRUE( get_reported_doomed_messages().empty() ) << get_reported_doomed_messages();
}

TEST_F(Part_Decomposition_Fixture, TwoSidedSideset_Two_Death_Criteria_Decompose_Throws)
{
  // Without parse-time validation, decomposition refuses the sideset.
  Block_Surface_Connectivity block_surface_info = addTwoSidedSideset();
  EXPECT_ANY_THROW( decompose(make_two_death_decompositions(false), block_surface_info) );
}

TEST_F(Part_Decomposition_Fixture, Subset_Touching_Two_Death_Criteria_Decompose_Throws)
{
  Block_Surface_Connectivity block_surface_info = addSidesetWithSubsetTouchingTwoBlocks();
  EXPECT_ANY_THROW( decompose(make_two_death_decompositions(false), block_surface_info) );
}

TEST_F(Part_Decomposition_Fixture, Subset_Shared_By_Two_Entries_With_Same_Phases)
{
  // setup_phases creates one decomposition entry per block, so a topology-split subset touching
  // block_1 and block_2 is touched by two entries.  With matching phases it is phased once.
  Block_Surface_Connectivity block_surface_info = addSidesetWithSubsetTouchingTwoBlocks();
  const DecompositionPackage decomps = make_two_death_decompositions(true);
  EXPECT_TRUE( validate(decomps, block_surface_info) );
  decompose(decomps, block_surface_info);
  EXPECT_EQ(0u, stk::get_doomed_count()); // the fixture resets the count before each test

  const stk::mesh::Part * subset = findPart("surface_tetra4_tri3_1");
  ASSERT_TRUE( subset != nullptr );
  EXPECT_TRUE( phase_support().is_decomposed(*subset) );
  EXPECT_TRUE( findSuperset("surface_1", subset) != nullptr );

  const stk::mesh::Part * deadSubset = findPart("surface_1_dead_tetra4_tri3");
  ASSERT_TRUE( deadSubset != nullptr );
  EXPECT_TRUE( phase_support().is_conformal(*deadSubset) );
  EXPECT_TRUE( findSuperset("surface_1_dead", deadSubset) != nullptr );
}

TEST_F(Part_Decomposition_Fixture, Shared_Surface_Common_Phase_Propagates_To_Undecomposed_Subset)
{
  // surface_1 touches block_1, block_2 (decomposed by the SAME criterion/phases) and block_3
  // (undecomposed).  The decomposed blocks share one criterion, so the whole sideset, including
  // the block_3 subset, is phased by it (no error).
  Block_Surface_Connectivity block_surface_info = addThreeSidedSideset();
  const DecompositionPackage decomps = make_one_death_decomposition_of_two_blocks();
  EXPECT_TRUE( validate(decomps, block_surface_info) );
  decompose(decomps, block_surface_info);
  EXPECT_EQ(0u, stk::get_doomed_count()); // the fixture resets the count before each test

  for (const std::string subsetName : {"surface_block_1_tri3_1", "surface_block_2_tri3_1", "surface_block_3_tri3_1"})
  {
    const stk::mesh::Part * subset = findPart(subsetName);
    ASSERT_TRUE( subset != nullptr ) << "Missing subset part " << subsetName;
    EXPECT_TRUE( phase_support().is_conformal(*subset) ) << subsetName << " should be phased (conforming)";
    EXPECT_TRUE( findSuperset("surface_1", subset) != nullptr ) << subsetName << " should be parented under surface_1";
  }

  // The superset is the conforming part of the alive phase shared by all of its subsets.
  const stk::mesh::Part * superset = findPart("surface_1");
  ASSERT_TRUE( superset != nullptr );
  EXPECT_TRUE( phase_support().is_conformal(*superset) );
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
  const stk::mesh::Part & block_1_A = phase_support().get_conforming_part(*block_1, get_phase_A());
  EXPECT_EQ( "block_1_A", block_1_A.name() );

  const stk::mesh::Part & block_1_B = phase_support().get_conforming_part(*block_1, get_phase_B());
  EXPECT_EQ( "block_1_B", block_1_B.name() );

  const stk::mesh::Part & block_1_C = phase_support().get_conforming_part(*block_1, get_phase_C());
  EXPECT_EQ( "block_1_C", block_1_C.name() );

  const stk::mesh::Part & block_1_D = phase_support().get_conforming_part(*block_1, get_phase_D());
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

static std::set<std::string> get_names_of_blocks_with_phase(const Phase_Support & phaseSupport, const unsigned phaseIndex)
{
  std::set<std::string> blockNames;
  for (const stk::mesh::Part * block : phaseSupport.get_blocks_decomposed_by_surface({phaseIndex}))
    blockNames.insert(block->name());
  return blockNames;
}

static PhaseInput make_phase_input_with_phases_A_B_AB(const PartnamePhasenameMap & blockPhasesByName)
{
  PhaseInput phaseInput;
  phaseInput.meshPhases = {NamedPhase("A"), NamedPhase("B"), NamedPhase("AB")};
  phaseInput.blockPhasesByName = blockPhasesByName;
  return phaseInput;
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_ExactNames)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_1", {"A", "B"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(0u, stk::get_doomed_count());
  EXPECT_EQ(std::set<std::string>({"block_1"}), get_names_of_blocks_with_phase(phase_support(), 0));
  EXPECT_EQ(std::set<std::string>({"block_1"}), get_names_of_blocks_with_phase(phase_support(), 1));
  EXPECT_TRUE(get_names_of_blocks_with_phase(phase_support(), 2).empty());
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_BlockPatternMatchesAllMatchingElementBlocks)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_*", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(0u, stk::get_doomed_count());
  EXPECT_EQ(std::set<std::string>({"block_1", "block_2", "block_3"}), get_names_of_blocks_with_phase(phase_support(), 0));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_BlockPatternDoesNotMatchSidesets)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"*_1", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(0u, stk::get_doomed_count());
  EXPECT_EQ(std::set<std::string>({"block_1"}), get_names_of_blocks_with_phase(phase_support(), 0));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_PhasePatternMatchesAllMatchingPhases)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_2", {"A*"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(0u, stk::get_doomed_count());
  EXPECT_EQ(std::set<std::string>({"block_2"}), get_names_of_blocks_with_phase(phase_support(), 0));
  EXPECT_TRUE(get_names_of_blocks_with_phase(phase_support(), 1).empty());
  EXPECT_EQ(std::set<std::string>({"block_2"}), get_names_of_blocks_with_phase(phase_support(), 2));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_UnmatchedBlockPattern_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"layer*", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("layer*"));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_MissingBlock_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_4", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("block_4"));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_SidesetNamedAsBlock_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"surface_1", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("surface_1"));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_MissingPhase_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_1", {"C"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("Failed to find phase C"));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_UnmatchedPhasePattern_ReportsDoomedForEachMatchingBlock)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_*", {"C*"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(3u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("matching pattern C*"));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_BlockPatternWithTwoWildcards_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block*_*", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("part name block*_*"));
  EXPECT_TRUE(get_names_of_blocks_with_phase(phase_support(), 0).empty());
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_BlockNameWithQuestionMark_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_?", {"A"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("part name block_?"));
  EXPECT_TRUE(get_names_of_blocks_with_phase(phase_support(), 0).empty());
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_SubdomainPatternWithTwoWildcards_ReportsDoomedOnce)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_*", {"*A*"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("subdomain name *A*"));
}

TEST_F(Part_Decomposition_Fixture, DetermineBlockPhases_SubdomainNameWithQuestionMark_ReportsDoomed)
{
  phase_support().set_phase_input(make_phase_input_with_phases_A_B_AB({{"block_1", {"A", "?"}}}));
  phase_support().determine_block_phases();

  EXPECT_EQ(1u, stk::get_doomed_count());
  EXPECT_NE(std::string::npos, get_reported_doomed_messages().find("subdomain name ?"));
  EXPECT_TRUE(get_names_of_blocks_with_phase(phase_support(), 0).empty());
}

} // namespace krino
