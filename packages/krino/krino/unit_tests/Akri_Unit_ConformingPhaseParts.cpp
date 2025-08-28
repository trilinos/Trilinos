// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_ConformingPhaseParts.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Unit_Part_Decomposition_Fixture.hpp>

namespace krino {

class ConformingPhasePartFixture : public Part_Decomposition_Fixture
{
public:
  ConformingPhasePartFixture() {}

  CDFEM_Support & cdfem_support() { return CDFEM_Support::get(meta_data()); }

  bool have_part(const stk::mesh::PartVector & parts, stk::mesh::Part * part)
  {
    return std::find(parts.begin(), parts.end(), part) != parts.end();
  }

  stk::mesh::Part & active_part() { return aux_meta().active_part(); }
  stk::mesh::Part & block_boundary_part() { return aux_meta().block_boundary_part(); }
  stk::mesh::Part & exposed_boundary_part() { return aux_meta().exposed_boundary_part(); }
  stk::mesh::Part & child_part() { return cdfem_support().get_child_part(); }
  stk::mesh::Part & parent_part() { return cdfem_support().get_parent_part(); }
  stk::topology side_topology() { return stk::topology::TRIANGLE_3; }
  stk::topology element_topology() { return stk::topology::TETRAHEDRON_4; }
  stk::mesh::Part & side_topology_part() { return meta_data().get_topology_root_part(side_topology()); }
  stk::mesh::Part & element_topology_part() { return meta_data().get_topology_root_part(element_topology()); }

  void expect_matching_parts(const stk::mesh::PartVector & goldParts, const stk::mesh::PartVector & parts, const std::string & desc)
  {
    for (auto * part : goldParts)
    {
      EXPECT_TRUE(have_part(parts, part)) << "In " << desc << ", Gold part not found in actual parts: " << part->name();
    }
    for (auto * part : parts)
    {
      EXPECT_TRUE(have_part(goldParts, part)) << "In " << desc << ", Actual part not found in gold parts: " << part->name();
    }
  }

  void expect_part_changes_for_side(const std::vector<std::string> & existingPartNames, const bool doesSideHaveActiveElements, const std::vector<std::string> & goldAddPartNames, const std::vector<std::string> & goldRemovePartNames)
  {
    stk::mesh::PartVector existingParts = get_parts(existingPartNames);
    stk::mesh::PartVector goldAddParts = get_parts(goldAddPartNames);
    stk::mesh::PartVector goldRemoveParts = get_parts(goldRemovePartNames);
    stk::mesh::PartVector addParts;
    stk::mesh::PartVector removeParts;
    determine_part_changes_for_side(phase_support(),
        meta_data().side_rank(),
        existingParts,
        doesSideHaveActiveElements,
        block_boundary_part(),
        child_part(),
        parent_part(),
        active_part(),
        addParts,
        removeParts);
    expect_matching_parts(goldAddParts, addParts, "add parts");
    expect_matching_parts(goldRemoveParts, removeParts, "remove parts");
  }

  void expect_conforming_part_changes(const std::vector<std::string> & existingPartNames, const stk::mesh::EntityRank entityRank, const PhaseTag & phase, const std::vector<std::string> & goldAddPartNames, const std::vector<std::string> & goldRemovePartNames)
  {
    stk::mesh::PartVector existingParts = get_parts(existingPartNames);
    stk::mesh::PartVector goldAddParts = get_parts(goldAddPartNames);
    stk::mesh::PartVector goldRemoveParts = get_parts(goldRemovePartNames);
    stk::mesh::PartVector addParts;
    stk::mesh::PartVector removeParts;
    append_conforming_part_changes_for_current_parts(phase_support(),
        existingParts,
        entityRank,
        phase,
        addParts,
        removeParts);
    expect_matching_parts(goldAddParts, addParts, "add parts");
    expect_matching_parts(goldRemoveParts, removeParts, "remove parts");
  }

  void expect_nonconforming_part_changes(const std::vector<std::string> & existingPartNames, const stk::mesh::EntityRank entityRank, const std::vector<std::string> & goldAddPartNames, const std::vector<std::string> & goldRemovePartNames)
  {
    stk::mesh::PartVector existingParts = get_parts(existingPartNames);
    stk::mesh::PartVector goldAddParts = get_parts(goldAddPartNames);
    stk::mesh::PartVector goldRemoveParts = get_parts(goldRemovePartNames);
    stk::mesh::PartVector addParts;
    stk::mesh::PartVector removeParts;
    append_nonconforming_part_changes_for_current_parts(phase_support(),
        existingParts,
        entityRank,
        child_part(),
        parent_part(),
        active_part(),
        addParts,
        removeParts);
    expect_matching_parts(goldAddParts, addParts, "add parts");
    expect_matching_parts(goldRemoveParts, removeParts, "remove parts");
  }

  void expect_child_creation_parts(const stk::topology topology, const std::vector<std::string> & parentPartNames, const stk::mesh::EntityRank /*entityRank*/, const PhaseTag & phase, const std::vector<std::string> & goldChildPartNames)
  {
    stk::mesh::PartVector parentParts = get_parts(parentPartNames);
    stk::mesh::PartVector goldChildParts = get_parts(goldChildPartNames);
    stk::mesh::PartVector attributeParts;
    stk::mesh::PartVector childParts;
    stk::mesh::insert(attributeParts, active_part());
    stk::mesh::insert(attributeParts, exposed_boundary_part());
    stk::mesh::insert(attributeParts, block_boundary_part());

    determine_child_conforming_parts(meta_data(),
        phase_support(),
        topology,
        parentParts,
        attributeParts,
        child_part(),
        active_part(),
        phase,
        childParts);
    expect_matching_parts(goldChildParts, childParts, "child parts");
  }
};

TEST_F(ConformingPhasePartFixture, OneLS_TwoSidedSideset_side_part_changes)
{
  Block_Surface_Connectivity blockSurfInfo = addTwoSidedSideset();
  performDecomposition({meta_data().get_part("block_1")}, blockSurfInfo, false, 2);

  expect_part_changes_for_side(
      {"block_1_A", "block_1_B"}, true,
      {"surface_block_1_A_B", "surface_block_1_B_A", active_part().name(), block_boundary_part().name()},
      {} );

  expect_part_changes_for_side(
      {"block_1_A", "surface_A_B", "surface_block_1_A_B", "surface_block_1_B_A"}, true,
      {active_part().name()},
      {"surface_A_B", "surface_block_1_A_B", "surface_block_1_B_A", block_boundary_part().name()} );

  expect_part_changes_for_side(
      {"block_1_nonconformal", "surface_1"}, false,
      {"surface_1_nonconformal"},
      {"surface_1", active_part().name(), block_boundary_part().name()} );

  expect_part_changes_for_side(
      {"block_2", "surface_1"}, true,
      {active_part().name()},
      {block_boundary_part().name()} );

  expect_part_changes_for_side(
      {"block_2", "surface_1"}, false,
      {},
      {active_part().name(), block_boundary_part().name()} );

  expect_part_changes_for_side(
        {"block_1_B", "block_1_C", "surface_A_B", "surface_block_1_A_B", "surface_block_1_B_A"}, true,
        {"surface_block_1_B_C", "surface_block_1_C_B", active_part().name(), block_boundary_part().name()},
        {"surface_A_B", "surface_block_1_A_B", "surface_block_1_B_A"} );
}

TEST_F(ConformingPhasePartFixture, OneLS_TwoSidedSideset_conforming_part_changes)
{
  Block_Surface_Connectivity blockSurfInfo = addTwoSidedSideset();
  performDecomposition({meta_data().get_part("block_1")}, blockSurfInfo, false);

  expect_conforming_part_changes(
      {"block_1", "block_2", "surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"}, meta_data().side_rank(), get_phase_A(),
      {"surface_1_A", "surface_1_A_block_1_A_tri3", "surface_1_A_block_2_tri3"},
      {"surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"} );

  expect_conforming_part_changes(
      {"block_1", "block_2", "surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"}, meta_data().side_rank(), get_phase_B(),
      {"surface_1_B", "surface_1_B_block_1_B_tri3", "surface_1_B_block_2_tri3"},
      {"surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"} );

  expect_conforming_part_changes(
      {"block_1"}, stk::topology::ELEMENT_RANK, get_phase_A(),
      {"block_1_A"},
      {"block_1"} );

  expect_conforming_part_changes(
      {"block_1"}, stk::topology::ELEMENT_RANK, get_phase_B(),
      {"block_1_B"},
      {"block_1"} );
}

TEST_F(ConformingPhasePartFixture, OneLS_TwoSidedSideset_nonconforming_part_changes)
{
  Block_Surface_Connectivity blockSurfInfo = addTwoSidedSideset();
  performDecomposition({meta_data().get_part("block_1")}, blockSurfInfo, false);

  expect_nonconforming_part_changes(
      {"block_1", "block_2", "surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"}, meta_data().side_rank(),
      {"surface_1_nonconformal", "surface_1_nonconformal_block_1_nonconformal_tri3", "surface_1_nonconformal_block_2_tri3"},
      {"surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1", active_part().name()} );

  expect_nonconforming_part_changes(
      {"block_1_A", "block_2", "surface_1_A", "surface_1_A_block_1_A_tri3", "surface_1_A_block_2_tri3"}, meta_data().side_rank(),
      {"surface_1_nonconformal", "surface_1_nonconformal_block_1_nonconformal_tri3", "surface_1_nonconformal_block_2_tri3"},
      {"surface_1_A", "surface_1_A_block_1_A_tri3", "surface_1_A_block_2_tri3", active_part().name()} );

  expect_nonconforming_part_changes(
      {"block_1_A"}, stk::topology::ELEMENT_RANK,
      {"block_1_nonconformal", parent_part().name()},
      {"block_1_A", active_part().name(), child_part().name()} );

  expect_nonconforming_part_changes(
      {"block_1"}, stk::topology::ELEMENT_RANK,
      {"block_1_nonconformal", parent_part().name()},
      {"block_1", active_part().name(), child_part().name()} );
}

TEST_F(ConformingPhasePartFixture, OneLS_TwoSidedSideset_child_creation_parts)
{
  Block_Surface_Connectivity blockSurfInfo = addTwoSidedSideset();
  performDecomposition({meta_data().get_part("block_1")}, blockSurfInfo, false);

  expect_child_creation_parts(side_topology(),
      {"block_1", "block_2", "surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"}, meta_data().side_rank(), get_phase_A(),
      {"surface_1_A", "surface_1_A_block_1_A_tri3", "surface_1_A_block_2_tri3", side_topology_part().name(), active_part().name()} );

  expect_child_creation_parts(side_topology(),
      {"block_1", "block_2", "surface_1", "surface_block_1_tri3_1", "surface_block_2_tri3_1"}, meta_data().side_rank(), get_phase_B(),
      {"surface_1_B", "surface_1_B_block_1_B_tri3", "surface_1_B_block_2_tri3", side_topology_part().name(), active_part().name()} );

  expect_child_creation_parts(element_topology(),
      {"block_1"}, stk::topology::ELEMENT_RANK, get_phase_A(),
      {"block_1_A", element_topology_part().name(), child_part().name(), active_part().name()} );

  expect_child_creation_parts(element_topology(),
      {"block_1"}, stk::topology::ELEMENT_RANK, get_phase_B(),
      {"block_1_B", element_topology_part().name(), child_part().name(), active_part().name()} );
}

}






