// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_UNIT_PHASE_SUPPORT_H_
#define AKRI_UNIT_PHASE_SUPPORT_H_

#include <gtest/gtest.h>

#include <Akri_Unit_Single_Element_Fixtures.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_AuxMetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <string>

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class MetaData; } }
namespace krino { class Block_Surface_Connectivity; }
namespace krino { class Interface_Name_Generator; }
namespace krino { class NamedPhase; }

namespace krino
{

struct PartNameIs
{
  PartNameIs(const std::string & name) : match_name(name) {}
  bool operator()(const stk::mesh::Part * part) {
    return part->name() == match_name;
  }
  std::string match_name;
};

class Part_Decomposition_Fixture : public ::testing::Test
{
public:
  Part_Decomposition_Fixture();
  virtual ~Part_Decomposition_Fixture();

  // Captures doomed (error) messages for the duration of the test, and resets the doomed count in
  // both SetUp and TearDown, so each test starts at zero and can check exact counts.
  void SetUp() override;
  void TearDown() override;
  static const std::string & get_reported_doomed_messages() { return theReportedDoomedMessages; }

  void performDecomposition(const std::vector<stk::mesh::Part *> & used_blocks, const Block_Surface_Connectivity & input_block_surface_info, bool cdfem_death, int num_ls = 1, bool one_ls_per_phase=false);

  // block_1 and block_2 each decomposed by their own CDFEM death entry (death_by_temp_1 and
  // death_by_temp_2).  With sameCriterion both use Surface_Identifier 0; otherwise block_2 uses 1.
  // Both criteria name their phases "" (alive) and "dead", so a sideset touching both blocks with
  // different criteria has different criteria for phases with the same name (an error).
  DecompositionPackage make_two_death_decompositions(bool sameCriterion);

  // block_1 and block_2 decomposed together by one CDFEM death criterion.
  DecompositionPackage make_one_death_decomposition_of_two_blocks();

  // Sets the block-surface connectivity and runs the parse-time validation, returning its result.
  bool validate(const DecompositionPackage & decomps, const Block_Surface_Connectivity & input_block_surface_info);
  // Sets the block-surface connectivity and decomposes, without validating first.
  void decompose(const DecompositionPackage & decomps, const Block_Surface_Connectivity & input_block_surface_info);

  stk::mesh::MetaData & meta_data() { return fixture.meta_data(); }
  AuxMetaData & aux_meta() { return myAuxMeta; }
  Phase_Support & phase_support() { return Phase_Support::get(meta_data()); }

  Block_Surface_Connectivity addOneSidedSideset();
  Block_Surface_Connectivity addTwoSidedSideset();
  Block_Surface_Connectivity addThreeSidedSideset();
  // surface_1 with one topology-split subset (surface_tetra4_tri3_1) touching both block_1 and block_2.
  Block_Surface_Connectivity addSidesetWithSubsetTouchingTwoBlocks();

protected:
  stk::mesh::Part & get_part(const std::string & partName);
  stk::mesh::PartVector get_parts(const std::vector<std::string> & partNames);
  stk::mesh::Part * findPart(const std::string & part_name);
  stk::mesh::Part * findSuperset(const std::string & superset_name, const stk::mesh::Part * const part);
  void assert_conformal_part_exists(const std::string & conformal_part_name, const std::string & nonconformal_part_name)
  {
    const stk::mesh::Part * conformal_part = findPart(conformal_part_name);
    ASSERT_TRUE( conformal_part != nullptr );
    EXPECT_TRUE( phase_support().is_conformal(*conformal_part) || phase_support().is_interface(*conformal_part) );
    if (phase_support().is_conformal(*conformal_part))
    {
      EXPECT_EQ(nonconformal_part_name, phase_support().get_nonconforming_part(*conformal_part).name() );
    }
  }
  PhaseTag & get_phase_A() { return myPhases[0].tag(); }
  PhaseTag & get_phase_B() { return myPhases[1].tag(); }
  PhaseTag & get_phase_C() { return myPhases[2].tag(); }
  PhaseTag & get_phase_D() { return myPhases[3].tag(); }

private:
  void set_ls_phases(int num_ls, bool one_ls_per_phase=false);
  void set_death_phases();
  // Alive/dead phase pair over the given surface id; GREATER_THAN death: dead -> +1,
  // alive -> -1. Alive name "" (alive conforming part == io part); dead name "dead".
  PhaseVec make_death_phases(unsigned surfaceId) const;

  static void capture_report(const char * message, int type);

  SimpleStkFixture fixture;
  AuxMetaData & myAuxMeta;
  PhaseVec myPhases;
  stk::REH myPreviousReportHandler = nullptr;
  static std::string theReportedDoomedMessages;
};

} // namespace krino

#endif /* AKRI_UNIT_PHASE_SUPPORT_H_ */
