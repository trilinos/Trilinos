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
#include <Akri_Unit_LogRedirecter.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_PhaseTag.hpp>

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

  void performDecomposition(const std::vector<stk::mesh::Part *> & used_blocks, const Block_Surface_Connectivity & input_block_surface_info, bool cdfem_death, int num_ls = 1, bool one_ls_per_phase=false);

  stk::mesh::MetaData & get_meta_data();

  Block_Surface_Connectivity addOneSidedSideset();
  Block_Surface_Connectivity addTwoSidedSideset();

protected:
  stk::mesh::Part * findPart(const std::string & part_name);
  stk::mesh::Part * findSuperset(const std::string & superset_name, const stk::mesh::Part * const part);
  void assert_conformal_part_exists(const std::string & conformal_part_name, const std::string & nonconformal_part_name)
  {
    const stk::mesh::Part * conformal_part = findPart(conformal_part_name);
    Phase_Support & phase_support = Phase_Support::get(get_meta_data());
    ASSERT_TRUE( conformal_part != NULL );
    EXPECT_TRUE( phase_support.is_conformal(conformal_part) );
    EXPECT_EQ(nonconformal_part_name, phase_support.find_nonconformal_part(*conformal_part)->name() );
  }

private:
  static PhaseVec ls_phases(int num_ls, bool one_ls_per_phase=false);
  static PhaseVec death_phases();

  static const Interface_Name_Generator & ls_name_generator();
  static const Interface_Name_Generator & death_name_generator();

  SimpleStkFixture fixture;
  LogRedirecter log;
};

} // namespace krino

#endif /* AKRI_UNIT_PHASE_SUPPORT_H_ */
