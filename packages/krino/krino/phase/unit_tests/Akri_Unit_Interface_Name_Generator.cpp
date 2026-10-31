// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Interface_Name_Generator.hpp>

#include <string>

namespace krino {

namespace {
const std::string BLOCK = "block_1";
const std::string PHASE = "phaseA";
}

//--------------------------------------------------------------------
// LS_Name_Generator
//--------------------------------------------------------------------

TEST(LS_Name_Generator, InterfaceName_CombinesPartAndPhase)
{
  const LS_Name_Generator gen;
  EXPECT_EQ("surface_block_1_phaseA", gen.interface_name(BLOCK, PHASE));
}

TEST(LS_Name_Generator, SupersetName_UsesPhaseOnly)
{
  const LS_Name_Generator gen;
  EXPECT_EQ("surface_phaseA", gen.interface_superset_name(PHASE));
}

TEST(LS_Name_Generator, SubsetName_MatchesInterfaceName)
{
  // For the LS convention the subset name and the interface name coincide.
  const LS_Name_Generator gen;
  EXPECT_EQ(gen.interface_name(BLOCK, PHASE), gen.interface_subset_name(BLOCK, PHASE));
}

TEST(LS_Name_Generator, UsableThroughBaseInterface)
{
  const LS_Name_Generator concrete;
  const Interface_Name_Generator & gen = concrete;
  EXPECT_EQ("surface_block_1_phaseA", gen.interface_name(BLOCK, PHASE));
  EXPECT_EQ("surface_phaseA", gen.interface_superset_name(PHASE));
  EXPECT_EQ("surface_block_1_phaseA", gen.interface_subset_name(BLOCK, PHASE));
}

TEST(LS_Name_Generator, HasPartPerPhase)
{
  const LS_Name_Generator gen;
  EXPECT_TRUE(gen.has_part_per_phase());
}

//--------------------------------------------------------------------
// Death_Name_Generator
//--------------------------------------------------------------------

TEST(Death_Name_Generator, InterfaceName_UsesDeathSpecInsteadOfPhase)
{
  const Death_Name_Generator gen("death_spec");
  // The phase name is ignored; the death spec name is used in its place.
  EXPECT_EQ("surface_block_1_death_spec", gen.interface_name(BLOCK, PHASE));
}

TEST(Death_Name_Generator, SupersetName_UsesDeathSpecOnly)
{
  const Death_Name_Generator gen("death_spec");
  EXPECT_EQ("surface_death_spec", gen.interface_superset_name(PHASE));
}

TEST(Death_Name_Generator, SubsetName_MatchesInterfaceName)
{
  const Death_Name_Generator gen("death_spec");
  EXPECT_EQ(gen.interface_name(BLOCK, PHASE), gen.interface_subset_name(BLOCK, PHASE));
}

TEST(Death_Name_Generator, PhaseNameDoesNotAffectAnyName)
{
  const Death_Name_Generator gen("death_spec");
  EXPECT_EQ(gen.interface_name(BLOCK, "phaseA"), gen.interface_name(BLOCK, "phaseB"));
  EXPECT_EQ(gen.interface_superset_name("phaseA"), gen.interface_superset_name("phaseB"));
  EXPECT_EQ(gen.interface_subset_name(BLOCK, "phaseA"), gen.interface_subset_name(BLOCK, "phaseB"));
}

TEST(Death_Name_Generator, DifferentSpecNamesProduceDifferentNames)
{
  const Death_Name_Generator genA("specA");
  const Death_Name_Generator genB("specB");
  EXPECT_NE(genA.interface_name(BLOCK, PHASE), genB.interface_name(BLOCK, PHASE));
}

TEST(Death_Name_Generator, UsableThroughBaseInterface)
{
  const Death_Name_Generator concrete("death_spec");
  const Interface_Name_Generator & gen = concrete;
  EXPECT_EQ("surface_block_1_death_spec", gen.interface_name(BLOCK, PHASE));
  EXPECT_EQ("surface_death_spec", gen.interface_superset_name(PHASE));
  EXPECT_EQ("surface_block_1_death_spec", gen.interface_subset_name(BLOCK, PHASE));
}

TEST(Death_Name_Generator, DoesNotHavePartPerPhase)
{
  const Death_Name_Generator gen("death_spec");
  EXPECT_FALSE(gen.has_part_per_phase());
}

} // namespace krino
