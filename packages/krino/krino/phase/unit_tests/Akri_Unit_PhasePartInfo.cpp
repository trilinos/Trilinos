// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_PhasePartInfo.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <memory>

namespace krino {

namespace {
PhaseTag make_phase(const unsigned lsId, const int lsSign)
{
  PhaseTag phase;
  phase.add(Surface_Identifier(lsId), lsSign);
  return phase;
}
}

class PhasePartInfoFixture : public ::testing::Test
{
protected:
  PhasePartInfoFixture()
    : bulk(stk::mesh::MeshBuilder(MPI_COMM_SELF).set_spatial_dimension(3).create()),
      meta(bulk->mesh_meta_data())
  {
  }

  unsigned declare_block(const std::string & name)
  {
    return meta.declare_part_with_topology(name, stk::topology::TET_4).mesh_meta_data_ordinal();
  }

  unsigned declare_part(const std::string & name)
  {
    return meta.declare_part(name).mesh_meta_data_ordinal();
  }

  // An ordinal that is guaranteed to be past the end of any part that has been set up.
  unsigned past_end_ordinal() const { return meta.get_parts().size() + 100; }

  std::unique_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData & meta;
  PhasePartInfo phasePartInfo;
};

TEST_F(PhasePartInfoFixture, DefaultConstructed_IsEmptyAndNothingClassified)
{
  EXPECT_TRUE(phasePartInfo.empty());

  const unsigned ord = declare_block("block_1");
  EXPECT_FALSE(phasePartInfo.is_decomposed(ord));
  EXPECT_FALSE(phasePartInfo.is_original_part(ord));
  EXPECT_FALSE(phasePartInfo.is_nonconforming_part(ord));
  EXPECT_FALSE(phasePartInfo.is_conforming_part(ord));
  EXPECT_FALSE(phasePartInfo.is_interface_part(ord));
  EXPECT_FALSE(phasePartInfo.is_interface_superset_part(ord));
}

TEST_F(PhasePartInfoFixture, PredicatesGuardOrdinalsPastEndOfVector)
{
  const unsigned ord = past_end_ordinal();
  EXPECT_FALSE(phasePartInfo.is_decomposed(ord));
  EXPECT_FALSE(phasePartInfo.is_original_part(ord));
  EXPECT_FALSE(phasePartInfo.is_nonconforming_part(ord));
  EXPECT_FALSE(phasePartInfo.is_conforming_part(ord));
  EXPECT_FALSE(phasePartInfo.is_interface_part(ord));
  EXPECT_FALSE(phasePartInfo.is_interface_superset_part(ord));
}

TEST_F(PhasePartInfoFixture, ConformingPhaseOfUnknownPart_ReturnsEmptyPhase)
{
  const unsigned ord = declare_block("block_1");
  EXPECT_TRUE(phasePartInfo.get_conforming_part_phase(ord).empty());
}

TEST_F(PhasePartInfoFixture, SetupNonconformingPart_ClassifiesOriginalAndNonconforming)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");

  phasePartInfo.setup_nonconforming_part(meta.get_parts().size(), nonconfPart, origPart);

  EXPECT_FALSE(phasePartInfo.empty());

  EXPECT_TRUE(phasePartInfo.is_original_part(origPart));
  EXPECT_TRUE(phasePartInfo.is_nonconforming_part(nonconfPart));
  EXPECT_TRUE(phasePartInfo.is_decomposed(origPart));
  EXPECT_TRUE(phasePartInfo.is_decomposed(nonconfPart));

  // Original resolves to itself, nonconforming resolves to itself.
  EXPECT_EQ(origPart, phasePartInfo.get_original_part(origPart));
  EXPECT_EQ(nonconfPart, phasePartInfo.get_nonconforming_part(nonconfPart));

  // Cross lookups walk the relationship.
  EXPECT_EQ(nonconfPart, phasePartInfo.get_nonconforming_part(origPart));
  EXPECT_EQ(origPart, phasePartInfo.get_original_part(nonconfPart));
}

TEST_F(PhasePartInfoFixture, SetupConformingPart_SetsUpFullChainAndPhase)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");
  const PhaseTag phase = make_phase(0, -1);

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, phase);

  // The nonconforming relationship is set up transitively.
  EXPECT_TRUE(phasePartInfo.is_original_part(origPart));
  EXPECT_TRUE(phasePartInfo.is_nonconforming_part(nonconfPart));
  EXPECT_TRUE(phasePartInfo.is_conforming_part(confPart));
  EXPECT_TRUE(phasePartInfo.is_decomposed(confPart));

  EXPECT_EQ(nonconfPart, phasePartInfo.get_nonconforming_part(confPart));
  EXPECT_EQ(origPart, phasePartInfo.get_original_part(confPart));

  EXPECT_EQ(phase, phasePartInfo.get_conforming_part_phase(confPart));
}

TEST_F(PhasePartInfoFixture, ConformingTakesPrecedenceOverOriginal_DeathCase)
{
  // In death, the conforming part is the same part as the original part; conforming wins.
  const unsigned origAndConfPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const PhaseTag phase = make_phase(0, -1);

  // First register it as an original part...
  phasePartInfo.setup_nonconforming_part(meta.get_parts().size(), nonconfPart, origAndConfPart);
  EXPECT_TRUE(phasePartInfo.is_original_part(origAndConfPart));

  // ...then as a conforming part on the same ordinal.
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), origAndConfPart, nonconfPart, origAndConfPart, phase);

  EXPECT_TRUE(phasePartInfo.is_conforming_part(origAndConfPart));
  EXPECT_FALSE(phasePartInfo.is_original_part(origAndConfPart));
  EXPECT_EQ(phase, phasePartInfo.get_conforming_part_phase(origAndConfPart));
  EXPECT_EQ(nonconfPart, phasePartInfo.get_nonconforming_part(origAndConfPart));
}

TEST_F(PhasePartInfoFixture, OriginalOverConforming_KeepsConformingClassification)
{
  // Reverse order of the death case: conforming already recorded, then original setup arrives.
  const unsigned origAndConfPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const PhaseTag phase = make_phase(0, -1);

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), origAndConfPart, nonconfPart, origAndConfPart, phase);
  EXPECT_TRUE(phasePartInfo.is_conforming_part(origAndConfPart));

  // Setting it up again as a nonconforming/original relationship must not clobber conforming.
  phasePartInfo.setup_nonconforming_part(meta.get_parts().size(), nonconfPart, origAndConfPart);

  EXPECT_TRUE(phasePartInfo.is_conforming_part(origAndConfPart));
  EXPECT_FALSE(phasePartInfo.is_original_part(origAndConfPart));
  EXPECT_EQ(phase, phasePartInfo.get_conforming_part_phase(origAndConfPart));
}

TEST_F(PhasePartInfoFixture, RepeatedSetup_IsIdempotent)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");
  const PhaseTag phase = make_phase(0, -1);

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, phase);
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, phase);
  phasePartInfo.setup_nonconforming_part(meta.get_parts().size(), nonconfPart, origPart);

  EXPECT_TRUE(phasePartInfo.is_original_part(origPart));
  EXPECT_TRUE(phasePartInfo.is_nonconforming_part(nonconfPart));
  EXPECT_TRUE(phasePartInfo.is_conforming_part(confPart));
  EXPECT_EQ(phase, phasePartInfo.get_conforming_part_phase(confPart));
  EXPECT_EQ(nonconfPart, phasePartInfo.get_nonconforming_part(confPart));
  EXPECT_EQ(origPart, phasePartInfo.get_original_part(confPart));
}

TEST_F(PhasePartInfoFixture, SetupInterfacePart_ClassifiesAndChasesPhases)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const unsigned interfacePart = declare_block("surface_block_1_A_B");
  const PhaseTag phaseA = make_phase(0, -1);
  const PhaseTag phaseB = make_phase(0, +1);

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, phaseA);
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, phaseB);
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfacePart, confPartA, confPartB);

  EXPECT_TRUE(phasePartInfo.is_interface_part(interfacePart));
  // An interface part is not one of the decomposed part types.
  EXPECT_FALSE(phasePartInfo.is_decomposed(interfacePart));

  EXPECT_EQ(confPartA, phasePartInfo.get_interface_part_touching_part(interfacePart));
  EXPECT_EQ(confPartB, phasePartInfo.get_interface_part_opposite_part(interfacePart));

  EXPECT_EQ(phaseA, phasePartInfo.get_interface_part_touching_phase(interfacePart));
  EXPECT_EQ(phaseB, phasePartInfo.get_interface_part_opposite_phase(interfacePart));

  // Nonconforming lookup on the interface resolves through the touching conforming part.
  EXPECT_EQ(nonconfPart, phasePartInfo.get_nonconforming_part(interfacePart));
}

TEST_F(PhasePartInfoFixture, SetupInterfaceSupersetPart_ClassifiesAndIsIdempotent)
{
  const unsigned supersetPart = declare_block("surface_interface_superset");

  phasePartInfo.setup_interface_superset_part(meta.get_parts().size(), supersetPart);
  EXPECT_TRUE(phasePartInfo.is_interface_superset_part(supersetPart));
  EXPECT_FALSE(phasePartInfo.is_decomposed(supersetPart));

  // Repeated setup passes the internal require and does not change classification.
  phasePartInfo.setup_interface_superset_part(meta.get_parts().size(), supersetPart);
  EXPECT_TRUE(phasePartInfo.is_interface_superset_part(supersetPart));
}

TEST_F(PhasePartInfoFixture, InterfaceSupersetPart_HasNoOriginalOrNonconformingPart)
{
  const unsigned supersetPart = declare_block("surface_interface_superset");
  phasePartInfo.setup_interface_superset_part(meta.get_parts().size(), supersetPart);

  EXPECT_ANY_THROW(phasePartInfo.get_original_part(supersetPart));
  EXPECT_ANY_THROW(phasePartInfo.get_nonconforming_part(supersetPart));
}

TEST_F(PhasePartInfoFixture, GetNonconformingPartOfUnknownPart_Throws)
{
  const unsigned ord = declare_block("block_1");
  EXPECT_ANY_THROW(phasePartInfo.get_nonconforming_part(ord));
}

TEST_F(PhasePartInfoFixture, ConflictingNonconformingSetup_Throws)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned otherOrigPart = declare_block("block_2");

  phasePartInfo.setup_nonconforming_part(meta.get_parts().size(), nonconfPart, origPart);

  // Re-setup of the same nonconforming part with a different original part must be rejected.
  EXPECT_ANY_THROW(phasePartInfo.setup_nonconforming_part(meta.get_parts().size(), nonconfPart, otherOrigPart));
}

TEST_F(PhasePartInfoFixture, ConflictingConformingPhase_Throws)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, make_phase(0, -1));

  // Re-setup of the same conforming part with a different phase must be rejected.
  EXPECT_ANY_THROW(phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, make_phase(0, +1)));
}

TEST_F(PhasePartInfoFixture, ConflictingInterfaceParts_Throws)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const unsigned confPartC = declare_block("block_1_C");
  const unsigned interfacePart = declare_block("surface_block_1_A_B");

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, make_phase(0, +1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartC, nonconfPart, origPart, make_phase(1, -1));
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfacePart, confPartA, confPartB);

  // Re-setup of the same interface part with different touching/opposite parts must be rejected.
  EXPECT_ANY_THROW(phasePartInfo.setup_interface_part(meta.get_parts().size(), interfacePart, confPartA, confPartC));
}

TEST_F(PhasePartInfoFixture, GetInterfacePartInfoOnNonInterfacePart_Throws)
{
  // This guard was STK_ThrowAssert (compiled out in release); it is now STK_ThrowRequire.
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, make_phase(0, -1));

  EXPECT_ANY_THROW(phasePartInfo.get_interface_part_touching_part(confPart));
}

TEST_F(PhasePartInfoFixture, SetupWithShrinkingNumAllParts_Throws)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned numParts = meta.get_parts().size();

  phasePartInfo.setup_nonconforming_part(numParts, nonconfPart, origPart);

  // A later call that claims fewer total parts would truncate already-classified
  // entries, so it must be rejected rather than silently discarding state.
  EXPECT_ANY_THROW(phasePartInfo.setup_nonconforming_part(numParts - 1, nonconfPart, origPart));
}

//--------------------------------------------------------------------
// Interface part lookup: sort_interface_parts / find_interface_part / get_interface_part
//--------------------------------------------------------------------

TEST_F(PhasePartInfoFixture, FindInterfacePartBeforeSort_Throws)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const unsigned interfacePart = declare_block("surface_block_1_A_B");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, make_phase(0, +1));
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfacePart, confPartA, confPartB);

  unsigned found = 0;
  EXPECT_ANY_THROW(phasePartInfo.find_interface_part(confPartA, confPartB, found));
}

TEST_F(PhasePartInfoFixture, SortAndFindInterfacePart_FindsByTouchingOppositeKey)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const unsigned confPartC = declare_block("block_1_C");
  const unsigned interfaceAB = declare_block("surface_block_1_A_B");
  const unsigned interfaceBC = declare_block("surface_block_1_B_C");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, make_phase(0, +1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartC, nonconfPart, origPart, make_phase(1, -1));

  // Set up out of key order to exercise the sort.
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfaceBC, confPartB, confPartC);
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfaceAB, confPartA, confPartB);

  phasePartInfo.sort_interface_parts();

  unsigned found = 0;
  EXPECT_TRUE(phasePartInfo.find_interface_part(confPartA, confPartB, found));
  EXPECT_EQ(interfaceAB, found);
  EXPECT_TRUE(phasePartInfo.find_interface_part(confPartB, confPartC, found));
  EXPECT_EQ(interfaceBC, found);

  // Both classification and the index mapping survive the sort.
  EXPECT_TRUE(phasePartInfo.is_interface_part(interfaceAB));
  EXPECT_EQ(confPartA, phasePartInfo.get_interface_part_touching_part(interfaceAB));
  EXPECT_EQ(confPartB, phasePartInfo.get_interface_part_opposite_part(interfaceAB));
}

TEST_F(PhasePartInfoFixture, FindInterfacePart_UnknownKeyReturnsFalse)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const unsigned interfaceAB = declare_block("surface_block_1_A_B");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, make_phase(0, +1));
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfaceAB, confPartA, confPartB);
  phasePartInfo.sort_interface_parts();

  unsigned found = 123;
  // Interface parts are directed: (B,A) was never set up.
  EXPECT_FALSE(phasePartInfo.find_interface_part(confPartB, confPartA, found));
  EXPECT_EQ(123u, found); // output left untouched on failure
}

TEST_F(PhasePartInfoFixture, GetInterfacePart_ReturnsPartOrThrowsWhenMissing)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const unsigned interfaceAB = declare_block("surface_block_1_A_B");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, make_phase(0, +1));
  phasePartInfo.setup_interface_part(meta.get_parts().size(), interfaceAB, confPartA, confPartB);
  phasePartInfo.sort_interface_parts();

  EXPECT_EQ(interfaceAB, phasePartInfo.get_interface_part(confPartA, confPartB));
  EXPECT_ANY_THROW(phasePartInfo.get_interface_part(confPartB, confPartA));
}

//--------------------------------------------------------------------
// Conforming part lookup: sort_conforming_parts / find_conforming_part / get_conforming_part
//--------------------------------------------------------------------

TEST_F(PhasePartInfoFixture, FindConformingPartBeforeSort_Throws)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, make_phase(0, -1));

  unsigned found = 0;
  EXPECT_ANY_THROW(phasePartInfo.find_conforming_part(origPart, make_phase(0, -1), found));
}

TEST_F(PhasePartInfoFixture, SortAndFindConformingPart_ExactPhaseMatch)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const PhaseTag phaseA = make_phase(0, -1);
  const PhaseTag phaseB = make_phase(0, +1);

  // Set up in reverse phase order to exercise the sort.
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, phaseB);
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, phaseA);
  phasePartInfo.sort_conforming_parts();

  unsigned found = 0;
  // Lookup works from the original, nonconforming, or a conforming part ordinal.
  EXPECT_TRUE(phasePartInfo.find_conforming_part(origPart, phaseA, found));
  EXPECT_EQ(confPartA, found);
  EXPECT_TRUE(phasePartInfo.find_conforming_part(nonconfPart, phaseB, found));
  EXPECT_EQ(confPartB, found);
  EXPECT_TRUE(phasePartInfo.find_conforming_part(confPartA, phaseB, found));
  EXPECT_EQ(confPartB, found);

  // Classification and the index mapping survive the sort.
  EXPECT_TRUE(phasePartInfo.is_conforming_part(confPartA));
  EXPECT_EQ(phaseA, phasePartInfo.get_conforming_part_phase(confPartA));
}

TEST_F(PhasePartInfoFixture, FindConformingPart_MatchesContainedPhase)
{
  // Conforming parts are keyed by single-level-set phases; a richer query phase
  // that contains one of them must resolve through the containment fallback.
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartNeg = declare_block("block_1_A");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartNeg, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.sort_conforming_parts();

  PhaseTag richQuery = make_phase(0, -1);
  richQuery.add(Surface_Identifier(1), -1); // (LS0,-1),(LS1,-1) contains (LS0,-1)

  unsigned found = 0;
  EXPECT_TRUE(phasePartInfo.find_conforming_part(origPart, richQuery, found));
  EXPECT_EQ(confPartNeg, found);
}

TEST_F(PhasePartInfoFixture, FindConformingPart_UndecomposedPartReturnsFalse)
{
  const unsigned plainPart = declare_block("block_1");
  phasePartInfo.sort_conforming_parts(); // sorting an empty set is fine

  unsigned found = 55;
  EXPECT_FALSE(phasePartInfo.find_conforming_part(plainPart, make_phase(0, -1), found));
  EXPECT_EQ(55u, found); // output left untouched on failure
}

TEST_F(PhasePartInfoFixture, FindConformingPart_NoMatchingPhaseReturnsFalse)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.sort_conforming_parts();

  unsigned found = 0;
  // Decomposed, but no conforming part exists for the opposite side.
  EXPECT_FALSE(phasePartInfo.find_conforming_part(origPart, make_phase(0, +1), found));
}

TEST_F(PhasePartInfoFixture, GetConformingPart_ReturnsPartOrThrowsWhenNotDecomposedByPhase)
{
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPart = declare_block("block_1_A");
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPart, nonconfPart, origPart, make_phase(0, -1));
  phasePartInfo.sort_conforming_parts();

  EXPECT_EQ(confPart, phasePartInfo.get_conforming_part(origPart, make_phase(0, -1)));
  EXPECT_ANY_THROW(phasePartInfo.get_conforming_part(origPart, make_phase(0, +1)));
}

TEST_F(PhasePartInfoFixture, SortConformingParts_DuplicateKeyToDifferentParts_Throws)
{
  // Two distinct conforming part ordinals claim the same (nonconforming, phase) key.
  const unsigned origPart = declare_block("block_1");
  const unsigned nonconfPart = declare_block("block_1_nonconformal");
  const unsigned confPartA = declare_block("block_1_A");
  const unsigned confPartB = declare_block("block_1_B");
  const PhaseTag phase = make_phase(0, -1);

  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartA, nonconfPart, origPart, phase);
  phasePartInfo.setup_conforming_part(meta.get_parts().size(), confPartB, nonconfPart, origPart, phase);

  EXPECT_ANY_THROW(phasePartInfo.sort_conforming_parts());
}

}
