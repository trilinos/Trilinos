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
#include <Akri_PhaseUtils.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <map>
#include <memory>
#include <vector>

namespace krino {

namespace {
constexpr int INSIDE = -1;
constexpr int OUTSIDE = +1;

Surface_Identifier ls(const unsigned id) { return Surface_Identifier(id); }

NamedPhase make_phase(const std::string & name, const Surface_Identifier surf, const int sign)
{
  NamedPhase phase(name);
  phase.tag().add(surf, sign);
  return phase;
}
}

//--------------------------------------------------------------------
// get_indices_of_phases_that_depend_on_surface
//--------------------------------------------------------------------

TEST(GetIndicesOfPhasesThatDependOnSurface, OneLSPerPhase_ReturnsEveryIndex)
{
  PhaseVec meshPhases;
  meshPhases.push_back(make_phase("A", ls(0), INSIDE));
  meshPhases.push_back(make_phase("B", ls(1), INSIDE));
  meshPhases.push_back(make_phase("C", ls(2), INSIDE));

  // With one level set per phase, every level set impacts every phase,
  // so the queried surface identifier is irrelevant.
  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(true, meshPhases, ls(7));

  const std::vector<unsigned> expected = {0u, 1u, 2u};
  EXPECT_EQ(expected, indices);
}

TEST(GetIndicesOfPhasesThatDependOnSurface, OneLSPerPhase_EmptyMesh_ReturnsEmpty)
{
  const PhaseVec meshPhases;
  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(true, meshPhases, ls(0));
  EXPECT_TRUE(indices.empty());
}

TEST(GetIndicesOfPhasesThatDependOnSurface, NotOneLSPerPhase_ReturnsOnlyDependentPhases)
{
  PhaseVec meshPhases;
  meshPhases.push_back(make_phase("neg0", ls(0), INSIDE));  // depends on surface 0
  meshPhases.push_back(make_phase("neg1", ls(1), INSIDE));  // depends on surface 1
  meshPhases.push_back(make_phase("pos0", ls(0), OUTSIDE)); // depends on surface 0

  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(false, meshPhases, ls(0));

  // Phases 0 and 2 reference surface 0 (either sign); phase 1 does not.
  const std::vector<unsigned> expected = {0u, 2u};
  EXPECT_EQ(expected, indices);
}

TEST(GetIndicesOfPhasesThatDependOnSurface, NotOneLSPerPhase_PositiveSideCounts)
{
  PhaseVec meshPhases;
  meshPhases.push_back(make_phase("pos5", ls(5), OUTSIDE));

  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(false, meshPhases, ls(5));

  const std::vector<unsigned> expected = {0u};
  EXPECT_EQ(expected, indices);
}

TEST(GetIndicesOfPhasesThatDependOnSurface, NotOneLSPerPhase_NoDependence_ReturnsEmpty)
{
  PhaseVec meshPhases;
  meshPhases.push_back(make_phase("neg0", ls(0), INSIDE));
  meshPhases.push_back(make_phase("neg1", ls(1), INSIDE));

  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(false, meshPhases, ls(9));

  EXPECT_TRUE(indices.empty());
}

TEST(GetIndicesOfPhasesThatDependOnSurface, NotOneLSPerPhase_MultiSurfacePhase)
{
  // A phase whose tag mentions several surfaces should be reported for each
  // of those surfaces, but not for an unrelated one.
  NamedPhase multi("multi");
  multi.tag().add(ls(0), INSIDE);
  multi.tag().add(ls(1), OUTSIDE);

  PhaseVec meshPhases;
  meshPhases.push_back(multi);

  EXPECT_EQ((std::vector<unsigned>{0u}),
      get_indices_of_phases_that_depend_on_surface(false, meshPhases, ls(0)));
  EXPECT_EQ((std::vector<unsigned>{0u}),
      get_indices_of_phases_that_depend_on_surface(false, meshPhases, ls(1)));
  EXPECT_TRUE(
      get_indices_of_phases_that_depend_on_surface(false, meshPhases, ls(2)).empty());
}

//--------------------------------------------------------------------
// build_phase_conversion_map_for_element_and_side_parts
//--------------------------------------------------------------------

namespace {
PhaseTag phase_tag(const unsigned lsId, const int sign)
{
  PhaseTag phase;
  phase.add(Surface_Identifier(lsId), sign);
  return phase;
}
}

class PhaseConversionFixture : public ::testing::Test
{
protected:
  PhaseConversionFixture()
    : bulk(stk::mesh::MeshBuilder(MPI_COMM_SELF).set_spatial_dimension(3).create()),
      meta(bulk->mesh_meta_data())
  {
  }

  // Element-rank (block) part; passes the ELEMENT_RANK filter.
  unsigned declare_block(const std::string & name)
  {
    return meta.declare_part_with_topology(name, stk::topology::TET_4).mesh_meta_data_ordinal();
  }

  // Side-rank (surface) part in 3D; passes the side_rank filter and is used for
  // interface parts.
  unsigned declare_surface(const std::string & name)
  {
    return meta.declare_part_with_topology(name, stk::topology::TRI_3).mesh_meta_data_ordinal();
  }

  // Rank-less part; should be skipped by the element/side rank filter.
  unsigned declare_unranked(const std::string & name)
  {
    return meta.declare_part(name).mesh_meta_data_ordinal();
  }

  void declare_superset(const unsigned supersetOrd, const unsigned subsetOrd)
  {
    meta.declare_part_subset(meta.get_part(supersetOrd), meta.get_part(subsetOrd));
  }

  unsigned numParts() const { return meta.get_parts().size(); }

  std::map<int,int> build_map(const PhaseTag & fromPhase, const PhaseTag & toPhase)
  {
    return build_phase_conversion_map_for_element_and_side_parts(meta, phasePartInfo, fromPhase, toPhase);
  }

  std::unique_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData & meta;
  PhasePartInfo phasePartInfo;
};

TEST_F(PhaseConversionFixture, ConformingPart_MappedToToPhaseConformingPart)
{
  const unsigned orig = declare_block("block_1");
  const unsigned nonconf = declare_block("block_1_nonconformal");
  const unsigned confA = declare_block("block_1_A");
  const unsigned confB = declare_block("block_1_B");
  const PhaseTag phaseA = phase_tag(0, INSIDE);
  const PhaseTag phaseB = phase_tag(0, OUTSIDE);

  phasePartInfo.setup_conforming_part(numParts(), confA, nonconf, orig, phaseA);
  phasePartInfo.setup_conforming_part(numParts(), confB, nonconf, orig, phaseB);
  phasePartInfo.sort_conforming_parts();

  const std::map<int,int> partMap = build_map(phaseA, phaseB);

  // Only the conforming part whose phase matches fromPhase is converted, and it
  // maps to the toPhase conforming part on the same nonconforming part.
  ASSERT_EQ(1u, partMap.count(confA));
  EXPECT_EQ(static_cast<int>(confB), partMap.at(confA));
  EXPECT_EQ(0u, partMap.count(confB));
}

TEST_F(PhaseConversionFixture, InterfacePart_TouchingSideConverted_RemapsToNewInterface)
{
  const unsigned orig = declare_block("block_1");
  const unsigned nonconf = declare_block("block_1_nonconformal");
  const unsigned confA = declare_block("block_1_A");
  const unsigned confB = declare_block("block_1_B");
  const unsigned confC = declare_block("block_1_C");
  const unsigned ifaceAB = declare_surface("surface_block_1_A_B");
  const unsigned ifaceCB = declare_surface("surface_block_1_C_B");
  const PhaseTag phaseA = phase_tag(0, INSIDE);
  const PhaseTag phaseB = phase_tag(0, OUTSIDE);
  const PhaseTag phaseC = phase_tag(1, INSIDE);

  phasePartInfo.setup_conforming_part(numParts(), confA, nonconf, orig, phaseA);
  phasePartInfo.setup_conforming_part(numParts(), confB, nonconf, orig, phaseB);
  phasePartInfo.setup_conforming_part(numParts(), confC, nonconf, orig, phaseC);
  phasePartInfo.setup_interface_part(numParts(), ifaceAB, confA, confB);
  phasePartInfo.setup_interface_part(numParts(), ifaceCB, confC, confB);
  phasePartInfo.sort_conforming_parts();
  phasePartInfo.sort_interface_parts();

  const std::map<int,int> partMap = build_map(phaseA, phaseC);

  // Converting the touching side A->C turns interface (A,B) into interface (C,B).
  ASSERT_EQ(1u, partMap.count(ifaceAB));
  EXPECT_EQ(static_cast<int>(ifaceCB), partMap.at(ifaceAB));

  // The conforming block A is also remapped to its C counterpart.
  ASSERT_EQ(1u, partMap.count(confA));
  EXPECT_EQ(static_cast<int>(confC), partMap.at(confA));
}

TEST_F(PhaseConversionFixture, InterfacePart_OppositeSideConverted_RemapsToNewInterface)
{
  const unsigned orig = declare_block("block_1");
  const unsigned nonconf = declare_block("block_1_nonconformal");
  const unsigned confA = declare_block("block_1_A");
  const unsigned confB = declare_block("block_1_B");
  const unsigned confC = declare_block("block_1_C");
  const unsigned ifaceAB = declare_surface("surface_block_1_A_B");
  const unsigned ifaceAC = declare_surface("surface_block_1_A_C");
  const PhaseTag phaseA = phase_tag(0, INSIDE);
  const PhaseTag phaseB = phase_tag(0, OUTSIDE);
  const PhaseTag phaseC = phase_tag(1, INSIDE);

  phasePartInfo.setup_conforming_part(numParts(), confA, nonconf, orig, phaseA);
  phasePartInfo.setup_conforming_part(numParts(), confB, nonconf, orig, phaseB);
  phasePartInfo.setup_conforming_part(numParts(), confC, nonconf, orig, phaseC);
  phasePartInfo.setup_interface_part(numParts(), ifaceAB, confA, confB);
  phasePartInfo.setup_interface_part(numParts(), ifaceAC, confA, confC);
  phasePartInfo.sort_conforming_parts();
  phasePartInfo.sort_interface_parts();

  // fromPhase = B matches the opposite side of (A,B); convert B->C gives (A,C).
  const std::map<int,int> partMap = build_map(phaseB, phaseC);

  ASSERT_EQ(1u, partMap.count(ifaceAB));
  EXPECT_EQ(static_cast<int>(ifaceAC), partMap.at(ifaceAB));
}

TEST_F(PhaseConversionFixture, InterfacePart_CollapsesToRemoval_AlsoRemovesInterfaceSuperset)
{
  const unsigned orig = declare_block("block_1");
  const unsigned nonconf = declare_block("block_1_nonconformal");
  const unsigned confA = declare_block("block_1_A");
  const unsigned confB = declare_block("block_1_B");
  const unsigned ifaceAB = declare_surface("surface_block_1_A_B");
  const unsigned supersetIface = declare_surface("surface_block_1");
  declare_superset(supersetIface, ifaceAB);
  const PhaseTag phaseA = phase_tag(0, INSIDE);
  const PhaseTag phaseB = phase_tag(0, OUTSIDE);

  phasePartInfo.setup_conforming_part(numParts(), confA, nonconf, orig, phaseA);
  phasePartInfo.setup_conforming_part(numParts(), confB, nonconf, orig, phaseB);
  phasePartInfo.setup_interface_part(numParts(), ifaceAB, confA, confB);
  phasePartInfo.setup_interface_superset_part(numParts(), supersetIface);
  phasePartInfo.sort_conforming_parts();
  phasePartInfo.sort_interface_parts();

  // Converting touching side A->B makes both sides B, so the interface collapses
  // and is flagged for removal (-1), along with its interface superset part.
  const std::map<int,int> partMap = build_map(phaseA, phaseB);

  ASSERT_EQ(1u, partMap.count(ifaceAB));
  EXPECT_EQ(-1, partMap.at(ifaceAB));
  ASSERT_EQ(1u, partMap.count(supersetIface));
  EXPECT_EQ(-1, partMap.at(supersetIface));
}

TEST_F(PhaseConversionFixture, InterfacePartMatchingNeitherSide_IsAbsentAndSupersetKept)
{
  // An interface part whose touching/opposite phases both differ from fromPhase
  // must not appear in the map at all, and its interface superset must not be
  // flagged for removal.  (Regression: operator[] used to default-insert a 0
  // entry for such parts and then treat it as a removal target.)
  const unsigned orig = declare_block("block_1");
  const unsigned nonconf = declare_block("block_1_nonconformal");
  const unsigned confA = declare_block("block_1_A");
  const unsigned confB = declare_block("block_1_B");
  const unsigned confC = declare_block("block_1_C");
  const unsigned ifaceBC = declare_surface("surface_block_1_B_C");
  const unsigned supersetIface = declare_surface("surface_block_1");
  declare_superset(supersetIface, ifaceBC);
  const PhaseTag phaseA = phase_tag(0, INSIDE);
  const PhaseTag phaseB = phase_tag(0, OUTSIDE);
  const PhaseTag phaseC = phase_tag(1, INSIDE);

  phasePartInfo.setup_conforming_part(numParts(), confA, nonconf, orig, phaseA);
  phasePartInfo.setup_conforming_part(numParts(), confB, nonconf, orig, phaseB);
  phasePartInfo.setup_conforming_part(numParts(), confC, nonconf, orig, phaseC);
  phasePartInfo.setup_interface_part(numParts(), ifaceBC, confB, confC);
  phasePartInfo.setup_interface_superset_part(numParts(), supersetIface);
  phasePartInfo.sort_conforming_parts();
  phasePartInfo.sort_interface_parts();

  // Convert phase A, but interface (B,C) touches neither A side.
  const std::map<int,int> partMap = build_map(phaseA, phaseB);

  EXPECT_EQ(0u, partMap.count(ifaceBC));
  EXPECT_EQ(0u, partMap.count(supersetIface));
}

TEST_F(PhaseConversionFixture, NonElementNonSidePartsAreIgnored)
{
  const unsigned orig = declare_block("block_1");
  const unsigned nonconf = declare_block("block_1_nonconformal");
  const unsigned confA = declare_block("block_1_A");
  const unsigned confB = declare_block("block_1_B");
  const PhaseTag phaseA = phase_tag(0, INSIDE);
  const PhaseTag phaseB = phase_tag(0, OUTSIDE);
  phasePartInfo.setup_conforming_part(numParts(), confA, nonconf, orig, phaseA);
  phasePartInfo.setup_conforming_part(numParts(), confB, nonconf, orig, phaseB);

  // A rank-less conforming part on a separate nonconforming part, with phase == fromPhase
  // but no toPhase counterpart. If the rank filter did not exclude it, the toPhase
  // conforming-part lookup would throw; the fact that it does not proves it is skipped.
  const unsigned orig2 = declare_block("block_2");
  const unsigned nonconf2 = declare_block("block_2_nonconformal");
  const unsigned confUnranked = declare_unranked("block_2_A_unranked");
  phasePartInfo.setup_conforming_part(numParts(), confUnranked, nonconf2, orig2, phaseA);
  phasePartInfo.sort_conforming_parts();

  const std::map<int,int> partMap = build_map(phaseA, phaseB);

  EXPECT_EQ(0u, partMap.count(confUnranked));
  EXPECT_EQ(1u, partMap.count(confA)); // the element-rank part is still converted
}

} // namespace krino
