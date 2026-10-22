// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Phase_Support.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <memory>
#include <vector>

#include <Akri_PhaseUtils.hpp>
namespace krino {

namespace {
Surface_Identifier ls(const unsigned id) { return Surface_Identifier(id); }

NamedPhase make_named_phase(const std::string & name, const std::vector<std::pair<unsigned,int>> & sides)
{
  PhaseTag tag;
  for (auto && side : sides)
    tag.add(ls(side.first), side.second);
  return NamedPhase(name, tag);
}
}

//--------------------------------------------------------------------
// Phase_Support::get_indices_of_phases_that_depend_on_surface  (static, pure)
//--------------------------------------------------------------------

TEST(Get_Indices_Of_Phases_That_Depend_On_Surface, OneLevelSetPerPhase_ReturnsAllPhases)
{
  // With one level set per phase, every level set is considered to impact every phase,
  // so all phase indices are returned regardless of the surface.
  PhaseVec phases;
  phases.push_back(make_named_phase("LS1", {{0,-1}}));
  phases.push_back(make_named_phase("LS2", {{1,-1}}));
  phases.push_back(make_named_phase("LS3", {{2,-1}}));

  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(true, phases, ls(0));

  EXPECT_EQ((std::vector<unsigned>{0u, 1u, 2u}), indices);
}

TEST(Get_Indices_Of_Phases_That_Depend_On_Surface, OneLevelSetPerPhase_IgnoresSurfaceIdentity)
{
  // The returned set does not depend on which surface is queried in the oneLSPerPhase case.
  PhaseVec phases;
  phases.push_back(make_named_phase("LS1", {{0,-1}}));
  phases.push_back(make_named_phase("LS2", {{1,-1}}));

  const std::vector<unsigned> forLs0 =
      get_indices_of_phases_that_depend_on_surface(true, phases, ls(0));
  const std::vector<unsigned> forUnrelated =
      get_indices_of_phases_that_depend_on_surface(true, phases, ls(99));

  EXPECT_EQ((std::vector<unsigned>{0u, 1u}), forLs0);
  EXPECT_EQ(forLs0, forUnrelated);
}

TEST(Get_Indices_Of_Phases_That_Depend_On_Surface, ExplicitDependence_MatchesEitherSign)
{
  // Not oneLSPerPhase: only phases whose tag contains (surf,-1) or (surf,+1) are returned.
  PhaseVec phases;
  phases.push_back(make_named_phase("A", {{0,-1}}));   // depends on surface 0 (negative side)
  phases.push_back(make_named_phase("B", {{0,+1}}));   // depends on surface 0 (positive side)
  phases.push_back(make_named_phase("C", {{1,-1}}));   // depends on surface 1 only

  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(false, phases, ls(0));

  EXPECT_EQ((std::vector<unsigned>{0u, 1u}), indices);
}

TEST(Get_Indices_Of_Phases_That_Depend_On_Surface, ExplicitDependence_NoneDependOnSurface)
{
  PhaseVec phases;
  phases.push_back(make_named_phase("A", {{0,-1}}));
  phases.push_back(make_named_phase("B", {{1,+1}}));

  const std::vector<unsigned> indices =
      get_indices_of_phases_that_depend_on_surface(false, phases, ls(2));

  EXPECT_TRUE(indices.empty());
}

TEST(Get_Indices_Of_Phases_That_Depend_On_Surface, ExplicitDependence_MultiSidedPhaseIsMatched)
{
  // A phase that contains several sides depends on the surface if any side references it.
  PhaseVec phases;
  phases.push_back(make_named_phase("A", {{0,-1},{1,-1}}));
  phases.push_back(make_named_phase("B", {{2,-1},{3,+1}}));

  const std::vector<unsigned> forLs1 =
      get_indices_of_phases_that_depend_on_surface(false, phases, ls(1));
  const std::vector<unsigned> forLs3 =
      get_indices_of_phases_that_depend_on_surface(false, phases, ls(3));

  EXPECT_EQ((std::vector<unsigned>{0u}), forLs1);
  EXPECT_EQ((std::vector<unsigned>{1u}), forLs3);
}

//--------------------------------------------------------------------
// DecompositionPackage
//--------------------------------------------------------------------

class DecompositionPackageFixture : public ::testing::Test
{
protected:
  DecompositionPackageFixture()
    : bulk(stk::mesh::MeshBuilder(MPI_COMM_SELF).set_spatial_dimension(3).create()),
      meta(bulk->mesh_meta_data())
  {
  }

  stk::mesh::Part & declare_block(const std::string & name)
  {
    return meta.declare_part_with_topology(name, stk::topology::TET_4);
  }

  PhaseVec make_phases() const
  {
    PhaseVec phases;
    phases.push_back(make_named_phase("A", {{0,-1}}));
    phases.push_back(make_named_phase("B", {{0,+1}}));
    return phases;
  }

  std::unique_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData & meta;
};

TEST_F(DecompositionPackageFixture, DefaultConstructed_IsEmpty)
{
  const DecompositionPackage decomps;
  EXPECT_EQ(0u, decomps.size());
}

TEST_F(DecompositionPackageFixture, AddLevelsetDecomposition_IncrementsSize)
{
  DecompositionPackage decomps;
  decomps.add_levelset_decomposition({&declare_block("block_1")}, make_phases());

  ASSERT_EQ(1u, decomps.size());
  EXPECT_EQ(1u, decomps.get_blocks_to_decompose(0).size());
  EXPECT_EQ(2u, decomps.get_phases(0).size());
  // LS decomposition uses a generator with a distinct part per phase.
  EXPECT_TRUE(decomps.get_interface_name_generator(0).has_part_per_phase());
}

TEST_F(DecompositionPackageFixture, AddDeathDecomposition_UsesGeneratorWithoutPartPerPhase)
{
  DecompositionPackage decomps;
  decomps.add_death_decomposition({&declare_block("block_1")}, make_phases(), "death_spec");

  ASSERT_EQ(1u, decomps.size());
  EXPECT_FALSE(decomps.get_interface_name_generator(0).has_part_per_phase());
}

TEST_F(DecompositionPackageFixture, EmptyBlocks_IsSkipped)
{
  DecompositionPackage decomps;
  decomps.add_levelset_decomposition({}, make_phases());
  EXPECT_EQ(0u, decomps.size());
}

TEST_F(DecompositionPackageFixture, EmptyPhases_IsSkipped)
{
  DecompositionPackage decomps;
  decomps.add_levelset_decomposition({&declare_block("block_1")}, PhaseVec{});
  EXPECT_EQ(0u, decomps.size());
}

TEST_F(DecompositionPackageFixture, EmptyDeathDecomposition_IsSkipped)
{
  DecompositionPackage decomps;
  decomps.add_death_decomposition({}, make_phases(), "death_spec");
  decomps.add_death_decomposition({&declare_block("block_1")}, PhaseVec{}, "death_spec");
  EXPECT_EQ(0u, decomps.size());
}

} // namespace krino
