// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Options_Parser.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Parser.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Unit tests for CDFEM_Options_Parser::parse.  The parser reads the optional
// cdfem_options map from a region node and mutates CDFEM_Support (and
// RefinementSupport / Phase_Support) attached to the MetaData.  It takes no
// timer and needs no mesh/MPI, so a bare MetaData fixture drives every branch.
// RuntimeDoomedAdHoc increments a monotonic process-global counter rather than
// throwing; error paths capture stk::get_doomed_count() before the call and
// assert it grew, happy paths assert it is unchanged.
class CDFEMOptionsParserTest : public ::testing::Test
{
protected:
  CDFEMOptionsParserTest()
    : meta(3)
  {
    // CDFEM_Support / Phase_Support construction reaches into AuxMetaData::get,
    // which requires AuxMetaData to have been registered on the MetaData first.
    AuxMetaData::create(meta);
  }

  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }

  stk::mesh::MetaData meta;
};

// -----------------------------------------------------------------------------
// No cdfem_options map -> no-op
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, NoCdfemOptionsMapIsNoOp)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
some_other_option: 1
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// cdfem_edge_degeneracy_handling
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, ParsesEdgeDegeneracyHandlingSnapToNode)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  cdfem_edge_degeneracy_handling: SNAP_TO_NODE
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(SNAP_TO_NODE, CDFEM_Support::get(meta).get_cdfem_edge_degeneracy_handling());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, ParsesEdgeDegeneracyHandlingSnapToInterface)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // Case-insensitive value; also triggers the default-edge-tol-when-snapping side effect.
  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  cdfem_edge_degeneracy_handling: snap_to_interface_when_quality_allows
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE,
            CDFEM_Support::get(meta).get_cdfem_edge_degeneracy_handling());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, InvalidEdgeDegeneracyHandlingDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  cdfem_edge_degeneracy_handling: NOT_A_HANDLING
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// max_edge_snap
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, ParsesMaxEdgeSnapInRange)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  max_edge_snap: 0.25
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_DOUBLE_EQ(0.25, CDFEM_Support::get(meta).get_max_edge_snap());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, MaxEdgeSnapBelowRangeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  max_edge_snap: -0.1
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, MaxEdgeSnapAboveRangeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  max_edge_snap: 1.5
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// cdfem_simplex_generation_method
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, ParsesSimplexGenerationMethod)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  cdfem_simplex_generation_method: CUT_QUADS_BY_LARGEST_ANGLE
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(CUT_QUADS_BY_LARGEST_ANGLE, CDFEM_Support::get(meta).get_simplex_generation_method());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, InvalidSimplexGenerationMethodDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  cdfem_simplex_generation_method: NOT_A_METHOD
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// Scalar setters that never doom
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, ParsesSharpFeatureAngleWhenSnapping)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  sharp_feature_angle_when_snapping: 90.0
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_DOUBLE_EQ(90.0, CDFEM_Support::get(meta).get_snapping_sharp_feature_angle_in_degrees());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, ParsesUseFacetsToPerformDecomposition)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  use_facets_to_perform_decomposition: true
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_TRUE(CDFEM_Support::get(meta).use_facets_instead_of_levelset_fields());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, ParsesNumberOfInitialDecompositionCycles)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  number_of_initial_decomposition_cycles: 3
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(3, CDFEM_Support::get(meta).get_num_initial_decomposition_cycles());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, ParsesPostCdfemRefinementLevels)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  post_cdfem_refinement_levels: 2
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(2, CDFEM_Support::get(meta).get_post_cdfem_refinement_levels());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// cdfem_edge_tolerance -- below-minimum clamps (logs) rather than dooms
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, TinyEdgeToleranceClampsWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // A tolerance below the 1e-12 minimum is clamped (with a log message), not doomed.
  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  cdfem_edge_tolerance: 1.0e-20
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// small_cluster_cleanup
// -----------------------------------------------------------------------------

TEST_F(CDFEMOptionsParserTest, SmallClusterCleanupNeitherMapNorSequenceDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // A bare scalar for small_cluster_cleanup is neither a map nor a sequence.
  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  small_cluster_cleanup: not_a_container
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CDFEMOptionsParserTest, SmallClusterCleanupSequenceMissingPhaseNamesThrows)
{
  // Sequence form, but the entry omits the required from/to/count keys.  Each
  // missing key dooms, yet the parser still emplaces an ("","",0) cleanup tuple
  // and calls Phase_Support::set_small_cluster_cleanup, which looks up the empty
  // phase name and STK_ThrowRequire's that the phase exists.  With no phases
  // registered that requirement fails, so this path surfaces as a throw.
  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  small_cluster_cleanup:
    - unrelated_key: 1
)yaml");

  EXPECT_ANY_THROW(CDFEM_Options_Parser::parse(region_node, meta));
}

TEST_F(CDFEMOptionsParserTest, SmallClusterCleanupMapClusterWrongOperandCountDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // Map form: each cluster in convert_all_but_the_largest_cluster_of_phase_from_to
  // must be a [from_phase, to_phase] pair.  A single-element list dooms.
  Parser::Node region_node = load_node(R"yaml(
cdfem_options:
  small_cluster_cleanup:
    convert_all_but_the_largest_cluster_of_phase_from_to:
      - [only_one_phase]
)yaml");

  CDFEM_Options_Parser::parse(region_node, meta);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

#else

TEST(CDFEMOptionsParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
