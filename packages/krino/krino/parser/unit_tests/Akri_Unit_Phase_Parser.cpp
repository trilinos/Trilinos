// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Phase_Parser.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>
#include <vector>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Unit tests for Phase_Parser::parse.  The parser reads "subdomains" and "parts"
// from a finite-element-model node and returns them as a PhaseInput, including the
// surface names its level set identifiers refer to.  It needs no MetaData / mesh / MPI,
// and keeps no state between calls.  RuntimeDoomedAdHoc
// increments a monotonic process-global counter rather than throwing; error paths
// capture stk::get_doomed_count() before the call and assert it grew.
class PhaseParserTest : public ::testing::Test
{
protected:
  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }

  // Unique model name per test so the name-keyed registries stay isolated.
  std::string model_name() const
  {
    const ::testing::TestInfo * info = ::testing::UnitTest::GetInstance()->current_test_info();
    return std::string("PhaseParserTest.") + info->name();
  }
};

TEST_F(PhaseParserTest, EmptyNodeYieldsNoPhases)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
unrelated: 1
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, fem);
  EXPECT_TRUE(phaseInput.meshPhases.empty());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, ParsesSubdomainsWithLevelSetRegions)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: solid
    level_set_regions:
      - negative_level_set: ls1
  - name: void
    level_set_regions:
      - positive_level_set: ls1
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, fem);
  SurfaceNames surfaceNames = phaseInput.surfaceNames;
  const PhaseVec & phases = phaseInput.meshPhases;

  ASSERT_EQ(2u, phases.size());
  EXPECT_EQ("solid", phases[0].name());
  EXPECT_EQ("void", phases[1].name());

  const Surface_Identifier ls1 = surfaceNames.get_identifier("ls1");
  EXPECT_TRUE(phases[0].tag().contain(ls1, -1));
  EXPECT_TRUE(phases[1].tag().contain(ls1, +1));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, MissingSubdomainNameDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - level_set_regions:
      - negative_level_set: ls1
)yaml");

  Phase_Parser::parse(fem_node, fem);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, ParsesSmallestLevelSetSyntax)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: phaseA
    smallest_level_set: lsA
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, fem);
  SurfaceNames surfaceNames = phaseInput.surfaceNames;
  const PhaseVec & phases = phaseInput.meshPhases;

  ASSERT_EQ(1u, phases.size());
  const Surface_Identifier lsA = surfaceNames.get_identifier("lsA");
  EXPECT_TRUE(phases[0].tag().contain(lsA, -1));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, CombiningSmallestLevelSetWithLevelSetRegionsDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: phaseA
    level_set_regions:
      - negative_level_set: ls1
    smallest_level_set: lsA
)yaml");

  Phase_Parser::parse(fem_node, fem);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, ParsesBoundingSurfaceSyntax)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: inside
    bounding_surface: surf1
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, fem);
  SurfaceNames surfaceNames = phaseInput.surfaceNames;
  const PhaseVec & phases = phaseInput.meshPhases;

  ASSERT_EQ(1u, phases.size());
  const Surface_Identifier surf1 = surfaceNames.get_identifier("surf1");
  EXPECT_TRUE(phases[0].tag().contain(surf1, -1));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, CombiningBoundingSurfaceWithLevelSetRegionsDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: phaseA
    level_set_regions:
      - negative_level_set: ls1
    bounding_surface: surf1
)yaml");

  Phase_Parser::parse(fem_node, fem);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, MoreThanOneElsewherePhaseDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  // Two subdomains with no levelset / bounding surface each yield an empty tag,
  // so both are "elsewhere" phases -- which is ambiguous and must doom.
  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: phaseA
  - name: phaseB
)yaml");

  Phase_Parser::parse(fem_node, fem);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, ParsesPartsWithSubdomains)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: solid
    level_set_regions:
      - negative_level_set: ls1
  - name: void
    level_set_regions:
      - positive_level_set: ls1
parts:
  - name: block_1
    subdomains:
      - solid
      - void
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, fem);
  const auto & blockPhases = phaseInput.blockPhasesByName;

  auto it = blockPhases.find("block_1");
  ASSERT_NE(blockPhases.end(), it);
  ASSERT_EQ(2u, it->second.size());
  EXPECT_EQ("solid", it->second[0]);
  EXPECT_EQ("void", it->second[1]);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, MissingPartNameDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const std::string fem = model_name();

  Parser::Node fem_node = load_node(R"yaml(
parts:
  - subdomains:
      - solid
)yaml");

  Phase_Parser::parse(fem_node, fem);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// between_consecutive_level_sets
// -----------------------------------------------------------------------------

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_IsRecordedNotExpanded)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: before
    level_set_regions:
      - negative_level_set: other
  - between_consecutive_level_sets:
      subdomain_prefix: layer
      level_set_prefix: LS
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  ASSERT_EQ(1u, phaseInput.meshPhases.size());
  ASSERT_EQ(1u, phaseInput.generatedSubdomains.size());
  EXPECT_EQ(GeneratedSubdomains::BETWEEN_CONSECUTIVE_LEVEL_SETS, phaseInput.generatedSubdomains[0].type);
  EXPECT_FALSE(phaseInput.oneLevelSetPerPhase);
  EXPECT_EQ("layer", phaseInput.generatedSubdomains[0].subdomainPrefix);
  EXPECT_EQ("LS", phaseInput.generatedSubdomains[0].levelSetPrefix);
  EXPECT_EQ(1u, phaseInput.generatedSubdomains[0].phaseIndex);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_ExpandsToLayersBetweenConsecutiveLevelSets)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: before
    level_set_regions:
      - negative_level_set: other
  - between_consecutive_level_sets:
      subdomain_prefix: layer
      level_set_prefix: LS
  - name: after
    level_set_regions:
      - positive_level_set: other
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  // Padded, out of order and mixed case, with a non-matching name.
  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS10", "other", "ls02", "LS00", "LS01", "LS03", "LS04", "LS05", "LS06", "LS07", "LS08", "LS09"});

  const PhaseVec & phases = phaseInput.meshPhases;
  ASSERT_EQ(12u, phases.size());
  EXPECT_EQ("before", phases[0].name());
  EXPECT_EQ("after", phases[11].name());
  EXPECT_TRUE(phaseInput.generatedSubdomains.empty());

  SurfaceNames surfaceNames = phaseInput.surfaceNames;
  const auto padded = [](const unsigned i) { return (i < 10) ? "0" + std::to_string(i) : std::to_string(i); };
  for (unsigned k=1; k<=10; ++k)
  {
    const NamedPhase & layer = phases[k];
    EXPECT_EQ("layer" + padded(k), layer.name());
    EXPECT_TRUE(layer.tag().contain(surfaceNames.get_identifier("LS" + padded(k-1)), +1)) << layer.name();
    EXPECT_TRUE(layer.tag().contain(surfaceNames.get_identifier("LS" + padded(k)), -1)) << layer.name();
  }
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_LayerIsNamedForUpperLevelSet)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - between_consecutive_level_sets:
      subdomain_prefix: layer
      level_set_prefix: LS
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS_top", "LS_bottom"});

  ASSERT_EQ(1u, phaseInput.meshPhases.size());
  EXPECT_EQ("layer_top", phaseInput.meshPhases[0].name());
  SurfaceNames surfaceNames = phaseInput.surfaceNames;
  EXPECT_TRUE(phaseInput.meshPhases[0].tag().contain(surfaceNames.get_identifier("LS_bottom"), +1));
  EXPECT_TRUE(phaseInput.meshPhases[0].tag().contain(surfaceNames.get_identifier("LS_top"), -1));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_ElsewherePhaseExcludesLayers)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - between_consecutive_level_sets:
      subdomain_prefix: layer
      level_set_prefix: LS
  - name: outside
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  ASSERT_EQ(1u, phaseInput.meshPhases.size());
  EXPECT_TRUE(phaseInput.meshPhases[0].tag().empty()); // elsewhere phase not set until expansion

  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS0", "LS1", "LS2"});

  ASSERT_EQ(3u, phaseInput.meshPhases.size());
  const NamedPhase & outside = phaseInput.meshPhases[2];
  EXPECT_EQ("outside", outside.name());
  SurfaceNames surfaceNames = phaseInput.surfaceNames;
  // Opposite sides of the layers' sides.
  EXPECT_TRUE(outside.tag().contain(surfaceNames.get_identifier("LS0"), -1));
  EXPECT_TRUE(outside.tag().contain(surfaceNames.get_identifier("LS2"), +1));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_FewerThanTwoMatchesDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - between_consecutive_level_sets:
      subdomain_prefix: layer
      level_set_prefix: LS
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS0", "other"});
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_DuplicateGeneratedNameDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: Layer1
    level_set_regions:
      - negative_level_set: other
  - between_consecutive_level_sets:
      subdomain_prefix: layer
      level_set_prefix: LS
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS0", "LS1"}); // generates layer1
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_MissingPrefixDooms)
{
  for (const std::string & yamlText : {
      std::string("subdomains:\n  - between_consecutive_level_sets:\n      level_set_prefix: LS\n"),
      std::string("subdomains:\n  - between_consecutive_level_sets:\n      subdomain_prefix: layer\n"),
      std::string("subdomains:\n  - between_consecutive_level_sets:\n      subdomain_prefix: ''\n      level_set_prefix: LS\n")})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    Phase_Parser::parse(load_node(yamlText), model_name());
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

TEST_F(PhaseParserTest, BetweenConsecutiveLevelSets_CombinedWithOtherKeywordsDooms)
{
  for (const std::string & otherLine : {
      std::string("    name: layers\n"),
      std::string("    smallest_level_set: ls1\n"),
      std::string("    bounding_surface: bs1\n"),
      std::string("    level_set_regions:\n      - negative_level_set: ls1\n")})
  {
    const std::string yamlText = "subdomains:\n  - between_consecutive_level_sets:\n      subdomain_prefix: layer\n      level_set_prefix: LS\n" + otherLine;
    const unsigned oldDoomedCount = stk::get_doomed_count();
    Phase_Parser::parse(load_node(yamlText), model_name());
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

TEST_F(PhaseParserTest, PerSmallestLevelSet_IsRecordedAndMakesModelOneLevelSetPerPhase)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - per_smallest_level_set:
      subdomain_prefix: p
      level_set_prefix: LS
)yaml");

  const PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  EXPECT_TRUE(phaseInput.meshPhases.empty());
  ASSERT_EQ(1u, phaseInput.generatedSubdomains.size());
  EXPECT_EQ(GeneratedSubdomains::PER_SMALLEST_LEVEL_SET, phaseInput.generatedSubdomains[0].type);
  EXPECT_EQ("p", phaseInput.generatedSubdomains[0].subdomainPrefix);
  EXPECT_EQ("LS", phaseInput.generatedSubdomains[0].levelSetPrefix);
  EXPECT_TRUE(phaseInput.oneLevelSetPerPhase);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, PerSmallestLevelSet_MatchesSmallestLevelSetSubdomains)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node shorthandNode = load_node(R"yaml(
subdomains:
  - name: e
  - per_smallest_level_set:
      subdomain_prefix: p
      level_set_prefix: LS
)yaml");
  Parser::Node explicitNode = load_node(R"yaml(
subdomains:
  - name: e
  - name: p001
    smallest_level_set: LS001
  - name: p007
    smallest_level_set: LS007
  - name: p100
    smallest_level_set: LS100
)yaml");

  PhaseInput shorthand = Phase_Parser::parse(shorthandNode, model_name());
  // Out of order and mixed case, with non-matching names.
  Phase_Parser::expand_generated_subdomains(shorthand, {"ls100", "other", "LS001", "LS007"});
  PhaseInput expected = Phase_Parser::parse(explicitNode, model_name() + ".explicit");

  EXPECT_TRUE(shorthand.oneLevelSetPerPhase);
  EXPECT_EQ(expected.oneLevelSetPerPhase, shorthand.oneLevelSetPerPhase);
  ASSERT_EQ(expected.meshPhases.size(), shorthand.meshPhases.size());

  // Compare by surface name, since each input numbers its surfaces in the order it first saw them.
  for (size_t i=0; i<expected.meshPhases.size(); ++i)
  {
    EXPECT_EQ(expected.meshPhases[i].name(), shorthand.meshPhases[i].name());
    for (const auto ls : {"LS001", "LS007", "LS100"})
      for (const int sign : {-1, +1})
        EXPECT_EQ(expected.meshPhases[i].tag().contain(expected.surfaceNames.get_identifier(ls), sign),
                  shorthand.meshPhases[i].tag().contain(shorthand.surfaceNames.get_identifier(ls), sign))
            << expected.meshPhases[i].name() << " " << ls << " " << sign;
  }
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, PerSmallestLevelSet_OneMatchIsEnough)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - per_smallest_level_set:
      subdomain_prefix: p
      level_set_prefix: LS
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS1", "other"});
  ASSERT_EQ(1u, phaseInput.meshPhases.size());
  EXPECT_EQ("p1", phaseInput.meshPhases[0].name());
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, PerSmallestLevelSet_NoMatchesDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - per_smallest_level_set:
      subdomain_prefix: p
      level_set_prefix: LS
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  Phase_Parser::expand_generated_subdomains(phaseInput, {"other"});
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(PhaseParserTest, PerSmallestLevelSet_MissingPrefixOrCombinedWithOtherKeywordsDooms)
{
  for (const std::string & yamlText : {
      std::string("subdomains:\n  - per_smallest_level_set:\n      level_set_prefix: LS\n"),
      std::string("subdomains:\n  - per_smallest_level_set:\n      subdomain_prefix: p\n"),
      std::string("subdomains:\n  - per_smallest_level_set:\n      subdomain_prefix: p\n      level_set_prefix: LS\n    name: particles\n"),
      std::string("subdomains:\n  - per_smallest_level_set:\n      subdomain_prefix: p\n      level_set_prefix: LS\n    smallest_level_set: LS1\n")})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    Phase_Parser::parse(load_node(yamlText), model_name());
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

// -----------------------------------------------------------------------------
// Flat layout of between_consecutive_level_sets and per_smallest_level_set: options beside the keyword
// -----------------------------------------------------------------------------

TEST_F(PhaseParserTest, GeneratedSubdomains_FlatLayoutExpandsLikeNestedLayout)
{
  for (const std::string keyword : {"between_consecutive_level_sets", "per_smallest_level_set"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();

    const std::string nestedYaml = "subdomains:\n  - name: e\n  - " + keyword + ":\n      subdomain_prefix: p\n      level_set_prefix: LS\n";
    const std::string flatYaml = "subdomains:\n  - name: e\n  - " + keyword + ":\n    subdomain_prefix: p\n    level_set_prefix: LS\n";
    const std::vector<std::string> surfaceNames{"LS2", "other", "LS0", "LS1"};

    PhaseInput nested = Phase_Parser::parse(load_node(nestedYaml), model_name() + "." + keyword + ".nested");
    PhaseInput flat = Phase_Parser::parse(load_node(flatYaml), model_name() + "." + keyword + ".flat");

    ASSERT_EQ(1u, flat.generatedSubdomains.size()) << keyword;
    EXPECT_EQ(nested.generatedSubdomains[0].type, flat.generatedSubdomains[0].type) << keyword;
    EXPECT_EQ("p", flat.generatedSubdomains[0].subdomainPrefix) << keyword;
    EXPECT_EQ("LS", flat.generatedSubdomains[0].levelSetPrefix) << keyword;
    EXPECT_EQ(nested.generatedSubdomains[0].phaseIndex, flat.generatedSubdomains[0].phaseIndex) << keyword;
    EXPECT_EQ(nested.oneLevelSetPerPhase, flat.oneLevelSetPerPhase) << keyword;

    Phase_Parser::expand_generated_subdomains(nested, surfaceNames);
    Phase_Parser::expand_generated_subdomains(flat, surfaceNames);

    ASSERT_EQ(nested.meshPhases.size(), flat.meshPhases.size()) << keyword;
    for (size_t i=0; i<nested.meshPhases.size(); ++i)
    {
      EXPECT_EQ(nested.meshPhases[i].name(), flat.meshPhases[i].name()) << keyword;
      for (const auto ls : {"LS0", "LS1", "LS2"})
        for (const int sign : {-1, +1})
          EXPECT_EQ(nested.meshPhases[i].tag().contain(nested.surfaceNames.get_identifier(ls), sign),
                    flat.meshPhases[i].tag().contain(flat.surfaceNames.get_identifier(ls), sign))
              << keyword << " " << nested.meshPhases[i].name() << " " << ls << " " << sign;
    }
    EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount) << keyword;
  }
}

TEST_F(PhaseParserTest, GeneratedSubdomains_FlatLayoutMissingPrefixOrCombinedWithOtherKeywordsDooms)
{
  for (const std::string keyword : {"between_consecutive_level_sets", "per_smallest_level_set"})
  {
    for (const std::string & entryLines : {
        std::string("    level_set_prefix: LS\n"),
        std::string("    subdomain_prefix: p\n"),
        std::string("    subdomain_prefix: p\n    level_set_prefix: LS\n    name: particles\n"),
        std::string("    subdomain_prefix: p\n    level_set_prefix: LS\n    smallest_level_set: LS1\n")})
    {
      const std::string yamlText = "subdomains:\n  - " + keyword + ":\n" + entryLines;
      const unsigned oldDoomedCount = stk::get_doomed_count();
      Phase_Parser::parse(load_node(yamlText), model_name());
      EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
    }
  }
}

TEST_F(PhaseParserTest, GeneratedSubdomains_ScalarOrSequenceValueDooms)
{
  for (const std::string keyword : {"between_consecutive_level_sets", "per_smallest_level_set"})
  {
    for (const std::string & value : {std::string(" yes"), std::string(" [p, LS]")})
    {
      const std::string yamlText = "subdomains:\n  - " + keyword + ":" + value + "\n    subdomain_prefix: p\n    level_set_prefix: LS\n";
      const unsigned oldDoomedCount = stk::get_doomed_count();
      Phase_Parser::parse(load_node(yamlText), model_name());
      EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
    }
  }
}

TEST_F(PhaseParserTest, ExpandWithoutLayeredSubdomains_LeavesPhasesUnchanged)
{
  Parser::Node fem_node = load_node(R"yaml(
subdomains:
  - name: solid
    level_set_regions:
      - negative_level_set: ls1
  - name: void
)yaml");

  PhaseInput phaseInput = Phase_Parser::parse(fem_node, model_name());
  const PhaseVec before = phaseInput.meshPhases;
  Phase_Parser::expand_generated_subdomains(phaseInput, {"LS0", "LS1"});

  ASSERT_EQ(before.size(), phaseInput.meshPhases.size());
  for (size_t i=0; i<before.size(); ++i)
  {
    EXPECT_EQ(before[i].name(), phaseInput.meshPhases[i].name());
    EXPECT_TRUE(before[i].tag() == phaseInput.meshPhases[i].tag());
  }
}

#else

TEST(PhaseParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
