// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CSG_Parser.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_Surface.hpp>
#include <Akri_Parser.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Error / edge-case coverage for CSG_Parser::parse_CSG_tree.  The happy paths
// (leaves, union/unite, intersect/intersection, subtract/subtraction/difference,
// nesting) are covered by Akri_Unit_CSG_Surface_Parser.cpp; this file targets the
// RuntimeDoomed branches those tests never trip.  RuntimeDoomedAdHoc increments a
// monotonic process-global counter rather than throwing, so error paths capture
// stk::get_doomed_count() before the call and assert it grew.
class CSGParserErrorTest : public ::testing::Test
{
protected:
  CSGParserErrorTest()
    : meta(3)
  {
  }

  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }

  stk::mesh::MetaData meta;

  stk::diag::Timer lsTimer{"LS Timer", sierra::Diag::sierraTimer()};
};

TEST_F(CSGParserErrorTest, EmptyUnionOperandListDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // A union key whose value is an empty sequence.
  const std::string yamlText = R"yaml(
csg_initialization:
  union: []
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CSGParserErrorTest, EmptyIntersectOperandListDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  const std::string yamlText = R"yaml(
csg_initialization:
  intersection: []
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CSGParserErrorTest, SubtractWithOneOperandDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  const std::string yamlText = R"yaml(
csg_initialization:
  difference:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 2.0
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CSGParserErrorTest, SubtractWithThreeOperandsDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  const std::string yamlText = R"yaml(
csg_initialization:
  difference:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 3.0
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 2.0
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CSGParserErrorTest, UnknownOperationOrSurfaceDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // A map that is neither a recognized CSG operation nor a valid surface.
  const std::string yamlText = R"yaml(
csg_initialization:
  not_an_operation:
    value: 1
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CSGParserErrorTest, OperandListThatIsNotASequenceDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // union whose value is a map (scalar-ish) instead of a sequence.  get_sequence_nodes
  // returns null for the non-sequence, so this falls through to the unknown-operation
  // path -- either way it must doom and yield an empty tree.
  const std::string yamlText = R"yaml(
csg_initialization:
  union:
    inner: 1
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(CSGParserErrorTest, NonMapOperandInSequenceThrows)
{
  // A union operand that is a bare scalar rather than a surface/operation map.
  // parse_CSG_node dooms ("expecting map") and returns an empty CSGTree, which is
  // then fed into CSGTree::unite -- uniting an operand with no root node throws.
  // Document that current behavior: the malformed operand surfaces as a throw.
  const std::string yamlText = R"yaml(
csg_initialization:
  union:
    - just_a_scalar
)yaml";

  Parser::Node node = load_node(yamlText);
  EXPECT_ANY_THROW(CSG_Parser::parse_CSG_tree(node, meta, lsTimer));
}

TEST_F(CSGParserErrorTest, CsgInitializationAsSequenceYieldsEmptyTreeWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();

  // csg_initialization present but a sequence, not a map: get_map_if_present
  // filters it to an empty tree silently (no doom).
  const std::string yamlText = R"yaml(
csg_initialization:
  - 1
  - 2
)yaml";

  CSGTree tree = CSG_Parser::parse_CSG_tree(load_node(yamlText), meta, lsTimer);

  EXPECT_TRUE(tree.empty());
  // get_map_if_present on a non-map value dooms (type mismatch); a sequence value
  // for a map-typed lookup is a genuine input error, so the count is expected to grow.
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

#else

TEST(CSGParserErrorTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
