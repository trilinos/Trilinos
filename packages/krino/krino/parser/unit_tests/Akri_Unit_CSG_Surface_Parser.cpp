/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Akri_CSG_Parser.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_Surface.hpp>
#include <Akri_Parser.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <memory>
#include <string>
#include <vector>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

class CSGParserTest : public ::testing::Test
{
protected:
  CSGParserTest()
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

TEST_F(CSGParserTest, MissingCsgInitializationReturnsEmptyTree)
{
  const std::string yamlText = R"yaml(
some_other_block:
  value: 1
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  EXPECT_TRUE(tree.empty());
}

TEST_F(CSGParserTest, ParsesSingleSphereAsLeaf)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  sphere:
  center: [0.0, 0.0, 0.0]
  radius: 1.0
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Leaf);
  EXPECT_NE(root->surface(), nullptr);

  std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();
  ASSERT_EQ(surfaces.size(), 1u);
  EXPECT_EQ(surfaces[0], root->surface());
}

TEST_F(CSGParserTest, ParsesSingleCuboidUsingExistingFlatPrimitiveSyntax)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  cuboid:
  center: [0.0, 0.0, 0.0]
  dimensions: [90.0, 60.0, 18.0]
  rotation: [0.0, 0.0, 0.0]
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Leaf);
  EXPECT_NE(root->surface(), nullptr);

  EXPECT_EQ(tree.gather_surfaces().size(), 1u);
}

TEST_F(CSGParserTest, ParsesUnionOfTwoSurfaces)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  union:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0

    - sphere:
      center: [1.0, 0.0, 0.0]
      radius: 2.0
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Union);

  const std::vector<std::shared_ptr<CSGTreeNode>> & operands =
      root->vector_operands();

  ASSERT_EQ(operands.size(), 2u);
  ASSERT_NE(operands[0], nullptr);
  ASSERT_NE(operands[1], nullptr);

  EXPECT_EQ(operands[0]->node_operator(), CSGOperation::Leaf);
  EXPECT_EQ(operands[1]->node_operator(), CSGOperation::Leaf);

  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesUniteAlias)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  unite:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0

    - sphere:
      center: [1.0, 0.0, 0.0]
      radius: 2.0
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());
  ASSERT_NE(tree.root(), nullptr);

  EXPECT_EQ(tree.root()->node_operator(), CSGOperation::Union);
  EXPECT_EQ(tree.root()->vector_operands().size(), 2u);
  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesNaryUnion)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  union:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0

    - sphere:
      center: [1.0, 0.0, 0.0]
      radius: 2.0

    - cylinder:
      p1: [0.0, 0.0, -1.0]
      p2: [0.0, 0.0,  1.0]
      radius: 0.25
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Union);
  EXPECT_EQ(root->vector_operands().size(), 3u);

  EXPECT_EQ(tree.gather_surfaces().size(), 3u);
}

TEST_F(CSGParserTest, ParsesIntersectionOfTwoSurfaces)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  intersection:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0

    - cuboid:
      center: [0.0, 0.0, 0.0]
      dimensions: [1.0, 1.0, 1.0]
      rotation: [0.0, 0.0, 0.0]
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Intersection);

  const auto & operands = root->vector_operands();
  ASSERT_EQ(operands.size(), 2u);

  EXPECT_EQ(operands[0]->node_operator(), CSGOperation::Leaf);
  EXPECT_EQ(operands[1]->node_operator(), CSGOperation::Leaf);

  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesIntersectAlias)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  intersect:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0

    - sphere:
      center: [0.5, 0.0, 0.0]
      radius: 1.0
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());
  ASSERT_NE(tree.root(), nullptr);

  EXPECT_EQ(tree.root()->node_operator(), CSGOperation::Intersection);
  EXPECT_EQ(tree.root()->vector_operands().size(), 2u);
  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesNaryIntersection)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  intersection:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 2.0

    - cuboid:
      center: [0.0, 0.0, 0.0]
      dimensions: [2.0, 2.0, 2.0]
      rotation: [0.0, 0.0, 0.0]

    - cylinder:
      p1: [0.0, 0.0, -2.0]
      p2: [0.0, 0.0,  2.0]
      radius: 0.75
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Intersection);
  EXPECT_EQ(root->vector_operands().size(), 3u);

  EXPECT_EQ(tree.gather_surfaces().size(), 3u);
}

TEST_F(CSGParserTest, ParsesDifferenceOfTwoSurfaces)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  difference:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 2.0

    - cylinder:
      p1: [0.0, 0.0, -3.0]
      p2: [0.0, 0.0,  3.0]
      radius: 0.5
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Difference);

  ASSERT_NE(root->left_operand(), nullptr);
  ASSERT_NE(root->right_operand(), nullptr);

  EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Leaf);
  EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Leaf);

  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesSubtractAlias)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  subtract:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 2.0

    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());
  ASSERT_NE(tree.root(), nullptr);

  EXPECT_EQ(tree.root()->node_operator(), CSGOperation::Difference);
  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesSubtractionAlias)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  subtraction:
    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 2.0

    - sphere:
      center: [0.0, 0.0, 0.0]
      radius: 1.0
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());
  ASSERT_NE(tree.root(), nullptr);

  EXPECT_EQ(tree.root()->node_operator(), CSGOperation::Difference);
  EXPECT_EQ(tree.gather_surfaces().size(), 2u);
}

TEST_F(CSGParserTest, ParsesNestedCsgTree)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  difference:
    - union:
        - sphere:
          center: [0.0, 0.0, 0.0]
          radius: 2.0

        - sphere:
          center: [1.0, 0.0, 0.0]
          radius: 2.0

        - cylinder:
          p1: [0.0, 0.0, -2.0]
          p2: [0.0, 0.0,  2.0]
          radius: 0.5

    - intersection:
        - cuboid:
          center: [0.0, 0.0, 0.0]
          dimensions: [1.0, 1.0, 1.0]
          rotation: [0.0, 0.0, 0.0]

        - sphere:
          center: [0.0, 0.0, 0.0]
          radius: 0.75
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  /*
   * Expected:
   *
   *   Difference
   *   /        \
   * Union    Intersection
   * / | \      /       \
   * S S C    Cuboid    S
   */

  EXPECT_EQ(root->node_operator(), CSGOperation::Difference);

  std::shared_ptr<CSGTreeNode> left = root->left_operand();
  std::shared_ptr<CSGTreeNode> right = root->right_operand();

  ASSERT_NE(left, nullptr);
  ASSERT_NE(right, nullptr);

  EXPECT_EQ(left->node_operator(), CSGOperation::Union);
  EXPECT_EQ(right->node_operator(), CSGOperation::Intersection);

  EXPECT_EQ(left->vector_operands().size(), 3u);
  EXPECT_EQ(right->vector_operands().size(), 2u);

  EXPECT_EQ(tree.gather_surfaces().size(), 5u);
}

TEST_F(CSGParserTest, ParsesLargeExistingStyleStructure)
{
  const std::string yamlText = R"yaml(
csg_initialization:
  difference:
    - union:
        - difference:
            - cuboid:
              center: [0.0, 0.0, 0.0]
              dimensions: [90.0, 60.0, 18.0]
              rotation: [0.0, 0.0, 0.0]

            - union:
                - cuboid:
                  center: [0.0, 0.0, 0.0]
                  dimensions: [58.0, 34.0, 24.0]
                  rotation: [0.0, 0.0, 0.0]

                - cylinder:
                  p1: [-36.0, -40.0, 0.0]
                  p2: [-36.0, 40.0, 0.0]
                  radius: 8.0

                - cylinder:
                  p1: [36.0, -40.0, 0.0]
                  p2: [36.0, 40.0, 0.0]
                  radius: 8.0

        - sphere:
          center: [0.0, 0.0, 0.0]
          radius: 14.0

        - cylinder:
          p1: [-39.0, 0.0, 0.0]
          p2: [39.0, 0.0, 0.0]
          radius: 4.0

    - union:
        - cylinder:
          p1: [-50.0, 0.0, 0.0]
          p2: [50.0, 0.0, 0.0]
          radius: 3.0

        - cylinder:
          p1: [0.0, -40.0, 0.0]
          p2: [0.0, 40.0, 0.0]
          radius: 2.5

        - cuboid:
          center: [0.0, 0.0, -18.0]
          dimensions: [40.0, 40.0, 12.0]
          rotation: [0.0, 0.0, 0.0]
)yaml";

  Parser::Node node = load_node(yamlText);

  CSGTree tree = CSG_Parser::parse_CSG_tree(node, meta, lsTimer);

  ASSERT_FALSE(tree.empty());

  std::shared_ptr<CSGTreeNode> root = tree.root();
  ASSERT_NE(root, nullptr);

  EXPECT_EQ(root->node_operator(), CSGOperation::Difference);

  ASSERT_NE(root->left_operand(), nullptr);
  ASSERT_NE(root->right_operand(), nullptr);

  EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Union);
  EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Union);

  /*
   * Count:
   * left union:
   *   difference:
   *     cuboid
   *     union: cuboid + cylinder + cylinder
   *   sphere
   *   cylinder
   *
   * right union:
   *   cylinder
   *   cylinder
   *   cuboid
   *
   * Total surfaces = 1 + 3 + 1 + 1 + 3 = 9
   */
  EXPECT_EQ(tree.gather_surfaces().size(), 9u);
}

#else

TEST(SurfaceParserTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
