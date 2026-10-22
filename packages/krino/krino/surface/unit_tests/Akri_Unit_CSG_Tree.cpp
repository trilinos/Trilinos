// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_CSG_Tree.hpp>

namespace krino {

static std::shared_ptr<Surface> makeSimpleSurface()
{
    return std::make_shared<Sphere>(stk::math::Vector3d(0,0,0), 1.);
}

static std::shared_ptr<Surface> makeSimpleSurface(int id)
{
    return makeSimpleSurface();
}

// -----------------------------------------------------------------------------
// Default construction
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, DefaultConstructorCreatesEmptyTree)
{
    CSGTree tree;

    EXPECT_TRUE(tree.empty());
    EXPECT_EQ(tree.root(), nullptr);
}

TEST(CSGTreeTest, RequireRootThrowsForEmptyTree)
{
    CSGTree tree;

    EXPECT_THROW(tree.requireRoot(), std::runtime_error);
}

// -----------------------------------------------------------------------------
// Construction from Surface
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, ConstructFromSurfaceCreatesLeafNode)
{
    auto surface = makeSimpleSurface();

    CSGTree tree(surface);

    ASSERT_FALSE(tree.empty());

    auto root = tree.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->surface(), surface);
    EXPECT_TRUE(root->vector_operands().empty());
}

TEST(CSGTreeTest, ConstructFromNullSurfaceThrows)
{
    std::shared_ptr<Surface> surface;

    EXPECT_THROW(CSGTree tree(surface), std::invalid_argument);
}

// -----------------------------------------------------------------------------
// Construction from CSGTreeNode
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, ConstructFromRootNodeStoresRoot)
{
    auto surface = makeSimpleSurface();
    auto node = std::make_shared<CSGTreeNode>(surface);

    CSGTree tree(node);

    EXPECT_FALSE(tree.empty());
    EXPECT_EQ(tree.root(), node);
    EXPECT_EQ(tree.requireRoot(), node);
}

TEST(CSGTreeTest, ConstructFromNullRootNodeThrows)
{
    std::shared_ptr<CSGTreeNode> root;

    EXPECT_THROW(CSGTree tree(root), std::invalid_argument);
}

// -----------------------------------------------------------------------------
// Static unite with CSGTree arguments
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, StaticUniteCreatesUnionNode)
{
    CSGTree a(makeSimpleSurface());
    CSGTree b(makeSimpleSurface());

    auto aRoot = a.root();
    auto bRoot = b.root();

    CSGTree result = CSGTree::unite(a, b);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Union);
    EXPECT_EQ(root->left_operand(), aRoot);
    EXPECT_EQ(root->right_operand(), bRoot);
}

TEST(CSGTreeTest, StaticUniteThrowsIfLeftTreeIsEmpty)
{
    CSGTree empty;
    CSGTree b(makeSimpleSurface());

    EXPECT_THROW(CSGTree::unite(empty, b), std::runtime_error);
}

TEST(CSGTreeTest, StaticUniteThrowsIfRightTreeIsEmpty)
{
    CSGTree a(makeSimpleSurface());
    CSGTree empty;

    EXPECT_THROW(CSGTree::unite(a, empty), std::runtime_error);
}

// -----------------------------------------------------------------------------
// Static unite with Surface arguments
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, StaticUniteSurfaceOverloadCreatesUnionOfLeaves)
{
    auto surfaceA = makeSimpleSurface();
    auto surfaceB = makeSimpleSurface();

    CSGTree result = CSGTree::unite(surfaceA, surfaceB);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Union);
    ASSERT_NE(root->left_operand(), nullptr);
    ASSERT_NE(root->right_operand(), nullptr);

    EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Leaf);

    EXPECT_EQ(root->left_operand()->surface(), surfaceA);
    EXPECT_EQ(root->right_operand()->surface(), surfaceB);
}

// -----------------------------------------------------------------------------
// Static intersect with CSGTree arguments
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, StaticIntersectCreatesIntersectionNode)
{
    CSGTree a(makeSimpleSurface());
    CSGTree b(makeSimpleSurface());

    auto aRoot = a.root();
    auto bRoot = b.root();

    CSGTree result = CSGTree::intersect(a, b);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Intersection);
    EXPECT_EQ(root->left_operand(), aRoot);
    EXPECT_EQ(root->right_operand(), bRoot);
}

TEST(CSGTreeTest, StaticIntersectThrowsIfLeftTreeIsEmpty)
{
    CSGTree empty;
    CSGTree b(makeSimpleSurface());

    EXPECT_THROW(CSGTree::intersect(empty, b), std::runtime_error);
}

TEST(CSGTreeTest, StaticIntersectThrowsIfRightTreeIsEmpty)
{
    CSGTree a(makeSimpleSurface());
    CSGTree empty;

    EXPECT_THROW(CSGTree::intersect(a, empty), std::runtime_error);
}

// -----------------------------------------------------------------------------
// Static intersect with Surface arguments
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, StaticIntersectSurfaceOverloadCreatesIntersectionOfLeaves)
{
    auto surfaceA = makeSimpleSurface();
    auto surfaceB = makeSimpleSurface();

    CSGTree result = CSGTree::intersect(surfaceA, surfaceB);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Intersection);
    ASSERT_NE(root->left_operand(), nullptr);
    ASSERT_NE(root->right_operand(), nullptr);

    EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Leaf);

    EXPECT_EQ(root->left_operand()->surface(), surfaceA);
    EXPECT_EQ(root->right_operand()->surface(), surfaceB);
}

// -----------------------------------------------------------------------------
// Static subtract with CSGTree arguments
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, StaticSubtractCreatesDifferenceNode)
{
    CSGTree a(makeSimpleSurface());
    CSGTree b(makeSimpleSurface());

    auto aRoot = a.root();
    auto bRoot = b.root();

    CSGTree result = CSGTree::subtract(a, b);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Difference);
    EXPECT_EQ(root->left_operand(), aRoot);
    EXPECT_EQ(root->right_operand(), bRoot);
}

TEST(CSGTreeTest, StaticSubtractThrowsIfLeftTreeIsEmpty)
{
    CSGTree empty;
    CSGTree b(makeSimpleSurface());

    EXPECT_THROW(CSGTree::subtract(empty, b), std::runtime_error);
}

TEST(CSGTreeTest, StaticSubtractThrowsIfRightTreeIsEmpty)
{
    CSGTree a(makeSimpleSurface());
    CSGTree empty;

    EXPECT_THROW(CSGTree::subtract(a, empty), std::runtime_error);
}

// -----------------------------------------------------------------------------
// Static subtract with Surface arguments
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, StaticSubtractSurfaceOverloadCreatesDifferenceOfLeaves)
{
    auto surfaceA = makeSimpleSurface();
    auto surfaceB = makeSimpleSurface();

    CSGTree result = CSGTree::subtract(surfaceA, surfaceB);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Difference);
    ASSERT_NE(root->left_operand(), nullptr);
    ASSERT_NE(root->right_operand(), nullptr);

    EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Leaf);

    EXPECT_EQ(root->left_operand()->surface(), surfaceA);
    EXPECT_EQ(root->right_operand()->surface(), surfaceB);
}

// -----------------------------------------------------------------------------
// setRoot
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, SetRootOnEmptyTreeCreatesLeafNode)
{
    CSGTree tree;
    auto surface = makeSimpleSurface();

    tree.setRoot(surface);

    ASSERT_FALSE(tree.empty());

    auto root = tree.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->surface(), surface);
    EXPECT_TRUE(root->vector_operands().empty());
}

TEST(CSGTreeTest, SetRootReplacesExistingRoot)
{
    auto surfaceA = makeSimpleSurface();
    auto surfaceB = makeSimpleSurface();

    CSGTree tree(surfaceA);

    auto oldRoot = tree.root();

    tree.setRoot(surfaceB);

    auto newRoot = tree.root();

    ASSERT_NE(newRoot, nullptr);
    EXPECT_NE(newRoot, oldRoot);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(newRoot->surface(), surfaceB);
}

TEST(CSGTreeTest, SetRootWithNullSurfaceThrows)
{
    CSGTree tree;
    std::shared_ptr<Surface> surface;

    EXPECT_THROW(tree.setRoot(surface), std::invalid_argument);
}

// -----------------------------------------------------------------------------
// makeUnionWith
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, MakeUnionWithTreeCreatesUnionUsingExistingRootAsLeftChild)
{
    CSGTree tree(makeSimpleSurface());
    CSGTree other(makeSimpleSurface());

    auto oldRoot = tree.root();
    auto otherRoot = other.root();

    tree.makeUnionWith(other);

    auto newRoot = tree.root();
    ASSERT_NE(newRoot, nullptr);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Union);
    EXPECT_EQ(newRoot->left_operand(), oldRoot);
    EXPECT_EQ(newRoot->right_operand(), otherRoot);
}

TEST(CSGTreeTest, MakeUnionWithSurfaceCreatesUnionWithSurfaceLeaf)
{
    CSGTree tree(makeSimpleSurface());
    auto otherSurface = makeSimpleSurface();

    auto oldRoot = tree.root();

    tree.makeUnionWith(otherSurface);

    auto newRoot = tree.root();
    ASSERT_NE(newRoot, nullptr);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Union);
    EXPECT_EQ(newRoot->left_operand(), oldRoot);

    ASSERT_NE(newRoot->right_operand(), nullptr);
    EXPECT_EQ(newRoot->right_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(newRoot->right_operand()->surface(), otherSurface);
}

TEST(CSGTreeTest, MakeUnionWithThrowsWhenThisTreeIsEmpty)
{
    CSGTree tree;
    CSGTree other(makeSimpleSurface());

    EXPECT_THROW(tree.makeUnionWith(other), std::runtime_error);
}

TEST(CSGTreeTest, MakeUnionWithThrowsWhenOtherTreeIsEmpty)
{
    CSGTree tree(makeSimpleSurface());
    CSGTree other;

    EXPECT_THROW(tree.makeUnionWith(other), std::runtime_error);
}

// -----------------------------------------------------------------------------
// makeIntersectionWith
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, MakeIntersectionWithTreeCreatesIntersectionUsingExistingRootAsLeftChild)
{
    CSGTree tree(makeSimpleSurface());
    CSGTree other(makeSimpleSurface());

    auto oldRoot = tree.root();
    auto otherRoot = other.root();

    tree.makeIntersectionWith(other);

    auto newRoot = tree.root();
    ASSERT_NE(newRoot, nullptr);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Intersection);
    EXPECT_EQ(newRoot->left_operand(), oldRoot);
    EXPECT_EQ(newRoot->right_operand(), otherRoot);
}

TEST(CSGTreeTest, MakeIntersectionWithSurfaceCreatesIntersectionWithSurfaceLeaf)
{
    CSGTree tree(makeSimpleSurface());
    auto otherSurface = makeSimpleSurface();

    auto oldRoot = tree.root();

    tree.makeIntersectionWith(otherSurface);

    auto newRoot = tree.root();
    ASSERT_NE(newRoot, nullptr);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Intersection);
    EXPECT_EQ(newRoot->left_operand(), oldRoot);

    ASSERT_NE(newRoot->right_operand(), nullptr);
    EXPECT_EQ(newRoot->right_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(newRoot->right_operand()->surface(), otherSurface);
}

TEST(CSGTreeTest, MakeIntersectionWithThrowsWhenThisTreeIsEmpty)
{
    CSGTree tree;
    CSGTree other(makeSimpleSurface());

    EXPECT_THROW(tree.makeIntersectionWith(other), std::runtime_error);
}

TEST(CSGTreeTest, MakeIntersectionWithThrowsWhenOtherTreeIsEmpty)
{
    CSGTree tree(makeSimpleSurface());
    CSGTree other;

    EXPECT_THROW(tree.makeIntersectionWith(other), std::runtime_error);
}

// -----------------------------------------------------------------------------
// subtractInPlace
// -----------------------------------------------------------------------------

TEST(CSGTreeTest, SubtractInPlaceWithTreeCreatesDifferenceUsingExistingRootAsLeftChild)
{
    CSGTree tree(makeSimpleSurface());
    CSGTree other(makeSimpleSurface());

    auto oldRoot = tree.root();
    auto otherRoot = other.root();

    tree.subtractInPlace(other);

    auto newRoot = tree.root();
    ASSERT_NE(newRoot, nullptr);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Difference);
    EXPECT_EQ(newRoot->left_operand(), oldRoot);
    EXPECT_EQ(newRoot->right_operand(), otherRoot);
}

TEST(CSGTreeTest, SubtractInPlaceWithSurfaceCreatesDifferenceWithSurfaceLeaf)
{
    CSGTree tree(makeSimpleSurface());
    auto otherSurface = makeSimpleSurface();

    auto oldRoot = tree.root();

    tree.subtractInPlace(otherSurface);

    auto newRoot = tree.root();
    ASSERT_NE(newRoot, nullptr);

    EXPECT_EQ(newRoot->node_operator(), CSGOperation::Difference);
    EXPECT_EQ(newRoot->left_operand(), oldRoot);

    ASSERT_NE(newRoot->right_operand(), nullptr);
    EXPECT_EQ(newRoot->right_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(newRoot->right_operand()->surface(), otherSurface);
}

TEST(CSGTreeTest, SubtractInPlaceThrowsWhenThisTreeIsEmpty)
{
    CSGTree tree;
    CSGTree other(makeSimpleSurface());

    EXPECT_THROW(tree.subtractInPlace(other), std::runtime_error);
}

TEST(CSGTreeTest, SubtractInPlaceThrowsWhenOtherTreeIsEmpty)
{
    CSGTree tree(makeSimpleSurface());
    CSGTree other;

    EXPECT_THROW(tree.subtractInPlace(other), std::runtime_error);
}

static bool containsSurface(
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const std::shared_ptr<Surface>& surface)
{
    return std::find(surfaces.begin(), surfaces.end(), surface) != surfaces.end();
}

static std::size_t countSurface(
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const std::shared_ptr<Surface>& surface)
{
    return static_cast<std::size_t>(
        std::count(surfaces.begin(), surfaces.end(), surface));
}

// -----------------------------------------------------------------------------
// gather_surfaces tests
// -----------------------------------------------------------------------------

TEST(CSGTreeGatherSurfacesTest, EmptyTreeReturnsEmptyVector)
{
    CSGTree tree;

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    EXPECT_TRUE(surfaces.empty());
}

TEST(CSGTreeGatherSurfacesTest, SinglePrimitiveReturnsOneSurface)
{
    auto surface = makeSimpleSurface(1);

    CSGTree tree(surface);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 1u);
    EXPECT_EQ(surfaces[0], surface);
}

TEST(CSGTreeGatherSurfacesTest, UnionOfTwoSurfacesReturnsBothSurfaces)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);

    CSGTree tree = CSGTree::unite(surfaceA, surfaceB);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 2u);

    EXPECT_TRUE(containsSurface(surfaces, surfaceA));
    EXPECT_TRUE(containsSurface(surfaces, surfaceB));
}

TEST(CSGTreeGatherSurfacesTest, IntersectionOfTwoSurfacesReturnsBothSurfaces)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);

    CSGTree tree = CSGTree::intersect(surfaceA, surfaceB);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 2u);

    EXPECT_TRUE(containsSurface(surfaces, surfaceA));
    EXPECT_TRUE(containsSurface(surfaces, surfaceB));
}

TEST(CSGTreeGatherSurfacesTest, DifferenceOfTwoSurfacesReturnsBothSurfaces)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);

    CSGTree tree = CSGTree::subtract(surfaceA, surfaceB);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 2u);

    EXPECT_TRUE(containsSurface(surfaces, surfaceA));
    EXPECT_TRUE(containsSurface(surfaces, surfaceB));
}

TEST(CSGTreeGatherSurfacesTest, ComplexTreeReturnsAllLeafSurfaces)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);
    auto surfaceC = makeSimpleSurface(3);
    auto surfaceD = makeSimpleSurface(4);

    CSGTree treeA(surfaceA);
    CSGTree treeB(surfaceB);
    CSGTree treeC(surfaceC);
    CSGTree treeD(surfaceD);

    /* Build:

            Difference
            /        \
        Union      Intersection
        /   \        /       \
       A     B      C         D

    */
    CSGTree left = CSGTree::unite(treeA, treeB);
    CSGTree right = CSGTree::intersect(treeC, treeD);
    CSGTree tree = CSGTree::subtract(left, right);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 4u);

    EXPECT_TRUE(containsSurface(surfaces, surfaceA));
    EXPECT_TRUE(containsSurface(surfaces, surfaceB));
    EXPECT_TRUE(containsSurface(surfaces, surfaceC));
    EXPECT_TRUE(containsSurface(surfaces, surfaceD));
}

TEST(CSGTreeGatherSurfacesTest, MutatingBuilderMethodsAreReflectedInGatherSurfaces)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);
    auto surfaceC = makeSimpleSurface(3);
    auto surfaceD = makeSimpleSurface(4);

    CSGTree tree(surfaceA);

    tree.makeUnionWith(surfaceB);
    tree.makeIntersectionWith(surfaceC);
    tree.subtractInPlace(surfaceD);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 4u);

    EXPECT_TRUE(containsSurface(surfaces, surfaceA));
    EXPECT_TRUE(containsSurface(surfaces, surfaceB));
    EXPECT_TRUE(containsSurface(surfaces, surfaceC));
    EXPECT_TRUE(containsSurface(surfaces, surfaceD));
}

TEST(CSGTreeGatherSurfacesTest, SetRootReplacesPreviouslyGatheredSurfaces)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);
    auto surfaceC = makeSimpleSurface(3);

    CSGTree tree = CSGTree::unite(surfaceA, surfaceB);

    std::vector<std::shared_ptr<Surface>> surfacesBefore = tree.gather_surfaces();

    ASSERT_EQ(surfacesBefore.size(), 2u);
    EXPECT_TRUE(containsSurface(surfacesBefore, surfaceA));
    EXPECT_TRUE(containsSurface(surfacesBefore, surfaceB));

    tree.setRoot(surfaceC);

    std::vector<std::shared_ptr<Surface>> surfacesAfter = tree.gather_surfaces();

    ASSERT_EQ(surfacesAfter.size(), 1u);
    EXPECT_EQ(surfacesAfter[0], surfaceC);

    EXPECT_FALSE(containsSurface(surfacesAfter, surfaceA));
    EXPECT_FALSE(containsSurface(surfacesAfter, surfaceB));
}

TEST(CSGTreeGatherSurfacesTest, DuplicateSurfacePointerAppearsMultipleTimes)
{
    auto sharedSurface = makeSimpleSurface(1);

    CSGTree tree = CSGTree::unite(sharedSurface, sharedSurface);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 2u);

    EXPECT_EQ(countSurface(surfaces, sharedSurface), 2u);
}

TEST(CSGTreeGatherSurfacesTest, TraversalOrderIsLeftThenRightForSimpleUnion)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);

    CSGTree tree = CSGTree::unite(surfaceA, surfaceB);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 2u);

    EXPECT_EQ(surfaces[0], surfaceA);
    EXPECT_EQ(surfaces[1], surfaceB);
}

TEST(CSGTreeGatherSurfacesTest, TraversalOrderIsDepthFirstLeftThenRight)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);
    auto surfaceC = makeSimpleSurface(3);
    auto surfaceD = makeSimpleSurface(4);

    CSGTree left = CSGTree::unite(surfaceA, surfaceB);
    CSGTree right = CSGTree::intersect(surfaceC, surfaceD);
    CSGTree tree = CSGTree::subtract(left, right);

    std::vector<std::shared_ptr<Surface>> surfaces = tree.gather_surfaces();

    ASSERT_EQ(surfaces.size(), 4u);

    EXPECT_EQ(surfaces[0], surfaceA);
    EXPECT_EQ(surfaces[1], surfaceB);
    EXPECT_EQ(surfaces[2], surfaceC);
    EXPECT_EQ(surfaces[3], surfaceD);
}

// -----------------------------------------------------------------------------
// n-ary unite tests
// -----------------------------------------------------------------------------

TEST(CSGTreeNaryTest, UniteWithEmptyVectorThrows)
{
    std::vector<CSGTree> operands;

    EXPECT_THROW(
        CSGTree::unite(operands),
        std::invalid_argument);
}

TEST(CSGTreeNaryTest, UniteWithOneOperandReturnsEquivalentTree)
{
    auto surfaceA = makeSimpleSurface(1);
    CSGTree treeA(surfaceA);

    std::vector<CSGTree> operands{treeA};

    CSGTree result = CSGTree::unite(operands);

    ASSERT_FALSE(result.empty());
    ASSERT_NE(result.root(), nullptr);

    EXPECT_EQ(result.root()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(result.root()->surface(), surfaceA);
    EXPECT_TRUE(result.root()->vector_operands().empty());
}

TEST(CSGTreeNaryTest, UniteWithTwoOperandsCreatesUnion)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);

    CSGTree treeA(surfaceA);
    CSGTree treeB(surfaceB);

    std::vector<CSGTree> operands{treeA, treeB};

    CSGTree result = CSGTree::unite(operands);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Union);

    ASSERT_NE(root->left_operand(), nullptr);
    ASSERT_NE(root->right_operand(), nullptr);

    EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->left_operand()->surface(), surfaceA);

    EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->right_operand()->surface(), surfaceB);
}

TEST(CSGTreeNaryTest, UniteWithThreeOperandsCreatesNaryUnionTree)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);
    auto surfaceC = makeSimpleSurface(3);

    CSGTree treeA(surfaceA);
    CSGTree treeB(surfaceB);
    CSGTree treeC(surfaceC);

    std::vector<CSGTree> operands{treeA, treeB, treeC};

    CSGTree result = CSGTree::unite(operands);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Union);

    const auto& rootOperands = root->vector_operands();

    ASSERT_EQ(rootOperands.size(), 3u);

    ASSERT_NE(rootOperands[0], nullptr);
    ASSERT_NE(rootOperands[1], nullptr);
    ASSERT_NE(rootOperands[2], nullptr);

    EXPECT_EQ(rootOperands[0]->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(rootOperands[0]->surface(), surfaceA);

    EXPECT_EQ(rootOperands[1]->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(rootOperands[1]->surface(), surfaceB);

    EXPECT_EQ(rootOperands[2]->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(rootOperands[2]->surface(), surfaceC);

    EXPECT_THROW(root->left_operand(), std::invalid_argument);
    EXPECT_THROW(root->right_operand(), std::invalid_argument);
}

// -----------------------------------------------------------------------------
// n-ary intersect tests
// -----------------------------------------------------------------------------

TEST(CSGTreeNaryTest, IntersectWithEmptyVectorThrows)
{
    std::vector<CSGTree> operands;

    EXPECT_THROW(
        CSGTree::intersect(operands),
        std::invalid_argument);
}

TEST(CSGTreeNaryTest, IntersectWithOneOperandReturnsEquivalentTree)
{
    auto surfaceA = makeSimpleSurface(1);
    CSGTree treeA(surfaceA);

    std::vector<CSGTree> operands{treeA};

    CSGTree result = CSGTree::intersect(operands);

    ASSERT_FALSE(result.empty());
    ASSERT_NE(result.root(), nullptr);

    EXPECT_EQ(result.root()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(result.root()->surface(), surfaceA);
    EXPECT_TRUE(result.root()->vector_operands().empty());
}

TEST(CSGTreeNaryTest, IntersectWithTwoOperandsCreatesIntersection)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);

    CSGTree treeA(surfaceA);
    CSGTree treeB(surfaceB);

    std::vector<CSGTree> operands{treeA, treeB};

    CSGTree result = CSGTree::intersect(operands);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Intersection);

    ASSERT_NE(root->left_operand(), nullptr);
    ASSERT_NE(root->right_operand(), nullptr);

    EXPECT_EQ(root->left_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->left_operand()->surface(), surfaceA);

    EXPECT_EQ(root->right_operand()->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(root->right_operand()->surface(), surfaceB);
}

TEST(CSGTreeNaryTest, IntersectWithThreeOperandsCreatesNaryIntersectionTree)
{
    auto surfaceA = makeSimpleSurface(1);
    auto surfaceB = makeSimpleSurface(2);
    auto surfaceC = makeSimpleSurface(3);

    CSGTree treeA(surfaceA);
    CSGTree treeB(surfaceB);
    CSGTree treeC(surfaceC);

    std::vector<CSGTree> operands{treeA, treeB, treeC};

    CSGTree result = CSGTree::intersect(operands);

    ASSERT_FALSE(result.empty());

    auto root = result.root();
    ASSERT_NE(root, nullptr);

    EXPECT_EQ(root->node_operator(), CSGOperation::Intersection);

    const auto& rootOperands = root->vector_operands();

    ASSERT_EQ(rootOperands.size(), 3u);

    ASSERT_NE(rootOperands[0], nullptr);
    ASSERT_NE(rootOperands[1], nullptr);
    ASSERT_NE(rootOperands[2], nullptr);

    EXPECT_EQ(rootOperands[0]->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(rootOperands[0]->surface(), surfaceA);

    EXPECT_EQ(rootOperands[1]->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(rootOperands[1]->surface(), surfaceB);

    EXPECT_EQ(rootOperands[2]->node_operator(), CSGOperation::Leaf);
    EXPECT_EQ(rootOperands[2]->surface(), surfaceC);

    EXPECT_THROW(root->left_operand(), std::invalid_argument);
    EXPECT_THROW(root->right_operand(), std::invalid_argument);
}


}

