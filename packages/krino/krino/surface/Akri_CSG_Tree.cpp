// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CSG_Tree.hpp>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace krino{

// ------------------------------------------------------------
// CSGTreeNode implementation
// ------------------------------------------------------------

CSGTreeNode::CSGTreeNode(std::shared_ptr<Surface> s)
    : myOp(CSGOperation::Leaf),
      mySurface{std::move(s)}
{
    if (!surface())
    {
        throw std::invalid_argument("CSG leaf cannot contain a null surface.");
    }
}

CSGTreeNode::CSGTreeNode(
    CSGOperation operation,
    std::shared_ptr<CSGTreeNode> lhs,
    std::shared_ptr<CSGTreeNode> rhs)
    : myOp(operation),
      myOperands{std::move(lhs), std::move(rhs)}
{
    if (myOp == CSGOperation::Leaf)
    {
        throw std::invalid_argument("Use the leaf constructor for Leaf nodes.");
    }

    if (!left_operand() || !right_operand())
    {
        throw std::invalid_argument("CSG operation nodes must be non-null.");
    }
}

CSGTreeNode::CSGTreeNode(
    CSGOperation operation,
    const std::vector<std::shared_ptr<CSGTreeNode>> & vecOperands)
: myOp(operation),
  myOperands(vecOperands)
{
  if (myOp != CSGOperation::Union && myOp != CSGOperation::Intersection)
  {
      throw std::invalid_argument(
          "Vector operand constructor is only valid for Union and Intersection nodes.");
  }
  if (myOperands.empty())
  {
      throw std::invalid_argument("Expects a non-empty vector.");
  }
  for (const auto & operand : myOperands)
  {
      if (!operand)
      {
          throw std::invalid_argument("CSG operation nodes must all be non-null.");
      }
  }
}

const std::shared_ptr<Surface> & CSGTreeNode::surface() const
{
  if (myOp != CSGOperation::Leaf)
  {
      throw std::invalid_argument("surface() method requires a Leaf node.");
  }

  return mySurface;
}

const std::shared_ptr<CSGTreeNode> & CSGTreeNode::left_operand() const
{
  if (myOperands.size() != 2)
  {
    throw std::invalid_argument("left_operand() method requires a node with exact two operands.");
  }
  return myOperands[0];
}

const std::shared_ptr<CSGTreeNode> & CSGTreeNode::right_operand() const
{
  if (myOperands.size() != 2)
  {
    throw std::invalid_argument("right_operand() method requires a node with exact two operands.");
  }
  return myOperands[1];
}

const std::vector<std::shared_ptr<CSGTreeNode>> & CSGTreeNode::vector_operands() const
{
  return myOperands;
}

// ------------------------------------------------------------
// CSGTree implementation
// ------------------------------------------------------------

CSGTree::CSGTree() = default;

CSGTree::CSGTree(std::shared_ptr<Surface> surface)
    : myRoot(std::make_shared<CSGTreeNode>(std::move(surface)))
{
}

CSGTree::CSGTree(std::shared_ptr<CSGTreeNode> root)
    : myRoot(std::move(root))
{
    if (!myRoot)
    {
        throw std::invalid_argument("CSGTree root cannot be null.");
    }
}

CSGTree CSGTree::unite(const CSGTree& a, const CSGTree& b)
{
    return CSGTree(
        std::make_shared<CSGTreeNode>(
            CSGOperation::Union,
            a.requireRoot(),
            b.requireRoot()));
}

static std::vector<std::shared_ptr<CSGTreeNode>> operands_from_trees(const std::vector<CSGTree> & trees)
{
  if (trees.empty())
  {
      throw std::invalid_argument(
          "operands_from_trees expects a non-empty vector.");
  }
  std::vector<std::shared_ptr<CSGTreeNode>> operands;
  operands.reserve(trees.size());
  for (const auto & tree : trees)
  {
      operands.emplace_back(tree.requireRoot());
  }
  return operands;
}

static std::vector<CSGTree> trees_from_surfaces(const std::vector<std::shared_ptr<Surface>> & surfaces)
{
if (surfaces.empty())
  {
      throw std::invalid_argument(
          "trees_from_surfaces expects a non-empty vector.");
  }
  std::vector<CSGTree> trees;
  trees.reserve(surfaces.size());
  for (const auto & surface : surfaces)
  {
    trees.emplace_back(CSGTree(surface));
  }
  return trees;
}

CSGTree CSGTree::unite(const std::vector<CSGTree> & trees)
{
  if (trees.size() == 1)
  {
      return CSGTree(trees[0].requireRoot());
  }

  return CSGTree(
      std::make_shared<CSGTreeNode>(
          CSGOperation::Union,
          operands_from_trees(trees)));
}

CSGTree CSGTree::unite(const std::vector<std::shared_ptr<Surface>> & surfaces)
{
  return unite(trees_from_surfaces(surfaces));
}

CSGTree CSGTree::intersect(const std::vector<CSGTree> & trees)
{
  if (trees.size() == 1)
  {
      return CSGTree(trees[0].requireRoot());
  }

  return CSGTree(
      std::make_shared<CSGTreeNode>(
          CSGOperation::Intersection,
          operands_from_trees(trees)));
}

CSGTree CSGTree::intersect(const std::vector<std::shared_ptr<Surface>> & surfaces)
{
  return intersect(trees_from_surfaces(surfaces));
}

CSGTree CSGTree::intersect(const CSGTree& a, const CSGTree& b)
{
    return CSGTree(
        std::make_shared<CSGTreeNode>(
            CSGOperation::Intersection,
            a.requireRoot(),
            b.requireRoot()));
}

CSGTree CSGTree::subtract(const CSGTree& a, const CSGTree& b)
{
    return CSGTree(
        std::make_shared<CSGTreeNode>(
            CSGOperation::Difference,
            a.requireRoot(),
            b.requireRoot()));
}

void CSGTree::setRoot(std::shared_ptr<Surface> surface)
{
  myRoot = std::make_shared<CSGTreeNode>(std::move(surface));
}

void CSGTree::makeUnionWith(const CSGTree& other)
{
  myRoot = std::make_shared<CSGTreeNode>(
        CSGOperation::Union,
        requireRoot(),
        other.requireRoot());
}

void CSGTree::makeIntersectionWith(const CSGTree& other)
{
  myRoot = std::make_shared<CSGTreeNode>(
        CSGOperation::Intersection,
        requireRoot(),
        other.requireRoot());
}

void CSGTree::subtractInPlace(const CSGTree& other)
{
  myRoot = std::make_shared<CSGTreeNode>(
        CSGOperation::Difference,
        requireRoot(),
        other.requireRoot());
}

bool CSGTree::empty() const
{
    return myRoot == nullptr;
}

void CSGTree::clear()
{
    myRoot.reset();
}

std::shared_ptr<CSGTreeNode> CSGTree::root() const
{
    return myRoot;
}

std::shared_ptr<CSGTreeNode> CSGTree::requireRoot() const
{
    if (!myRoot)
    {
        throw std::runtime_error("CSGTree has no root node.");
    }

    return myRoot;
}

std::vector<std::shared_ptr<Surface>> CSGTree::gather_surfaces() const
{
    std::vector<std::shared_ptr<Surface>> surfaces;

    gatherSurfacesRecursive(myRoot, surfaces);

    return surfaces;
}

void CSGTree::gatherSurfacesRecursive(
    const std::shared_ptr<CSGTreeNode>& node,
    std::vector<std::shared_ptr<Surface>>& surfaces)
{
    if (!node)
    {
        return;
    }

    if (node->node_operator() == CSGOperation::Leaf)
    {
        surfaces.push_back(node->surface());
    }

    for (auto & operand : node->vector_operands())
    {
        gatherSurfacesRecursive(operand, surfaces);
    }
}

std::string CSGTree::operation_name(CSGOperation op)
{
    switch (op)
    {
        case CSGOperation::Leaf:
            return "Leaf";

        case CSGOperation::Union:
            return "Union";

        case CSGOperation::Intersection:
            return "Intersection";

        case CSGOperation::Difference:
            return "Difference";

        default:
            return "Unknown";
    }
}

}
