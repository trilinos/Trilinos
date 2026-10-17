// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Filtered_CSG_Tree.hpp>
#include <Akri_CSG_Tree.hpp>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace krino{

namespace {

using SurfacePtrSet = std::unordered_set<const Surface*>;

bool is_nearby_surface(
    const std::shared_ptr<Surface>& surface,
    const SurfacePtrSet& nearbySurfacePtrs)
{
  return surface && nearbySurfacePtrs.count(surface.get()) > 0;
}

std::shared_ptr<CSGTreeNode> filter_node(
    const std::shared_ptr<CSGTreeNode>& node,
    const SurfacePtrSet& nearbySurfacePtrs)
{
  if (!node)
  {
    return nullptr;
  }

  switch (node->node_operator())
  {
    case CSGOperation::Leaf:
    {
      if (is_nearby_surface(node->surface(), nearbySurfacePtrs))
      {
        return std::make_shared<CSGTreeNode>(node->surface());
      }

      return nullptr;
    }

    case CSGOperation::Union:
    case CSGOperation::Intersection:
    {
      std::vector<std::shared_ptr<CSGTreeNode>> filteredOperands;

      for (const auto& operand : node->vector_operands())
      {
        std::shared_ptr<CSGTreeNode> filteredOperand =
            filter_node(operand, nearbySurfacePtrs);

        if (filteredOperand)
        {
          filteredOperands.push_back(filteredOperand);
        }
      }

      if (filteredOperands.empty())
      {
        return nullptr;
      }

      if (filteredOperands.size() == 1)
      {
        return filteredOperands[0];
      }

      return std::make_shared<CSGTreeNode>(
          node->node_operator(),
          filteredOperands);
    }

    case CSGOperation::Difference:
    {
      std::shared_ptr<CSGTreeNode> filteredLeft =
          filter_node(node->left_operand(), nearbySurfacePtrs);

      std::shared_ptr<CSGTreeNode> filteredRight =
          filter_node(node->right_operand(), nearbySurfacePtrs);

      if (!filteredLeft)
      {
        return nullptr;
      }

      if (!filteredRight)
      {
        return filteredLeft;
      }

      return std::make_shared<CSGTreeNode>(
          CSGOperation::Difference,
          filteredLeft,
          filteredRight);
    }

    default:
    {
      throw std::runtime_error("Unknown CSG operation in filter_tree.");
    }
  }
}

SurfacePtrSet make_surface_ptr_set(
    const std::vector<std::shared_ptr<Surface>>& surfaces)
{
  SurfacePtrSet surfacePtrs;
  surfacePtrs.reserve(surfaces.size());

  for (const auto& surface : surfaces)
  {
    if (surface)
    {
      surfacePtrs.insert(surface.get());
    }
  }

  return surfacePtrs;
}

} // namespace

CSGTree filtered_CSG_tree(
    const CSGTree& tree,
    const std::vector<std::shared_ptr<Surface>>& nearbySurfaces)
{
  if (tree.empty() || nearbySurfaces.empty())
  {
    return CSGTree();
  }

  const SurfacePtrSet nearbySurfacePtrs =
      make_surface_ptr_set(nearbySurfaces);

  std::shared_ptr<CSGTreeNode> filteredRoot =
      filter_node(tree.root(), nearbySurfacePtrs);

  if (!filteredRoot)
  {
    return CSGTree();
  }

  return CSGTree(filteredRoot);
}

}


