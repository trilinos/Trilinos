/*
 * Akri_OrientedSideNodes.cpp
 *
 *  Created on: Aug 4, 2025
 *      Author: drnoble
 */

#include <Akri_OrientedSideNodes.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

static bool determine_polarity_for_negative_side_of_interface(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side)
{
  const stk::topology sideTopology = mesh.bucket(side).topology();
  const unsigned numSideElems = mesh.num_elements(side);
  const stk::mesh::Entity * sideElems = mesh.begin_elements(side);
  const stk::mesh::Permutation * sideElemPermutatons = mesh.begin_permutations(side, stk::topology::ELEMENT_RANK);

  for (unsigned iElem = 0; iElem < numSideElems; ++iElem)
    if (negativeSideElementSelector(mesh.bucket(sideElems[iElem])))
      return sideTopology.is_positive_polarity(sideElemPermutatons[iElem]);

  STK_ThrowRequireMsg(false, "determine_polarity_for_negative_side_of_interface has no selected element.");
  return false;
}

std::array<stk::mesh::Entity,3> get_oriented_triangle_side_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side)
{
  const stk::mesh::Entity* sideNodes = mesh.begin_nodes(side);
  const bool polarity = determine_polarity_for_negative_side_of_interface(mesh, negativeSideElementSelector, side);

  if (polarity)
    return {{sideNodes[0], sideNodes[1], sideNodes[2]}};
  return {{sideNodes[0], sideNodes[2], sideNodes[1]}};
}

std::array<stk::mesh::Entity,2> get_oriented_line_side_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side)
{
  const stk::mesh::Entity* sideNodes = mesh.begin_nodes(side);
  const bool polarity = determine_polarity_for_negative_side_of_interface(mesh, negativeSideElementSelector, side);

  if (polarity)
    return {{sideNodes[0], sideNodes[1]}};
  return {{sideNodes[1], sideNodes[0]}};
}

}

