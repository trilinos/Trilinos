// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDMesh.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_SubElementChildNodeAncestry.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace krino {

SubElementChildNodeAncestry::SubElementChildNodeAncestry( stk::CommBuffer & b )
{
  size_t ancestrySize = 0;
  b.unpack(ancestrySize);
  STK_ThrowAssert(ancestrySize > 0);

  myAncestry.reserve(ancestrySize);

  size_t numParents = 0;
  std::vector<stk::mesh::EntityId> parentIDs;
  std::vector<double> parentWeights;

  for(unsigned index=0; index<ancestrySize; ++index)
  {
    b.unpack(numParents);
    parentIDs.resize(numParents);
    parentWeights.resize(numParents);
    for (size_t i=0; i<numParents; ++i)
    {
      b.unpack(parentIDs[i]);
      b.unpack(parentWeights[i]);
    }
    myAncestry.emplace_back(parentIDs, parentWeights);
  }
}

void
SubElementChildNodeAncestry::build_ancestry(const SubElementNode * node)
{
  if (node->is_mesh_node())
  {
    myAncestry.emplace_back(std::vector<stk::mesh::EntityId>{node->entityId()}, std::vector<double>{0.});
    return;
  }

  const NodeVec & parents = node->get_parents();
  const auto & parentWeights = node->get_parent_weights();
  std::vector<stk::mesh::EntityId> parentIDs;
  parentIDs.reserve(parents.size());

  for (auto && parent : node->get_parents())
  {
    if (parent->is_mesh_node())
      parentIDs.push_back(parent->entityId());
    else
      parentIDs.push_back(0);
  }

  myAncestry.emplace_back(parentIDs, parentWeights);

  for (auto && parent : node->get_parents())
    if (!parent->is_mesh_node())
      build_ancestry(parent);
}

void
SubElementChildNodeAncestry::get_parent_node_keys(std::vector<stk::mesh::EntityKey> & parentNodeKeys) const
{
  parentNodeKeys.clear();
  for (auto&& cut : myAncestry)
    for (auto && parentID : cut.myParentIDs)
      if (parentID != 0)
        parentNodeKeys.emplace_back(stk::topology::NODE_RANK, parentID);
}

bool
SubElementChildNodeAncestry::is_shared(const stk::mesh::BulkData & mesh, const SubElementNode * node)
{
  if (node->is_mesh_node())
  {
    return mesh.bucket(node->entity()).shared();
  }
  for (auto && parent : node->get_parents())
  {
    if (!is_shared(mesh, parent))
    {
      return false;
    }
  }
  return true;
}

void
SubElementChildNodeAncestry::pack_into_buffer(stk::CommBuffer & b) const
{
  STK_ThrowAssert(!myAncestry.empty());
  const size_t ancestrySize = myAncestry.size();
  b.pack(ancestrySize);
  for (auto&& cut : myAncestry)
  {
    const size_t numParents = cut.myParentIDs.size();
    b.pack(numParents);
    for (size_t i=0; i<numParents; ++i)
    {
      b.pack(cut.myParentIDs[i]);
      b.pack(cut.myParentWeights[i]);
    }
  }
}

const SubElementNode *
SubElementChildNodeAncestry::find_subelement_node(CDMesh & mesh) const
{
  unsigned ancestryIndex = 0;
  return find_subelement_node(mesh, ancestryIndex);
}

const SubElementNode *
SubElementChildNodeAncestry::find_subelement_node(CDMesh & mesh, unsigned & ancestryIndex) const
{
  STK_ThrowAssert(ancestryIndex < myAncestry.size());
  const Cut & cut = myAncestry[ancestryIndex];
  const size_t numParents = cut.myParentIDs.size();
  std::vector<const SubElementNode*> parents(numParents, nullptr);

  for (size_t i=0; i<numParents; ++i)
  {
    if (cut.myParentIDs[i] == 0)
      parents[i] = find_subelement_node(mesh, ++ancestryIndex);
    else
      parents[i] = mesh.get_mesh_node(cut.myParentIDs[i]);

    if (parents[i] == nullptr)
      return nullptr;
  }

  if (numParents == 1)
  {
    return parents[0];
  }

  return SubElementNode::common_child(parents);
}

void
SubElementChildNodeAncestry::build_missing_child_nodes(CDMesh & mesh) const
{
  unsigned ancestry_index = 0;
  build_missing_child_nodes(mesh, ancestry_index);
}

std::vector<SubElementChildNodeAncestry>
SubElementChildNodeAncestry::get_constrained_node_ancestries(const std::unordered_map<stk::mesh::EntityId, std::vector<stk::mesh::EntityId> > & constrainedNodeMap) const
{
  std::set<stk::mesh::EntityId> parentNodeIDs;
  for (auto&& cut : myAncestry)
    for (auto && parentID : cut.myParentIDs)
      if (parentID != 0)
        parentNodeIDs.insert(parentID);

  std::vector<SubElementChildNodeAncestry> constrainedNodeAncestries;
  constrainedNodeAncestries.push_back(*this);
  for (auto&& parentNodeID : parentNodeIDs)
  {
    std::vector<SubElementChildNodeAncestry> toBeConverted;
    toBeConverted.swap(constrainedNodeAncestries);
    auto it = constrainedNodeMap.find(parentNodeID);
    if (it == constrainedNodeMap.end()) break;
    constrainedNodeAncestries.reserve(it->second.size() * toBeConverted.size());
    for (auto&& constrainedParentNodeID : it->second)
    {
      for (auto&& unconverted : toBeConverted)
      {
        SubElementChildNodeAncestry converted = unconverted;
        for (auto&& cut : converted.myAncestry)
        {
          for (auto && cutParentNodeID : cut.myParentIDs)
            if (cutParentNodeID == parentNodeID)
              cutParentNodeID = constrainedParentNodeID;
        }
        constrainedNodeAncestries.push_back(converted);
      }
    }
  }
  return constrainedNodeAncestries;
}

const SubElementNode *
SubElementChildNodeAncestry::build_missing_child_nodes(CDMesh & mesh, unsigned & ancestryIndex) const
{
  STK_ThrowAssert(ancestryIndex < myAncestry.size());
  const Cut & cut = myAncestry[ancestryIndex];
  const size_t numParents = cut.myParentIDs.size();
  std::vector<const SubElementNode*> parents(numParents, nullptr);

  for (size_t i=0; i<numParents; ++i)
  {
    if (cut.myParentIDs[i] == 0)
      parents[i] = find_subelement_node(mesh, ++ancestryIndex);
    else
      parents[i] = mesh.get_mesh_node(cut.myParentIDs[i]);

    if (parents[i] == nullptr)
      return nullptr;
  }

  NodeSet ancestors;
  for (auto && parent : parents)
    parent->get_ancestors(ancestors);

  std::vector<stk::mesh::Entity> ancestorNodes;
  ancestorNodes.reserve(ancestors.size());
  for (auto&& ancestor : ancestors)
  {
    STK_ThrowAssert(ancestor->entity_is_valid(mesh.stk_bulk()));
    ancestorNodes.push_back(ancestor->entity());
  }

  std::vector<stk::mesh::Entity> elems;
  stk::mesh::get_entities_through_relations(mesh.stk_bulk(), ancestorNodes, stk::topology::ELEMENT_RANK, elems);

  const Mesh_Element * owner = nullptr;
  for(auto && elem : elems)
  {
    owner = mesh.find_mesh_element(mesh.stk_bulk().identifier(elem));
    if(owner != nullptr) break;
  }
  if(!owner)
  {
    // It is possible in parallel for both nodes of a parent edge to be shared with a processor,
    // but that processor does not own any of the elements connected to that edge. In that case
    // we do not need to build the missing edge node.
    return nullptr;
  }

  if (parents.size() == 2)
    return mesh.create_edge_node(owner, parents[0], parents[1], cut.myParentWeights[1]);

  return mesh.create_child_internal_or_face_node(owner, parents, cut.myParentWeights);
}

}
