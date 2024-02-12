#include <Akri_ChildNodeStencil.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>

namespace krino {

static void build_node_stencil(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & childNodeSelector,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    stk::mesh::Entity node,
    const double selfWeight,
    std::vector<stk::mesh::Entity> & parentNodes,
    std::vector<double> & parentWeights)
{
  if (!childNodeSelector(mesh.bucket(node)))
  {
    parentNodes.push_back(node);
    parentWeights.push_back(selfWeight);
    return;
  }

  const std::vector<std::pair<stk::mesh::Entity, double>> nodeParentsAndWeights = get_child_node_parents_and_weights(mesh, parentIdsField, parentWtsField, node);
  for (auto & [parent, parentWeight] : nodeParentsAndWeights)
    build_node_stencil(mesh, childNodeSelector, parentIdsField, parentWtsField, parent, selfWeight*parentWeight, parentNodes, parentWeights);
}

static void build_child_node_stencil(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & childNodeSelector,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    const stk::mesh::Entity childNode,
    std::vector<stk::mesh::Entity> & workParentNodes,
    std::vector<double> & workParentWeights)
{
  workParentNodes.clear();
  workParentWeights.clear();
  build_node_stencil(mesh, childNodeSelector, parentIdsField, parentWtsField, childNode, 1.0, workParentNodes, workParentWeights);
}

void fill_child_node_stencils(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    std::vector<ChildNodeStencil> & childNodeStencils)
{
  std::vector<stk::mesh::Entity> workParentNodes;
  std::vector<double> workParentWeights;

  const stk::mesh::Selector childNodeSelector = childNodePart;
  const stk::mesh::Selector ownedOrSharedChildNodeSelector = childNodePart & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part());
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedOrSharedChildNodeSelector))
  {
    for(const auto childNode : *bucketPtr)
    {
      build_child_node_stencil(mesh, childNodeSelector, parentIdsField, parentWtsField, childNode, workParentNodes, workParentWeights);
      childNodeStencils.emplace_back(childNode, workParentNodes, workParentWeights);
    }
  }
}

}


