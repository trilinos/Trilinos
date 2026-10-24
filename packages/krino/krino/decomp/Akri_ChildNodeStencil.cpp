#include <Akri_ChildNodeStencil.hpp>

#include <iterator>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>

namespace krino {

void insert_parent_into_stencil(const stk::mesh::Entity node,
    const double wt,
    std::vector<stk::mesh::Entity> & parentNodes,
    std::vector<double> & parentWeights)
{
  auto result = std::find(parentNodes.begin(), parentNodes.end(), node);
  if(result == parentNodes.end())
  {
    parentNodes.push_back(node);
    parentWeights.push_back(wt);
  }
  else
  {
    auto idx = std::distance(parentNodes.begin(), result);
    parentWeights[idx] += wt;
  }
}

static void build_node_stencil_wrt_background_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    stk::mesh::Entity node,
    const double selfWeight,
    std::vector<stk::mesh::Entity> & parentNodes,
    std::vector<double> & parentWeights)
{
  if (mesh.bucket(node).member(childNodePart))
  {
    const std::vector<std::pair<stk::mesh::Entity, double>> nodeParentsAndWeights = get_child_node_parents_and_weights(mesh, parentIdsField, parentWtsField, node);
    for (auto & [parent, parentWeight] : nodeParentsAndWeights)
      build_node_stencil_wrt_background_nodes(mesh, childNodePart, parentIdsField, parentWtsField, parent, selfWeight*parentWeight, parentNodes, parentWeights);
  }
  else
  {
    insert_parent_into_stencil(node, selfWeight, parentNodes, parentWeights);
  }
}

static void build_child_node_stencil_wrt_background_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    const stk::mesh::Entity childNode,
    std::vector<stk::mesh::Entity> & workParentNodes,
    std::vector<double> & workParentWeights)
{
  workParentNodes.clear();
  workParentWeights.clear();
  build_node_stencil_wrt_background_nodes(mesh, childNodePart, parentIdsField, parentWtsField, childNode, 1.0, workParentNodes, workParentWeights);
}

static void fill_child_node_stencils_for_selected_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const stk::mesh::Selector & interfaceChildNodeSelector,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    std::vector<ChildNodeStencil> & childNodeStencils)
{
  std::vector<stk::mesh::Entity> workParentNodes;
  std::vector<double> workParentWeights;

  const stk::mesh::Selector ownedOrSharedInterfaceChildNodeSelector = interfaceChildNodeSelector & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part());
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedOrSharedInterfaceChildNodeSelector))
  {
    for(const auto childNode : *bucketPtr)
    {
      build_child_node_stencil_wrt_background_nodes(mesh, childNodePart, parentIdsField, parentWtsField, childNode, workParentNodes, workParentWeights);
      childNodeStencils.emplace_back(childNode, workParentNodes, workParentWeights);
    }
  }
}

void fill_child_node_stencils_for_interface(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const stk::mesh::Selector & interfaceSelector,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    std::vector<ChildNodeStencil> & childNodeStencils)
{
  const stk::mesh::Selector interfaceChildNodeSelector = interfaceSelector & childNodePart;
  fill_child_node_stencils_for_selected_nodes(mesh, childNodePart, interfaceChildNodeSelector, parentIdsField, parentWtsField, childNodeStencils);
}

void fill_child_node_stencils_for_all_interfaces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    std::vector<ChildNodeStencil> & childNodeStencils)
{
  const stk::mesh::Selector interfaceChildNodeSelector = childNodePart;
  fill_child_node_stencils_for_selected_nodes(mesh, childNodePart, interfaceChildNodeSelector, parentIdsField, parentWtsField, childNodeStencils);
}

}


