#include "Akri_NodeRefiner.hpp"

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "Akri_ChildNodeCreator.hpp"
#include "Akri_Edge.hpp"
#include "Akri_EntityIdPool.hpp"
#include "Akri_FieldRef.hpp"
#include "Akri_TransitionElementEdgeMarker.hpp"

namespace krino {

size_t NodeRefiner::get_num_edges_to_refine() const
{
  return myRefinedEdgesToChildNodes.size();
}

void NodeRefiner::assign_refined_edge_node_parent_ids(const stk::mesh::BulkData & mesh, FieldRef refinedEdgeNodeParentIdsField) const
{
  for (auto && myRefinedEdgeToChildNodes : myRefinedEdgesToChildNodes)
  {
    const auto & edgeNodes = get_edge_nodes(myRefinedEdgeToChildNodes.first);
    const stk::mesh::Entity refinedEdgeNode = myRefinedEdgeToChildNodes.second;
    auto * refinedEdgeNodeParentIds = field_data<uint64_t>(refinedEdgeNodeParentIdsField, refinedEdgeNode);
    refinedEdgeNodeParentIds[0] = mesh.identifier(edgeNodes[0]);
    refinedEdgeNodeParentIds[1] = mesh.identifier(edgeNodes[1]);
  }
}

void NodeRefiner::create_refined_edge_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedEdgeNodeParts, FieldRef refinedEdgeNodeParentIdsField)
{
  const size_t numEdgesToRefine = myRefinedEdgesToChildNodes.size();
  std::vector<std::array<stk::mesh::Entity,2>> edgesParentNodes;
  edgesParentNodes.reserve(numEdgesToRefine);
  std::vector<ChildNodeRequest> childNodeRequests;
  childNodeRequests.reserve(numEdgesToRefine);
  for (auto && myRefinedEdgeToChildNodes : myRefinedEdgesToChildNodes)
  {
    edgesParentNodes.emplace_back(get_edge_nodes(myRefinedEdgeToChildNodes.first));
    std::array<stk::mesh::Entity,2> & edgeParentNodes = edgesParentNodes.back();
    childNodeRequests.emplace_back(std::vector<stk::mesh::Entity*>{&(edgeParentNodes[0]), &(edgeParentNodes[1])}, &(myRefinedEdgeToChildNodes.second));
  }

  auto generate_new_ids = [&](stk::topology::rank_t entityRank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId>& requestedIds)
  {
     EntityIdPool::generate_new_ids(mesh, entityRank, numIdsNeeded, requestedIds, myAssert32Bit, myForce64Bit);
  };

  batch_create_child_nodes(mesh, childNodeRequests, refinedEdgeNodeParts, generate_new_ids);

  assign_refined_edge_node_parent_ids(mesh, refinedEdgeNodeParentIdsField);
}

stk::mesh::Entity NodeRefiner::get_edge_child_node(const Edge edge) const
{
  const auto iter = myRefinedEdgesToChildNodes.find(edge);
  if (iter != myRefinedEdgesToChildNodes.end())
    return iter->second;
  return stk::mesh::Entity();
}

static void prolong_edge_node(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,2> & parentNodes, const stk::mesh::Entity childNode)
{
  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();
  for ( auto && stkField : allFields )
  {
    const FieldRef field(stkField);

    if( field.entity_rank() == stk::topology::NODE_RANK && field.type_is<double>() )
    {
      const auto * parentNodeData0 = field_data<double>(field, parentNodes[0]);
      const auto * parentNodeData1 = field_data<double>(field, parentNodes[1]);
      if (nullptr != parentNodeData0 && nullptr != parentNodeData1)
      {
        auto * childNodeData = field_data<double>(field, childNode);
        if (nullptr != childNodeData)
        {
          const unsigned fieldLength = field.length();
          for (unsigned i=0; i<fieldLength; ++i)
          {
            childNodeData[i] = 0.5*(parentNodeData0[i] + parentNodeData1[i]);
          }
        }
      }
    }
  }
}

void NodeRefiner::prolong_refined_edge_nodes(const stk::mesh::BulkData & mesh) const
{
  for (auto && myRefinedEdgeToChildNodes : myRefinedEdgesToChildNodes)
  {
    prolong_edge_node(mesh, get_edge_nodes(myRefinedEdgeToChildNodes.first), myRefinedEdgeToChildNodes.second);
  }

  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();
  const std::vector<const stk::mesh::FieldBase *> constFields(allFields.begin(), allFields.end());
  stk::mesh::communicate_field_data(mesh, constFields);
}

static
void pack_shared_edges_to_refine(const std::vector<std::pair<std::array<stk::mesh::EntityKey,2>,std::vector<int>>> & edgesNodeKeysAndSharingProcs,
    stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    //for (const auto & [edgeNodeKeys, sharingProcs] : edgesNodeKeysAndSharingProcs)
    for (const auto & entry : edgesNodeKeysAndSharingProcs)
    {
      const auto & edgeNodeKeys = std::get<0>(entry);
      const auto & sharingProcs = std::get<1>(entry);

      for (auto&& procId : sharingProcs)
      {
        if (procId != commSparse.parallel_rank())
        {
          commSparse.send_buffer(procId).pack(edgeNodeKeys[0]);
          commSparse.send_buffer(procId).pack(edgeNodeKeys[1]);
        }
      }
    }
  });
}

static
void unpack_shared_edges_to_refine(const stk::mesh::BulkData & mesh,
    typename NodeRefiner::RefinedEdgeMap & edgesToRefine,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      std::array<stk::mesh::EntityKey,2> edgeNodeKeys;
      commSparse.recv_buffer(procId).unpack(edgeNodeKeys[0]);
      commSparse.recv_buffer(procId).unpack(edgeNodeKeys[1]);
      stk::mesh::Entity node0 = mesh.get_entity(edgeNodeKeys[0]);
      if (mesh.is_valid(node0))
      {
        stk::mesh::Entity node1 = mesh.get_entity(edgeNodeKeys[1]);
        if (mesh.is_valid(node1))
          edgesToRefine.emplace(edge_from_edge_nodes(mesh, node0, node1), stk::mesh::Entity::InvalidEntity);
      }
    }
  });
}

static void fill_procs_that_own_elements_using_edge(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,2> & edgeNodes, std::vector<stk::mesh::Entity> & workspaceForEdgeElems, std::vector<int> & procsThatOwnElementsUsingEdge)
{
  stk::mesh::get_entities_through_relations(mesh, {edgeNodes[0], edgeNodes[1]}, stk::topology::ELEMENT_RANK, workspaceForEdgeElems);
  procsThatOwnElementsUsingEdge.clear();
  for (auto edgeElem : workspaceForEdgeElems)
    procsThatOwnElementsUsingEdge.push_back(mesh.parallel_owner_rank(edgeElem));
  stk::util::sort_and_unique(procsThatOwnElementsUsingEdge);
}

static
std::vector<std::pair<std::array<stk::mesh::EntityKey,2>,std::vector<int>>> get_shared_edges_and_sharing_procs(const stk::mesh::BulkData & mesh, const typename NodeRefiner::RefinedEdgeMap & edgesToRefine)
{
  ThrowAssert(mesh.is_automatic_aura_on()); // NOTE: Uses AURA to determine which procs have elements that use this edge and therefore should create the child node
  std::vector<int> sharingProcs;
  std::vector<std::pair<std::array<stk::mesh::EntityKey,2>,std::vector<int>>> edgesNodeKeysAndSharingProcs;
  std::vector<stk::mesh::Entity> edgeElems;
  for (auto && edgeToRefine : edgesToRefine)
  {
    const Edge edge = edgeToRefine.first;
    const std::array<stk::mesh::Entity,2> edgeNodes = get_edge_nodes(edge);
    if (mesh.bucket(edgeNodes[0]).shared() && mesh.bucket(edgeNodes[1]).shared())
    {
      fill_procs_that_own_elements_using_edge(mesh, edgeNodes, edgeElems, sharingProcs);
      edgesNodeKeysAndSharingProcs.emplace_back(std::array<stk::mesh::EntityKey,2>{mesh.entity_key(edgeNodes[0]), mesh.entity_key(edgeNodes[1])}, sharingProcs);
    }
  }
  return edgesNodeKeysAndSharingProcs;
}

void NodeRefiner::sync_shared_edges_from_other_procs_to_refine(const stk::mesh::BulkData & mesh)
{
  if (mesh.parallel_size() < 2) return;

  const auto edgesNodeKeysAndSharingProcs = get_shared_edges_and_sharing_procs(mesh, myRefinedEdgesToChildNodes);

  stk::CommSparse commSparse(mesh.parallel());
  pack_shared_edges_to_refine(edgesNodeKeysAndSharingProcs, commSparse);
  unpack_shared_edges_to_refine(mesh, myRefinedEdgesToChildNodes, commSparse);
}

void NodeRefiner::clear_edges_to_refine()
{
  myRefinedEdgesToChildNodes.clear();
}

bool NodeRefiner::mark_edge_for_refinement(const Edge & edge)
{
  auto result = myRefinedEdgesToChildNodes.emplace(edge, stk::mesh::Entity::InvalidEntity);
  return result.second;
}


bool NodeRefiner::mark_already_refined_edge(const Edge & edge, const stk::mesh::Entity refinedEdgeNode)
{
  auto result = myRefinedEdgesToChildNodes.emplace(edge, refinedEdgeNode);
  return result.second;
}

bool NodeRefiner::is_edge_marked_for_refinement(const Edge & edge) const
{
  return myRefinedEdgesToChildNodes.end() != myRefinedEdgesToChildNodes.find(edge);
}

}
