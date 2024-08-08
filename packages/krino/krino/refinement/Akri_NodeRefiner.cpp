#include <Akri_MeshHelpers.hpp>
#include "Akri_NodeRefiner.hpp"

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
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

void NodeRefiner::create_refined_element_centroid_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedElemCentroidNodeParts)
{
  if(stk::is_true_on_all_procs(mesh.parallel(), myRefinedElementsToChildNodes.empty()))
    return;

  std::sort(myRefinedElementsToChildNodes.begin(), myRefinedElementsToChildNodes.end());

  const size_t numElementsToRefine = myRefinedElementsToChildNodes.size();

  std::vector<stk::mesh::EntityId> newNodeIds;
  EntityIdPool::generate_new_ids(mesh, stk::topology::NODE_RANK, numElementsToRefine, newNodeIds, myAssert32Bit, myForce64Bit);

  for (size_t i=0; i<numElementsToRefine; ++i)
  {
    myRefinedElementsToChildNodes[i].second = mesh.declare_node(newNodeIds[i], refinedElemCentroidNodeParts);
  }
}

void NodeRefiner::assign_refined_quad_face_node_parent_ids(const stk::mesh::BulkData & mesh, FieldRef refinedQuadFaceNodeParentIdsField) const
{
  for (auto && refinedQuadFaceToChildNodes : myRefinedQuadFacesToChildNodes)
  {
    const auto & quadFaceNodes = get_quad_face_nodes_sorted_by_id(refinedQuadFaceToChildNodes.first);
    const stk::mesh::Entity refinedQuadFaceNode = refinedQuadFaceToChildNodes.second;
    auto * refinedQuadFaceNodeParentIds = field_data<uint64_t>(refinedQuadFaceNodeParentIdsField, refinedQuadFaceNode);
    refinedQuadFaceNodeParentIds[0] = mesh.identifier(quadFaceNodes[0]);
    refinedQuadFaceNodeParentIds[1] = mesh.identifier(quadFaceNodes[1]);
    refinedQuadFaceNodeParentIds[2] = mesh.identifier(quadFaceNodes[2]);
    refinedQuadFaceNodeParentIds[3] = mesh.identifier(quadFaceNodes[3]);
  }
}

void NodeRefiner::create_refined_quad_face_nodes(stk::mesh::BulkData & mesh, const stk::mesh::PartVector & refinedQuadFaceNodeParts, FieldRef refinedQuadFaceNodeParentIdsField)
{
  if(stk::is_true_on_all_procs(mesh.parallel(), myRefinedQuadFacesToChildNodes.empty()))
    return;

  const size_t numQuadFacesToRefine = myRefinedQuadFacesToChildNodes.size();
  std::vector<std::array<stk::mesh::Entity,4>> quadFacesParentNodes;
  quadFacesParentNodes.reserve(numQuadFacesToRefine);
  std::vector<ChildNodeRequest> childNodeRequests;
  childNodeRequests.reserve(numQuadFacesToRefine);
  for (auto && refinedQuadFaceToChildNodes : myRefinedQuadFacesToChildNodes)
  {
    quadFacesParentNodes.emplace_back(get_quad_face_nodes_sorted_by_id(refinedQuadFaceToChildNodes.first));
    auto & quadFaceParentNodes = quadFacesParentNodes.back();
    childNodeRequests.emplace_back(std::vector<stk::mesh::Entity*>{&(quadFaceParentNodes[0]), &(quadFaceParentNodes[1]), &(quadFaceParentNodes[2]), &(quadFaceParentNodes[3])}, &(refinedQuadFaceToChildNodes.second));
  }

  auto generate_new_ids = [&](stk::topology::rank_t entityRank, size_t numIdsNeeded, std::vector<stk::mesh::EntityId>& requestedIds)
  {
     EntityIdPool::generate_new_ids(mesh, entityRank, numIdsNeeded, requestedIds, myAssert32Bit, myForce64Bit);
  };

  batch_create_child_nodes(mesh, childNodeRequests, refinedQuadFaceNodeParts, generate_new_ids);

  assign_refined_quad_face_node_parent_ids(mesh, refinedQuadFaceNodeParentIdsField);
}

stk::mesh::Entity NodeRefiner::get_edge_child_node(const Edge edge) const
{
  const auto iter = myRefinedEdgesToChildNodes.find(edge);
  if (iter != myRefinedEdgesToChildNodes.end())
    return iter->second;
  return stk::mesh::Entity();
}

stk::mesh::Entity NodeRefiner::get_element_centroid_child_node(const stk::mesh::Entity elem) const
{
  auto lb_cmp = [](const std::pair<stk::mesh::Entity,stk::mesh::Entity> & elemAndNode, stk::mesh::Entity searchElem) { return elemAndNode.first < searchElem; };
  const auto iter = std::lower_bound(myRefinedElementsToChildNodes.begin(), myRefinedElementsToChildNodes.end(), elem, lb_cmp);
  if (iter != myRefinedElementsToChildNodes.end())
    return iter->second;
  return stk::mesh::Entity();
}

stk::mesh::Entity NodeRefiner::get_element_child_face_node(const QuadFace quadFace) const
{
  const auto iter = myRefinedQuadFacesToChildNodes.find(quadFace);
  if (iter != myRefinedQuadFacesToChildNodes.end())
    return iter->second;
  return stk::mesh::Entity();
}

std::vector<stk::mesh::Entity> NodeRefiner::get_element_child_face_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity elem) const
{
  const stk::topology elemTopology = mesh.bucket(elem).topology();
  STK_ThrowAssert(elemTopology == stk::topology::HEXAHEDRON_8);
  const stk::mesh::Entity* elemNodes = mesh.begin_nodes(elem);

  constexpr unsigned numSides = 6;
  constexpr unsigned nodesPerSide = 4;
  std::array<stk::mesh::Entity,nodesPerSide> elemSideNodes;

  std::vector<stk::mesh::Entity> childFaceNodes;
  childFaceNodes.reserve(numSides);
  for (unsigned iSide=0; iSide<numSides; ++iSide)
  {
    elemTopology.side_nodes(elemNodes, iSide, elemSideNodes.data());
    const QuadFace quadFace = quad_face_from_unordered_nodes(mesh, elemSideNodes);
    childFaceNodes.push_back(get_element_child_face_node(quadFace));
  }

  return childFaceNodes;
}

template<typename NODECONTAINER>
static void prolong_child_node(const stk::mesh::BulkData & mesh, const NODECONTAINER & parentNodes, const stk::mesh::Entity childNode)
{
  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();
  for ( auto && stkField : allFields )
  {
    const FieldRef field(stkField);

    if( field.entity_rank() == stk::topology::NODE_RANK && field.type_is<double>() )
    {
      auto * childNodeData = field_data<double>(field, childNode);

      if (nullptr != childNodeData)
      {
        const unsigned fieldLength = field.length();
        for (unsigned i=0; i<fieldLength; ++i)
          childNodeData[i] = 0.;

        for (auto parentNode : parentNodes)
        {
          const auto * parentNodeData = field_data<double>(field, parentNode);
          STK_ThrowAssertMsg(parentNodeData, "Child centroid node has field " << field.name() << " but the parent node " << mesh.identifier(parentNode) << " does not.");
          for (unsigned i=0; i<fieldLength; ++i)
            childNodeData[i] += parentNodeData[i]/parentNodes.size();
        }
      }
    }
  }
}

void NodeRefiner::prolong_refined_nodes(const stk::mesh::BulkData & mesh) const
{
  for (auto && refinedEdgeToChildNodes : myRefinedEdgesToChildNodes)
  {
    prolong_child_node(mesh, get_edge_nodes(refinedEdgeToChildNodes.first), refinedEdgeToChildNodes.second);
  }

  std::vector<stk::mesh::Entity> elemNodes;
  for (auto & refinedElementToChildNodes : myRefinedElementsToChildNodes)
  {
    const stk::mesh::Entity elem = refinedElementToChildNodes.first;
    elemNodes.assign(mesh.begin_nodes(elem), mesh.end_nodes(elem));
    prolong_child_node(mesh, elemNodes, refinedElementToChildNodes.second);
  }

  for (auto & refinedQuadFacesToChildNodes : myRefinedQuadFacesToChildNodes)
  {
    prolong_child_node(mesh, get_quad_face_nodes_sorted_by_id(refinedQuadFacesToChildNodes.first), refinedQuadFacesToChildNodes.second);
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
    for (const auto & [edgeNodeKeys, sharingProcs] : edgesNodeKeysAndSharingProcs)
    {
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

static bool have_locally_owned_element_using_edge_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node0, const stk::mesh::Entity node1, std::vector<stk::mesh::Entity> & workspaceForEdgeElems)
{
  stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{node0, node1}, stk::topology::ELEMENT_RANK, workspaceForEdgeElems);
  for (auto & elem : workspaceForEdgeElems)
    if (mesh.bucket(elem).owned())
      return true;
  return false;
}

static
void unpack_shared_edges_to_refine_if_there_is_a_locally_owned_element_using_edge(const stk::mesh::BulkData & mesh,
    typename NodeRefiner::RefinedEdgeMap & edgesToRefine,
    stk::CommSparse &commSparse)
{
  std::vector<stk::mesh::Entity> workspaceForEdgeElems;
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
        if (mesh.is_valid(node1) && have_locally_owned_element_using_edge_nodes(mesh, node0, node1, workspaceForEdgeElems))
          edgesToRefine.emplace(edge_from_edge_nodes(mesh, node0, node1), stk::mesh::Entity::InvalidEntity);
      }
    }
  });
}

static void fill_procs_that_might_own_elements_using_edge(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,2> & edgeNodes, std::vector<int> & procsThatMightOwnElementsUsingEdge)
{
  mesh.shared_procs_intersection({edgeNodes[0], edgeNodes[1]}, procsThatMightOwnElementsUsingEdge);
}

static
std::vector<std::pair<std::array<stk::mesh::EntityKey,2>,std::vector<int>>> get_shared_edges_and_sharing_procs(const stk::mesh::BulkData & mesh, const typename NodeRefiner::RefinedEdgeMap & edgesToRefine)
{
  std::vector<int> sharingProcs;
  std::vector<std::pair<std::array<stk::mesh::EntityKey,2>,std::vector<int>>> edgesNodeKeysAndSharingProcs;

  for (auto && edgeToRefine : edgesToRefine)
  {
    const Edge edge = edgeToRefine.first;
    const std::array<stk::mesh::Entity,2> edgeNodes = get_edge_nodes(edge);
    if (mesh.bucket(edgeNodes[0]).shared() && mesh.bucket(edgeNodes[1]).shared())
    {
      fill_procs_that_might_own_elements_using_edge(mesh, edgeNodes, sharingProcs);
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
  unpack_shared_edges_to_refine_if_there_is_a_locally_owned_element_using_edge(mesh, myRefinedEdgesToChildNodes, commSparse);
}

void NodeRefiner::clear_entities_to_refine()
{
  myRefinedEdgesToChildNodes.clear();
  myRefinedElementsToChildNodes.clear();
  myRefinedQuadFacesToChildNodes.clear();
}

bool NodeRefiner::mark_edge_for_refinement(const Edge & edge)
{
  auto result = myRefinedEdgesToChildNodes.emplace(edge, stk::mesh::Entity::InvalidEntity);
  return result.second;
}

bool NodeRefiner::mark_quad_face_for_refinement(const QuadFace & quadFace)
{
  auto result = myRefinedQuadFacesToChildNodes.emplace(quadFace, stk::mesh::Entity::InvalidEntity);
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

void NodeRefiner::mark_element_with_child_centroid_node_for_refinement(const stk::mesh::Entity & elem)
{
  myRefinedElementsToChildNodes.emplace_back(elem, stk::mesh::Entity::InvalidEntity);
}


}
