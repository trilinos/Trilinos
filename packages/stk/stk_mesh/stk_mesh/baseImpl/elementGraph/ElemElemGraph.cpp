#include "ElemElemGraph.hpp"
#include "ElemFilter.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "ElemGraphShellConnections.hpp"
#include "BulkDataIdMapper.hpp"

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/DeletedElementInfo.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk { namespace mesh {

void convert_keys_to_entities(stk::mesh::BulkData& bulkData,
                              const std::vector<stk::mesh::EntityKey>& keys,
                              stk::mesh::EntityVector& entities)
{
    entities.resize(keys.size());
    for(size_t i=0; i<keys.size(); ++i) {
        entities[i] = bulkData.get_entity(keys[i]);
    }
}

void SharedSidesCommunication::pack_shared_side_nodes_of_elements(stk::CommSparse &comm) const
{
    impl::ElemSideProcVector::const_iterator iter = m_elementSidesToSend.begin();
    impl::ElemSideProcVector::const_iterator end = m_elementSidesToSend.end();
    stk::mesh::EntityVector side_nodes;
    std::vector<stk::mesh::EntityKey> side_node_entity_keys;
    for(; iter!= end; ++iter)
    {
        stk::mesh::Entity elem = iter->elemSidePair.entity;
        unsigned side_index    = iter->elemSidePair.side_id;
        int sharing_proc       = iter->proc;
        stk::mesh::EntityId element_id     = m_bulkData.identifier(elem);
        stk::topology topology = m_bulkData.bucket(elem).topology();

        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(element_id);
        comm.send_buffer(sharing_proc).pack<stk::topology>(topology);
        comm.send_buffer(sharing_proc).pack<unsigned>(side_index);

        const stk::mesh::Entity* elemNodes = m_bulkData.begin_nodes(elem);
        impl:: fill_element_side_nodes_from_topology(topology, elemNodes, side_index, side_nodes);
        side_node_entity_keys.resize(side_nodes.size());
        for(size_t i=0; i<side_nodes.size(); ++i)
        {
            side_node_entity_keys[i] = m_bulkData.entity_key(side_nodes[i]);
        }
        stk::pack_vector_to_proc(comm, side_node_entity_keys, sharing_proc);
    }
}

void SharedSidesCommunication::unpack_side_data_map(SideNodeToReceivedElementDataMap& elementSidesReceived)
{
    std::vector<stk::mesh::EntityKey> node_keys;
    stk::mesh::impl::ParallelElementData elementData;
    stk::mesh::EntityVector sideNodes;
    for(int proc_id = 0; proc_id < m_bulkData.parallel_size(); ++proc_id)
    {
        if(proc_id != m_bulkData.parallel_rank())
        {
            while(m_commSparse.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityId elementIdentifier;
                stk::topology topology;
                unsigned side_index = 0;
                m_commSparse.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(elementIdentifier);
                m_commSparse.recv_buffer(proc_id).unpack<stk::topology>(topology);
                m_commSparse.recv_buffer(proc_id).unpack<unsigned>(side_index);

                elementData.set_proc_rank(proc_id);
                elementData.set_element_identifier(elementIdentifier);
                elementData.set_element_topology(topology);
                elementData.set_element_side_index(side_index);

                stk::unpack_vector_from_proc(m_commSparse, node_keys, proc_id);
                convert_keys_to_entities(m_bulkData, node_keys, sideNodes);
                elementData.set_element_side_nodes(sideNodes);
                std::sort(sideNodes.begin(), sideNodes.end());
                elementSidesReceived[sideNodes].push_back(elementData);
            }
        }
    }
}

void SharedSidesCommunication::unpack_side_data_vec(impl::ParallelElementDataVector& elementSidesReceived)
{
    std::vector<stk::mesh::EntityKey> node_keys;
    stk::mesh::impl::ParallelElementData elementData;
    stk::mesh::EntityVector sideNodes;
    for(int proc_id = 0; proc_id < m_bulkData.parallel_size(); ++proc_id)
    {
        if(proc_id != m_bulkData.parallel_rank())
        {
            while(m_commSparse.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityId elementIdentifier;
                stk::topology topology;
                unsigned side_index = 0;
                m_commSparse.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(elementIdentifier);
                m_commSparse.recv_buffer(proc_id).unpack<stk::topology>(topology);
                m_commSparse.recv_buffer(proc_id).unpack<unsigned>(side_index);

                elementData.set_proc_rank(proc_id);
                elementData.set_element_identifier(elementIdentifier);
                elementData.set_element_topology(topology);
                elementData.set_element_side_index(side_index);

                stk::unpack_vector_from_proc(m_commSparse, node_keys, proc_id);
                convert_keys_to_entities(m_bulkData, node_keys, sideNodes);
                elementData.set_element_side_nodes(sideNodes);
                elementSidesReceived.push_back(elementData);
            }
        }
    }
}



void SharedSidesCommunication::communicate_element_sides()
{
  stk::CommSparse& commSparse = m_commSparse;
  stk::pack_and_communicate(commSparse, [this,&commSparse](){pack_shared_side_nodes_of_elements(commSparse);} );
}

ElemElemGraph::ElemElemGraph(stk::mesh::BulkData& bulkData) :
        m_bulk_data(bulkData),
        m_parallelInfoForGraphEdges(bulkData.parallel_rank()),
        m_modCycleWhenGraphModified(0),
        m_any_shell_elements_exist(false),
        m_idMapper(),
        m_sideConnector(m_bulk_data, m_graph, m_coincidentGraph, m_idMapper),
        m_sideNodeConnector(m_bulk_data, m_graph, m_coincidentGraph, m_parallelInfoForGraphEdges, m_idMapper)
{
    fill_from_mesh();
}

void ElemElemGraph::fill_from_mesh()
{
    clear_data_members();

    impl::ElemSideProcVector elementSideIdsToSend;

    int numElems = size_data_members();
    if (numElems > 0) {
        m_any_shell_elements_exist = impl::fill_topologies(*this, m_idMapper, m_element_topologies);
        elementSideIdsToSend = impl::gather_element_side_ids_to_send(m_bulk_data);
        std::sort(elementSideIdsToSend.begin(), elementSideIdsToSend.end());
    }

    m_edgesToAdd.clear();

    m_parallelInfoForGraphEdges.clear();
    impl::ParallelElementDataVector elementSidesReceived;
    communicate_shared_sides(elementSideIdsToSend, elementSidesReceived);
    fill_parallel_graph(elementSideIdsToSend, elementSidesReceived);

    if (numElems > 0) {
        fill_graph();
    }

    stk::util::sort_and_unique(m_edgesToAdd, GraphEdgeLessByElem1());
    m_graph.replace_sorted_edges(m_edgesToAdd);

    GraphInfo graphInfo(m_graph, m_parallelInfoForGraphEdges, m_element_topologies);
    remove_graph_edges_blocked_by_shell(graphInfo);

    m_graph.delete_sorted_edges(m_edgesToDelete);

    {
      std::vector<GraphEdge> emptyDelete;
      m_edgesToDelete.swap(emptyDelete);
    }

    m_modCycleWhenGraphModified = m_bulk_data.synchronized_count();
}

size_t ElemElemGraph::num_parallel_edges() const
{
  size_t numParallelEdges = 0;
  for (size_t i = 0; i < m_graph.get_num_elements_in_graph(); ++i) {
    const GraphEdgesForElement graphEdges = m_graph.get_edges_for_element(i);
    for (const GraphEdge& graphEdge : graphEdges) {
      if (!graphEdge.is_elem2_local()) {
        ++numParallelEdges;
      }
    }
  }
  return numParallelEdges;
}

ElemElemGraph::~ElemElemGraph() {}

size_t ElemElemGraph::get_num_connected_elems(stk::mesh::Entity local_element) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return m_graph.get_num_edges_for_element(local_id);
}

bool ElemElemGraph::is_connected_elem_locally_owned(stk::mesh::Entity localElement, size_t indexConnElement) const
{
    impl::LocalId local_id = get_local_element_id(localElement);
    return m_graph.get_edge_for_element(local_id, indexConnElement).elem2() >= 0;
}

impl::ElementViaSidePair ElemElemGraph::get_connected_element_and_via_side(stk::mesh::Entity localElement, size_t indexConnElement) const
{
    impl::LocalId local_id = get_local_element_id(localElement);
    const GraphEdge & graphEdge = m_graph.get_edge_for_element(local_id, indexConnElement);
    impl::LocalId other_element_id = graphEdge.elem2();
    return {m_idMapper.local_to_entity(other_element_id),graphEdge.side1()};
}

impl::IdViaSidePair ElemElemGraph::get_connected_remote_id_and_via_side(stk::mesh::Entity localElement, size_t indexConnElement) const
{
    STK_ThrowRequireWithSierraHelpMsg(!is_connected_elem_locally_owned(localElement, indexConnElement));
    impl::LocalId local_id = get_local_element_id(localElement);
    const GraphEdge & graphEdge = m_graph.get_edge_for_element(local_id, indexConnElement);
    return {m_parallelInfoForGraphEdges.convert_negative_local_id_to_remote_global_id(graphEdge.elem2()), graphEdge.side1()};
}

int ElemElemGraph::get_connected_elements_side(stk::mesh::Entity localElement, size_t indexConnElement) const
{
    impl::LocalId local_id = get_local_element_id(localElement);
    const GraphEdge & graphEdge = m_graph.get_edge_for_element(local_id, indexConnElement);
    return graphEdge.side2();
}

int ElemElemGraph::get_owning_proc_id_of_remote_element(stk::mesh::Entity localElement, size_t indexConnElement) const
{
    impl::LocalId local_id = get_local_element_id(localElement);
    const GraphEdge & graphEdge = m_graph.get_edge_for_element(local_id, indexConnElement);
    const impl::ParallelInfo& parallelInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    return parallelInfo.get_proc_rank_of_neighbor();
}

bool ElemElemGraph::is_connected_to_other_element_via_side_ordinal(stk::mesh::Entity element, int sideOrdinal) const
{
    impl::LocalId elementLocalId = get_local_element_id(element);
    std::vector<GraphEdge> graphEdges = m_graph.get_edges_for_element_side(elementLocalId, sideOrdinal);
    return !graphEdges.empty();
}

size_t ElemElemGraph::num_edges() const
{
    return m_graph.get_num_edges();
}

impl::ParallelInfo& ElemElemGraph::get_parallel_edge_info(stk::mesh::Entity element, int side1, stk::mesh::EntityId remote_id, int side2)
{
    impl::LocalId this_elem_local_id = get_local_element_id(element);
    impl::LocalId elem2 = m_parallelInfoForGraphEdges.convert_remote_global_id_to_negative_local_id(remote_id);
    return m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(GraphEdge(this_elem_local_id, side1, elem2, side2));
}

const impl::ParallelInfo& ElemElemGraph::get_const_parallel_edge_info(stk::mesh::Entity element, int side1, stk::mesh::EntityId remote_id, int side2) const
{
    impl::LocalId this_elem_local_id = get_local_element_id(element);
    impl::LocalId elem2 = m_parallelInfoForGraphEdges.convert_remote_global_id_to_negative_local_id(remote_id);
    return m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(GraphEdge(this_elem_local_id, side1, elem2, side2));
}

impl::LocalId ElemElemGraph::get_local_element_id(stk::mesh::Entity local_element, bool require_valid_id) const
{
    impl::LocalId local_id = m_idMapper.entity_to_local(local_element);
    if (require_valid_id)
    {
        STK_ThrowRequireWithSierraHelpMsg(local_id != impl::INVALID_LOCAL_ID);
    }
    return local_id;
}

int ElemElemGraph::size_data_members()
{
    unsigned numElems = impl::get_num_local_elems(m_bulk_data);
    m_graph.set_num_local_elements(numElems);
    m_element_topologies.resize(numElems);

    m_idMapper.initialize(m_bulk_data);

    return numElems;
}

void ElemElemGraph::clear_data_members()
{
    m_graph.clear();
    m_edgesToAdd.clear();
    m_edgesToDelete.clear();
    m_idMapper.clear();
    m_element_topologies.clear();
    m_any_shell_elements_exist = false;
    m_parallelInfoForGraphEdges.clear();
    m_modCycleWhenGraphModified = m_bulk_data.synchronized_count();
    m_deleted_element_local_id_pool.clear();
    m_deleted_elem_pool.clear();
    m_coincidentGraph.clear();
}

void ElemElemGraph::fill_elements_attached_to_local_nodes(const stk::mesh::EntityVector& sideNodesOfReceivedElement,
                                                                                  stk::mesh::EntityId elementId,
                                                                                  stk::topology elementTopology,
                                                                                  impl::SerialElementDataVector& connectedElementDataVector) const
{
    connectedElementDataVector.clear();
    impl::get_elements_with_larger_ids_connected_via_sidenodes<impl::SerialElementData>(m_bulk_data, elementId, elementTopology, m_idMapper, sideNodesOfReceivedElement,
                                                                                        m_scratchEntityVector, connectedElementDataVector);
}

void ElemElemGraph::get_elements_attached_to_remote_nodes(const stk::mesh::EntityVector& sideNodesOfReceivedElement, stk::mesh::EntityId elementId, stk::topology elementTopology,
                                                                                     impl::ParallelElementDataVector& connectedElementDataVector) const
{
    connectedElementDataVector.clear();
    impl::get_elements_connected_via_sidenodes<impl::ParallelElementData>(m_bulk_data, elementId, elementTopology, m_idMapper,
                                                                                 sideNodesOfReceivedElement, m_scratchEntityVector, connectedElementDataVector);
}

bool is_coincident(const stk::mesh::impl::TopologyChecker& topoChecker,
                   const unsigned sideIndex,
                   const stk::mesh::Permutation& elemPerm,
                   const unsigned otherSideIndex,
                   const stk::mesh::Permutation& otherElemPerm)
{
  if (topoChecker.are_both_shells()) {
    return true;
  }
  if (topoChecker.are_both_not_shells()) {
    const bool elemIsPos = topoChecker.localTopology.side_topology(sideIndex).is_positive_polarity(elemPerm);
    const bool otherElemIsPos = topoChecker.remoteTopology.side_topology(otherSideIndex).is_positive_polarity(otherElemPerm);
    return elemIsPos == otherElemIsPos;
  }
  return false;
}

void ElemElemGraph::insert_edge_between_elements(impl::LocalId local_elem_id,
                                                 int side_index,
                                                 const impl::SerialElementData& otherElem,
                                                 std::vector<stk::mesh::GraphEdge>& graphEdges) const
{
    graphEdges.emplace_back(local_elem_id, side_index, otherElem.get_element_local_id(), otherElem.get_element_side_index());
}

void ElemElemGraph::add_local_graph_edges_for_elem(const stk::mesh::MeshIndex &meshIndex,
                                                   impl::LocalId local_elem_id,
                                                   std::vector<stk::mesh::GraphEdge> &graphEdges,
                                                   std::vector<stk::mesh::GraphEdge> &coincidentGraphEdges,
                                                   stk::mesh::EntityVector& side_nodes,
                                                   impl::SerialElementDataVector& connectedElementDataVector,
                                                   bool only_consider_upper_symmetry) const
{
    stk::mesh::Entity element = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
    stk::topology elemTopology = meshIndex.bucket->topology();
    int num_sides = elemTopology.num_sides();
    graphEdges.clear();
    coincidentGraphEdges.clear();
    stk::mesh::EntityId elemGlobalId = m_bulk_data.identifier(element);
    if(!only_consider_upper_symmetry) elemGlobalId = 0;
    const stk::mesh::Entity* elemNodes = m_bulk_data.begin_nodes(element);
    for(int side_index=0; side_index<num_sides; ++side_index)
    {
        impl::fill_element_side_nodes_from_topology(elemTopology, elemNodes, side_index, side_nodes);
        fill_elements_attached_to_local_nodes(side_nodes, elemGlobalId, elemTopology, connectedElementDataVector);
        const stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
        for (const impl::SerialElementData & otherElem: connectedElementDataVector)
        {
            if(local_elem_id != otherElem.get_element_local_id())
            {
                stk::mesh::impl::TopologyChecker topoChecker{elemTopology, otherElem.get_element_topology()};
                
                const stk::mesh::Permutation otherPerm = otherElem.get_permutation();
                const bool isCoincidentConnection = is_coincident(topoChecker, side_index, perm, otherElem.get_element_side_index(), otherPerm);
                if(isCoincidentConnection)
                    insert_edge_between_elements(local_elem_id, side_index, otherElem, coincidentGraphEdges);
                else
                    insert_edge_between_elements(local_elem_id, side_index, otherElem, graphEdges);
            }
        }
    }
    stk::util::sort_and_unique(graphEdges, GraphEdgeLessByElem2());
}

class GraphEdgeAdder
{
public:
  GraphEdgeAdder(std::vector<GraphEdge>& edges) : m_edges(edges) {}
  void add_edge(const GraphEdge& edge) { m_edges.emplace_back(edge); }
private:
  std::vector<GraphEdge>& m_edges;
};

void ElemElemGraph::fill_graph()
{
    const stk::mesh::BucketVector& elemBuckets = m_bulk_data.get_buckets(stk::topology::ELEM_RANK, m_bulk_data.mesh_meta_data().locally_owned_part());
    stk::mesh::EntityVector side_nodes;
    impl::SerialElementDataVector connectedElementDataVector;
    const bool only_consider_upper_symmetry = true;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            Entity elemToAdd = bucket[j];
            impl::LocalId local_elem_id = get_local_element_id(elemToAdd);
            add_local_graph_edges_for_elem(m_bulk_data.mesh_index(elemToAdd), local_elem_id, m_newGraphEdges, m_coincidentGraphEdges, side_nodes, connectedElementDataVector, only_consider_upper_symmetry);
            GraphEdgeAdder graphEdgeAdder(m_edgesToAdd);
            add_new_local_edges_to_graph(graphEdgeAdder, m_newGraphEdges);
            add_new_local_edges_to_graph(m_coincidentGraph, m_coincidentGraphEdges);
        }
    }

    stk::util::sort_and_unique(m_edgesToAdd, GraphEdgeLessByElem1());
}

void ElemElemGraph::communicate_shared_sides(impl::ElemSideProcVector& elementSidesToSend, impl::ParallelElementDataVector& elementSidesReceived)
{
     SharedSidesCommunication sharedSidesCommunication(m_bulk_data, elementSidesToSend);
     sharedSidesCommunication.communicate_element_sides();
     sharedSidesCommunication.unpack_side_data_vec(elementSidesReceived);
}

void add_shared_edge(const impl::ParallelElementData& elementDataOtherProc, stk::mesh::EntityId elem_id, const unsigned side_index,
        const stk::mesh::EntityVector& localElemSideNodes, stk::topology elem_topology, std::vector<impl::SharedEdgeInfo> &newlySharedEdges)
{
    impl::SharedEdgeInfo sharedEdgeInfo;
    impl::GraphEdgeProc graphEdgeProc(elem_id, side_index, elementDataOtherProc.get_element_identifier(), elementDataOtherProc.get_element_side_index(),
            elementDataOtherProc.get_proc_rank_of_neighbor());
    sharedEdgeInfo.set_graph_edge_proc(graphEdgeProc);

    sharedEdgeInfo.m_sharedNodes = localElemSideNodes;
    sharedEdgeInfo.m_remoteElementTopology = elem_topology;
    newlySharedEdges.push_back(sharedEdgeInfo);
}

std::string ElemElemGraph::print_edge(const GraphEdge& graphEdge)
{
  std::ostringstream os;
  Entity localElem = get_entity(graphEdge.elem1());

  stk::topology localTopology, otherTopology;
  if(m_bulk_data.is_valid(localElem)) {
    localTopology = m_bulk_data.bucket(localElem).topology();
  }

  if(!impl::is_local_element(graphEdge.elem2())) {
    auto iter = m_parallelInfoForGraphEdges.get_parallel_info_iterator_for_graph_edge(graphEdge);
    if(iter != m_parallelInfoForGraphEdges.get_parallel_graph_info().end()) {
      otherTopology = iter->second.m_remote_element_topology;
    }
    os << "ParallelGraphEdge{" << m_bulk_data.identifier(localElem) << "|" << graphEdge.side1() << "|" << localTopology << " , "
                               << -graphEdge.elem2()                << "|" << graphEdge.side2() << "|" << otherTopology << "}";
  } else {
    Entity otherElem = get_entity(graphEdge.elem2());
    if(m_bulk_data.is_valid(otherElem)) {
      otherTopology = m_bulk_data.bucket(otherElem).topology();
    }
    os << "LocalGraphEdge{"    << m_bulk_data.identifier(localElem) << "|" << graphEdge.side1() << "|" << localTopology << " , "
                               << m_bulk_data.identifier(otherElem) << "|" << graphEdge.side2() << "|" << otherTopology << "}";
  }

  return os.str();
}

void ElemElemGraph::create_parallel_graph_edge(const impl::ParallelElementData &localElementData,
                                               const stk::mesh::impl::ParallelElementData &remoteElementData,
                                               impl::ElemSideProcVector & elementSidesToSend,
                                               std::vector<impl::SharedEdgeInfo> &newlySharedEdges,
                                               impl::ParallelGraphInfo& newParallelGraphInfo)
{
    impl::LocalId local_elem_id = localElementData.get_element_local_id();
    impl::LocalId negativeRemoteElemId = -1 * static_cast<impl::LocalId>(remoteElementData.get_element_identifier());
    stk::mesh::GraphEdge graphEdge(local_elem_id, localElementData.get_element_side_index(), negativeRemoteElemId, remoteElementData.get_element_side_index());

    impl::ParallelInfo parInfo(remoteElementData.get_proc_rank_of_neighbor(),
                               localElementData.get_permutation(),
                                remoteElementData.get_element_topology());

    bool parallelEdgeAlreadyExists = m_parallelInfoForGraphEdges.find_parallel_info_for_graph_edge(graphEdge);

    if(parallelEdgeAlreadyExists) {
      return;
    }

    m_edgesToAdd.push_back(graphEdge);
    newParallelGraphInfo.push_back(std::make_pair(graphEdge, parInfo));

    stk::mesh::Entity localElem = get_entity(local_elem_id);
    STK_ThrowRequireWithSierraHelpMsg(m_bulk_data.is_valid(localElem));

    const bool did_this_proc_send_info_about_this_side_to_other_proc = std::binary_search(elementSidesToSend.begin(), elementSidesToSend.end(), impl::EntitySidePair(localElem, localElementData.get_element_side_index()));
    // communicate back to originating processor
    if(!did_this_proc_send_info_about_this_side_to_other_proc && !parallelEdgeAlreadyExists)
    {
        // for cases of mesh modification, only proc with "created" element will communicate
        // new side was sent to this proc (to fish for possible face connections)
        // if fish is caught, need to send info back to original proc so they both will create edge in graph
        add_shared_edge(remoteElementData, m_bulk_data.identifier(localElem), localElementData.get_element_side_index(), localElementData.get_side_nodes(), m_bulk_data.bucket(localElem).topology(), newlySharedEdges);
    }
}

void ElemElemGraph::fill_parallel_graph(impl::ElemSideProcVector & elementSidesToSend, impl::ParallelElementDataVector& elementSidesReceived)
{
    std::vector<impl::SharedEdgeInfo> newlySharedEdges;
    impl::ParallelElementDataVector localElementsAttachedToReceivedNodes;
    impl::ParallelGraphInfo newParallelGraphInfo;
    for(const stk::mesh::impl::ParallelElementData &remoteElementData : elementSidesReceived)
    {
         get_elements_attached_to_remote_nodes(remoteElementData.get_side_nodes(), remoteElementData.get_element_identifier(),
                                               remoteElementData.get_element_topology(), localElementsAttachedToReceivedNodes);
        for(const impl::ParallelElementData &localElemAttachedToNodes : localElementsAttachedToReceivedNodes)
            create_parallel_graph_edge(localElemAttachedToNodes, remoteElementData, elementSidesToSend, newlySharedEdges, newParallelGraphInfo);
    }

    {
      impl::ElemSideProcVector().swap(elementSidesToSend);
      impl::ParallelElementDataVector().swap(elementSidesReceived);
    }

    stk::util::sort_and_unique(newParallelGraphInfo, GraphEdgeLessByElem2());
    m_parallelInfoForGraphEdges.insert_sorted_edges(newParallelGraphInfo);

    std::vector<impl::SharedEdgeInfo> receivedSharedEdges;
    communicate_remote_edges_for_pre_existing_graph_nodes(newlySharedEdges, receivedSharedEdges);

    for (const impl::SharedEdgeInfo &receivedSharedEdge : receivedSharedEdges)
        connect_remote_element_to_existing_graph(receivedSharedEdge);
}

void ElemElemGraph::write_graph(std::ostream& out, const std::string preamble) const
{
    std::ostringstream os;
    os << preamble;
    os << "Graph for processor " << m_bulk_data.parallel_rank() << std::endl;
    for(size_t localId=0;localId<this->size();++localId)
    {
        stk::mesh::Entity e = m_idMapper.local_to_entity(localId);
        os << "Element " << m_bulk_data.identifier(e) << " has connections: ";
        for(const stk::mesh::GraphEdge &graphEdge : m_graph.get_edges_for_element(localId))
            write_graph_edge(os, graphEdge);
        os << "Coincident Elements:  ";

        const std::vector<GraphEdge> &coincidentEdges = get_coincident_edges_for_element(localId);

        for(const stk::mesh::GraphEdge &graphEdge : coincidentEdges)
            write_graph_edge(os, graphEdge);

        if(coincidentEdges.empty())
            os << "none";
        os << std::endl;
    }
    os << "Parallel Info for Graph Edges:  " << std::endl;
    const impl::ParallelGraphInfo& parallelGraphInfo = const_cast<ParallelInfoForGraphEdges&>(m_parallelInfoForGraphEdges).get_parallel_graph_info();
    for (const impl::ParallelGraphInfo::value_type& iter : parallelGraphInfo)
    {
        write_graph_edge(os, iter.first);
        os << " parallel info = " << iter.second << std::endl;
    }
    out << os.str();
}


void ElemElemGraph::write_graph_edge(std::ostringstream& os,
                                     const stk::mesh::GraphEdge& graphEdge) const
{
    stk::mesh::EntityId elem2EntityId;
    std::string remoteMarker = "";
    if(graphEdge.elem2() < 0)
    {
        elem2EntityId = -1 * graphEdge.elem2();
        remoteMarker = "-";
    }
    else
    {
        elem2EntityId = m_bulk_data.identifier(m_idMapper.local_to_entity(graphEdge.elem2()));
    }
    stk::mesh::Entity elem1 = m_idMapper.local_to_entity(graphEdge.elem1());
    stk::mesh::EntityId elem1EntityId = m_bulk_data.identifier(elem1);
    os << "("                 << elem1EntityId << ", " << graphEdge.side1() << ")->";
    os << "(" << remoteMarker << elem2EntityId << ", " << graphEdge.side2() << ") ";
}

void ElemElemGraph::communicate_remote_edges_for_pre_existing_graph_nodes(const std::vector<impl::SharedEdgeInfo> &newlySharedEdges,
                                                                          std::vector<impl::SharedEdgeInfo> &receivedSharedEdges)
{
    stk::CommSparse comm(m_bulk_data.parallel());

    impl::pack_newly_shared_remote_edges(comm, m_bulk_data, newlySharedEdges);
    comm.allocate_buffers();

    impl::pack_newly_shared_remote_edges(comm, m_bulk_data, newlySharedEdges);
    comm.communicate();

    for(int proc_id=0; proc_id<m_bulk_data.parallel_size(); ++proc_id)
    {
        if (proc_id != m_bulk_data.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                impl::SharedEdgeInfo remoteEdgeInfo;

                stk::mesh::EntityId remoteElementId, localElementId;
                int remoteSide = -1, localSide = -1;

                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(remoteElementId);
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(localElementId);
                comm.recv_buffer(proc_id).unpack<int>(remoteSide);
                comm.recv_buffer(proc_id).unpack<int>(localSide);

                impl::GraphEdgeProc graphEdgeProc(localElementId, localSide, remoteElementId, remoteSide, proc_id);
                remoteEdgeInfo.set_graph_edge_proc(graphEdgeProc);

                comm.recv_buffer(proc_id).unpack<stk::topology>(remoteEdgeInfo.m_remoteElementTopology);
                unsigned numNodes = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(numNodes);
                stk::mesh::EntityVector nodes;
                for(unsigned i=0; i<numNodes; ++i)
                {
                    stk::mesh::EntityKey key;
                    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(key);
                    nodes.push_back(m_bulk_data.get_entity(key));
                }
                remoteEdgeInfo.m_sharedNodes = nodes;
                receivedSharedEdges.push_back(remoteEdgeInfo);
            }
        }
    }
}

void ElemElemGraph::connect_remote_element_to_existing_graph( const impl::SharedEdgeInfo &receivedSharedEdge)
{
    const stk::mesh::EntityVector &sideNodes = receivedSharedEdge.m_sharedNodes;

    unsigned side_index = receivedSharedEdge.get_local_element_side_index();
    stk::mesh::Entity localElem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, receivedSharedEdge.get_local_element_global_id());
    stk::mesh::EntityVector localElemSideNodes;
    stk::topology elemTopo = m_bulk_data.bucket(localElem).topology();
    const stk::mesh::Entity* elemNodes = m_bulk_data.begin_nodes(localElem);
    impl::fill_element_side_nodes_from_topology(elemTopo, elemNodes, side_index, localElemSideNodes);

    stk::EquivalentPermutation permutationIfConnected = stk::mesh::side_equivalent(m_bulk_data, localElem, side_index, sideNodes.data());
    STK_ThrowRequireWithSierraHelpMsg(permutationIfConnected.is_equivalent);

    impl::LocalId local_elem_id = get_local_element_id(localElem);
    impl::LocalId negSgnRemoteElemId = -1 * static_cast<impl::LocalId>(receivedSharedEdge.get_remote_element_global_id());

    GraphEdge graphEdge(local_elem_id, side_index, negSgnRemoteElemId, receivedSharedEdge.get_remote_element_side_index());
    m_edgesToAdd.push_back(graphEdge);

    impl::ParallelInfo parInfo(receivedSharedEdge.get_remote_processor_rank(),
                                permutationIfConnected.permutation_number,
                                receivedSharedEdge.m_remoteElementTopology);

    m_parallelInfoForGraphEdges.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
}

stk::topology ElemElemGraph::get_topology_of_remote_element(const GraphEdge &graphEdge)
{
    impl::ParallelInfo &parallel_edge_info = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    return parallel_edge_info.m_remote_element_topology;
}

stk::topology ElemElemGraph::get_topology_of_connected_element(const GraphEdge &graphEdge)
{
    if(graphEdge.elem2()<0)
    {
        return this->get_topology_of_remote_element(graphEdge);
    }
    return m_element_topologies[graphEdge.elem2()];
}

stk::mesh::SideIdChooser ElemElemGraph::get_side_id_chooser()
{
    return stk::mesh::SideIdChooser(m_bulk_data, m_idMapper, m_graph, m_coincidentGraph);
}

impl::LocalId ElemElemGraph::get_new_local_element_id_from_pool()
{
    impl::LocalId new_local_id;
    if (!m_deleted_element_local_id_pool.empty())
    {
        new_local_id = m_deleted_element_local_id_pool.back();
        m_deleted_element_local_id_pool.pop_back();
    }
    else
    {
        new_local_id = m_graph.get_num_elements_in_graph();
        m_graph.add_new_element();
        m_idMapper.make_space_for_local_id(new_local_id);
        m_element_topologies.push_back(stk::topology::INVALID_TOPOLOGY);
    }
    return new_local_id;
}

bool ElemElemGraph::is_valid_graph_element(stk::mesh::Entity localElement) const
{
    bool value = false;
//    if (m_bulk_data.is_valid(local_element))
    {
        impl::LocalId max_elem_id = static_cast<impl::LocalId>(m_graph.get_num_elements_in_graph());
        if (m_idMapper.does_entity_have_local_id(localElement))
        {
            impl::LocalId elem_id = m_idMapper.entity_to_local(localElement);
            value = elem_id >= 0 && elem_id < max_elem_id;
        }
    }
    return value;
}

void ElemElemGraph::pack_deleted_element_comm(stk::CommSparse &comm,
                                              const std::vector<impl::DeletedElementData> &local_elem_and_remote_connected_elem)
{
    for(size_t i=0;i<local_elem_and_remote_connected_elem.size();++i)
    {
        int remote_proc = local_elem_and_remote_connected_elem[i].m_remoteProc;

        impl::LocalId elem_to_delete_local_id = local_elem_and_remote_connected_elem[i].m_deletedElement;
        stk::mesh::Entity elem_to_delete = m_idMapper.local_to_entity(elem_to_delete_local_id);
        stk::mesh::EntityId elem_to_delete_global_id = m_bulk_data.identifier(elem_to_delete);
        stk::mesh::EntityId connected_elem_global_id = local_elem_and_remote_connected_elem[i].m_remoteElement;

        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(elem_to_delete_global_id);
        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(connected_elem_global_id);
    }
}

void ElemElemGraph::pack_shell_connectivity(stk::CommSparse & comm,
                                            const std::vector<impl::ShellConnectivityData> & shellConnectivityList,
                                            const std::vector<impl::ElementSidePair> &deletedShells)
{
    for(size_t i=0; i<shellConnectivityList.size(); i++)
    {
        const impl::ShellConnectivityData & shellConnectivityData = shellConnectivityList[i];
        if (shellConnectivityData.m_farElementIsRemote) {
            impl::LocalId deletedElementId = deletedShells[i].first;
            int deletedElementSide = deletedShells[i].second;
            GraphEdge graphEdge(deletedElementId, deletedElementSide, m_parallelInfoForGraphEdges.convert_remote_global_id_to_negative_local_id(shellConnectivityData.m_farElementId), shellConnectivityData.m_farElementSide);
            const impl::ParallelInfo& parallelInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
            int remote_proc = parallelInfo.get_proc_rank_of_neighbor();

            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(shellConnectivityData.m_nearElementId);
            comm.send_buffer(remote_proc).pack<int>(shellConnectivityData.m_nearElementSide);
            comm.send_buffer(remote_proc).pack<int>(shellConnectivityData.m_nearElementProc);
            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(shellConnectivityData.m_shellElementId);
            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(shellConnectivityData.m_farElementId);
            comm.send_buffer(remote_proc).pack<int>(shellConnectivityData.m_farElementSide);
        }
    }
}

bool ElemElemGraph::are_connectivities_for_same_graph_edge(const impl::ShellConnectivityData& a,
                                                           const impl::ShellConnectivityData& b)
{
    return a.m_nearElementId == b.m_nearElementId &&
           a.m_nearElementSide == b.m_nearElementSide &&
           a.m_farElementId == b.m_farElementId &&
           a.m_farElementSide == b.m_farElementSide;
}

bool ElemElemGraph::did_already_delete_a_shell_between_these_elements(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                                                      const impl::ShellConnectivityData& shellConnectivityData)
{
    for(const impl::ShellConnectivityData& shellConn : shellConnectivityList)
        if(are_connectivities_for_same_graph_edge(shellConn, shellConnectivityData))
            return true;
    return false;
}

void ElemElemGraph::collect_local_shell_connectivity_data(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                                          std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                                          std::vector<impl::ElementSidePair> &deletedShells)
{
    for (const stk::mesh::impl::DeletedElementInfo &deletedElementInfo : elements_to_delete) {
        stk::mesh::Entity elem_to_delete = deletedElementInfo.entity;

        if (deletedElementInfo.isShell) {
            STK_ThrowRequireWithSierraHelpMsg(is_valid_graph_element(elem_to_delete));
            impl::LocalId shell_to_delete_id = get_local_element_id( elem_to_delete);
            size_t num_connected_elems = m_graph.get_num_edges_for_element(shell_to_delete_id);
            for (int near_elem_index = num_connected_elems - 1; near_elem_index >= 0; --near_elem_index) {
                impl::ShellConnectivityData shellConnectivityData;

                const GraphEdge & graphEdge = m_graph.get_edge_for_element(shell_to_delete_id, near_elem_index);
                impl::LocalId near_elem_id = graphEdge.elem2();

                if (near_elem_id >= 0) {
                    stk::mesh::Entity nearElem = m_idMapper.local_to_entity(near_elem_id);
                    shellConnectivityData.m_nearElementSide = graphEdge.side2();
                    shellConnectivityData.m_nearElementId = m_bulk_data.identifier(nearElem);
                    shellConnectivityData.m_nearElementProc = m_bulk_data.parallel_rank();
                }
                else {
                    shellConnectivityData.m_nearElementId = -near_elem_id;

                    impl::ParallelInfo& parallelInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
                    shellConnectivityData.m_nearElementSide = graphEdge.side2();
                    shellConnectivityData.m_nearElementProc = parallelInfo.get_proc_rank_of_neighbor();
                }

                shellConnectivityData.m_shellElementId = deletedElementInfo.identifier;

                for(size_t i=0; i<num_connected_elems; ++i) {
                    const GraphEdge & nestedGraphEdge = m_graph.get_edge_for_element(shell_to_delete_id, i);
                    impl::LocalId graphElemId = nestedGraphEdge.elem2();
                    if (graphElemId != near_elem_id) {
                        if (graphElemId < 0) {
                            // Remote Element
                            shellConnectivityData.m_farElementId = -graphElemId;
                            shellConnectivityData.m_farElementIsRemote = true;

                            impl::ParallelInfo& parallelInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(nestedGraphEdge);
                            shellConnectivityData.m_farElementSide = nestedGraphEdge.side2();
                            shellConnectivityData.m_farElementProc = parallelInfo.get_proc_rank_of_neighbor();
                        }
                        else {
                            // Local Element
                            stk::mesh::Entity remoteElement = m_idMapper.local_to_entity(graphElemId);
                            shellConnectivityData.m_farElementId = m_bulk_data.identifier(remoteElement);
                            shellConnectivityData.m_farElementSide = nestedGraphEdge.side2();
                            shellConnectivityData.m_farElementIsRemote = false;
                            shellConnectivityData.m_farElementProc = m_bulk_data.parallel_rank();
                        }

                        if(!did_already_delete_a_shell_between_these_elements(shellConnectivityList, shellConnectivityData))
                        {
                            deletedShells.emplace_back(nestedGraphEdge.elem1(),nestedGraphEdge.side1());
                            shellConnectivityList.push_back(shellConnectivityData);
                        }
                        break;
                    }
                }
            }
        }
    }
}

bool ElemElemGraph::communicate_if_shell_connectivity(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                                   const std::vector<impl::ElementSidePair> &deletedShells) {

    unsigned localNumShells = shellConnectivityList.size();
    unsigned globalMaxShells = localNumShells;
    if (m_bulk_data.parallel_size() > 1) {
      stk::all_reduce_max(m_bulk_data.parallel(), &localNumShells, &globalMaxShells, 1);
    }
    if (globalMaxShells == 0) {
        return false;
    }

    stk::CommSparse shellComm(m_bulk_data.parallel());
    pack_shell_connectivity(shellComm, shellConnectivityList, deletedShells);
    shellComm.allocate_buffers();
    pack_shell_connectivity(shellComm, shellConnectivityList, deletedShells);
    shellComm.communicate();

    for (int proc = 0; proc < m_bulk_data.parallel_size(); ++proc) {
        while (shellComm.recv_buffer(proc).remaining()) {
            impl::ShellConnectivityData shellConnectivityData;
            shellComm.recv_buffer(proc).unpack<stk::mesh::EntityId>(shellConnectivityData.m_farElementId); // Flip remote and local
            shellComm.recv_buffer(proc).unpack<int>(shellConnectivityData.m_farElementSide);
            shellComm.recv_buffer(proc).unpack<int>(shellConnectivityData.m_farElementProc);
            shellComm.recv_buffer(proc).unpack<stk::mesh::EntityId>(shellConnectivityData.m_shellElementId);
            shellComm.recv_buffer(proc).unpack<stk::mesh::EntityId>(shellConnectivityData.m_nearElementId); // Flip remote and local
            shellComm.recv_buffer(proc).unpack<int>(shellConnectivityData.m_nearElementSide);
            shellConnectivityData.m_nearElementProc = m_bulk_data.parallel_rank();

            shellConnectivityList.push_back(shellConnectivityData);
        }
    }

    return globalMaxShells > 0;
}

void ElemElemGraph::delete_local_connections_and_collect_remote(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                                                std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem)
{
    std::vector<GraphEdge> pllEdgesToDelete;
    for (const stk::mesh::impl::DeletedElementInfo &deletedElementInfo : elements_to_delete) {
        stk::mesh::Entity elem_to_delete = deletedElementInfo.entity;
        impl::LocalId elem_to_delete_id = get_local_element_id(elem_to_delete);

        size_t num_connected_elems = m_graph.get_num_edges_for_element(elem_to_delete_id);
        for (int conn_elem_index = num_connected_elems - 1; conn_elem_index >= 0; --conn_elem_index) {
            const GraphEdge & graphEdge = m_graph.get_edge_for_element(elem_to_delete_id, conn_elem_index);
            impl::LocalId connected_elem_id = graphEdge.elem2();

            bool local_connection = connected_elem_id >= 0;
            if (local_connection) {
                m_edgesToDelete.push_back(create_symmetric_edge(graphEdge));
            } else {
                impl::ParallelInfo& parallelInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
                int remote_proc = parallelInfo.get_proc_rank_of_neighbor();

                stk::mesh::EntityId connected_global_id = -connected_elem_id;
                impl::DeletedElementData deletedElementData;
                deletedElementData.m_deletedElement = elem_to_delete_id;
                deletedElementData.m_remoteElement = connected_global_id;
                deletedElementData.m_remoteProc = remote_proc;
                local_elem_and_remote_connected_elem.push_back(deletedElementData);
                pllEdgesToDelete.push_back(graphEdge);
            }
        }
        for(const stk::mesh::GraphEdge &graphEdge : m_coincidentGraph.get_edges_for_element(elem_to_delete_id))
        {
            m_coincidentGraph.delete_edge(create_symmetric_edge(graphEdge));
        }
        m_coincidentGraph.delete_all_edges(elem_to_delete_id);

        GraphEdgesForElement graphEdgesForElement = m_graph.get_edges_for_element(elem_to_delete_id);
        for(const GraphEdge& edge : graphEdgesForElement) {
          m_edgesToDelete.push_back(edge);
        }
    }
    m_parallelInfoForGraphEdges.erase_edges(pllEdgesToDelete);
}

void ElemElemGraph::communicate_remote_connections_to_delete(const std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem,
                                                             std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges)
{
    stk::CommSparse comm(m_bulk_data.parallel());
    pack_deleted_element_comm(comm, local_elem_and_remote_connected_elem);
    comm.allocate_buffers();
    pack_deleted_element_comm(comm, local_elem_and_remote_connected_elem);
    comm.communicate();

    for (int proc = 0; proc < m_bulk_data.parallel_size(); ++proc) {
        while (comm.recv_buffer(proc).remaining()) {
            stk::mesh::EntityId deleted_elem_global_id;
            stk::mesh::EntityId connected_elem_global_id;
            comm.recv_buffer(proc).unpack<stk::mesh::EntityId>(deleted_elem_global_id);
            comm.recv_buffer(proc).unpack<stk::mesh::EntityId>(connected_elem_global_id);
            remote_edges.emplace_back(deleted_elem_global_id, connected_elem_global_id);
        }
    }
}

void ElemElemGraph::clear_deleted_element_connections(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete)
{
    for (const stk::mesh::impl::DeletedElementInfo &deletedElementInfo : elements_to_delete)
    {
        stk::mesh::Entity elem = deletedElementInfo.entity;
        impl::LocalId elem_to_delete_id = m_idMapper.entity_to_local(elem);
        m_deleted_element_local_id_pool.push_back(elem_to_delete_id);
        m_element_topologies[elem_to_delete_id] = stk::topology::INVALID_TOPOLOGY;
        m_idMapper.clear_local_id_for_elem(elem);
    }
}

void ElemElemGraph::delete_remote_connection(stk::mesh::Entity connected_elem,
                                             stk::mesh::EntityId deleted_elem_global_id,
                                             std::vector<GraphEdge>& pllEdgesToDelete)
{
  impl::LocalId connected_elem_id = m_idMapper.entity_to_local(connected_elem);
  size_t num_conn_elem = m_graph.get_num_edges_for_element(connected_elem_id);
  bool found_deleted_elem = false;
  for (size_t conn_elem_index = 0; conn_elem_index < num_conn_elem; ++conn_elem_index) {
      const GraphEdge & graphEdge = m_graph.get_edge_for_element(connected_elem_id, conn_elem_index);
      if (graphEdge.elem2() == static_cast<int64_t>(-deleted_elem_global_id)) {
          pllEdgesToDelete.push_back(graphEdge);
          m_edgesToDelete.push_back(graphEdge);
          found_deleted_elem = true;
          break;
      }
  }
  STK_ThrowRequireWithSierraHelpMsg(found_deleted_elem);
}

void ElemElemGraph::delete_remote_connections(const std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges)
{
    std::vector<GraphEdge> pllEdgesToDelete;
    for (const std::pair<stk::mesh::EntityId, stk::mesh::EntityId> & edge : remote_edges) {
        stk::mesh::EntityId deleted_elem_global_id = edge.first;
        stk::mesh::EntityId connected_elem_global_id = edge.second;
        stk::mesh::Entity connected_elem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, connected_elem_global_id);

        if (!is_valid_graph_element(connected_elem)) {
            continue;
        }

        delete_remote_connection(connected_elem, deleted_elem_global_id, pllEdgesToDelete);
    }
    m_parallelInfoForGraphEdges.erase_edges(pllEdgesToDelete);
}

bool ElemElemGraph::is_connected_to_shell_on_side(stk::mesh::impl::LocalId localId, int side)
{
    for(const stk::mesh::GraphEdge& graphEdge : m_graph.get_edges_for_element(localId))
        if(graphEdge.side1() == side && get_topology_of_connected_element(graphEdge).is_shell())
            return true;
    return false;
}

void ElemElemGraph::reconnect_volume_elements_across_deleted_shells(std::vector<impl::ShellConnectivityData> & shellConnectivityList)
{
    // Prune the shell connectivity list for entries that we do not need to act on
    std::vector<impl::ShellConnectivityData>::iterator it = shellConnectivityList.begin();
    while (it != shellConnectivityList.end()) {
        const impl::ShellConnectivityData& shellConnectivityData = *it;
        if (shellConnectivityData.m_nearElementProc != m_bulk_data.parallel_rank()) {
            it = shellConnectivityList.erase(it);
        } else {
            ++it;
        }
    }

    impl::ElemSideProcVector shellNeighborsToReconnect;
    shellNeighborsToReconnect.reserve(shellConnectivityList.size());
    for (const impl::ShellConnectivityData& data : shellConnectivityList) {
        stk::mesh::Entity localElement = m_bulk_data.get_entity(stk::topology::ELEM_RANK, data.m_nearElementId);
        if (m_bulk_data.is_valid(localElement))
        {
            stk::mesh::impl::LocalId localElemLocalId = m_idMapper.entity_to_local(localElement);
            if(!is_connected_to_shell_on_side(localElemLocalId, data.m_nearElementSide))
            {
                if (data.m_nearElementProc == m_bulk_data.parallel_rank() && data.m_farElementProc == m_bulk_data.parallel_rank()) {
                    impl::LocalId localElementId = get_local_element_id(localElement);
                    stk::mesh::Entity remoteElement = m_bulk_data.get_entity(stk::topology::ELEM_RANK, data.m_farElementId);
                    if (m_bulk_data.is_valid(remoteElement)) {
                        impl::LocalId remoteElementId = get_local_element_id(remoteElement);
                        m_edgesToAdd.push_back(stk::mesh::GraphEdge(localElementId, data.m_nearElementSide, remoteElementId, data.m_farElementSide));
                    }
                }
                else {
                    shellNeighborsToReconnect.push_back(impl::ElemSideProc(localElement, data.m_nearElementSide, data.m_farElementProc));
                }
            }
        }
    }

    std::sort(shellNeighborsToReconnect.begin(), shellNeighborsToReconnect.end());
    impl::ParallelElementDataVector elementSidesReceived;
    communicate_shared_sides(shellNeighborsToReconnect, elementSidesReceived);
    fill_parallel_graph(shellNeighborsToReconnect, elementSidesReceived);
    m_modCycleWhenGraphModified = m_bulk_data.synchronized_count();
}

stk::mesh::impl::DeletedElementInfoVector ElemElemGraph::filter_delete_elements_argument(const stk::mesh::impl::DeletedElementInfoVector& elements_to_delete_argument) const
{
    stk::mesh::impl::DeletedElementInfoVector filtered_elements;
    filtered_elements.reserve(elements_to_delete_argument.size());
    for(stk::mesh::impl::DeletedElementInfo deletedElemInfo : elements_to_delete_argument) {
        if (is_valid_graph_element(deletedElemInfo.entity)) {
            filtered_elements.push_back(deletedElemInfo);
        }
    }
    return filtered_elements;
}

void ElemElemGraph::delete_elements(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete_argument)
{
    m_edgesToAdd.clear();
    m_edgesToDelete.clear();
    stk::mesh::impl::DeletedElementInfoVector elements_to_delete = filter_delete_elements_argument(elements_to_delete_argument);

    std::vector<impl::ShellConnectivityData> shellConnectivityList;
    std::vector<impl::ElementSidePair> deletedNonCoincidentShells;
    collect_local_shell_connectivity_data(elements_to_delete, shellConnectivityList, deletedNonCoincidentShells);

    const bool shellsExist = communicate_if_shell_connectivity(shellConnectivityList, deletedNonCoincidentShells);

    std::vector<impl::DeletedElementData> local_elem_and_remote_connected_elem;
    delete_local_connections_and_collect_remote(elements_to_delete, local_elem_and_remote_connected_elem);

    std::vector< std::pair< stk::mesh::EntityId, stk::mesh::EntityId > > remote_edges;
    communicate_remote_connections_to_delete(local_elem_and_remote_connected_elem, remote_edges);

    clear_deleted_element_connections(elements_to_delete);

    delete_remote_connections(remote_edges);

    stk::util::sort_and_unique(m_edgesToDelete, GraphEdgeLessByElem1());
    m_graph.delete_sorted_edges(m_edgesToDelete);
    {
      std::vector<GraphEdge> emptyDelete;
      m_edgesToDelete.swap(emptyDelete);
    }

    if (shellsExist) {
        reconnect_volume_elements_across_deleted_shells(shellConnectivityList);
    }

    stk::util::sort_and_unique(m_edgesToAdd, GraphEdgeLessByElem1());
    m_graph.add_sorted_edges(m_edgesToAdd);
    {
      std::vector<GraphEdge> emptyAdd;
      m_edgesToAdd.swap(emptyAdd);
    }

    m_modCycleWhenGraphModified = m_bulk_data.synchronized_count();
}

template <typename GraphType>
void ElemElemGraph::add_new_local_edges_to_graph(GraphType &graph, const std::vector<stk::mesh::GraphEdge> &newGraphEdges)
{
    for (const stk::mesh::GraphEdge& graphEdge : newGraphEdges)
    {
        stk::mesh::Entity neighbor = m_idMapper.local_to_entity(graphEdge.elem2());
        if (is_valid_graph_element(neighbor)) {
          graph.add_edge(graphEdge);
          graph.add_edge(create_symmetric_edge(graphEdge));
        }
    }
}

void ElemElemGraph::add_vertex(impl::LocalId newElemLocalId, stk::mesh::Entity elem)
{
    m_idMapper.add_new_entity_with_local_id(elem, newElemLocalId);

    stk::topology elem_topology = m_bulk_data.bucket(elem).topology();
    m_element_topologies[newElemLocalId] = elem_topology;
    if (elem_topology.is_shell()) {
        m_any_shell_elements_exist = true;
    }
}

std::pair<stk::mesh::EntityVector,stk::mesh::EntityVector> ElemElemGraph::filter_add_elements_arguments(const stk::mesh::EntityVector& allUnfilteredElementsNotAlreadyInGraph) const
{
  stk::mesh::EntityVector allElementsNotAlreadyInGraph;
  stk::mesh::EntityVector shellElementsNotAlreadyInGraph;
  allElementsNotAlreadyInGraph.reserve(allUnfilteredElementsNotAlreadyInGraph.size());
  for(stk::mesh::Entity element : allUnfilteredElementsNotAlreadyInGraph) {
    if(m_bulk_data.is_valid(element) && m_bulk_data.bucket(element).owned()) {
      STK_ThrowRequire(!is_valid_graph_element(element));
      allElementsNotAlreadyInGraph.push_back(element);
      if (m_bulk_data.bucket(element).topology().is_shell()) {
        shellElementsNotAlreadyInGraph.push_back(element);
      }
    }
  }
  return std::make_pair(allElementsNotAlreadyInGraph,shellElementsNotAlreadyInGraph);
}

void ElemElemGraph::add_elements_locally(const stk::mesh::EntityVector& allElementsNotAlreadyInGraph)
{
    m_idMapper.make_space_for_new_elements(allElementsNotAlreadyInGraph);
    stk::mesh::EntityVector side_nodes;
    impl::SerialElementDataVector connectedElementDataVector;
    const bool only_consider_upper_symmetry = false;
    for(stk::mesh::Entity newElem : allElementsNotAlreadyInGraph)
    {
        impl::LocalId newElemLocalId = get_new_local_element_id_from_pool();
        add_vertex(newElemLocalId, newElem);
        add_local_graph_edges_for_elem(m_bulk_data.mesh_index(newElem), newElemLocalId, m_newGraphEdges, m_coincidentGraphEdges, side_nodes, connectedElementDataVector, only_consider_upper_symmetry);
        GraphEdgeAdder graphEdgeAdder(m_edgesToAdd);
        add_new_local_edges_to_graph(graphEdgeAdder, m_newGraphEdges);
        add_new_local_edges_to_graph(m_coincidentGraph, m_coincidentGraphEdges);
    }
}

void ElemElemGraph::add_elements(const stk::mesh::EntityVector &allUnfilteredElementsNotAlreadyInGraph)
{
  stk::mesh::EntityVector allElementsNotAlreadyInGraph;
  stk::mesh::EntityVector addedShellsVector;
  std::tie(allElementsNotAlreadyInGraph,addedShellsVector) = filter_add_elements_arguments(allUnfilteredElementsNotAlreadyInGraph);

  add_elements_locally(allElementsNotAlreadyInGraph);

  impl::ElemSideProcVector only_added_elements = impl::build_element_side_ids_to_proc_map(m_bulk_data, allElementsNotAlreadyInGraph);
  std::sort(allElementsNotAlreadyInGraph.begin(), allElementsNotAlreadyInGraph.end());

  std::vector<int> sharingProcs;
  stk::mesh::EntityVector elemNodes;
  for(stk::mesh::Entity addedShell : addedShellsVector) {
    elemNodes.assign(m_bulk_data.begin_nodes(addedShell), m_bulk_data.end_nodes(addedShell));
    m_bulk_data.shared_procs_intersection(elemNodes, sharingProcs);

    if (!sharingProcs.empty()) {
      stk::mesh::EntityVector solidElems = impl::gather_solid_elements_connected_to_shell(m_bulk_data, addedShell); 

      for(stk::mesh::Entity elem : solidElems) {
        stk::mesh::EntityVector::iterator elem_iter = std::lower_bound(allElementsNotAlreadyInGraph.begin(), allElementsNotAlreadyInGraph.end(), elem);
        if (elem_iter == allElementsNotAlreadyInGraph.end()) {
          unsigned side = stk::mesh::get_ordinal_and_permutation(m_bulk_data, elem, m_bulk_data.mesh_meta_data().side_rank(), elemNodes).first;
          for(int p : sharingProcs) {
            only_added_elements.push_back(impl::ElemSideProc(elem, side, p));
          }
        }
      }
    }
  }

  std::sort(only_added_elements.begin(), only_added_elements.end());

  bool anyShellsExist = m_any_shell_elements_exist;
  if (m_bulk_data.parallel_size() > 1) {
    anyShellsExist = stk::is_true_on_any_proc(m_bulk_data.parallel(), m_any_shell_elements_exist);
  }

  impl::ParallelElementDataVector elementSidesReceived;
  communicate_shared_sides(only_added_elements, elementSidesReceived);
  fill_parallel_graph(only_added_elements, elementSidesReceived);

  stk::util::sort_and_unique(m_edgesToAdd, GraphEdgeLessByElem1());
  m_graph.add_sorted_edges(m_edgesToAdd);

  {
    std::vector<GraphEdge> emptyAdd;
    m_edgesToAdd.swap(emptyAdd);
  }

  if (anyShellsExist) {
    std::sort(addedShellsVector.begin(), addedShellsVector.end());
    STK_ThrowRequire(m_edgesToDelete.empty());
    GraphInfo graphInfo(m_graph, m_parallelInfoForGraphEdges, m_element_topologies);
    remove_graph_edges_blocked_by_shell(graphInfo);

    stk::CommSparse comm(m_bulk_data.parallel());
    for (int phase = 0; phase < 2; ++phase) {
      pack_remote_edge_across_shell(comm, addedShellsVector, phase);
      if (0 == phase) {
          comm.allocate_buffers();
      }
      if (1 == phase) {
          comm.communicate();
      }
    }

    unpack_remote_edge_across_shell(comm);
  }

  stk::util::sort_and_unique(m_edgesToDelete, GraphEdgeLessByElem1());
  m_graph.delete_sorted_edges(m_edgesToDelete);

  {
    std::vector<GraphEdge> emptyDelete;
    m_edgesToDelete.swap(emptyDelete);
  }

  m_modCycleWhenGraphModified = m_bulk_data.synchronized_count();
}

void ElemElemGraph::break_remote_shell_connectivity_and_pack(stk::CommSparse &comm, impl::LocalId leftId, impl::LocalId rightId, int phase)
{
    size_t index = 0;
    size_t numConnected = m_graph.get_num_edges_for_element(leftId);
    while(index < numConnected)
    {
        const GraphEdge & graphEdge = m_graph.get_edge_for_element(leftId, index);
        if(graphEdge.elem2() == rightId)
        {
            stk::mesh::Entity localElem = m_idMapper.local_to_entity(leftId);
            int sharingProc = get_owning_proc_id_of_remote_element(localElem, index);
            if (phase == 1)
            {
                m_edgesToDelete.push_back(graphEdge);
            }
            else
            {
                ++index;
            }

            stk::mesh::EntityId localElemId = m_bulk_data.identifier(localElem);

            comm.send_buffer(sharingProc).pack<stk::mesh::EntityId>(localElemId);
            comm.send_buffer(sharingProc).pack<stk::mesh::EntityId>(-rightId);
            break;
        }
        else
        {
            ++index;
        }
    }
}

void ElemElemGraph::pack_both_remote_shell_connectivity(stk::CommSparse &comm, impl::LocalId shellId, impl::LocalId leftId, impl::LocalId rightId)
{
    size_t index = 0;
    size_t num_connected = m_graph.get_num_edges_for_element(shellId);
    while(index < num_connected)
    {
        const GraphEdge & graphEdge = m_graph.get_edge_for_element(shellId, index);
        if (graphEdge.elem2() == rightId)
        {
            stk::mesh::Entity shell = m_idMapper.local_to_entity(shellId);
            int sharingProc = get_owning_proc_id_of_remote_element(shell, index);

            comm.send_buffer(sharingProc).pack<stk::mesh::EntityId>(-leftId);
            comm.send_buffer(sharingProc).pack<stk::mesh::EntityId>(-rightId);
        }
        else if (graphEdge.elem2() == leftId)
        {
            stk::mesh::Entity shell = m_idMapper.local_to_entity(shellId);
            int sharingProc = get_owning_proc_id_of_remote_element(shell, index);

            comm.send_buffer(sharingProc).pack<stk::mesh::EntityId>(-rightId);
            comm.send_buffer(sharingProc).pack<stk::mesh::EntityId>(-leftId);
        }
        ++index;
    }
}

void ElemElemGraph::pack_remote_edge_across_shell(stk::CommSparse &comm, stk::mesh::EntityVector &addedShells, int phase)
{
    for(stk::mesh::Entity &shell : addedShells)
    {
        impl::LocalId shellId = m_idMapper.entity_to_local(shell);
        impl::LocalId leftId = impl::INVALID_LOCAL_ID;
        impl::LocalId rightId = impl::INVALID_LOCAL_ID;
        size_t numConnected = m_graph.get_num_edges_for_element(shellId);
        for(size_t i=0; i<numConnected; i++)
        {
            const GraphEdge & graphEdge = m_graph.get_edge_for_element(shellId, i);
            bool connectedElemIsShell = get_topology_of_connected_element(graphEdge).is_shell();
            if (leftId == impl::INVALID_LOCAL_ID && !connectedElemIsShell)
            {
                leftId = graphEdge.elem2();
                continue;
            }
            if (rightId == impl::INVALID_LOCAL_ID && !connectedElemIsShell)
            {
                rightId = graphEdge.elem2();
                continue;
            }
        }
        bool isLeftRemote = leftId < 0;
        bool isRightRemote = rightId < 0;

        if (leftId == impl::INVALID_LOCAL_ID || rightId == impl::INVALID_LOCAL_ID)
        {
            continue;
        }

        if (!isLeftRemote && !isRightRemote)
        {
            continue;
        }

        if (isLeftRemote && isRightRemote) {
            pack_both_remote_shell_connectivity(comm, shellId, leftId, rightId);
        }
        else if (!isLeftRemote) {
            break_remote_shell_connectivity_and_pack(comm, leftId, rightId, phase);
        }
        else if (!isRightRemote) {
            break_remote_shell_connectivity_and_pack(comm, rightId, leftId, phase);
        }
    }
}

void ElemElemGraph::unpack_remote_edge_across_shell(stk::CommSparse &comm)
{
    for(int i=0;i<m_bulk_data.parallel_size();++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            stk::mesh::EntityId localElemId;
            stk::mesh::EntityId remoteElemId;
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(remoteElemId);
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(localElemId);
            stk::mesh::Entity localElem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, localElemId);
            impl::LocalId localId = m_idMapper.entity_to_local(localElem);

            size_t index = 0;
            size_t numConnected = m_graph.get_num_edges_for_element(localId);
            while(index <numConnected)
            {
                const GraphEdge & graphEdge = m_graph.get_edge_for_element(localId, index);
                if(graphEdge.elem2() == (-1 * static_cast<impl::LocalId>(remoteElemId)))
                {
                    m_edgesToDelete.push_back(graphEdge);
                    break;
                }
                else
                {
                    ++index;
                }
            }
        }
    }
}

stk::mesh::Entity ElemElemGraph::add_side_to_mesh(const stk::mesh::impl::ElementSidePair& sidePair, const stk::mesh::PartVector& skinParts)
{
    stk::mesh::Entity element = m_idMapper.local_to_entity(sidePair.first);
    int sideOrd = sidePair.second;
    return m_bulk_data.declare_element_side(element, sideOrd, skinParts);
}

void add_shared_side_to_element(stk::mesh::BulkData& bulkData,
                                const stk::mesh::GraphEdge& graphEdge,
                                const impl::ParallelInfo& parallel_edge_info,
                                stk::mesh::Entity local_element,
                                const stk::mesh::PartVector& parts_for_creating_side,
                                std::vector<stk::mesh::sharing_info> &shared_modified)
{
    int sideOrd = graphEdge.side1();
    STK_ThrowRequireWithSierraHelpMsg(sideOrd != -1);

    stk::mesh::Entity side = bulkData.declare_element_side(local_element,
                                                           sideOrd,
                                                           parts_for_creating_side);

    int other_proc = parallel_edge_info.get_proc_rank_of_neighbor();
    int owning_proc = bulkData.state(side)==Created ? std::min(other_proc, bulkData.parallel_rank()) : bulkData.parallel_owner_rank(side);
    shared_modified.push_back(stk::mesh::sharing_info(side, other_proc, owning_proc));
}

struct LocalEdge
{
    stk::mesh::Entity m_element;
    stk::mesh::Entity m_other_element;
    LocalEdge(stk::mesh::Entity element, stk::mesh::Entity other_element) :
        m_element(element), m_other_element(other_element)
    {}

};

struct RemoteEdge
{
    const GraphEdge& m_graphEdge;
    stk::mesh::impl::ParallelInfo m_parallel_edge_info;
    RemoteEdge(const GraphEdge& graphEdge, stk::mesh::impl::ParallelInfo parallel_edge_info) :
        m_graphEdge(graphEdge), m_parallel_edge_info(parallel_edge_info)
    {}

};

const stk::mesh::BulkData& ElemElemGraph::get_mesh() const
{
    return m_bulk_data;
}

void ElemElemGraph::create_side_entities(const std::vector<int> &exposedSides,
                                         impl::LocalId localId,
                                         const stk::mesh::PartVector& skinParts,
                                         std::vector<stk::mesh::sharing_info> &sharedModified)
{
    stk::mesh::Entity element = m_idMapper.local_to_entity(localId);
    for(size_t i=0;i<exposedSides.size();++i)
    {
        for(const GraphEdge & graphEdge : m_graph.get_edges_for_element(localId))
            add_side_for_remote_edge(graphEdge, exposedSides[i], element, skinParts, sharedModified);

        for(const GraphEdge & graphEdge : get_coincident_edges_for_element(localId))
            add_side_for_remote_edge(graphEdge, exposedSides[i], element, skinParts, sharedModified);

        add_side_to_mesh({localId, exposedSides[i]}, skinParts);
    }
}

void ElemElemGraph::add_side_for_remote_edge(const GraphEdge & graphEdge,
                                                            int elemSide,
                                                            stk::mesh::Entity element,
                                                            const stk::mesh::PartVector& skin_parts,
                                                            std::vector<stk::mesh::sharing_info> &shared_modified)
{
    if(graphEdge.side1() == elemSide)
    {
        if(!impl::is_local_element(graphEdge.elem2()))
        {
            impl::ParallelInfo &parallel_edge_info = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
            add_shared_side_to_element(m_bulk_data, graphEdge, parallel_edge_info, element, skin_parts, shared_modified);
        }
    }
}

unsigned ElemElemGraph::get_max_num_sides_per_element() const
{
    return 6;
}

bool ElemElemGraph::is_valid_graph_edge(const GraphEdge &graphEdge) const
{
  if(!impl::is_local_element(graphEdge.elem1())) {
    return false;
  }

  if(graphEdge.side1() < 0 || graphEdge.side2() < 0) {
    return false;
  }

  if((size_t)graphEdge.elem1() >= m_idMapper.size()) {
    return false;
  }

  return true;
}

namespace impl {

impl::ElemSideProcVector gather_element_side_ids_to_send(const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector elements_to_communicate;
    std::vector<bool> visitedElem(bulkData.get_size_of_entity_index_space(), false);
    const stk::mesh::BucketVector& shared_node_buckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().globally_shared_part());
    for(size_t i=0; i<shared_node_buckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *shared_node_buckets[i];
        for(stk::mesh::Entity node : bucket)
        {
            const stk::mesh::Entity* elements = bulkData.begin_elements(node);
            const unsigned num_elements = bulkData.num_elements(node);
            for(unsigned elemIdx=0; elemIdx<num_elements; ++elemIdx)
            {
                if (!visitedElem[elements[elemIdx].local_offset()]) {
                    visitedElem[elements[elemIdx].local_offset()] = true;
                    if (bulkData.bucket(elements[elemIdx]).owned()) {
                        elements_to_communicate.push_back(elements[elemIdx]);
                    }
                }
            }
        }
    }

    return impl::build_element_side_ids_to_proc_map(bulkData, elements_to_communicate);
}

stk::mesh::EntityVector gather_solid_elements_connected_to_shell(const stk::mesh::BulkData& bulkData, stk::mesh::Entity shellElement)
{
  stk::mesh::EntityVector solidElementsConnectedToShell;

  stk::topology elemTopo = bulkData.bucket(shellElement).topology();
  if(bulkData.is_valid(shellElement) && elemTopo.is_shell()) {

    const stk::mesh::Entity* elemNodes = bulkData.begin_nodes(shellElement);
    stk::mesh::EntityVector connectedElements;
    find_entities_these_nodes_have_in_common_and(bulkData, stk::topology::ELEMENT_RANK, elemTopo.num_nodes(), elemNodes, connectedElements,
      [&](const Entity& entity){
        return !bulkData.bucket(entity).topology().is_shell();
      }
    );

    for (stk::mesh::Entity connectedElement : connectedElements) {
      if (connectedElement != shellElement && bulkData.bucket(connectedElement).owned()) {
        solidElementsConnectedToShell.push_back(connectedElement);
      }
    }
  }

  return solidElementsConnectedToShell;
}

}//namespace impl
}} // end namespaces stk mesh

