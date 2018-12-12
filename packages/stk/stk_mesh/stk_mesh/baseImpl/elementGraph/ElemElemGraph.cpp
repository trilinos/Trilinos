#include "ElemElemGraph.hpp"
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
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk { namespace mesh {

stk::mesh::EntityVector convert_keys_to_entities(stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKey>& keys)
{
    stk::mesh::EntityVector entities(keys.size());
    for(size_t i=0;i<keys.size();++i)
    {
        entities[i] = bulkData.get_entity(keys[i]);
    }
    return entities;
}

void SharedSidesCommunication::pack_shared_side_nodes_of_elements(stk::CommSparse &comm) const
{
    impl::ElemSideToProcAndFaceId::const_iterator iter = m_elementSidesToSend.begin();
    impl::ElemSideToProcAndFaceId::const_iterator end = m_elementSidesToSend.end();
    stk::mesh::EntityVector side_nodes;
    for(; iter!= end; ++iter)
    {
        stk::mesh::Entity elem = iter->first.entity;
        unsigned side_index    = iter->first.side_id;
        int sharing_proc       = iter->second.proc;
        stk::mesh::EntityId element_id     = m_bulkData.identifier(elem);
        stk::topology topology = m_bulkData.bucket(elem).topology();

        comm.send_buffer(sharing_proc).pack<stk::mesh::EntityId>(element_id);
        comm.send_buffer(sharing_proc).pack<stk::topology>(topology);
        comm.send_buffer(sharing_proc).pack<unsigned>(side_index);

        impl:: fill_element_side_nodes_from_topology(m_bulkData, elem, side_index, side_nodes);
        std::vector<stk::mesh::EntityKey> side_node_entity_keys(side_nodes.size());
        for(size_t i=0; i<side_nodes.size(); ++i)
        {
            side_node_entity_keys[i] = m_bulkData.entity_key(side_nodes[i]);
        }
        stk::pack_vector_to_proc(comm, side_node_entity_keys, sharing_proc);
    }
}

SideNodeToReceivedElementDataMap SharedSidesCommunication::unpack_side_data(stk::CommSparse comm) const
{
    SideNodeToReceivedElementDataMap element_side_data_received;
    std::vector<stk::mesh::EntityKey> node_keys;
    stk::mesh::impl::ParallelElementData elementData;
    stk::mesh::EntityVector sortedSideNodes;
    for(int proc_id = 0; proc_id < m_bulkData.parallel_size(); ++proc_id)
    {
        if(proc_id != m_bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityId elementIdentifier;
                stk::topology topology;
                unsigned side_index = 0;
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(elementIdentifier);
                comm.recv_buffer(proc_id).unpack<stk::topology>(topology);
                comm.recv_buffer(proc_id).unpack<unsigned>(side_index);

                elementData.set_proc_rank(proc_id);
                elementData.set_element_identifier(elementIdentifier);
                elementData.set_element_topology(topology);
                elementData.set_element_side_index(side_index);

                stk::unpack_vector_from_proc(comm, node_keys, proc_id);
                elementData.set_element_side_nodes(convert_keys_to_entities(m_bulkData, node_keys));
                sortedSideNodes = elementData.get_side_nodes();
                std::sort(sortedSideNodes.begin(), sortedSideNodes.end());
                element_side_data_received[sortedSideNodes].push_back(elementData);
            }
        }
    }
    return element_side_data_received;
}



void SharedSidesCommunication::communicate_element_sides()
{
    stk::CommSparse comm(m_bulkData.parallel());
    stk::pack_and_communicate(comm, [this,&comm](){pack_shared_side_nodes_of_elements(comm);} );
    m_elementSidesReceived = unpack_side_data(comm);
}

ElemElemGraph::ElemElemGraph(stk::mesh::BulkData& bulkData) :
        m_bulk_data(bulkData),
        m_parallelInfoForGraphEdges(bulkData.parallel_rank()),
        m_idMapper()
{
    fill_from_mesh();
}

void ElemElemGraph::fill_from_mesh()
{
    clear_data_members();

    impl::ElemSideToProcAndFaceId elementSideIdsToSend;

    int numElems = size_data_members();
    if (numElems > 0)
    {
        impl::fill_topologies(m_bulk_data, m_idMapper, m_element_topologies);
        fill_graph();
        elementSideIdsToSend = impl::gather_element_side_ids_to_send(m_bulk_data);
    }

    m_parallelInfoForGraphEdges.clear();
    SideNodeToReceivedElementDataMap elementSidesReceived = communicate_shared_sides(elementSideIdsToSend);
    fill_parallel_graph(elementSideIdsToSend, elementSidesReceived);

    GraphInfo graphInfo(m_graph, m_parallelInfoForGraphEdges, m_element_topologies);
    remove_graph_edges_blocked_by_shell(graphInfo);

    update_number_of_parallel_edges();
}

void ElemElemGraph::update_number_of_parallel_edges()
{
    // Figure out the real number of sides that can be generated, using
    // the rule of one side entity per element side (with multiple
    // coincident elements each being connected to the same side)
    m_num_parallel_edges = 0;
    for(size_t i = 0; i < m_graph.get_num_elements_in_graph(); ++i)
    {
        size_t numConnectedElems = m_graph.get_num_edges_for_element(i);
        for(size_t j = 0; j < numConnectedElems; ++j)
        {
            const GraphEdge & graphEdge = m_graph.get_edge_for_element(i, j);
            if(graphEdge.elem2() < 0)
            {
                m_num_parallel_edges++;
            }
        }
    }
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
    ThrowRequireWithSierraHelpMsg(!is_connected_elem_locally_owned(localElement, indexConnElement));
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
        ThrowRequireWithSierraHelpMsg(local_id != impl::INVALID_LOCAL_ID);
    }
    return local_id;
}

int ElemElemGraph::size_data_members()
{
    unsigned numElems = impl::get_num_local_elems(m_bulk_data);
    m_graph.set_num_local_elements(numElems);
    m_element_topologies.resize(numElems);

    m_idMapper.initialize(m_bulk_data);

    m_num_parallel_edges = 0;
    return numElems;
}

void ElemElemGraph::clear_data_members()
{
    m_graph.clear();
    m_idMapper.clear();
    m_element_topologies.clear();
    m_num_parallel_edges = 0;
    m_parallelInfoForGraphEdges.clear();
    m_deleted_element_local_id_pool.clear();
    m_deleted_elem_pool.clear();
    m_coincidentGraph.clear();
}

void ElemElemGraph::fill_elements_attached_to_local_nodes(const stk::mesh::EntityVector& sideNodesOfReceivedElement,
                                                                                  stk::mesh::EntityId elementId,
                                                                                  stk::topology elementTopology,
                                                                                  stk::mesh::EntityVector& scratchEntityVector,
                                                                                  impl::SerialElementDataVector& connectedElementDataVector) const
{
    connectedElementDataVector.clear();
    impl::get_elements_with_larger_ids_connected_via_sidenodes<impl::SerialElementData>(m_bulk_data, elementId, elementTopology, m_idMapper, sideNodesOfReceivedElement,
                                                                                        scratchEntityVector, connectedElementDataVector);
}

void ElemElemGraph::get_elements_attached_to_remote_nodes(const stk::mesh::EntityVector& sideNodesOfReceivedElement, stk::mesh::EntityId elementId, stk::topology elementTopology,
                                                                                     stk::mesh::EntityVector& scratchEntityVector,
                                                                                     impl::ParallelElementDataVector& connectedElementDataVector) const
{
    connectedElementDataVector.clear();
    impl::get_elements_connected_via_sidenodes<impl::ParallelElementData>(m_bulk_data, elementId, elementTopology, m_idMapper,
                                                                                 sideNodesOfReceivedElement, scratchEntityVector, connectedElementDataVector);
}


void ElemElemGraph::insert_edge_between_elements(impl::LocalId local_elem_id,
                                                 int side_index,
                                                 const impl::SerialElementData& otherElem,
                                                 std::vector<stk::mesh::GraphEdge>& graphEdges) const
{
    graphEdges.push_back(stk::mesh::GraphEdge(local_elem_id, side_index, otherElem.get_element_local_id(), otherElem.get_element_side_index()));
}

void ElemElemGraph::add_local_graph_edges_for_elem(const stk::mesh::MeshIndex &meshIndex,
                                                   impl::LocalId local_elem_id,
                                                   std::vector<stk::mesh::GraphEdge> &graphEdges,
                                                   std::vector<stk::mesh::GraphEdge> &coincidentGraphEdges,
                                                   stk::mesh::EntityVector& scratchEntityVector,
                                                   stk::mesh::EntityVector& side_nodes,
                                                   impl::SerialElementDataVector& connectedElementDataVector,
                                                   bool only_consider_upper_symmetry) const
{
    stk::mesh::Entity element = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
    int num_sides = meshIndex.bucket->topology().num_sides();
    graphEdges.clear();
    coincidentGraphEdges.clear();
    stk::mesh::EntityId elemGlobalId = m_bulk_data.identifier(element);
    if(!only_consider_upper_symmetry) elemGlobalId = 0;
    for(int side_index=0; side_index<num_sides; ++side_index)
    {
        impl::fill_element_side_nodes_from_topology(m_bulk_data, element, side_index, side_nodes);
        fill_elements_attached_to_local_nodes(side_nodes, elemGlobalId, meshIndex.bucket->topology(), scratchEntityVector, connectedElementDataVector);
        for (const impl::SerialElementData & otherElem: connectedElementDataVector)
        {
            if(local_elem_id != otherElem.get_element_local_id())
            {
                if(impl::is_coincident_connection(m_bulk_data, element, side_nodes, side_index, otherElem.get_element_topology(), otherElem.get_side_nodes()))
                    insert_edge_between_elements(local_elem_id, side_index, otherElem, coincidentGraphEdges);
                else
                    insert_edge_between_elements(local_elem_id, side_index, otherElem, graphEdges);
            }
        }
    }
    stk::util::sort_and_unique(graphEdges, GraphEdgeLessByElem2());
}

template <typename GraphType>
void add_edges_to_graph(const std::vector<stk::mesh::GraphEdge> &newGraphEdges, GraphType &graph)
{
    for(const stk::mesh::GraphEdge &graphEdge : newGraphEdges)
    {
        graph.add_edge(graphEdge);
        graph.add_edge(create_symmetric_edge(graphEdge));
    }
}

void ElemElemGraph::fill_graph()
{
    const stk::mesh::BucketVector& elemBuckets = m_bulk_data.get_buckets(stk::topology::ELEM_RANK, m_bulk_data.mesh_meta_data().locally_owned_part());
    std::vector<stk::mesh::GraphEdge> newGraphEdges;
    std::vector<stk::mesh::GraphEdge> coincidentGraphEdges;
    stk::mesh::EntityVector scratchEntityVector;
    stk::mesh::EntityVector side_nodes;
    impl::SerialElementDataVector connectedElementDataVector;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            impl::LocalId local_elem_id = get_local_element_id(bucket[j]);
            add_local_graph_edges_for_elem(stk::mesh::MeshIndex(elemBuckets[i],j), local_elem_id, newGraphEdges, coincidentGraphEdges, scratchEntityVector, side_nodes, connectedElementDataVector);
            add_edges_to_graph(newGraphEdges, m_graph);
            add_edges_to_graph(coincidentGraphEdges, m_coincidentGraph);
        }
    }
}

SideNodeToReceivedElementDataMap ElemElemGraph::communicate_shared_sides(impl::ElemSideToProcAndFaceId& elementSidesToSend)
{
     impl::fill_suggested_side_ids(m_bulk_data, elementSidesToSend);
     SharedSidesCommunication sharedSidesCommunication(m_bulk_data, elementSidesToSend);
     SideNodeToReceivedElementDataMap elementSidesReceived = sharedSidesCommunication.get_received_element_sides();
     return elementSidesReceived;
}

void add_shared_edge(const impl::ParallelElementData& elementDataOtherProc, stk::mesh::EntityId elem_id, const unsigned side_index,
        const stk::mesh::EntityVector localElemSideNodes, stk::topology elem_topology, std::vector<impl::SharedEdgeInfo> &newlySharedEdges)
{
    impl::SharedEdgeInfo sharedEdgeInfo;
    impl::GraphEdgeProc graphEdgeProc(elem_id, side_index, elementDataOtherProc.get_element_identifier(), elementDataOtherProc.get_element_side_index(),
            elementDataOtherProc.get_proc_rank_of_neighbor());
    sharedEdgeInfo.set_graph_edge_proc(graphEdgeProc);

    sharedEdgeInfo.m_sharedNodes = localElemSideNodes;
    sharedEdgeInfo.m_remoteElementTopology = elem_topology;
    newlySharedEdges.push_back(sharedEdgeInfo);
}

void ElemElemGraph::create_parallel_graph_edge(const impl::ParallelElementData &localElementData,
                                               const stk::mesh::impl::ParallelElementData &remoteElementData,
                                               impl::ElemSideToProcAndFaceId & elementSidesToSend,
                                               std::vector<impl::SharedEdgeInfo> &newlySharedEdges)
{
    impl::LocalId local_elem_id = localElementData.get_element_local_id();
    impl::LocalId negativeRemoteElemId = -1 * static_cast<impl::LocalId>(remoteElementData.get_element_identifier());
    stk::mesh::GraphEdge graphEdge(local_elem_id, localElementData.get_element_side_index(), negativeRemoteElemId, remoteElementData.get_element_side_index());
    m_graph.add_edge(graphEdge);

    stk::mesh::Entity localElem = get_entity_from_local_id(local_elem_id);
    stk::mesh::impl::ElemSideToProcAndFaceId::const_iterator iter_elem_side_comm = elementSidesToSend.find(impl::EntitySidePair(localElem, localElementData.get_element_side_index()));
    bool did_this_proc_send_info_about_this_side_to_other_proc = iter_elem_side_comm != elementSidesToSend.end();
    // communicate back to originating processor
    if(!did_this_proc_send_info_about_this_side_to_other_proc)
    {
        // for cases of mesh modification, only proc with "created" elementl will communicate
        // new side was sent to this proc (to fish for possible face connections)
        // if fish is caught, need to send info back to original proc so they both will create edge in graph
        add_shared_edge(remoteElementData, m_bulk_data.identifier(localElem), localElementData.get_element_side_index(), localElementData.get_side_nodes(), m_bulk_data.bucket(localElem).topology(), newlySharedEdges);
    }

    impl::ParallelInfo parInfo(remoteElementData.get_proc_rank_of_neighbor(),
                               localElementData.get_permutation(),
                                remoteElementData.get_element_topology());

    m_parallelInfoForGraphEdges.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
}

void ElemElemGraph::fill_parallel_graph(impl::ElemSideToProcAndFaceId & elementSidesToSend, SideNodeToReceivedElementDataMap&  elementSidesReceived)
{
    std::vector<impl::SharedEdgeInfo> newlySharedEdges;
    stk::mesh::EntityVector scratchEntityVector;
    impl::ParallelElementDataVector localElementsAttachedToReceivedNodes;
    for (SideNodeToReceivedElementDataMap::value_type & receivedElementDataForNodes: elementSidesReceived)
    {
        for(const stk::mesh::impl::ParallelElementData &remoteElementData : receivedElementDataForNodes.second)
        {
             get_elements_attached_to_remote_nodes(remoteElementData.get_side_nodes(), remoteElementData.get_element_identifier(),
                                                   remoteElementData.get_element_topology(), scratchEntityVector, localElementsAttachedToReceivedNodes);
            for(const impl::ParallelElementData &localElemAttachedToNodes : localElementsAttachedToReceivedNodes)
                create_parallel_graph_edge(localElemAttachedToNodes, remoteElementData, elementSidesToSend, newlySharedEdges);
        }
    }

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
    impl::fill_element_side_nodes_from_topology(m_bulk_data, localElem, side_index, localElemSideNodes);

    std::pair<bool,unsigned> permutationIfConnected = stk::mesh::side_equivalent(m_bulk_data, localElem, side_index, sideNodes.data());
    ThrowRequireWithSierraHelpMsg(permutationIfConnected.first);

    impl::LocalId local_elem_id = get_local_element_id(localElem);
    impl::LocalId negSgnRemoteElemId = -1 * static_cast<impl::LocalId>(receivedSharedEdge.get_remote_element_global_id());

    GraphEdge graphEdge(local_elem_id, side_index, negSgnRemoteElemId, receivedSharedEdge.get_remote_element_side_index());
    m_graph.add_edge(graphEdge);

    impl::ParallelInfo parInfo(receivedSharedEdge.get_remote_processor_rank(),
                                permutationIfConnected.second,
                                receivedSharedEdge.m_remoteElementTopology);

    m_parallelInfoForGraphEdges.insert_parallel_info_for_graph_edge(graphEdge, parInfo);
}

stk::mesh::EntityId find_side_id(const stk::mesh::impl::ElemSideToProcAndFaceId& elemSideDataSent,
        const stk::mesh::Entity localElem, const unsigned side_index, const int elem_proc_id)
{
    stk::mesh::EntityId chosen_side_id = stk::mesh::InvalidEntityId;
    const auto iterRange = elemSideDataSent.equal_range(impl::EntitySidePair(localElem, side_index));
    for (impl::ElemSideToProcAndFaceId::const_iterator iter = iterRange.first; iter != iterRange.second; ++iter)
    {
        bool is_received_element_data_info_from_proc_in_elem_side_comm = iter->second.proc == elem_proc_id;
        //bool is_received_element_data_info_have_same_side_index_as_in_elem_side_comm = iter->first.side_id == elemData.m_localSideIndex;

        if (is_received_element_data_info_from_proc_in_elem_side_comm)
        {
            chosen_side_id = iter->second.side_id;
            break;
        }
    }
    ThrowRequireWithSierraHelpMsg(chosen_side_id!=stk::mesh::InvalidEntityId);
    return chosen_side_id;
}

stk::topology ElemElemGraph::get_topology_of_remote_element(const GraphEdge &graphEdge)
{
    impl::ParallelInfo &parallel_edge_info = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
    return parallel_edge_info.m_remote_element_toplogy;
}

stk::topology ElemElemGraph::get_topology_of_connected_element(const GraphEdge &graphEdge)
{
    if(graphEdge.elem2()<0)
    {
        return this->get_topology_of_remote_element(graphEdge);
    }
    return m_element_topologies[graphEdge.elem2()];
}

void report_error_with_invalid_ordinal(std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm, const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& side_nodes_vec,
        stk::mesh::Entity element_with_perm_0, stk::mesh::Entity element_with_perm_4)
{
    if(ord_and_perm.first == stk::mesh::INVALID_CONNECTIVITY_ORDINAL)
    {
        std::ostringstream os;
        os << "Proc: " << bulkData.parallel_rank() << std::endl;
        os << "this element: " << bulkData.identifier(element_with_perm_0) << std::endl;
        os << "other element: " << bulkData.identifier(element_with_perm_4) << std::endl;
        os << "Nodes: ";

        for(stk::mesh::Entity side_node : side_nodes_vec)
        {
            os << bulkData.identifier(side_node) << " ";
        }

        os << std::endl;
        std::cerr << os.str();
    }

    ThrowRequireMsg(ord_and_perm.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL, "yikes!");
    ThrowRequireMsg(ord_and_perm.second != stk::mesh::INVALID_PERMUTATION, "yikes!");
}

void ensure_fresh_modifiable_state(stk::mesh::BulkData& bulkData)
{
    if(bulkData.in_modifiable_state())
    {
        bulkData.modification_end();
    }
    bulkData.modification_begin();
}

class RemoteDeathBoundary
{
public:
    RemoteDeathBoundary(stk::mesh::BulkData& bulkData, ElemElemGraph& elementGraph,
        const stk::mesh::EntityVector& killedElements, const stk::mesh::PartVector& parts_for_creating_side, stk::mesh::Part& active, const stk::mesh::PartVector* boundary_mesh_parts) :
            m_bulkData(bulkData), m_elementGraph(elementGraph), m_killedElements(killedElements), m_parts_for_creating_side(parts_for_creating_side), m_active(active),
            m_boundary_mesh_parts(boundary_mesh_parts), m_topology_modified(false)
    {}
    ~RemoteDeathBoundary(){}

    void update_death_boundary_for_remotely_killed_elements(std::vector<stk::mesh::sharing_info> &shared_modified,
                                                            stk::mesh::EntityVector& deletedEntities,
                                                            stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector)
    {
        std::vector<impl::GraphEdgeProc> remote_edges = get_remote_edges();

        for(impl::GraphEdgeProc& re : remote_edges)
        {
            stk::mesh::EntityId local_id = re.get_local_element_global_id();
            int local_side = re.get_local_element_side_index();
            stk::mesh::EntityId remote_id = re.get_remote_element_global_id();
            int remote_side = re.get_remote_element_side_index();

            stk::mesh::Entity element = m_bulkData.get_entity(stk::topology::ELEM_RANK, local_id);

            impl::ParallelInfo &parallel_edge_info = m_elementGraph.get_parallel_edge_info(element, local_side, remote_id, remote_side);
            remoteActiveSelector[-remote_id] = false;

            m_topology_modified = true;

            bool create_side = m_bulkData.bucket(element).member(m_active);
            if(create_side==true)
            {
                impl::add_side_into_exposed_boundary(m_bulkData,
                                                     parallel_edge_info,
                                                     element,
                                                     local_side,
                                                     remote_id,
                                                     m_parts_for_creating_side,
                                                     shared_modified,
                                                     remoteActiveSelector,
                                                     m_boundary_mesh_parts);
            }
            else
            {
                impl::remove_side_from_death_boundary(m_bulkData, element, m_active, deletedEntities, local_side);
            }
        }
    }

    void set_topology_is_modified()
    {
        m_topology_modified = true;
    }

    bool get_topology_modification_status() const
    {
        return m_topology_modified;
    }

private:

    std::vector<impl::GraphEdgeProc> get_remote_edges() const
    {
        std::vector<impl::GraphEdgeProc> elements_to_comm = get_elements_to_communicate();
        return impl::communicate_killed_entities(m_bulkData.parallel(), elements_to_comm);
    }

    std::vector<impl::GraphEdgeProc> get_elements_to_communicate() const
    {
        std::vector<impl::GraphEdgeProc> elements_to_comm;

        for(stk::mesh::Entity this_element :m_killedElements)
        {
            for(size_t j=0;j<m_elementGraph.get_num_connected_elems(this_element);++j)
            {
                if(impl::does_element_have_side(m_bulkData, this_element) && !m_elementGraph.is_connected_elem_locally_owned(this_element, j))
                {
                    impl::IdViaSidePair idViaSidePair = m_elementGraph.get_connected_remote_id_and_via_side(this_element,j);
                    stk::mesh::EntityId other_element_id = idViaSidePair.id;
                    int side1 = idViaSidePair.side;
                    int side2 = m_elementGraph.get_connected_elements_side(this_element,j);
                    int other_proc = m_elementGraph.get_owning_proc_id_of_remote_element(this_element, j);
                    elements_to_comm.push_back(impl::GraphEdgeProc(m_bulkData.identifier(this_element), side1, other_element_id, side2, other_proc));
                }
            }
        }

        return elements_to_comm;
    }

    stk::mesh::BulkData& m_bulkData;
    ElemElemGraph& m_elementGraph;
    const stk::mesh::EntityVector& m_killedElements;
    const stk::mesh::PartVector& m_parts_for_creating_side;
    stk::mesh::Part& m_active;
    const stk::mesh::PartVector* m_boundary_mesh_parts;
    bool m_topology_modified;
};


bool process_killed_elements(stk::mesh::BulkData& bulkData,
                             ElemElemGraph& elementGraph,
                             const stk::mesh::EntityVector& killedElements,
                             stk::mesh::Part& active,
                             stk::mesh::impl::ParallelSelectedInfo &remoteActiveSelector,
                             const stk::mesh::PartVector& parts_for_creating_side,
                             const stk::mesh::PartVector* boundary_mesh_parts)
{
    ensure_fresh_modifiable_state(bulkData);
    impl::create_sides_created_during_death_part(bulkData.mesh_meta_data());

    std::vector<stk::mesh::sharing_info> shared_modified;
    stk::mesh::EntityVector deletedEntities;

    RemoteDeathBoundary remote_death_boundary(bulkData, elementGraph, killedElements, parts_for_creating_side, active, boundary_mesh_parts);
    remote_death_boundary.update_death_boundary_for_remotely_killed_elements(shared_modified, deletedEntities, remoteActiveSelector);

    std::vector<impl::ElementSidePair> element_side_pairs;
    element_side_pairs.reserve(impl::get_element_side_multiplier() * killedElements.size());

    for(size_t k = 0; k < killedElements.size(); ++k)
    {
        stk::mesh::Entity this_element = killedElements[k];

        for(size_t j = 0; j < elementGraph.get_num_connected_elems(this_element); ++j)
        {
            if(impl::does_element_have_side(bulkData, this_element))
            {
                remote_death_boundary.set_topology_is_modified();
                if(elementGraph.is_connected_elem_locally_owned(this_element, j))
                {
                    impl::ElementViaSidePair other_element_via_side = elementGraph.get_connected_element_and_via_side(this_element, j);
                    stk::mesh::Entity other_element = other_element_via_side.element;
                    if(impl::does_element_have_side(bulkData, other_element_via_side.element))
                    {
                        int side_id = other_element_via_side.side;
                        ThrowRequireWithSierraHelpMsg(side_id != -1);

                        bool is_other_element_alive = bulkData.bucket(other_element).member(active);
                        if(is_other_element_alive)
                        {
                            stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulkData, this_element, side_id);

                            if(bulkData.is_valid(side))
                            {
                                if(bulkData.bucket(side).owned())
                                {
                                    stk::mesh::ConstPartVector parts = impl::get_stk_parts_for_moving_parts_into_death_boundary(boundary_mesh_parts);
                                    bulkData.change_entity_parts(side, parts);
                                }
                            }
                            else
                            {
                                stk::mesh::PartVector parts = impl::get_parts_for_creating_side(bulkData, parts_for_creating_side, other_element, side_id);

//                                stk::mesh::EntityId side_global_id = elementGraph.get_available_side_id();
//
//                                stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
//                                ThrowRequireWithSierraHelpMsg(!impl::is_id_already_in_use_locally(bulkData, side_rank, side_global_id));

                                // switch elements
                                stk::mesh::Entity element_with_perm_0 = other_element;
                                stk::mesh::Entity element_with_perm_4 = this_element;

                                int side_id_needed = elementGraph.get_connected_elements_side(this_element, j);

                                ThrowRequireMsg(side_id_needed >= 0, "ERROR: proc " << bulkData.parallel_rank() << " found side_id_needed=" << side_id_needed
                                                << " between elem " << bulkData.identifier(element_with_perm_0)<< " and " << bulkData.identifier(element_with_perm_4)
                                                << " in elem-elem-graph");

                                side = bulkData.declare_element_side(element_with_perm_0, side_id_needed, parts);
//                                side = bulkData.declare_element_side(side_global_id, element_with_perm_0, side_id_needed, parts);
//
//                                const stk::mesh::Entity* side_nodes = bulkData.begin_nodes(side);
//                                unsigned num_side_nodes = bulkData.num_nodes(side);
//                                stk::mesh::EntityVector side_nodes_vec(side_nodes, side_nodes + num_side_nodes);
//
//                                std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm =
//                                        stk::mesh::get_ordinal_and_permutation(bulkData, element_with_perm_4, side_rank, side_nodes_vec);
//
//                                report_error_with_invalid_ordinal(ord_and_perm, bulkData, side_nodes_vec, element_with_perm_0, element_with_perm_4);
//
//                                bulkData.declare_relation(element_with_perm_4, side, ord_and_perm.first, ord_and_perm.second);
                            }
                        }
                        else
                        {
                            impl::remove_side_from_death_boundary(bulkData, this_element, active, deletedEntities, side_id);
                        }
                    }
                }
                else
                {
                    impl::IdViaSidePair remote_id_side_pair = elementGraph.get_connected_remote_id_and_via_side(this_element, j);
                    stk::mesh::EntityId remote_id = remote_id_side_pair.id;
                    int remote_side = elementGraph.get_connected_elements_side(this_element, j);
                    impl::ParallelInfo &parallel_edge_info = elementGraph.get_parallel_edge_info(this_element, remote_id_side_pair.side, remote_id, remote_side);
                    bool other_element_active = remoteActiveSelector[-remote_id];
                    bool create_side = other_element_active;

                    if(create_side)
                    {
                        impl::add_side_into_exposed_boundary(bulkData, parallel_edge_info, this_element, remote_id_side_pair.side, remote_id, parts_for_creating_side,
                                shared_modified, remoteActiveSelector, boundary_mesh_parts);
                    }
                    else
                    {
                        int side_id = remote_id_side_pair.side;
                        ThrowRequireWithSierraHelpMsg(side_id != -1);
                        impl::remove_side_from_death_boundary(bulkData, this_element, active, deletedEntities, side_id);
                    }
                }
            }
        }
    }
    stk::mesh::impl::delete_entities_and_upward_relations(bulkData, deletedEntities);
    bulkData.make_mesh_parallel_consistent_after_element_death(shared_modified, deletedEntities, elementGraph, killedElements, active);
    return remote_death_boundary.get_topology_modification_status();
}

stk::mesh::SideConnector ElemElemGraph::get_side_connector()
{
    return stk::mesh::SideConnector(m_bulk_data, m_graph, m_coincidentGraph, m_idMapper);
}

stk::mesh::SideNodeConnector ElemElemGraph::get_side_node_connector()
{
    return stk::mesh::SideNodeConnector(m_bulk_data, m_graph, m_coincidentGraph, m_parallelInfoForGraphEdges, m_idMapper);
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
            ThrowRequireWithSierraHelpMsg(is_valid_graph_element(elem_to_delete));
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

void ElemElemGraph::communicate_shell_connectivity(std::vector<impl::ShellConnectivityData>& shellConnectivityList,
                                                   const std::vector<impl::ElementSidePair> &deletedShells) {

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
}

void ElemElemGraph::delete_local_connections_and_collect_remote(const stk::mesh::impl::DeletedElementInfoVector &elements_to_delete,
                                                                std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem)
{
    for (const stk::mesh::impl::DeletedElementInfo &deletedElementInfo : elements_to_delete) {
        stk::mesh::Entity elem_to_delete = deletedElementInfo.entity;
        impl::LocalId elem_to_delete_id = get_local_element_id(elem_to_delete);

        size_t num_connected_elems = m_graph.get_num_edges_for_element(elem_to_delete_id);
        for (int conn_elem_index = num_connected_elems - 1; conn_elem_index >= 0; --conn_elem_index) {
            const GraphEdge & graphEdge = m_graph.get_edge_for_element(elem_to_delete_id, conn_elem_index);
            impl::LocalId connected_elem_id = graphEdge.elem2();

            bool local_connection = connected_elem_id >= 0;
            if (local_connection) {
                m_graph.delete_edge(create_symmetric_edge(graphEdge));
            } else {
                impl::ParallelInfo& parallelInfo = m_parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
                int remote_proc = parallelInfo.get_proc_rank_of_neighbor();

                stk::mesh::EntityId connected_global_id = -connected_elem_id;
                impl::DeletedElementData deletedElementData;
                deletedElementData.m_deletedElement = elem_to_delete_id;
                deletedElementData.m_remoteElement = connected_global_id;
                deletedElementData.m_remoteProc = remote_proc;
                local_elem_and_remote_connected_elem.push_back(deletedElementData);
                m_parallelInfoForGraphEdges.erase_parallel_info_for_graph_edge(graphEdge);
            }
        }
        for(const stk::mesh::GraphEdge &graphEdge : m_coincidentGraph.get_edges_for_element(elem_to_delete_id))
        {
            m_coincidentGraph.delete_edge(create_symmetric_edge(graphEdge));
        }

        m_graph.delete_all_edges(elem_to_delete_id);
        m_coincidentGraph.delete_all_edges(elem_to_delete_id);
    }
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

void ElemElemGraph::delete_remote_connections(const std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges)
{
    for (const std::pair<stk::mesh::EntityId, stk::mesh::EntityId> & edge : remote_edges) {
        stk::mesh::EntityId deleted_elem_global_id = edge.first;
        stk::mesh::EntityId connected_elem_global_id = edge.second;
        stk::mesh::Entity connected_elem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, connected_elem_global_id);
        if (!is_valid_graph_element(connected_elem)) {
            continue;
        }
        impl::LocalId connected_elem_id = m_idMapper.entity_to_local(connected_elem);
        size_t num_conn_elem = m_graph.get_num_edges_for_element(connected_elem_id);
        bool found_deleted_elem = false;
        for (size_t conn_elem_index = 0; conn_elem_index < num_conn_elem; ++conn_elem_index) {
            const GraphEdge & graphEdge = m_graph.get_edge_for_element(connected_elem_id, conn_elem_index);
            if (graphEdge.elem2() == static_cast<int64_t>(-deleted_elem_global_id)) {
                m_parallelInfoForGraphEdges.erase_parallel_info_for_graph_edge(graphEdge);
                m_graph.delete_edge_from_graph(connected_elem_id, conn_elem_index);
                found_deleted_elem = true;
                break;
            }
        }
        ThrowRequireWithSierraHelpMsg(found_deleted_elem);
    }
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

    impl::ElemSideToProcAndFaceId shellNeighborsToReconnect;
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
                        m_graph.add_edge(stk::mesh::GraphEdge(localElementId, data.m_nearElementSide, remoteElementId, data.m_farElementSide));
                    }
                }
                else {
                    shellNeighborsToReconnect.insert(
                            std::pair<impl::EntitySidePair, impl::ProcFaceIdPair>(impl::EntitySidePair(localElement, data.m_nearElementSide),
                                    impl::ProcFaceIdPair(data.m_farElementProc, 0)));
                }
            }
        }
    }

    SideNodeToReceivedElementDataMap elementSidesReceived = communicate_shared_sides(shellNeighborsToReconnect);
    fill_parallel_graph(shellNeighborsToReconnect, elementSidesReceived);
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
    stk::mesh::impl::DeletedElementInfoVector elements_to_delete = filter_delete_elements_argument(elements_to_delete_argument);

    std::vector<impl::ShellConnectivityData> shellConnectivityList;
    std::vector<impl::ElementSidePair> deletedNonCoincidentShells;
    collect_local_shell_connectivity_data(elements_to_delete, shellConnectivityList, deletedNonCoincidentShells);

    communicate_shell_connectivity(shellConnectivityList, deletedNonCoincidentShells);

    std::vector<impl::DeletedElementData> local_elem_and_remote_connected_elem;
    delete_local_connections_and_collect_remote(elements_to_delete, local_elem_and_remote_connected_elem);

    std::vector< std::pair< stk::mesh::EntityId, stk::mesh::EntityId > > remote_edges;
    communicate_remote_connections_to_delete(local_elem_and_remote_connected_elem, remote_edges);

    clear_deleted_element_connections(elements_to_delete);

    delete_remote_connections(remote_edges);

    unsigned localNumShells = shellConnectivityList.size();
    unsigned globalMaxShells = 0;
    stk::all_reduce_max(m_bulk_data.parallel(), &localNumShells, &globalMaxShells, 1);
    if (globalMaxShells > 0) {
        reconnect_volume_elements_across_deleted_shells(shellConnectivityList);
    }

    update_number_of_parallel_edges();
}

template <typename GraphType>
void add_both_edges_between_local_elements(const GraphEdge& graphEdge, GraphType &graph)
{
    graph.add_edge(graphEdge);
    graph.add_edge(create_symmetric_edge(graphEdge));
}

template <typename GraphType>
void ElemElemGraph::add_new_local_edges_to_graph(GraphType &graph, const std::vector<stk::mesh::GraphEdge> &newGraphEdges)
{
    for (const stk::mesh::GraphEdge& graphEdge : newGraphEdges)
    {
        stk::mesh::Entity neighbor = m_idMapper.local_to_entity(graphEdge.elem2());
        if (is_valid_graph_element(neighbor))
            add_both_edges_between_local_elements(graphEdge, graph);
    }
}

void ElemElemGraph::add_local_edges(stk::mesh::Entity elem_to_add, impl::LocalId new_elem_id,
                                    stk::mesh::EntityVector& scratchEntityVector,
                                    stk::mesh::EntityVector& side_nodes,
                                    impl::SerialElementDataVector& connectedElementDataVector)
{
    std::vector<stk::mesh::GraphEdge> newGraphEdges;
    std::vector<stk::mesh::GraphEdge> coincidentGraphEdges;
    bool only_consider_upper_symmetry = false;
    add_local_graph_edges_for_elem(m_bulk_data.mesh_index(elem_to_add), new_elem_id, newGraphEdges, coincidentGraphEdges, scratchEntityVector, side_nodes, connectedElementDataVector, only_consider_upper_symmetry);
    add_new_local_edges_to_graph(m_graph, newGraphEdges);
    add_new_local_edges_to_graph(m_coincidentGraph, coincidentGraphEdges);
}

void ElemElemGraph::add_vertex(impl::LocalId newElemLocalId, stk::mesh::Entity elem)
{
    m_idMapper.add_new_entity_with_local_id(elem, newElemLocalId);

    stk::topology elem_topology = m_bulk_data.bucket(elem).topology();
    m_element_topologies[newElemLocalId] = elem_topology;
}

stk::mesh::EntityVector ElemElemGraph::filter_add_elements_arguments(const stk::mesh::EntityVector& allUnfilteredElementsNotAlreadyInGraph) const
{
    stk::mesh::EntityVector allElementsNotAlreadyInGraph;
    allElementsNotAlreadyInGraph.reserve(allUnfilteredElementsNotAlreadyInGraph.size());
    for(stk::mesh::Entity element : allUnfilteredElementsNotAlreadyInGraph)
    {
        if(m_bulk_data.is_valid(element) && m_bulk_data.bucket(element).owned())
        {
            ThrowRequire(!is_valid_graph_element(element));
            allElementsNotAlreadyInGraph.push_back(element);
        }
    }
    return allElementsNotAlreadyInGraph;
}

void ElemElemGraph::add_elements_locally(const stk::mesh::EntityVector& allElementsNotAlreadyInGraph)
{
    m_idMapper.make_space_for_new_elements(allElementsNotAlreadyInGraph);
    stk::mesh::EntityVector scratchEntityVector;
    stk::mesh::EntityVector side_nodes;
    impl::SerialElementDataVector connectedElementDataVector;
    for(stk::mesh::Entity newElem : allElementsNotAlreadyInGraph)
    {
        impl::LocalId newElemLocalId = get_new_local_element_id_from_pool();
        add_vertex(newElemLocalId, newElem);
        add_local_edges(newElem, newElemLocalId, scratchEntityVector, side_nodes, connectedElementDataVector);
    }
}

void ElemElemGraph::add_elements(const stk::mesh::EntityVector &allUnfilteredElementsNotAlreadyInGraph)
{
    stk::mesh::EntityVector allElementsNotAlreadyInGraph = filter_add_elements_arguments(allUnfilteredElementsNotAlreadyInGraph);

    const size_t numEdgesBefore = num_edges();
    add_elements_locally(allElementsNotAlreadyInGraph);
    size_t numNewSideIdsNeeded = num_edges() - numEdgesBefore;
    numNewSideIdsNeeded += get_max_num_sides_per_element() * allUnfilteredElementsNotAlreadyInGraph.size();

    impl::ElemSideToProcAndFaceId elem_side_comm = impl::gather_element_side_ids_to_send(m_bulk_data);

    stk::mesh::EntityVector allElementsNotAlreadyInGraph_copy = allElementsNotAlreadyInGraph;
    std::sort(allElementsNotAlreadyInGraph_copy.begin(), allElementsNotAlreadyInGraph_copy.end());

    std::set< stk::mesh::Entity > addedShells;
    impl::ElemSideToProcAndFaceId only_added_elements;
    for(impl::ElemSideToProcAndFaceId::value_type &elemSideToProcAndFaceId : elem_side_comm)
    {
        stk::mesh::Entity element = elemSideToProcAndFaceId.first.entity;
        stk::mesh::EntityVector::iterator elem_iter = std::lower_bound(allElementsNotAlreadyInGraph_copy.begin(), allElementsNotAlreadyInGraph_copy.end(), element);
        if(elem_iter!=allElementsNotAlreadyInGraph_copy.end() && *elem_iter==element)
        {
            only_added_elements.insert(elemSideToProcAndFaceId);
            if (m_bulk_data.bucket(element).topology().is_shell())
            {
                addedShells.insert(element);
            }
        }
    }

    SideNodeToReceivedElementDataMap elementSidesReceived = communicate_shared_sides(only_added_elements);
    fill_parallel_graph(only_added_elements, elementSidesReceived);

    GraphInfo graphInfo(m_graph, m_parallelInfoForGraphEdges, m_element_topologies);
    remove_graph_edges_blocked_by_shell(graphInfo);

    stk::mesh::EntityVector addedShellsVector;
    for (auto &shell : addedShells)
    {
        addedShellsVector.push_back(shell);
    }

    stk::CommSparse comm(m_bulk_data.parallel());
    for (int phase = 0; phase < 2; ++phase)
    {
        pack_remote_edge_across_shell(comm, addedShellsVector, phase);
        if (0 == phase)
        {
            comm.allocate_buffers();
        }
        if (1 == phase)
        {
            comm.communicate();
        }
    }

    unpack_remote_edge_across_shell(comm);
    update_number_of_parallel_edges();
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
                m_graph.delete_edge_from_graph(leftId, index);
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
            if (leftId == impl::INVALID_LOCAL_ID)
            {
                leftId = graphEdge.elem2();
                continue;
            }
            if (rightId == impl::INVALID_LOCAL_ID)
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
                    m_graph.delete_edge_from_graph(localId, index);
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
    ThrowRequireWithSierraHelpMsg(sideOrd != -1);

    stk::mesh::Entity side = bulkData.declare_element_side(local_element,
                                                           sideOrd,
                                                           parts_for_creating_side);

    int other_proc = parallel_edge_info.get_proc_rank_of_neighbor();
    int owning_proc = std::min(other_proc, bulkData.parallel_rank());
    if(bulkData.state(side) != stk::mesh::Created)
        owning_proc = bulkData.parallel_owner_rank(side);
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

// SideSetEntry
// extract_skinned_sideset
// SideSetEntryLess
// create_side_entities_given_sideset

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

}} // end namespaces stk mesh

