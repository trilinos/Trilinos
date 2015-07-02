#include "ElemElemGraph.hpp"
#include "ElemElemGraphImpl.hpp"

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace stk { namespace mesh {

const impl::LocalId ElemElemGraph::INVALID_LOCAL_ID = std::numeric_limits<impl::LocalId>::max();

ElemElemGraph::ElemElemGraph(stk::mesh::BulkData& bulkData, const stk::mesh::Part &part) : m_bulk_data(bulkData), m_part(part)
{
    int numElems = size_data_members();

    impl::ElemSideToProcAndFaceId elem_side_comm;
    if (numElems > 0) {
        impl::fill_local_ids_and_fill_element_entities_and_topologies(m_bulk_data, m_local_id_to_element_entity, m_entity_to_local_id, m_element_topologies);
        fill_graph();

        elem_side_comm = impl::get_element_side_ids_to_communicate(bulkData);
    }

    size_t num_side_ids_needed = elem_side_comm.size();

    num_side_ids_needed += m_num_edges;

    bulkData.generate_new_ids(bulkData.mesh_meta_data().side_rank(), num_side_ids_needed, m_suggested_side_ids);

    fill_parallel_graph(elem_side_comm, part);

    update_number_of_parallel_edges();

}

void ElemElemGraph::update_number_of_parallel_edges()
{
    // Figure out the real number of sides that can be generated, using
    // the rule of one side entity per element side (with multiple
    // coincident elements each being connected to the same side)
    m_num_parallel_edges = 0;
    for(size_t i = 0; i < m_elem_graph.size(); ++i)
    {
        std::set<int> uniqueRemoteOrdinals;
        const std::vector<impl::LocalId> & localElement = m_elem_graph[i];
        for(size_t j = 0; j < localElement.size(); ++j)
        {
            if(localElement[j] < 0)
            {
                // Connected to remote element through this side
                uniqueRemoteOrdinals.insert(m_via_sides[i][j]);
            }
        }
        m_num_parallel_edges += uniqueRemoteOrdinals.size();
    }
}

std::vector<stk::mesh::EntityId> ElemElemGraph::get_suggested_side_ids() const
{
    std::vector<stk::mesh::EntityId> available_ids;
    available_ids.assign(m_suggested_side_ids.begin()+m_num_ids_used, m_suggested_side_ids.end());
    return available_ids;
}

void ElemElemGraph::set_num_side_ids_used(size_t num_used)
{
    m_num_ids_used += num_used;
    ThrowRequireMsg(m_suggested_side_ids.size() >= m_num_ids_used, "Program error (exceeded available ids). Contact sierra-help@sandia.gov for support.");
}

ElemElemGraph::~ElemElemGraph() {}

size_t ElemElemGraph::get_num_connected_elems(stk::mesh::Entity local_element) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return m_elem_graph[local_id].size();
}

bool ElemElemGraph::is_connected_elem_locally_owned(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return m_elem_graph[local_id][index_conn_elem] >= 0;
}

int ElemElemGraph::get_side_id_to_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    return m_via_sides[local_id][index_conn_elem];
}

stk::mesh::Entity ElemElemGraph::get_connected_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    impl::LocalId other_element_id = m_elem_graph[local_id][index_conn_elem];
    return m_local_id_to_element_entity[other_element_id];
}

stk::mesh::EntityId ElemElemGraph::get_entity_id_of_remote_element(stk::mesh::Entity local_element, size_t index_conn_elem) const
{
    ThrowRequireMsg(!is_connected_elem_locally_owned(local_element, index_conn_elem) , "Program error. Contact sierra-help@sandia.gov for support.");
    impl::LocalId local_id = get_local_element_id(local_element);
    stk::mesh::EntityId id = -m_elem_graph[local_id][index_conn_elem];
    return id;
}

int ElemElemGraph::get_owning_proc_id_of_remote_element(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
{
    impl::LocalId local_id = get_local_element_id(local_element);
    impl::ParallelGraphInfo::const_iterator iter = m_parallel_graph_info.find(std::make_pair(local_id, other_element_id));
    ThrowRequireMsg( iter != m_parallel_graph_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
    int other_proc = iter->second.m_other_proc;
    return other_proc;
}

int ElemElemGraph::get_side_from_element1_to_remote_element2(stk::mesh::Entity local_element, stk::mesh::EntityId other_element_id) const
{
    impl::LocalId remote_element_local_id = -other_element_id;
    impl::LocalId element1_local_id = get_local_element_id(local_element);

    int side = -1;
    const std::vector<impl::LocalId>& conn_elements = m_elem_graph[element1_local_id];

    std::vector<impl::LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), remote_element_local_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = m_via_sides[element1_local_id][index];
    }
    return side;
}

int ElemElemGraph::get_side_from_element1_to_locally_owned_element2(stk::mesh::Entity local_element, stk::mesh::Entity other_element) const
{
    impl::LocalId other_element_id = get_local_element_id(other_element);
    impl::LocalId element1_local_id = get_local_element_id(local_element);

    int side = -1;
    const std::vector<impl::LocalId>& conn_elements = m_elem_graph[element1_local_id];

    std::vector<impl::LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), other_element_id);
    if ( iter != conn_elements.end() )
    {
        int64_t index = iter - conn_elements.begin();
        side = m_via_sides[element1_local_id][index];
    }
    return side;
}

bool ElemElemGraph::is_connected_to_other_element_via_side_ordinal(stk::mesh::Entity element, int sideOrdinal) const
{
    impl::LocalId elementLocalId = get_local_element_id(element);
    auto iterOfSide = std::find(m_via_sides[elementLocalId].begin(), m_via_sides[elementLocalId].end(), sideOrdinal);
    return iterOfSide != m_via_sides[elementLocalId].end();
}

size_t ElemElemGraph::num_edges() const
{
    return m_num_edges;
}

impl::parallel_info& ElemElemGraph::get_parallel_edge_info(stk::mesh::Entity element, stk::mesh::EntityId remote_id)
{
    impl::LocalId this_elem_local_id = get_local_element_id(element);

    impl::ParallelGraphInfo::iterator iter = m_parallel_graph_info.find(std::make_pair(this_elem_local_id, remote_id));
    ThrowRequireMsg( iter != m_parallel_graph_info.end(), "ERROR: Proc " << m_bulk_data.parallel_rank() << " failed to find parallel graph info for <"
                     <<m_bulk_data.identifier(element)<<","<<remote_id<<">");
    return iter->second;
}

impl::LocalId ElemElemGraph::get_local_element_id(stk::mesh::Entity local_element, bool require_valid_id) const
{
    ThrowRequireMsg(m_entity_to_local_id.size() > local_element.local_offset(),"Program error. Contact sierra-help@sandia.gov for support.");
    impl::LocalId local_id = m_entity_to_local_id[local_element.local_offset()];
    if (require_valid_id)
    {
        ThrowRequireMsg(local_id != INVALID_LOCAL_ID, "Program error. Contact sierra-help@sandia.gov for support.");
    }
    return local_id;
}

int ElemElemGraph::size_data_members()
{
    std::vector<unsigned> counts(stk::topology::NUM_RANKS,0);
    stk::mesh::count_entities(m_bulk_data.mesh_meta_data().locally_owned_part(), m_bulk_data, counts);
    int numElems = counts[stk::topology::ELEM_RANK];

    m_elem_graph.resize(numElems);
    m_via_sides.resize(numElems);
    m_local_id_to_element_entity.resize(numElems, Entity());
    m_entity_to_local_id.resize(m_bulk_data.m_entity_keys.size(), INVALID_LOCAL_ID);
    m_element_topologies.resize(numElems);
    m_num_edges = 0;
    m_num_parallel_edges = 0;
    m_local_id_in_pool.resize(numElems, false);
    m_num_ids_used = 0;

    return numElems;
}

void ElemElemGraph::ensure_space_in_entity_to_local_id(size_t max_index)
{
    size_t needed = max_index + 1;

    if (m_elem_graph.size() < needed)
    {
        m_entity_to_local_id.resize(needed, INVALID_LOCAL_ID);
    }
}

void ElemElemGraph::get_element_side_pairs(const stk::mesh::MeshIndex &meshIndex, impl::LocalId local_elem_id, std::vector<impl::ElementSidePair> &elem_side_pairs) const
{
    stk::mesh::EntityVector side_nodes;
    stk::mesh::EntityVector connected_elements;
    stk::mesh::Entity element = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
    stk::topology topology = meshIndex.bucket->topology();
    int num_sides = topology.num_sides();
    const stk::mesh::Entity* elem_nodes = meshIndex.bucket->begin_nodes(meshIndex.bucket_ordinal);
    elem_side_pairs.clear();
    for(int side_index=0; side_index<num_sides; ++side_index)
    {
        unsigned num_side_nodes = topology.side_topology(side_index).num_nodes();
        side_nodes.resize(num_side_nodes);
        topology.side_nodes(elem_nodes, side_index, side_nodes.begin());
        connected_elements.clear();
        impl::ConnectedElementDataVector connectedElementDataVector;
        impl::find_locally_owned_elements_these_nodes_have_in_common(m_bulk_data, num_side_nodes, side_nodes.data(), connected_elements);
        add_local_elements_to_connected_list(connected_elements, side_nodes, connectedElementDataVector);
        impl::filter_for_candidate_elements_to_connect(m_bulk_data, element, side_index, connectedElementDataVector);

        for (const impl::ConnectedElementData & elemData: connectedElementDataVector)
        {
            if (local_elem_id != elemData.m_elementId)
            {
                elem_side_pairs.push_back(std::make_pair(elemData.m_elementId,side_index));
            }
        }
    }
    std::sort(elem_side_pairs.begin(), elem_side_pairs.end());
    std::vector<impl::ElementSidePair>::iterator new_end = std::unique(elem_side_pairs.begin(), elem_side_pairs.end());
    elem_side_pairs.resize(new_end - elem_side_pairs.begin());
}

void ElemElemGraph::fill_graph()
{
    const stk::mesh::BucketVector& elemBuckets = m_bulk_data.get_buckets(stk::topology::ELEM_RANK, m_bulk_data.mesh_meta_data().locally_owned_part());
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        std::vector<impl::ElementSidePair> elem_side_pairs;

        for(size_t j=0; j<bucket.size(); ++j)
        {
            impl::LocalId local_elem_id = get_local_element_id(bucket[j]);
            get_element_side_pairs(stk::mesh::MeshIndex(elemBuckets[i],j), local_elem_id, elem_side_pairs);
            for(size_t index=0; index<elem_side_pairs.size(); ++index)
            {
                m_elem_graph[local_elem_id].push_back(elem_side_pairs[index].first);
                m_via_sides[local_elem_id].push_back(elem_side_pairs[index].second);
                ++m_num_edges;
            }
        }
    }
}

void ElemElemGraph::fill_parallel_graph(impl::ElemSideToProcAndFaceId& elem_side_comm,
                                        const stk::mesh::Part &part)
{
    stk::mesh::EntityVector elements_to_ignore;
    fill_parallel_graph(elem_side_comm, elements_to_ignore);
}

void ElemElemGraph::fill_parallel_graph(impl::ElemSideToProcAndFaceId& elem_side_comm, const stk::mesh::EntityVector & elements_to_ignore)
{
    stk::CommSparse comm(m_bulk_data.parallel());

    impl::pack_shared_side_nodes_of_elements(comm, m_bulk_data, elem_side_comm, this->get_suggested_side_ids(), m_part);
    comm.allocate_buffers();

    size_t num_edge_ids_used = impl::pack_shared_side_nodes_of_elements(comm, m_bulk_data, elem_side_comm, this->get_suggested_side_ids(), m_part);
    this->set_num_side_ids_used(num_edge_ids_used);
    comm.communicate();

    stk::mesh::impl::ConnectedElementDataVector communicatedElementDataVector;
    std::map<EntityVector, stk::mesh::impl::ConnectedElementDataVector> communicatedElementDataMap;

    for(int proc_id=0; proc_id<m_bulk_data.parallel_size(); ++proc_id)
    {
        if (proc_id != m_bulk_data.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                stk::mesh::impl::ConnectedElementData elementData;
                elementData.m_procId = proc_id;
                comm.recv_buffer(proc_id).unpack<impl::LocalId>(elementData.m_elementId);
                comm.recv_buffer(proc_id).unpack<stk::topology>(elementData.m_elementTopology);
                comm.recv_buffer(proc_id).unpack<unsigned>(elementData.m_sideIndex);
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(elementData.m_suggestedFaceId);
                comm.recv_buffer(proc_id).unpack<bool>(elementData.m_isInPart);

                unsigned num_side_nodes = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(num_side_nodes);
                elementData.m_sideNodes.resize(num_side_nodes);
                for(unsigned i=0; i<num_side_nodes; ++i)
                {
                    stk::mesh::EntityKey key;
                    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(key);
                    elementData.m_sideNodes[i] = m_bulk_data.get_entity(key);
                }

                stk::mesh::EntityVector sortedSideNodes = elementData.m_sideNodes;
                std::sort(sortedSideNodes.begin(), sortedSideNodes.end());
                communicatedElementDataMap[sortedSideNodes].push_back(elementData);
            }
        }
    }

    std::vector<impl::SharedEdgeInfo> newlySharedEdges;

    for (auto & elementData: communicatedElementDataMap)
    {
        add_possibly_connected_elements_to_graph_using_side_nodes(elem_side_comm, elementData.second, elements_to_ignore, &newlySharedEdges);
    }

    std::vector<impl::SharedEdgeInfo> receivedSharedEdges;
    communicate_remote_edges_for_pre_existing_graph_nodes(newlySharedEdges, receivedSharedEdges);

    for (const auto &receivedSharedEdge : receivedSharedEdges)
    {
        connect_remote_element_to_existing_graph( receivedSharedEdge);
    }
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
                remoteEdgeInfo.m_procId = proc_id;

                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(remoteEdgeInfo.m_remoteElementId);
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(remoteEdgeInfo.m_locaElementlId);
                comm.recv_buffer(proc_id).unpack<unsigned>(remoteEdgeInfo.m_sideIndex);
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(remoteEdgeInfo.m_chosenSideId);
                comm.recv_buffer(proc_id).unpack<bool>(remoteEdgeInfo.m_isInPart);
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
    // If received element is a shell, it was pre-existing.  Thus, we don't have to unhook the local element graph
    // node from any remote element that has a face that shell is coincident with.

    const stk::mesh::EntityVector &sideNodes = receivedSharedEdge.m_sharedNodes;  // Pick a convenient set of side nodes (same for all elems)
    unsigned num_side_nodes = sideNodes.size();


    stk::mesh::Entity localElem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, receivedSharedEdge.m_locaElementlId);
    stk::topology localElemTopology = m_bulk_data.bucket(localElem).topology();

    const stk::mesh::Entity* localElemNodes = m_bulk_data.begin_nodes(localElem);
    unsigned localElemNumSides = localElemTopology.num_sides();
    for(unsigned side_index=0; side_index<localElemNumSides; ++side_index)
    {
        unsigned num_nodes_this_side = localElemTopology.side_topology(side_index).num_nodes();
        if (num_nodes_this_side != num_side_nodes)
        {
            continue;
        }
        stk::mesh::EntityVector localElemSideNodes(num_nodes_this_side);
        localElemTopology.side_nodes(localElemNodes, side_index, localElemSideNodes.begin());

        std::pair<bool,unsigned> result = localElemTopology.side_topology(side_index).equivalent(localElemSideNodes, sideNodes);
        bool same_type_of_element = localElemTopology == receivedSharedEdge.m_remoteElementTopology;
        bool negative_permutation = result.second < localElemTopology.side_topology(side_index).num_positive_permutations();
        if ((result.first) && (same_type_of_element == negative_permutation))
        {
            impl::LocalId local_elem_id = get_local_element_id(localElem);
            impl::LocalId negSgnRemoteElemId = -1*receivedSharedEdge.m_remoteElementId;

            m_elem_graph[local_elem_id].push_back(negSgnRemoteElemId);
            m_via_sides[local_elem_id].push_back(side_index);
            ++m_num_edges;
            stk::mesh::EntityId chosen_side_id = receivedSharedEdge.m_chosenSideId;

            m_parallel_graph_info.insert(std::make_pair(std::make_pair(local_elem_id, receivedSharedEdge.m_remoteElementId),
                                                        impl::parallel_info(receivedSharedEdge.m_procId, receivedSharedEdge.m_sideIndex, result.second, chosen_side_id, receivedSharedEdge.m_isInPart)));
            break;
        }
    }

    // impl::break_volume_element_connections_across_shells(localElementsConnectedToRemoteShell, m_elem_graph, m_via_sides);
}

void ElemElemGraph::filter_for_elements_in_graph(stk::mesh::EntityVector &localElements)
{
    stk::mesh::EntityVector filteredElements;
    for (stk::mesh::Entity &element : localElements)
    {
        if (is_valid_graph_element(element))
        {
            filteredElements.push_back(element);
        }
    }
    localElements.swap(filteredElements);
}

void ElemElemGraph::add_possibly_connected_elements_to_graph_using_side_nodes( const stk::mesh::impl::ElemSideToProcAndFaceId& elemSideComm,
                                                                               stk::mesh::impl::ConnectedElementDataVector & communicatedElementDataVector,
                                                                               const stk::mesh::EntityVector & elements_to_ignore,
                                                                               std::vector<impl::SharedEdgeInfo> *newlySharedEdges)
{
    if (communicatedElementDataVector.empty()) {
        return;  // Nothing to connect to
    }

    const stk::mesh::EntityVector sideNodes = communicatedElementDataVector[0].m_sideNodes;  // Pick a convenient set of side nodes (same for all elems)
    unsigned num_side_nodes = sideNodes.size();
    std::set<EntityId> localElementsConnectedToRemoteShell;
    stk::mesh::EntityVector localElements;
    impl::find_locally_owned_elements_these_nodes_have_in_common(m_bulk_data, num_side_nodes, sideNodes.data(), localElements);

    // Trim down the list of local elements, ignoring elements in the provided list (since, e.g., these elements may be deleted soon)
    stk::mesh::EntityVector localElementsToConsider(localElements.size());
    stk::mesh::EntityVector sortedElementsToIgnore(elements_to_ignore);
    std::sort(sortedElementsToIgnore.begin(), sortedElementsToIgnore.end());  // localElements should already be sorted
    stk::mesh::EntityVector::iterator resultIter = std::set_difference(localElements.begin(), localElements.end(), sortedElementsToIgnore.begin(), sortedElementsToIgnore.end(), localElementsToConsider.begin());
    localElementsToConsider.resize(resultIter - localElementsToConsider.begin());

    add_local_elements_to_connected_list(localElementsToConsider, sideNodes, communicatedElementDataVector);

    for (stk::mesh::Entity localElem: localElementsToConsider) {
        stk::topology localElemTopology = m_bulk_data.bucket(localElem).topology();

        const stk::mesh::Entity* localElemNodes = m_bulk_data.begin_nodes(localElem);
        unsigned localElemNumSides = localElemTopology.num_sides();
        for(unsigned side_index=0; side_index<localElemNumSides; ++side_index)
        {
            unsigned num_nodes_this_side = localElemTopology.side_topology(side_index).num_nodes();
            if (num_nodes_this_side == num_side_nodes)
            {
                stk::mesh::EntityVector localElemSideNodes(num_nodes_this_side);
                localElemTopology.side_nodes(localElemNodes, side_index, localElemSideNodes.begin());

                impl::ConnectedElementDataVector filteredCommunicatedElementData = communicatedElementDataVector;  // Modify copy for this local element
                impl::filter_for_candidate_elements_to_connect(m_bulk_data, localElem, side_index, filteredCommunicatedElementData);

                for (const impl::ConnectedElementData & elemData : filteredCommunicatedElementData) {
                    const bool isRemoteElement = (elemData.m_procId != m_bulk_data.parallel_rank());
                    if (isRemoteElement)
                    {
                        if (elemData.m_elementTopology.is_shell()) {
                            localElementsConnectedToRemoteShell.insert(get_local_element_id(localElem));  // Store connection to error-check later
                        }
                        impl::LocalId local_elem_id = get_local_element_id(localElem);
                        impl::LocalId negSgnRemoteElemId = -1*elemData.m_elementId;

                        m_elem_graph[local_elem_id].push_back(negSgnRemoteElemId);
                        m_via_sides[local_elem_id].push_back(side_index);
                        ++m_num_edges;
                        stk::mesh::EntityId chosen_side_id = 0;
                        bool foundChosenSide = false;
                        if(m_bulk_data.identifier(localElem) < static_cast<stk::mesh::EntityId>(elemData.m_elementId))
                        {
                            const auto iterRange = elemSideComm.equal_range(impl::EntitySidePair(localElem, side_index));
                            for (impl::ElemSideToProcAndFaceId::const_iterator iter = iterRange.first; iter != iterRange.second; ++iter) {
                                if ((iter->second.proc != elemData.m_procId) && (iter->first.side_id != elemData.m_sideIndex)) {
                                    continue;
                                }
                                chosen_side_id = iter->second.side_id;
                                foundChosenSide = true;
                                break;
                            }
                        }
                        if (!foundChosenSide)
                        {
                            chosen_side_id = elemData.m_suggestedFaceId;
                        }

                        std::pair<bool,unsigned> result = localElemTopology.side_topology(side_index).equivalent(localElemSideNodes, elemData.m_sideNodes);
                        m_parallel_graph_info.insert(std::make_pair(std::make_pair(local_elem_id, elemData.m_elementId),
                                                                    impl::parallel_info(elemData.m_procId, elemData.m_sideIndex, result.second, chosen_side_id, elemData.m_isInPart)));

                        stk::mesh::impl::EntitySidePair key(localElem, side_index);

                        if ((nullptr != newlySharedEdges) && (elemSideComm.find(key) == elemSideComm.end()))
                        {
                            impl::SharedEdgeInfo sharedEdgeInfo;
                            sharedEdgeInfo.m_locaElementlId = m_bulk_data.identifier(localElem);
                            sharedEdgeInfo.m_remoteElementId = elemData.m_elementId;
                            sharedEdgeInfo.m_procId = elemData.m_procId;
                            sharedEdgeInfo.m_sideIndex = side_index;
                            sharedEdgeInfo.m_sharedNodes = localElemSideNodes;
                            sharedEdgeInfo.m_chosenSideId = chosen_side_id;
                            sharedEdgeInfo.m_isInPart = elemData.m_isInPart;
                            sharedEdgeInfo.m_remoteElementTopology = m_bulk_data.bucket(localElem).topology();

                            newlySharedEdges->push_back(sharedEdgeInfo);
                        }
                    }
                }
            }
        }
    }

    break_volume_element_connections_across_shells(localElementsConnectedToRemoteShell);
}

void ElemElemGraph::break_volume_element_connections_across_shells(const std::set<stk::mesh::EntityId> & localElementsConnectedToRemoteShell)
{
    // Fix the case where the serial graph connected two volume elements together before
    // it was known that there was a remote shell wedged between them (the "sandwich" conundrum).
    // Also, cover the case where the mesh is modified after the graph is created to
    // add a shell between existing volume elements.
    //
    if (localElementsConnectedToRemoteShell.size() > 1) {
        for (impl::LocalId localElemId: localElementsConnectedToRemoteShell) {
            std::vector<impl::LocalId>::iterator it = m_elem_graph[localElemId].begin();
            while (it != m_elem_graph[localElemId].end()) {
                const impl::LocalId connectedElemId = *it;
                if (localElementsConnectedToRemoteShell.find(connectedElemId) != localElementsConnectedToRemoteShell.end()) {
                    const int offset = (it - m_elem_graph[localElemId].begin());
                    it = m_elem_graph[localElemId].erase(it);
                    m_via_sides[localElemId].erase(m_via_sides[localElemId].begin() + offset);
                    --m_num_edges;
                }
                else {
                    ++it;
                }
            }
        }
    }
}

bool perform_element_death(stk::mesh::BulkData& bulkData, ElemElemGraph& elementGraph, const stk::mesh::EntityVector& killedElements, stk::mesh::Part& active,
        const stk::mesh::PartVector& boundary_mesh_parts)
{
    bool topology_modified = false;

    const std::vector<stk::mesh::EntityId> requestedIds = elementGraph.get_suggested_side_ids();
    size_t id_counter = 0;

    std::vector<stk::mesh::sharing_info> shared_modified;
    stk::mesh::EntityVector deletedEntities;
    stk::mesh::EntityVector facesWithNodesToBeMarkedInactive;

    std::vector<impl::graphEdgeProc> elements_to_comm = impl::get_elements_to_communicate(bulkData, killedElements, elementGraph);
    std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> > remote_edges;

    impl::communicate_killed_entities(bulkData, elements_to_comm, remote_edges);

    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

    stk::mesh::Part& sides_created_during_death = bulkData.mesh_meta_data().declare_part("sides_created_during_death",
            side_rank, true);

    if(bulkData.in_modifiable_state())
    {
        bulkData.modification_end();
    }
    bulkData.modification_begin();

    for(size_t re = 0; re < remote_edges.size(); ++re)
    {
        stk::mesh::EntityId local_id = remote_edges[re].first;
        stk::mesh::EntityId remote_id = remote_edges[re].second;

        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, local_id);
        bool create_side = true;
        if(!bulkData.bucket(element).member(active))
        {
            create_side = false;
        }

        stk::mesh::PartVector add_parts_for_shared_sides = boundary_mesh_parts;
        {
            const stk::mesh::PartVector & supersets = bulkData.bucket(element).supersets();
            for (size_t part_i=0 ; part_i<supersets.size() ; ++part_i)
            {
                if(!stk::mesh::is_auto_declared_part(*supersets[part_i]))
                {
                    add_parts_for_shared_sides.push_back(supersets[part_i]);
                }
            }
        }

        impl::parallel_info &parallel_edge_info = elementGraph.get_parallel_edge_info(element, remote_id);
        parallel_edge_info.m_in_part = false;

        // Process sides where element on another processor was deactivated
        topology_modified = impl::create_or_delete_shared_side(bulkData, parallel_edge_info, elementGraph, element, remote_id, create_side, add_parts_for_shared_sides,
                active, shared_modified, deletedEntities, facesWithNodesToBeMarkedInactive, id_counter, requestedIds[id_counter], sides_created_during_death) || topology_modified;
    }

    std::vector<impl::ElementSidePair> element_side_pairs;
    element_side_pairs.reserve(impl::get_element_side_multiplier() * killedElements.size());

    for(size_t k = 0; k < killedElements.size(); ++k)
    {
        stk::mesh::PartVector add_parts = boundary_mesh_parts;
        stk::mesh::Entity this_elem_entity = killedElements[k];

        for(size_t j = 0; j < elementGraph.get_num_connected_elems(this_elem_entity); ++j)
        {
            if(elementGraph.is_connected_elem_locally_owned(this_elem_entity, j))
            {
                // Process a side between locally owned elements
                stk::mesh::Entity other_element = elementGraph.get_connected_element(this_elem_entity, j);
                int side_id = elementGraph.get_side_id_to_connected_element(this_elem_entity, j);
                ThrowRequireMsg(side_id != -1, "Program error. Please contact sierra-help@sandia.gov for support.");

                bool is_other_element_alive = bulkData.bucket(other_element).member(active);
                stk::topology side_top = bulkData.bucket(this_elem_entity).topology().side_topology(side_id);
                if(is_other_element_alive)
                {
                    // create or delete a side with a particular id

                    stk::mesh::PartVector parts = add_parts;
                    parts.push_back(&bulkData.mesh_meta_data().get_topology_root_part(side_top));

                    std::string msg = "Program error. Please contact sierra-help@sandia.gov for support.";

                    stk::mesh::Entity side = stk::mesh::impl::get_side_for_element(bulkData, this_elem_entity, side_id);
                    topology_modified = true;
                    if(bulkData.is_valid(side))
                    {
                        if(bulkData.bucket(side).owned())
                        {
                            bulkData.change_entity_parts(side, parts, stk::mesh::PartVector());
                        }
                    }
                    else
                    {
                        stk::mesh::EntityId side_global_id = requestedIds[id_counter];
                        ++id_counter;
                        ThrowRequireMsg(!impl::is_id_already_in_use_locally(bulkData, side_rank, side_global_id), msg);
                        parts.push_back(&sides_created_during_death);
                        {
                            const stk::mesh::PartVector & supersets = bulkData.bucket(other_element).supersets();
                            for (size_t part_i=0 ; part_i<supersets.size() ; ++part_i)
                            {
                                if(!stk::mesh::is_auto_declared_part(*supersets[part_i]))
                                {
                                    parts.push_back(supersets[part_i]);
                                }
                            }
                        }

                        // switch elements
                        stk::mesh::Entity element_with_perm_0 = other_element;
                        stk::mesh::Entity element_with_perm_4 = this_elem_entity;

                        int side_id_needed = elementGraph.get_side_from_element1_to_locally_owned_element2(element_with_perm_0,
                                element_with_perm_4);
                        ThrowRequireMsg(side_id_needed >= 0, "ERROR: proc " << bulkData.parallel_rank() << " found side_id_needed=" << side_id_needed
                                        << " between elem " << bulkData.identifier(element_with_perm_0)<< " and " << bulkData.identifier(element_with_perm_4)
                                        << " in elem-elem-graph");

                        side = stk::mesh::declare_element_side(bulkData, side_global_id, element_with_perm_0, side_id_needed, parts);

                        const stk::mesh::Entity* side_nodes = bulkData.begin_nodes(side);
                        unsigned num_side_nodes = bulkData.num_nodes(side);
                        stk::mesh::EntityVector side_nodes_vec(side_nodes, side_nodes + num_side_nodes);

                        std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ord_and_perm =
                                stk::mesh::get_ordinal_and_permutation(bulkData, element_with_perm_4, side_rank, side_nodes_vec);

                        if(ord_and_perm.first == stk::mesh::INVALID_CONNECTIVITY_ORDINAL)
                        {
                            std::ostringstream yo;
                            yo << "Proc: " << bulkData.parallel_rank() << std::endl;
                            yo << "this element: " << bulkData.identifier(element_with_perm_0) << std::endl;
                            yo << "other element: " << bulkData.identifier(element_with_perm_4) << std::endl;
                            yo << "Nodes: ";

                            for(stk::mesh::Entity side_node : side_nodes_vec)
                            {
                                yo << bulkData.identifier(side_node) << " "; // nodes of elem 960 (other_elem: 1, this_elem: 401)
                            }

                            yo << std::endl;
                            std::cerr << yo.str();
                        }

                        ThrowRequireMsg(ord_and_perm.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL, "yikes!");
                        ThrowRequireMsg(ord_and_perm.second != stk::mesh::INVALID_PERMUTATION, "yikes!");

                        bulkData.declare_relation(element_with_perm_4, side, ord_and_perm.first, ord_and_perm.second);
                    }
                }
                else
                {
                    stk::mesh::Entity side = stk::mesh::impl::get_side_for_element(bulkData, this_elem_entity, side_id);
                    facesWithNodesToBeMarkedInactive.push_back(side);
                    if(bulkData.is_valid(side) && bulkData.bucket(side).member(sides_created_during_death))
                    {
                        deletedEntities.push_back(side);
                        topology_modified = true;
                    }
                    else if(bulkData.is_valid(side) && bulkData.bucket(side).owned())
                    {
                        bulkData.change_entity_parts(side, {}, {&active});
                    }
                }
            }
            else
            {
                // Process a side where we deactivated the element on this graph edge.
                // If the element on the remote processor was also deactivated, we may have already processed this side.
                // TODO:  Determine if true and optimize

                // create or delete a side with a particular id

                stk::mesh::EntityId remote_id = elementGraph.get_entity_id_of_remote_element(this_elem_entity, j);

                impl::parallel_info &parallel_edge_info = elementGraph.get_parallel_edge_info(this_elem_entity, remote_id);
                bool other_element_active = parallel_edge_info.m_in_part;
                bool create_side = false;
                if(other_element_active)
                {
                    create_side = true;
                }

                {
                    add_parts = boundary_mesh_parts;
                    const stk::mesh::PartVector & supersets = bulkData.bucket(this_elem_entity).supersets();
                    for (size_t part_i=0 ; part_i<supersets.size() ; ++part_i)
                    {
                        if(!stk::mesh::is_auto_declared_part(*supersets[part_i]))
                        {
                            add_parts.push_back(supersets[part_i]);
                        }
                    }
                }

                topology_modified = impl::create_or_delete_shared_side(bulkData, parallel_edge_info, elementGraph, this_elem_entity, remote_id, create_side, add_parts,
                        active, shared_modified, deletedEntities, facesWithNodesToBeMarkedInactive,  id_counter, requestedIds[id_counter], sides_created_during_death) || topology_modified;
            }
        }
    }

    ThrowRequireMsg(id_counter==0 || id_counter<requestedIds.size(), "Program error. Please contact sierra-help@sandia.gov for support.");
    elementGraph.set_num_side_ids_used(id_counter);
    stk::mesh::impl::delete_entities_and_upward_relations(bulkData, deletedEntities);
    bulkData.modification_end_for_face_creation_and_deletion(shared_modified, deletedEntities, elementGraph, killedElements, active);
    return topology_modified;
}

void ElemElemGraph::add_local_elements_to_connected_list(const stk::mesh::EntityVector & connected_elements,
                                                         const stk::mesh::EntityVector & sideNodes,
                                                         impl::ConnectedElementDataVector & connectedElementDataVector) const
{
    for (const stk::mesh::Entity & connectedElem: connected_elements) {
        impl::ConnectedElementData elemData;
        stk::mesh::EntityVector connectedSideNodes;
        stk::mesh::OrdinalAndPermutation connectedOrdAndPerm = stk::mesh::get_ordinal_and_permutation(m_bulk_data, connectedElem, m_bulk_data.mesh_meta_data().side_rank(), sideNodes);
        const stk::mesh::Bucket & connectedBucket = m_bulk_data.bucket(connectedElem);
        const stk::mesh::Entity* connectedElemNodes = m_bulk_data.begin_nodes(connectedElem);

        elemData.m_procId = m_bulk_data.parallel_rank();
        elemData.m_elementId = get_local_element_id(connectedElem, false);
        if (elemData.m_elementId == INVALID_LOCAL_ID)
        {
            continue;
        }
        elemData.m_elementTopology = connectedBucket.topology();
        elemData.m_sideIndex = connectedOrdAndPerm.first;
        elemData.m_suggestedFaceId = 0;
        elemData.m_sideNodes.resize(sideNodes.size());
        elemData.m_elementTopology.side_nodes(connectedElemNodes, elemData.m_sideIndex, elemData.m_sideNodes.begin());
        connectedElementDataVector.push_back(elemData);

        if (elemData.m_elementTopology.is_shell()) {
            // Also add the back-side face if this is a shell
            elemData.m_sideIndex = (elemData.m_sideIndex == 0) ? 1 : 0;
            elemData.m_elementTopology.side_nodes(connectedElemNodes, elemData.m_sideIndex, elemData.m_sideNodes.begin());
            connectedElementDataVector.push_back(elemData);
        }
    }
}


void ElemElemGraph::pack_remote_connected_element(impl::LocalId elem_local_id, stk::mesh::EntityId connected_global_id,
                                                  stk::CommBuffer &buff, std::vector<moved_parallel_graph_info> &moved_graph_info_vector,
                                                  int destination_proc, int phase)
{
    std::pair<impl::LocalId, stk::mesh::EntityId> key(elem_local_id, connected_global_id);
    auto iter = m_parallel_graph_info.find(key);
    ThrowRequire(iter != m_parallel_graph_info.end());
    impl::parallel_info &p_info = iter->second;
    buff.pack<int>(p_info.m_other_proc);
    buff.pack<int>(p_info.m_other_side_ord);
    buff.pack<int>(p_info.m_permutation);
    buff.pack<bool>(p_info.m_in_part);
    buff.pack<stk::mesh::EntityId>(p_info.m_chosen_side_id);

    if (phase == 0 && p_info.m_other_proc != destination_proc)
    {
        stk::mesh::Entity elem = m_local_id_to_element_entity[elem_local_id];
        stk::mesh::EntityId elem_global_id = m_bulk_data.identifier(elem);
        moved_graph_info_vector.push_back(moved_parallel_graph_info(p_info.m_other_proc, connected_global_id, elem_global_id, destination_proc));
    }

    if (phase == 1)
    {
        m_parallel_graph_info.erase(iter);
    }
}

void ElemElemGraph::pack_local_connected_element(impl::LocalId local_id, int side_id, stk::CommBuffer &buff,
                                                 stk::mesh::EntityId suggested_face_id,
                                                 stk::mesh::Part *active_part)
{
    stk::mesh::Entity connected_element = m_local_id_to_element_entity[local_id];
    buff.pack<int>(m_bulk_data.parallel_rank());
    buff.pack<int>(side_id);
    stk::topology elem_topology = m_bulk_data.bucket(connected_element).topology();
    stk::topology side_topology = elem_topology.side_topology(side_id);
    std::vector<stk::mesh::Entity> side_nodes(side_topology.num_nodes());
    const stk::mesh::Entity *elem_nodes = m_bulk_data.begin_nodes(connected_element);
    stk::mesh::EntityRank side_rank = m_bulk_data.mesh_meta_data().side_rank();
    elem_topology.side_nodes(elem_nodes, side_id, side_nodes.begin());

    stk::mesh::OrdinalAndPermutation ordperm = get_ordinal_and_permutation(m_bulk_data, connected_element, side_rank, side_nodes);
    ThrowRequireMsg(ordperm.second != stk::mesh::INVALID_PERMUTATION, "Invalid permutation for connected_element");

    buff.pack<int>(ordperm.second);
    bool in_part = true;
    if (active_part != NULL)
    {
        in_part = m_bulk_data.bucket(connected_element).member(*active_part);
    }
    buff.pack<bool>(in_part);
    buff.pack<stk::mesh::EntityId>(suggested_face_id);
    buff.pack<size_t>(side_nodes.size());
    for(size_t i=0; i<side_nodes.size(); ++i)
    {
        buff.pack<stk::mesh::EntityId>(m_bulk_data.identifier(side_nodes[i]));
    }
}

void ElemElemGraph::unpack_and_store_connected_element(stk::CommBuffer &buf, impl::LocalId recvd_elem_local_id,
                                                       stk::mesh::EntityId recvd_elem_global_id)
{
    bool local_connection_on_source_proc;
    buf.unpack<bool>(local_connection_on_source_proc);
    stk::mesh::EntityId connected_elem_global_id;
    buf.unpack<stk::mesh::EntityId>(connected_elem_global_id);
    int side_from_recvd_elem_to_connected_elem;
    buf.unpack<int>(side_from_recvd_elem_to_connected_elem);
    stk::mesh::Entity connected_elem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, connected_elem_global_id);
    int other_proc, side_from_connected_elem_to_recvd_elem, other_permutation;
    bool in_part;
    stk::mesh::EntityId chosen_face_id;
    buf.unpack<int>(other_proc);
    buf.unpack<int>(side_from_connected_elem_to_recvd_elem);
    buf.unpack<int>(other_permutation);
    buf.unpack<bool>(in_part);
    buf.unpack<stk::mesh::EntityId>(chosen_face_id);
    impl::parallel_info p_info(other_proc, side_from_connected_elem_to_recvd_elem, other_permutation, chosen_face_id, in_part);
    size_t num_side_nodes = 0;
    std::vector<stk::mesh::EntityId> side_node_ids;
    stk::mesh::EntityVector side_nodes;
    if (local_connection_on_source_proc)
    {
        buf.unpack<size_t>(num_side_nodes);
        side_node_ids.resize(num_side_nodes);
        side_nodes.resize(num_side_nodes);
        for(size_t i=0; i<num_side_nodes; ++i)
        {
            buf.unpack<stk::mesh::EntityId>(side_node_ids[i]);
            side_nodes[i] = m_bulk_data.get_entity(stk::topology::NODE_RANK, side_node_ids[i]);
        }
    }

    if (p_info.m_other_proc == m_bulk_data.parallel_rank())
    {
        impl::LocalId connected_elem_local_id = get_local_element_id(connected_elem);
        m_elem_graph[recvd_elem_local_id].push_back(connected_elem_local_id);

        std::vector<impl::LocalId>& other_connected_elements = m_elem_graph[connected_elem_local_id];
        impl::LocalId former_remote_id = -recvd_elem_global_id;
        auto iter_found = std::find(other_connected_elements.begin(), other_connected_elements.end(), former_remote_id);
        ThrowRequireMsg(iter_found != other_connected_elements.end(), "Failed to find formerly-remote connected id in elem-elem-graph.");
        size_t index = iter_found - other_connected_elements.begin();
        other_connected_elements[index] = recvd_elem_local_id;

        std::pair<impl::LocalId, stk::mesh::EntityId> key(connected_elem_local_id, recvd_elem_global_id);
        auto iter = m_parallel_graph_info.find(key);
        if (iter != m_parallel_graph_info.end())
        {
            m_parallel_graph_info.erase(iter);
        }

        m_num_parallel_edges--;
    }
    else
    {
        stk::mesh::EntityRank side_rank = m_bulk_data.mesh_meta_data().side_rank();
        stk::mesh::Entity recvd_elem = m_local_id_to_element_entity[recvd_elem_local_id];
        stk::mesh::OrdinalAndPermutation ordperm = get_ordinal_and_permutation(m_bulk_data, recvd_elem, side_rank, side_nodes);
        p_info.m_permutation = ordperm.second;

        m_elem_graph[recvd_elem_local_id].push_back(-connected_elem_global_id);
        std::pair<impl::LocalId, stk::mesh::EntityId> recvd_elem_key(recvd_elem_local_id, connected_elem_global_id);
        m_parallel_graph_info.insert(std::make_pair(recvd_elem_key, p_info));

        m_num_parallel_edges++;
    }
    m_num_edges++;
    m_via_sides[recvd_elem_local_id].push_back(side_from_recvd_elem_to_connected_elem);
}

void ElemElemGraph::communicate_moved_graph_info(std::vector <moved_parallel_graph_info> &moved_graph_info_vector)
{
    stk::CommSparse comm2(m_bulk_data.parallel());
    for(int phase=0; phase <2; ++phase) {
        for (size_t i=0; i<moved_graph_info_vector.size(); ++i)
        {
            moved_parallel_graph_info &info = moved_graph_info_vector[i];
            stk::CommBuffer &buf = comm2.send_buffer(info.proc_to_tell);

            buf.pack<stk::mesh::EntityId>(info.elem_id);
            buf.pack<stk::mesh::EntityId>(info.moved_elem_id);
            buf.pack<int>(info.destination_proc);
        }
        if (phase == 0)
        {
            comm2.allocate_buffers();
        }
        else
        {
            comm2.communicate();
        }
    }

    for(int p = 0; p < m_bulk_data.parallel_size(); ++p)
    {
        stk::CommBuffer & buf = comm2.recv_buffer(p);
        while(buf.remaining())
        {
            stk::mesh::EntityId elem_id, moved_elem_id;
            int destination_proc;
            buf.unpack<stk::mesh::EntityId>(elem_id);
            buf.unpack<stk::mesh::EntityId>(moved_elem_id);
            buf.unpack<int>(destination_proc);

            stk::mesh::Entity local_elem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, elem_id);
            impl::LocalId elem_local_id = get_local_element_id(local_elem);
            std::pair<impl::LocalId, stk::mesh::EntityId> key(elem_local_id, moved_elem_id);
            auto iter = m_parallel_graph_info.find(key);
            if (iter != m_parallel_graph_info.end())
            {
                iter->second.m_other_proc = destination_proc;
            }
        }
    }
}

impl::LocalId ElemElemGraph::create_new_local_id(stk::mesh::Entity new_elem)
{
    if (m_entity_to_local_id.size() > new_elem.local_offset() && m_entity_to_local_id[new_elem.local_offset()] != INVALID_LOCAL_ID)
    {
        return m_entity_to_local_id[new_elem.local_offset()];
    }

    impl::LocalId new_local_id = m_elem_graph.size();
    if (m_deleted_element_local_id_pool.size() > 0)
    {
        new_local_id = m_deleted_element_local_id_pool.back();
        m_local_id_in_pool[new_local_id] = false;
        m_deleted_element_local_id_pool.pop_back();
    }
    else
    {
        if (m_local_id_to_element_entity.size() <= static_cast<size_t> (new_local_id))
        {
            m_local_id_to_element_entity.resize(new_local_id+1);
        }
        if (m_entity_to_local_id.size() <= new_elem.local_offset())
        {
            m_entity_to_local_id.resize(new_elem.local_offset()+1);
        }
        m_elem_graph.push_back(std::vector<impl::LocalId>());
        m_via_sides.push_back(std::vector<int>());
    }
    m_local_id_to_element_entity[new_local_id] = new_elem;
    m_entity_to_local_id[new_elem.local_offset()] = new_local_id;

    return new_local_id;
}

void ElemElemGraph::change_entity_owner(const stk::mesh::EntityProcVec &elem_proc_pairs_to_move, impl::ParallelGraphInfo &new_parallel_graph_entries, stk::mesh::Part *active_part)
{
    std::vector <moved_parallel_graph_info> moved_graph_info_vector;

    stk::CommSparse comm(m_bulk_data.parallel());
    for(int phase=0; phase <2; ++phase) {
        for (size_t i=0; i<elem_proc_pairs_to_move.size(); i++)
        {
            stk::mesh::Entity elem_to_send = elem_proc_pairs_to_move[i].first;
            int destination_proc = elem_proc_pairs_to_move[i].second;
            stk::mesh::EntityId elem_global_id = m_bulk_data.identifier(elem_to_send);
            stk::CommBuffer &buff = comm.send_buffer(destination_proc);

            buff.pack<stk::mesh::EntityId>(elem_global_id);
            impl::LocalId elem_local_id = m_entity_to_local_id[elem_to_send.local_offset()];
            std::vector <impl::LocalId> &connected_elements = m_elem_graph[elem_local_id];
            size_t num_connected_elements = connected_elements.size();
            buff.pack<size_t>(num_connected_elements);
            for (size_t k=0; k<num_connected_elements; k++)
            {
                impl::LocalId local_id = connected_elements[k];
                stk::mesh::EntityId connected_global_id;
                bool local_connection = local_id >= 0;

                if (local_connection)
                {
                    stk::mesh::Entity connected_element = m_local_id_to_element_entity[local_id];
                    connected_global_id = m_bulk_data.identifier(connected_element);
                }
                else
                {
                    connected_global_id = -local_id;
                }

                buff.pack<bool>(local_connection);
                buff.pack<stk::mesh::EntityId>(connected_global_id);
                buff.pack<int>(m_via_sides[elem_local_id][k]);

                if (local_connection)
                {
                    stk::mesh::Entity connected_element = m_local_id_to_element_entity[local_id];
                    int side_id = get_side_from_element1_to_locally_owned_element2(connected_element, elem_to_send);
                    std::pair<impl::LocalId,stk::mesh::EntityId> key(local_id,elem_global_id);
                    auto iter_p_info = new_parallel_graph_entries.find(key);
                    ThrowRequireMsg(iter_p_info != new_parallel_graph_entries.end(), "ERROR, failed to find entry in new_parallel_graph_entries"
                                    << " for elem_to_send="<<elem_global_id<<", local-connected="<<connected_global_id);
                    stk::mesh::EntityId suggested_face_id = iter_p_info->second.m_chosen_side_id;
                    pack_local_connected_element(local_id, side_id, buff, suggested_face_id, active_part);
                    if (phase == 1) {
                        std::vector<impl::LocalId>& other_connected_elements = m_elem_graph[local_id];
                        auto iter = std::find(other_connected_elements.begin(), other_connected_elements.end(), elem_local_id);
                        ThrowRequireMsg(iter != other_connected_elements.end(), "Failed to find connected element");
                        size_t index = iter - other_connected_elements.begin();
                        other_connected_elements[index] = -elem_global_id;
                    }
                }
                else
                {
                    pack_remote_connected_element(elem_local_id, connected_global_id, buff, moved_graph_info_vector, destination_proc, phase);
                }
            }

            if (phase == 1)
            {
                for (size_t j=0; j<connected_elements.size(); j++)
                {
                    if (connected_elements[j]>0)
                    {
                        std::vector <impl::LocalId> &elements_connected = m_elem_graph[connected_elements[j]];
                        auto iter = std::find(elements_connected.begin(), elements_connected.end(), elem_local_id);
                        ThrowRequireMsg(iter != elements_connected.end(), "Program error. Please contact sierra-help@sandia.gov for support.");
                        int index = iter - elements_connected.begin();
                        elements_connected[index] = -elem_global_id;
                    }
                }

                for(size_t k=0; k<m_elem_graph[elem_local_id].size(); ++k)
                {
                    if(m_elem_graph[elem_local_id][k] < 0)
                        m_num_parallel_edges--;
                    else
                        m_num_parallel_edges++;
                }

                m_num_edges -= m_elem_graph[elem_local_id].size();

                m_elem_graph[elem_local_id].clear();
                m_via_sides[elem_local_id].clear();
                m_deleted_element_local_id_pool.push_back(elem_local_id);
                m_local_id_in_pool[elem_local_id] = true;
            }
        }

        if (phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    for(int p = 0; p < m_bulk_data.parallel_size(); ++p)
    {
        stk::CommBuffer & buf = comm.recv_buffer(p);
        while(buf.remaining())
        {
            stk::mesh::EntityId recvd_elem_global_id;
            buf.unpack<stk::mesh::EntityId>(recvd_elem_global_id);
            stk::mesh::Entity recvd_elem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, recvd_elem_global_id);
            impl::LocalId recvd_elem_local_id = create_new_local_id(recvd_elem);

            size_t num_connected_elements;
            buf.unpack<size_t>(num_connected_elements);
            for (size_t k=0; k<num_connected_elements; k++)
            {
                unpack_and_store_connected_element(buf, recvd_elem_local_id, recvd_elem_global_id);
            }
        }
    }
    auto iter_begin = new_parallel_graph_entries.begin();
    auto iter_end   = new_parallel_graph_entries.end();
    for (auto iter=iter_begin; iter!=iter_end; ++iter)
    {
        m_parallel_graph_info.insert(std::make_pair(iter->first, iter->second));
    }

    communicate_moved_graph_info(moved_graph_info_vector);
}

impl::LocalId ElemElemGraph::get_new_local_element_id_from_pool()
{
    impl::LocalId new_local_id;
    if (!m_deleted_element_local_id_pool.empty())
    {
        new_local_id = m_deleted_element_local_id_pool.back();
        m_deleted_element_local_id_pool.pop_back();
        m_local_id_in_pool[new_local_id] = false;
    }
    else
    {
        new_local_id = m_elem_graph.size();
        std::vector<impl::LocalId> new_element_connectivity;
        std::vector<int> new_element_via_sides;
        m_elem_graph.push_back(new_element_connectivity);
        m_via_sides.push_back(new_element_via_sides);
        m_local_id_in_pool.push_back(false);
        m_local_id_to_element_entity.push_back(Entity());
        m_element_topologies.push_back(stk::topology::INVALID_TOPOLOGY);
    }
    return new_local_id;
}

bool ElemElemGraph::is_valid_graph_element(stk::mesh::Entity local_element)
{
    bool value = false;
    if (m_bulk_data.is_valid(local_element))
    {
        impl::LocalId max_elem_id = static_cast<impl::LocalId>(m_elem_graph.size());
        impl::LocalId elem_id = get_local_element_id(local_element, false);
        value = elem_id >= 0 && elem_id < max_elem_id && !m_local_id_in_pool[elem_id];
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
        stk::mesh::Entity elem_to_delete = m_local_id_to_element_entity[elem_to_delete_local_id];
        stk::mesh::EntityId elem_to_delete_global_id = m_bulk_data.identifier(elem_to_delete);
        stk::mesh::EntityId connected_elem_global_id = local_elem_and_remote_connected_elem[i].m_remoteElement;

        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(elem_to_delete_global_id);
        comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(connected_elem_global_id);
    }
}

void ElemElemGraph::pack_shell_connectivity(stk::CommSparse & comm,
                                            const std::vector<impl::ShellConnectivityData> & shellConnectivityList)
{
    for(const impl::ShellConnectivityData & shellConnectivityData: shellConnectivityList)
    {
        if (shellConnectivityData.m_farElementIsRemote) {
            stk::mesh::Entity deletedElement = m_bulk_data.get_entity(stk::topology::ELEM_RANK, shellConnectivityData.m_shellElementId);
            impl::LocalId deletedElementId = get_local_element_id(deletedElement);
            auto iter = m_parallel_graph_info.find( std::make_pair(deletedElementId, shellConnectivityData.m_farElementId) );
            ThrowRequire(iter != m_parallel_graph_info.end());
            int remote_proc = iter->second.m_other_proc;

            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(shellConnectivityData.m_nearElementId);
            comm.send_buffer(remote_proc).pack<int>(shellConnectivityData.m_nearElementProc);
            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(shellConnectivityData.m_shellElementId);
            comm.send_buffer(remote_proc).pack<stk::mesh::EntityId>(shellConnectivityData.m_farElementId);
        }
    }
}

void ElemElemGraph::collect_local_shell_connectivity_data(const stk::mesh::EntityVector& elements_to_delete,
                                                          std::vector<impl::ShellConnectivityData>& shellConnectivityList)
{
    for (const stk::mesh::Entity &elem_to_delete : elements_to_delete) {
        ThrowRequireMsg(m_bulk_data.is_valid(elem_to_delete), "Do not delete entities before removing from ElemElemGraph. Contact sierra-help@sandia.gov for support.");
        ThrowRequireMsg(is_valid_graph_element(elem_to_delete), "Program error. Not valid graph element. Contact sierra-help@sandia.gov for support.");

        if (!m_bulk_data.bucket(elem_to_delete).owned()) {
            continue;
        }

        stk::topology elem_to_delete_topology = m_bulk_data.bucket(elem_to_delete).topology();
        if (elem_to_delete_topology.is_shell()) {
            impl::LocalId shell_to_delete_id = get_local_element_id( elem_to_delete);
            size_t num_connected_elems = get_num_connected_elems(elem_to_delete);
            for (int near_elem_index = num_connected_elems - 1; near_elem_index >= 0; --near_elem_index) {
                impl::ShellConnectivityData shellConnectivityData;

                impl::LocalId near_elem_id = m_elem_graph[shell_to_delete_id][near_elem_index];
                if (near_elem_id >= 0) {
                    std::vector<impl::LocalId>::iterator pos_of_shell_in_near_elem = std::find(m_elem_graph[near_elem_id].begin(), m_elem_graph[near_elem_id].end(), shell_to_delete_id);
                    ThrowRequire(pos_of_shell_in_near_elem != m_elem_graph[near_elem_id].end());
                    int index_of_shell_in_near_elem = pos_of_shell_in_near_elem - m_elem_graph[near_elem_id].begin();
                    stk::mesh::Entity nearElem = m_local_id_to_element_entity[near_elem_id];
                    shellConnectivityData.m_nearElementSide = m_via_sides[near_elem_id][index_of_shell_in_near_elem];
                    shellConnectivityData.m_nearElementId = m_bulk_data.identifier(nearElem);
                    shellConnectivityData.m_nearElementProc = m_bulk_data.parallel_rank();
                }
                else {
                    shellConnectivityData.m_nearElementSide = -1; // Near element is remote, so let it figure out the side
                    shellConnectivityData.m_nearElementId = -near_elem_id;

                    auto iter = m_parallel_graph_info.find( std::make_pair(shell_to_delete_id, shellConnectivityData.m_nearElementId));
                    ThrowRequire(iter != m_parallel_graph_info.end());
                    shellConnectivityData.m_nearElementProc = iter->second.m_other_proc;
                }

                shellConnectivityData.m_shellElementId = m_bulk_data.identifier(elem_to_delete);

                std::vector<impl::LocalId> & graphElemIds = m_elem_graph[shell_to_delete_id];
                ThrowAssertMsg(graphElemIds.size() <= 2, "Coincident shells detected.  Please call (505)844-3041 for support.");
                for (impl::LocalId graphElemId : graphElemIds) {
                    if (graphElemId != near_elem_id) {
                        if (graphElemId < 0) {
                            // Remote Element
                            shellConnectivityData.m_farElementId = -graphElemId;
                            shellConnectivityData.m_farElementIsRemote = true;

                            auto iter = m_parallel_graph_info.find( std::make_pair(shell_to_delete_id, shellConnectivityData.m_farElementId));
                            ThrowRequire(iter != m_parallel_graph_info.end());
                            shellConnectivityData.m_farElementProc = iter->second.m_other_proc;
                        }
                        else {
                            // Local Element
                            stk::mesh::Entity remoteElement = m_local_id_to_element_entity[graphElemId];
                            shellConnectivityData.m_farElementId = m_bulk_data.identifier(remoteElement);
                            shellConnectivityData.m_farElementIsRemote = false;
                            shellConnectivityData.m_farElementProc = m_bulk_data.parallel_rank();
                        }
                        // Only add to list if there *is* something connected on the back-side of the shell
                        shellConnectivityList.push_back(shellConnectivityData);
                        break;
                    }
                }
            }
        }
    }
}

void ElemElemGraph::communicate_shell_connectivity(std::vector<impl::ShellConnectivityData>& shellConnectivityList) {

    stk::CommSparse shellComm(m_bulk_data.parallel());
    pack_shell_connectivity(shellComm, shellConnectivityList);
    shellComm.allocate_buffers();
    pack_shell_connectivity(shellComm, shellConnectivityList);
    shellComm.communicate();

    for (int proc = 0; proc < m_bulk_data.parallel_size(); ++proc) {
        while (shellComm.recv_buffer(proc).remaining()) {
            impl::ShellConnectivityData shellConnectivityData;
            shellComm.recv_buffer(proc).unpack<stk::mesh::EntityId>(shellConnectivityData.m_farElementId); // Flip remote and local
            shellComm.recv_buffer(proc).unpack<int>(shellConnectivityData.m_farElementProc);
            shellComm.recv_buffer(proc).unpack<stk::mesh::EntityId>(shellConnectivityData.m_shellElementId);
            shellComm.recv_buffer(proc).unpack<stk::mesh::EntityId>(shellConnectivityData.m_nearElementId); // Flip remote and local
            shellConnectivityData.m_nearElementProc = m_bulk_data.parallel_rank();

            stk::mesh::Entity nearElement = m_bulk_data.get_entity(stk::topology::ELEM_RANK, shellConnectivityData.m_nearElementId);
            impl::LocalId nearElementId = get_local_element_id(nearElement);

            std::vector<impl::LocalId>& graphElemIds = m_elem_graph[nearElementId];
            for (size_t i = 0; i < graphElemIds.size(); ++i) {
                impl::LocalId shellElementId = -shellConnectivityData.m_shellElementId;
                if (graphElemIds[i] == shellElementId) {
                    shellConnectivityData.m_nearElementSide = m_via_sides[nearElementId][i];
                    break;
                }
            }
            shellConnectivityList.push_back(shellConnectivityData);
        }
    }
}

void ElemElemGraph::delete_local_connections_and_collect_remote(const stk::mesh::EntityVector& elements_to_delete,
                                                                std::vector<impl::DeletedElementData>& local_elem_and_remote_connected_elem)
{
    for (const stk::mesh::Entity &elem_to_delete : elements_to_delete) {
        impl::LocalId elem_to_delete_id = get_local_element_id(elem_to_delete);

        if (!m_bulk_data.bucket(elem_to_delete).owned()) {
            continue;
        }

        size_t num_connected_elems = get_num_connected_elems(elem_to_delete);
        for (int conn_elem_index = num_connected_elems - 1; conn_elem_index >= 0; --conn_elem_index) {
            impl::LocalId connected_elem_id = m_elem_graph[elem_to_delete_id][conn_elem_index];
            bool local_connection = connected_elem_id >= 0;
            if (local_connection) {
                m_via_sides[elem_to_delete_id].erase(m_via_sides[elem_to_delete_id].begin() + conn_elem_index);
                --m_num_edges;

                std::vector<impl::LocalId>::iterator pos_of_elem1_in_elem2 = std::find(m_elem_graph[connected_elem_id].begin(), m_elem_graph[connected_elem_id].end(), elem_to_delete_id);
                ThrowRequire(pos_of_elem1_in_elem2 != m_elem_graph[connected_elem_id].end());
                int index_of_elem1_in_elem2 = pos_of_elem1_in_elem2 - m_elem_graph[connected_elem_id].begin();

                m_via_sides[connected_elem_id].erase(m_via_sides[connected_elem_id].begin() + index_of_elem1_in_elem2);
                --m_num_edges;

            } else {
                stk::mesh::EntityId connected_global_id = -connected_elem_id;
                auto iter = m_parallel_graph_info.find(std::make_pair(elem_to_delete_id, connected_global_id));
                ThrowRequire(iter != m_parallel_graph_info.end());
                int remote_proc = iter->second.m_other_proc;

                impl::DeletedElementData deletedElementData;
                deletedElementData.m_deletedElement = elem_to_delete_id;
                deletedElementData.m_remoteElement = connected_global_id;
                deletedElementData.m_remoteProc = remote_proc;
                local_elem_and_remote_connected_elem.push_back(deletedElementData);

                m_via_sides[elem_to_delete_id].erase(m_via_sides[elem_to_delete_id].begin() + conn_elem_index);
                --m_num_edges;
                m_parallel_graph_info.erase(std::make_pair(elem_to_delete_id, connected_global_id));
            }
        }

        ThrowRequire(m_via_sides[elem_to_delete_id].empty());
        for (size_t id = 0; id < m_elem_graph.size(); ++id) {
            if (id != static_cast<size_t>(elem_to_delete_id) && !m_local_id_in_pool[id]) {
                std::vector<impl::LocalId>::iterator pos_of_deleted_elem_in_current_id = std::find(m_elem_graph[id].begin(), m_elem_graph[id].end(), elem_to_delete_id);
                if (pos_of_deleted_elem_in_current_id != m_elem_graph[id].end()) {
                    int index_of_deleted_elem_in_current_id = pos_of_deleted_elem_in_current_id - m_elem_graph[id].begin();
                    m_elem_graph[id].erase(m_elem_graph[id].begin() + index_of_deleted_elem_in_current_id);
                }
            }
        }
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
            remote_edges.push_back(std::make_pair(deleted_elem_global_id, connected_elem_global_id));
        }
    }
}

void ElemElemGraph::clear_deleted_element_connections(const stk::mesh::EntityVector& elements_to_delete)
{
    for (auto & elem : elements_to_delete) {
        impl::LocalId elem_to_delete_id = m_entity_to_local_id[elem.local_offset()];
        m_elem_graph[elem_to_delete_id].clear();
        ThrowAssertMsg(m_via_sides[elem_to_delete_id].empty(), "Unable to delete element from graph.  Contact sierra-help@sandia.gov for help.");
        m_deleted_element_local_id_pool.push_back(elem_to_delete_id);
        m_local_id_in_pool[elem_to_delete_id] = true;
        m_element_topologies[elem_to_delete_id] = stk::topology::INVALID_TOPOLOGY;
        m_entity_to_local_id[elem.local_offset()] = INVALID_LOCAL_ID;
        m_local_id_to_element_entity[elem_to_delete_id] = stk::mesh::Entity::InvalidEntity;
    }
}

void ElemElemGraph::delete_remote_connections(const std::vector<std::pair<stk::mesh::EntityId, stk::mesh::EntityId> >& remote_edges)
{
    for (auto & edge : remote_edges) {
        stk::mesh::EntityId deleted_elem_global_id = edge.first;
        stk::mesh::EntityId connected_elem_global_id = edge.second;
        stk::mesh::Entity connected_elem = m_bulk_data.get_entity(stk::topology::ELEM_RANK, connected_elem_global_id);
        impl::LocalId connected_elem_id = m_entity_to_local_id[connected_elem.local_offset()];
        if (!is_valid_graph_element(connected_elem)) {
            continue;
        }

        size_t num_conn_elem = get_num_connected_elems(connected_elem);
        bool found_deleted_elem = false;
        for (size_t conn_elem_index = 0; conn_elem_index < num_conn_elem; ++conn_elem_index) {
            if (m_elem_graph[connected_elem_id][conn_elem_index] == static_cast<int64_t>(-deleted_elem_global_id)) {
                m_elem_graph[connected_elem_id].erase(m_elem_graph[connected_elem_id].begin() + conn_elem_index);
                m_via_sides[connected_elem_id].erase(m_via_sides[connected_elem_id].begin() + conn_elem_index);
                --m_num_edges;
                found_deleted_elem = true;
                m_parallel_graph_info.erase(std::make_pair(connected_elem_id, deleted_elem_global_id));
                break;
            }
        }
        ThrowRequireMsg(found_deleted_elem, "Error. Contact sierra-help@sandia.gov for support.");
    }
}

void ElemElemGraph::reconnect_volume_elements_across_deleted_shells(std::vector<impl::ShellConnectivityData> & shellConnectivityList,
                                                                    const stk::mesh::EntityVector& elements_to_delete)
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
        if (data.m_nearElementProc == m_bulk_data.parallel_rank() && data.m_farElementProc == m_bulk_data.parallel_rank()) {
            impl::LocalId localElementId = get_local_element_id(localElement);
            stk::mesh::Entity remoteElement = m_bulk_data.get_entity(stk::topology::ELEM_RANK, data.m_farElementId);
            impl::LocalId remoteElementId;
            remoteElementId = get_local_element_id(remoteElement);
            m_elem_graph[localElementId].push_back(remoteElementId);
            m_via_sides[localElementId].push_back(data.m_nearElementSide);
            ++m_num_edges;
        }
        else {
            shellNeighborsToReconnect.insert(
                    std::pair<impl::EntitySidePair, impl::ProcFaceIdPair>(impl::EntitySidePair(localElement, data.m_nearElementSide),
                            impl::ProcFaceIdPair(data.m_farElementProc, 0)));
        }
    }

    this->generate_additional_ids_collective(m_num_edges + shellConnectivityList.size());
    fill_parallel_graph(shellNeighborsToReconnect, elements_to_delete);

    update_number_of_parallel_edges();
}

void ElemElemGraph::delete_elements_from_graph(const stk::mesh::EntityVector &elements_to_delete)
{
    std::vector<impl::ShellConnectivityData> shellConnectivityList;
    collect_local_shell_connectivity_data(elements_to_delete, shellConnectivityList);

    communicate_shell_connectivity(shellConnectivityList);

    std::vector<impl::DeletedElementData> local_elem_and_remote_connected_elem;
    delete_local_connections_and_collect_remote(elements_to_delete, local_elem_and_remote_connected_elem);

    std::vector< std::pair< stk::mesh::EntityId, stk::mesh::EntityId > > remote_edges;
    communicate_remote_connections_to_delete(local_elem_and_remote_connected_elem, remote_edges);

    clear_deleted_element_connections(elements_to_delete);

    delete_remote_connections(remote_edges);

    reconnect_volume_elements_across_deleted_shells(shellConnectivityList, elements_to_delete);
}

stk::mesh::ConnectivityOrdinal ElemElemGraph::get_neighboring_side_ordinal(const stk::mesh::BulkData &mesh,
                                                                           stk::mesh::Entity currentElem,
                                                                           stk::mesh::ConnectivityOrdinal currentOrdinal,
                                                                           stk::mesh::Entity neighborElem)
{
    stk::topology currentElemTopology = mesh.bucket(currentElem).topology();
    stk::topology currentFaceTopology = currentElemTopology.face_topology(currentOrdinal);
    const stk::mesh::Entity* currentElemNodes = mesh.begin_nodes(currentElem);
    stk::mesh::EntityVector currentElemSideNodes(currentFaceTopology.num_nodes());
    currentElemTopology.side_nodes(currentElemNodes, currentOrdinal, currentElemSideNodes.begin());

    stk::topology neighborTopology = mesh.bucket(neighborElem).topology();
    stk::mesh::EntityVector neighborSideNodes;

    bool foundNeighborOrdinal = false;
    unsigned neighborOrdinal = 0;
    for (; neighborOrdinal < neighborTopology.num_faces(); ++neighborOrdinal)
    {
        stk::topology neighborFaceTopology = neighborTopology.face_topology(neighborOrdinal);
        neighborSideNodes.resize(neighborFaceTopology.num_nodes());
        const stk::mesh::Entity* neighborNodes = mesh.begin_nodes(neighborElem);
        neighborTopology.side_nodes(neighborNodes, neighborOrdinal, neighborSideNodes.begin());
        std::pair<bool,unsigned> result = neighborFaceTopology.equivalent(currentElemSideNodes, neighborSideNodes);

        if (result.first && result.second >= neighborFaceTopology.num_positive_permutations())
        {
            foundNeighborOrdinal = true;
            break;
        }
    }
    ThrowRequireMsg(foundNeighborOrdinal, "Error: neighborElem is not a true neighbor of currentElem.");
    return static_cast<stk::mesh::ConnectivityOrdinal>(neighborOrdinal);
}

size_t ElemElemGraph::find_max_local_offset_in_neighborhood(stk::mesh::Entity element)
{
    stk::mesh::EntityVector side_nodes;
    stk::mesh::EntityVector connected_elements;
    stk::mesh::Bucket &elem_bucket = m_bulk_data.bucket(element);
    stk::topology topology = elem_bucket.topology();
    size_t num_sides = topology.num_sides();
    const stk::mesh::Entity* elem_nodes = elem_bucket.begin_nodes(m_bulk_data.bucket_ordinal(element));
    size_t max_local_offset = 0;
    size_t current_offset = element.local_offset();
    if (current_offset > max_local_offset)
    {
        max_local_offset = current_offset;
    }
    for (size_t side_index = 0; side_index < num_sides; ++side_index)
    {
        unsigned num_side_nodes = topology.side_topology(side_index).num_nodes();
        side_nodes.resize(num_side_nodes);
        topology.side_nodes(elem_nodes, side_index, side_nodes.begin());
        connected_elements.clear();
        impl::ConnectedElementDataVector connectedElementDataVector;
        impl::find_locally_owned_elements_these_nodes_have_in_common(m_bulk_data, num_side_nodes, side_nodes.data(), connected_elements);

        for (stk::mesh::Entity & connected_element : connected_elements)
        {
            current_offset = connected_element.local_offset();
            if (current_offset > max_local_offset)
            {
                max_local_offset = current_offset;
            }
        }
    }
    return max_local_offset;
}

void ElemElemGraph::generate_additional_ids_collective(size_t num_additional_ids_needed)
{
    std::vector<stk::mesh::EntityId> new_ids;
    m_bulk_data.generate_new_ids_given_reserved_ids(m_bulk_data.mesh_meta_data().side_rank(), num_additional_ids_needed, m_suggested_side_ids, new_ids);
    m_suggested_side_ids.insert(m_suggested_side_ids.end(), new_ids.begin(), new_ids.end());
}

void ElemElemGraph::add_elements_to_graph(const stk::mesh::EntityVector &elements_to_add)
{
    size_t max_offset = 0;
    for (const stk::mesh::Entity & element_to_add : elements_to_add)
    {
        size_t local_max = find_max_local_offset_in_neighborhood(element_to_add);
        if (local_max > max_offset)
        {
            max_offset = local_max;
        }
    }
    ensure_space_in_entity_to_local_id(max_offset);

    std::vector<impl::ElementSidePair> elem_side_pairs;
    size_t num_local_edges_needed = 0;

    for(unsigned i=0; i<elements_to_add.size(); ++i)
    {
        std::set<EntityId> localElementsConnectedToNewShell;
        stk::mesh::Entity elem_to_add = elements_to_add[i];
        if (!m_bulk_data.bucket(elem_to_add).owned())
        {
            continue;
        }
        ThrowRequire(!is_valid_graph_element(elem_to_add));
        impl::LocalId new_elem_id = get_new_local_element_id_from_pool();
        m_local_id_to_element_entity[new_elem_id] = elem_to_add;
        m_entity_to_local_id[elem_to_add.local_offset()] = new_elem_id;
        stk::topology elem_topology = m_bulk_data.bucket(elem_to_add).topology();
        m_element_topologies[new_elem_id] = elem_topology;
        get_element_side_pairs(m_bulk_data.mesh_index(elem_to_add), new_elem_id, elem_side_pairs);
        for(size_t index=0; index<elem_side_pairs.size(); ++index)
        {
            stk::mesh::Entity neighbor = m_local_id_to_element_entity[elem_side_pairs[index].first];
            if (is_valid_graph_element(neighbor))
            {
                m_elem_graph[new_elem_id].push_back(elem_side_pairs[index].first);
                m_via_sides[new_elem_id].push_back(elem_side_pairs[index].second);
                ++m_num_edges;
                impl::LocalId neighbor_id = m_entity_to_local_id[neighbor.local_offset()];
                stk::mesh::ConnectivityOrdinal currentOrdinal = static_cast<stk::mesh::ConnectivityOrdinal>(elem_side_pairs[index].second);
                stk::mesh::ConnectivityOrdinal neighborOrdinal = get_neighboring_side_ordinal(m_bulk_data, elem_to_add, currentOrdinal, neighbor);
                m_elem_graph[neighbor_id].push_back(new_elem_id);
                m_via_sides[neighbor_id].push_back(neighborOrdinal);
                ++m_num_edges;
                num_local_edges_needed+=2;
                if (elem_topology.is_shell()) {
                    localElementsConnectedToNewShell.insert(neighbor_id);
                }
            }
        }
        break_volume_element_connections_across_shells(localElementsConnectedToNewShell);
    }

    impl::ElemSideToProcAndFaceId elem_side_comm = impl::get_element_side_ids_to_communicate(m_bulk_data);

    size_t num_additional_parallel_edges = elem_side_comm.size() - m_num_parallel_edges;
    size_t num_additional_side_ids_needed =  num_additional_parallel_edges + num_local_edges_needed;
    this->generate_additional_ids_collective(num_additional_side_ids_needed);

    stk::mesh::EntityVector elements_to_add_copy = elements_to_add;
    std::sort(elements_to_add_copy.begin(), elements_to_add_copy.end());

    impl::ElemSideToProcAndFaceId only_added_elements;
    impl::ElemSideToProcAndFaceId::iterator iter = elem_side_comm.begin();
    for(;iter!=elem_side_comm.end();++iter)
    {
        stk::mesh::Entity element = iter->first.entity;
        stk::mesh::EntityVector::iterator elem_iter = std::lower_bound(elements_to_add_copy.begin(), elements_to_add_copy.end(), element);
        if(elem_iter!=elements_to_add_copy.end() && *elem_iter==element)
        {
            only_added_elements.insert(*iter);
        }
    }

    fill_parallel_graph(only_added_elements, m_part);

    update_number_of_parallel_edges();
}

void change_entity_owner(stk::mesh::BulkData &bulkData, stk::mesh::ElemElemGraph &elem_graph, std::vector< std::pair< stk::mesh::Entity, int > > &elem_proc_pairs_to_move, stk::mesh::Part *active_part)
{
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
    impl::ParallelGraphInfo new_parallel_graph_entries;

    const std::vector<stk::mesh::EntityId> suggested_face_id_vector = elem_graph.get_suggested_side_ids();
    size_t num_suggested_face_ids_used = 0;

    for (size_t i=0; i<elem_proc_pairs_to_move.size(); i++)
    {
        stk::mesh::Entity elem_to_send = elem_proc_pairs_to_move[i].first;
        int destination_proc = elem_proc_pairs_to_move[i].second;
        stk::mesh::EntityId elem_global_id = bulkData.identifier(elem_to_send);

        size_t num_connected_elements = elem_graph.get_num_connected_elems(elem_to_send);
        stk::topology elem_topology = bulkData.bucket(elem_to_send).topology();
        const stk::mesh::Entity *elem_nodes = bulkData.begin_nodes(elem_to_send);
        for (size_t k=0; k<num_connected_elements; k++)
        {
            if (elem_graph.is_connected_elem_locally_owned(elem_to_send, k))
            {
                int side_id = elem_graph.get_side_id_to_connected_element(elem_to_send, k);
                stk::mesh::Entity connected_element = elem_graph.get_connected_element(elem_to_send, k);
                impl::LocalId local_id = elem_graph.get_local_element_id(connected_element);
                stk::topology side_topology = elem_topology.side_topology(side_id);
                std::vector<stk::mesh::Entity> side_nodes(side_topology.num_nodes());

                elem_topology.side_nodes(elem_nodes, side_id, side_nodes.begin());
                stk::mesh::OrdinalAndPermutation ordperm = get_ordinal_and_permutation(bulkData, elem_to_send, side_rank, side_nodes);

                std::pair<impl::LocalId, stk::mesh::EntityId> key(local_id, elem_global_id);
                stk::mesh::EntityId face_id = suggested_face_id_vector[num_suggested_face_ids_used];
                num_suggested_face_ids_used++;
                bool inActivePart = true;
                if(active_part != NULL)
                {
                    inActivePart = bulkData.bucket(connected_element).member(*active_part);
                }
                impl::parallel_info p_info(destination_proc, side_id, ordperm.second, face_id, inActivePart);
                new_parallel_graph_entries.insert(std::make_pair(key, p_info));
            }
        }
    }

    elem_graph.set_num_side_ids_used(num_suggested_face_ids_used);

    bulkData.change_entity_owner(elem_proc_pairs_to_move);

    elem_graph.change_entity_owner(elem_proc_pairs_to_move, new_parallel_graph_entries);

}

}} // end namespaces stk mesh

