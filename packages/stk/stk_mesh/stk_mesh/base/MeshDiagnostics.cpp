#include <stk_mesh/base/MeshDiagnostics.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include <stk_mesh/base/ExodusTranslator.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemFilter.hpp"
#include "stk_mesh/baseImpl/elementGraph/MeshDiagnosticObserver.hpp"
#include "stk_mesh/baseImpl/EquivalentEntityBlocks.hpp"
#include "stk_util/parallel/DistributedIndex.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemGraphCoincidentElems.hpp"
#include <stddef.h>                     // for size_t, nullptr
#include <map>
#include <string>

namespace stk { namespace mesh {

void fill_single_split_coincident_connection(const stk::mesh::BulkData &bulk, stk::mesh::Entity localElem, const impl::ParallelElementData & connectedElemParallelData,
                                               SplitCoincidentInfo & splitCoincidents)
{
    splitCoincidents[bulk.identifier(localElem)].push_back({connectedElemParallelData.get_element_identifier(), connectedElemParallelData.get_proc_rank_of_neighbor()});
}

void fill_split_coincident_connections_for_elem(const stk::mesh::BulkData &bulk, stk::mesh::Entity localElem, const stk::mesh::EntityVector& localElemSideNodes, const impl::ParallelElementData &localElementData, const impl::ParallelElementDataVector &elementsConnectedToThisElementSide,
                                                  SplitCoincidentInfo & splitCoincidents)
{
    for(const impl::ParallelElementData & connectedElemParallelData : elementsConnectedToThisElementSide)
        if(connectedElemParallelData.is_parallel_edge())
            if(impl::is_coincident_connection(bulk, localElem, localElemSideNodes, localElementData.get_element_side_index(), connectedElemParallelData.get_element_topology(), connectedElemParallelData.get_side_nodes()))
                fill_single_split_coincident_connection(bulk, localElem, connectedElemParallelData, splitCoincidents);
}

void fill_split_coincident_connections(const stk::mesh::BulkData & bulk, const impl::ElementLocalIdMapper & localMapper,
                                         const impl::ParallelElementDataVector &localElementsAttachedToReceivedNodes,
                                         stk::mesh::impl::ParallelElementDataVector & remoteElementsConnectedToSideNodes,
                                         SplitCoincidentInfo & splitCoincidents)
{
    const size_t numReceivedNodes = remoteElementsConnectedToSideNodes[0].get_side_nodes().size();
    for (const impl::ParallelElementData &localElementData: localElementsAttachedToReceivedNodes)
    {
        stk::mesh::Entity localElem = localMapper.local_to_entity(localElementData.get_element_local_id());
        if (localElementData.get_side_nodes().size() == numReceivedNodes)
            fill_split_coincident_connections_for_elem(bulk, localElem, localElementData.get_side_nodes(), localElementData, remoteElementsConnectedToSideNodes, splitCoincidents);
    }
}

SideNodeToReceivedElementDataMap get_element_sides_from_other_procs(stk::mesh::BulkData & bulkData)
{
    impl::ElemSideProcVector elementSideIdsToSend = impl::gather_element_side_ids_to_send(bulkData);
    SharedSidesCommunication sharedSidesCommunication(bulkData, elementSideIdsToSend);
    sharedSidesCommunication.communicate_element_sides();
    SideNodeToReceivedElementDataMap elementSidesReceived;
    sharedSidesCommunication.unpack_side_data_map(elementSidesReceived);
    return elementSidesReceived;
}

SplitCoincidentInfo get_split_coincident_elements_from_received_element_sides(stk::mesh::BulkData& bulkData, const impl::ElementLocalIdMapper & localIdMapper, SideNodeToReceivedElementDataMap & elementSidesReceived)
{
    SplitCoincidentInfo splitCoincidents;
    stk::mesh::EntityVector scratchNodeVector;
    impl::ParallelElementDataVector localElementsAttachedToReceivedNodes;
    for (SideNodeToReceivedElementDataMap::value_type & receivedElementData: elementSidesReceived)
    {
        localElementsAttachedToReceivedNodes.clear();
        stk::mesh::impl::ParallelElementDataVector &parallelElementDatas = receivedElementData.second;
        impl::get_elements_connected_via_sidenodes<impl::ParallelElementData>(bulkData,
                                                                              parallelElementDatas[0].get_element_identifier(),
                                                                              parallelElementDatas[0].get_element_topology(),
                                                                              localIdMapper,
                                                                              parallelElementDatas[0].get_side_nodes(),
                                                                              scratchNodeVector, localElementsAttachedToReceivedNodes);
        fill_split_coincident_connections(bulkData, localIdMapper, localElementsAttachedToReceivedNodes, parallelElementDatas, splitCoincidents);
    }
    return splitCoincidents;
}

SplitCoincidentInfo get_split_coincident_elements(stk::mesh::BulkData& bulkData)
{
    SideNodeToReceivedElementDataMap elementSidesReceived = get_element_sides_from_other_procs(bulkData);

    impl::ElementLocalIdMapper localIdMapper;
    localIdMapper.initialize(bulkData);
    return get_split_coincident_elements_from_received_element_sides(bulkData, localIdMapper, elementSidesReceived);
}

std::vector<std::string> get_messages_for_split_coincident_elements(const stk::mesh::BulkData& bulkData, const stk::mesh::SplitCoincidentInfo & splitCoincidentElements)
{
    std::vector<std::string> errorList;
    std::ostringstream out;
    for(const auto& item : splitCoincidentElements) {
        out.str(std::string());
        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK,item.first);
        const stk::mesh::PartVector& elementParts = bulkData.bucket(element).supersets();
        std::string blockNames;
        blockNames = "{";
        for (const stk::mesh::Part* part : elementParts) {
            if (stk::mesh::is_element_block(*part)) {
                blockNames += " " + part->name();
            }
        }
        blockNames += " }";
        for(size_t i=0;i<item.second.size();++i)
            out << "ERROR: [" << bulkData.parallel_rank() << "] Element " << item.first << " (" << bulkData.bucket(element).topology() << ") in blocks " << blockNames << " is coincident with element " << item.second[i].first << " on processor " << item.second[i].second << std::endl;
        errorList.push_back(out.str());
    }
    return errorList;
}

void throw_if_any_proc_has_false(MPI_Comm comm, bool is_all_ok_locally)
{
    bool is_all_ok_globally = stk::is_true_on_all_procs(comm, is_all_ok_locally);
    STK_ThrowRequireMsg(is_all_ok_globally, "Mesh diagnostics failed.");
}

stk::mesh::Selector get_owned_or_shared_selector(const stk::mesh::BulkData & bulkData)
{
    return bulkData.mesh_meta_data().locally_owned_part() | bulkData.mesh_meta_data().globally_shared_part();
}

stk::parallel::DistributedIndex::KeyTypeVector get_all_local_keys(const stk::mesh::BulkData & bulkData)
{
    stk::parallel::DistributedIndex::KeyTypeVector localKeys;
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK;rank < bulkData.mesh_meta_data().entity_rank_count();++rank)
    {
        stk::mesh::EntityVector entities;
        stk::mesh::get_selected_entities(get_owned_or_shared_selector(bulkData), bulkData.buckets(rank), entities);
        for(stk::mesh::Entity entity: entities)
            localKeys.push_back(bulkData.entity_key(entity));
    }
    return localKeys;
}

void add_keys_to_distributed_index(const stk::mesh::BulkData & bulkData, stk::parallel::DistributedIndex & distributedIndex)
{
    stk::parallel::DistributedIndex::KeyTypeVector localKeys = get_all_local_keys(bulkData);

    stk::parallel::DistributedIndex::KeyTypeVector::const_iterator begin = localKeys.begin();
    stk::parallel::DistributedIndex::KeyTypeVector::const_iterator end = localKeys.end();
    distributedIndex.update_keys( begin, end );
}

std::vector<stk::mesh::EntityKeyProc> get_non_unique_keys(const stk::mesh::BulkData& bulkData, const stk::parallel::DistributedIndex& distributedIndex,
        const stk::parallel::DistributedIndex::KeyTypeVector& localKeys)
{
    stk::parallel::DistributedIndex::KeyProcVector sharedKeyProcs;
    distributedIndex.query_to_usage(localKeys, sharedKeyProcs);

    std::vector<stk::mesh::EntityKeyProc> badKeys;
    for (const stk::parallel::DistributedIndex::KeyProc& sharedKeyProc : sharedKeyProcs)
    {
        stk::mesh::EntityKey key( static_cast<stk::mesh::EntityKey::entity_key_t>(sharedKeyProc.first) );
        if ( bulkData.parallel_rank() != sharedKeyProc.second )
        {
            if(!bulkData.in_shared(key, sharedKeyProc.second))
                badKeys.push_back({key, sharedKeyProc.second});
        }
    }
    return badKeys;
}

std::string get_topology(stk::topology topology)
{
    if(topology==stk::topology::INVALID_TOPOLOGY)
        return " ";
    return " (" + topology.name() + ") ";
}

std::vector<stk::mesh::EntityKeyProc> get_non_unique_key_procs(const stk::mesh::BulkData& bulkData)
{
    stk::parallel::DistributedIndex distributedIndex( bulkData.parallel(), stk::mesh::impl::convert_entity_keys_to_spans(bulkData.mesh_meta_data()));
    add_keys_to_distributed_index(bulkData, distributedIndex);
    stk::parallel::DistributedIndex::KeyTypeVector localKeys = get_all_local_keys(bulkData);
    return get_non_unique_keys(bulkData, distributedIndex, localKeys);
}

std::vector<std::string> get_non_unique_key_messages(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKeyProc> &badKeyProcs)
{
    std::vector<std::string> errorList;
    std::ostringstream os;
    for(const stk::mesh::EntityKeyProc& keyProc : badKeyProcs)
    {
        os.str(std::string());
        stk::mesh::Entity entity = bulkData.get_entity(keyProc.first);
        os << "ERROR: [" << bulkData.parallel_rank() << "] Key " << keyProc.first <<
                get_topology(bulkData.bucket(entity).topology()) << "is also present (inappropriately) on processor " <<
                keyProc.second << "." << std::endl;
        errorList.push_back(os.str());
    }
    return errorList;
}


std::vector<stk::mesh::Entity> get_orphaned_owned_sides(const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()), sides);
    std::vector<stk::mesh::Entity> badSides;
    for(stk::mesh::Entity side : sides)
    {
        unsigned num_elements = bulkData.num_elements(side);
        const stk::mesh::Entity* elements = bulkData.begin_elements(side);
        size_t num_owned_elements = 0;
        for(unsigned i=0;i<num_elements;++i)
        {
            if(bulkData.bucket(elements[i]).owned())
                num_owned_elements++;
        }
        if(num_owned_elements == 0)
            badSides.push_back(side);
    }
    return badSides;
}

std::vector<stk::mesh::Entity> get_local_solo_sides(const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()), sides);
    std::vector<stk::mesh::Entity> localSoloSides;
    for(stk::mesh::Entity side : sides)
    {
        unsigned num_elements = bulkData.num_elements(side);
        if (0 == num_elements)
        {
          stk::mesh::EntityVector connectedElements;
          stk::mesh::impl::find_entities_these_nodes_have_in_common(bulkData, stk::topology::ELEM_RANK, bulkData.num_nodes(side), bulkData.begin_nodes(side), connectedElements);
          if (connectedElements.empty())
          {
            localSoloSides.push_back(side);
          }
        }
    }
    return localSoloSides;
}

void pack_side_node_keys(const std::vector<stk::mesh::Entity>& orphanedSides,
                         const stk::mesh::BulkData& bulkData,
                         stk::CommSparse &comm)
{
    std::vector<int> sharingProcs;
    std::vector<stk::mesh::Entity> nodesVec;
    std::vector<stk::mesh::EntityKey> nodeKeys;
    for(stk::mesh::Entity side : orphanedSides)
    {
        const stk::mesh::Entity* nodes = bulkData.begin_nodes(side);
        unsigned numNodes = bulkData.num_nodes(side);
        nodesVec.resize(numNodes);
        nodeKeys.resize(numNodes);
        for(unsigned i = 0; i < numNodes; ++i) {
            nodesVec[i] = nodes[i];
            nodeKeys[i] = bulkData.entity_key(nodes[i]);
        }
        bulkData.shared_procs_intersection(nodesVec, sharingProcs);
        stk::mesh::EntityKey sideKey = bulkData.entity_key(side);
        for(int proc : sharingProcs)
        {
            comm.send_buffer(proc).pack<stk::mesh::EntityKey>(sideKey);
            stk::pack_vector_to_proc(comm, nodeKeys, proc);
        }
    }
}

void unpack_side_nodes_and_check_for_attached_locally_owned_elements(const stk::mesh::BulkData& bulkData,
                                                       stk::CommSparse &comm,
                                                       std::map< std::pair<stk::mesh::EntityKey, int>, bool > &sideKeyProcMap)
{
    for(int proc_id=0; proc_id<bulkData.parallel_size(); ++proc_id)
    {
        if (proc_id != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityKey sideKey;

                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(sideKey);

                std::vector<stk::mesh::EntityKey> nodeKeys;
                stk::unpack_vector_from_proc(comm, nodeKeys, proc_id);

                std::vector<stk::mesh::Entity> nodes(nodeKeys.size());
                for(unsigned i =0; i<nodeKeys.size(); ++i)
                    nodes[i] = bulkData.get_entity(nodeKeys[i]);

                stk::mesh::EntityVector connectedElements;
                stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(
                    bulkData, nodes.size(), nodes.data(), connectedElements);

                sideKeyProcMap[std::make_pair(sideKey, proc_id)] = !connectedElements.empty();
            }
        }
    }
}

void unpack_side_nodes_and_check_for_elements(const stk::mesh::BulkData& bulkData,
                                                       stk::CommSparse &comm,
                                                       std::map< std::pair<stk::mesh::EntityKey, int>, bool > &sideKeyProcMap)
{
    for(int proc_id=0; proc_id<bulkData.parallel_size(); ++proc_id)
    {
        if (proc_id != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityKey sideKey;

                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(sideKey);

                std::vector<stk::mesh::EntityKey> nodeKeys;
                stk::unpack_vector_from_proc(comm, nodeKeys, proc_id);

                std::vector<stk::mesh::Entity> nodes(nodeKeys.size());
                for(unsigned i =0; i<nodeKeys.size(); ++i)
                    nodes[i] = bulkData.get_entity(nodeKeys[i]);

                stk::mesh::EntityVector connectedElements;
                stk::mesh::impl::find_entities_these_nodes_have_in_common(
                    bulkData, stk::topology::ELEM_RANK, nodes.size(), nodes.data(), connectedElements);

                sideKeyProcMap[std::make_pair(sideKey, proc_id)] = !connectedElements.empty();
            }
        }
    }
}

void pack_side_key_and_response(const std::map<std::pair<stk::mesh::EntityKey, int>, bool>& sideKeyProcMap, stk::CommSparse& comm_to)
{
    for(const std::map<std::pair<stk::mesh::EntityKey, int>, bool>::value_type& data : sideKeyProcMap)
    {
        const std::pair<stk::mesh::EntityKey, int>& keyProcPair = data.first;
        int proc = keyProcPair.second;
        comm_to.send_buffer(proc).pack<stk::mesh::EntityKey>(keyProcPair.first);
        comm_to.send_buffer(proc).pack<unsigned>(data.second);
    }
}

void unpack_side_key_and_response(const stk::mesh::BulkData& bulkData,
                     stk::CommSparse& comm_to,
                     std::map<stk::mesh::EntityKey, bool>& sideKeyMap)
{
    for(int proc_id = 0; proc_id < bulkData.parallel_size(); ++proc_id)
    {
        if(proc_id != bulkData.parallel_rank())
        {
            while(comm_to.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityKey sideKey;
                unsigned sideHasRemoteElement;
                comm_to.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(sideKey);
                comm_to.recv_buffer(proc_id).unpack<unsigned>(sideHasRemoteElement);
                if(0 != sideHasRemoteElement)
                    sideKeyMap[sideKey] = true;
            }
        }
    }
}

void exchange_side_connection_info_for_orphan_sides(const stk::mesh::BulkData& bulkData,
                                   const std::vector<stk::mesh::Entity> &orphanedSides,
                                   std::map< std::pair<stk::mesh::EntityKey, int>, bool > &sideKeyProcMap)
{
    stk::CommSparse comm(bulkData.parallel());
    pack_side_node_keys(orphanedSides, bulkData, comm);
    comm.allocate_buffers();

    pack_side_node_keys(orphanedSides, bulkData, comm);
    comm.communicate();

    unpack_side_nodes_and_check_for_attached_locally_owned_elements(bulkData, comm, sideKeyProcMap);
}

void exchange_side_connection_info_for_solo_sides(const stk::mesh::BulkData& bulkData,
                                   const std::vector<stk::mesh::Entity> &orphanedSides,
                                   std::map< std::pair<stk::mesh::EntityKey, int>, bool > &sideKeyProcMap)
{
    stk::CommSparse comm(bulkData.parallel());
    pack_side_node_keys(orphanedSides, bulkData, comm);
    comm.allocate_buffers();

    pack_side_node_keys(orphanedSides, bulkData, comm);
    comm.communicate();

    unpack_side_nodes_and_check_for_elements(bulkData, comm, sideKeyProcMap);
}

void populate_side_key_map(const stk::mesh::BulkData& bulkData,
                           const std::map< std::pair<stk::mesh::EntityKey, int>, bool > &sideKeyProcMap,
                           std::map<stk::mesh::EntityKey, bool> &sideKeyMap)
{
    stk::CommSparse comm(bulkData.parallel());
    pack_side_key_and_response(sideKeyProcMap, comm);
    comm.allocate_buffers();

    pack_side_key_and_response(sideKeyProcMap, comm);
    comm.communicate();

    unpack_side_key_and_response(bulkData, comm, sideKeyMap);
}

void communicate_side_nodes_and_check_for_attached_elements(const stk::mesh::BulkData& bulkData,
                                                            const std::vector<stk::mesh::Entity> &orphanedSides,
                                                            std::map<stk::mesh::EntityKey, bool> &sideKeyMap)
{
    std::map< std::pair<stk::mesh::EntityKey, int>, bool > sideKeyProcMap;
    exchange_side_connection_info_for_orphan_sides(bulkData, orphanedSides, sideKeyProcMap);
    populate_side_key_map(bulkData, sideKeyProcMap, sideKeyMap);
}

void communicate_local_solo_side_nodes_and_check_for_elements(const stk::mesh::BulkData& bulkData,
                                                            const std::vector<stk::mesh::Entity> &orphanedSides,
                                                            std::map<stk::mesh::EntityKey, bool> &sideKeyMap)
{
    std::map< std::pair<stk::mesh::EntityKey, int>, bool > sideKeyProcMap;
    exchange_side_connection_info_for_solo_sides(bulkData, orphanedSides, sideKeyProcMap);
    populate_side_key_map(bulkData, sideKeyProcMap, sideKeyMap);
}

std::vector<stk::mesh::Entity> get_orphaned_sides_with_attached_element_on_different_proc(const stk::mesh::BulkData& bulkData)
{
    std::vector<stk::mesh::Entity> orphanedSides = get_orphaned_owned_sides(bulkData);
    std::map<stk::mesh::EntityKey, bool> sideKeyMap;

    for(stk::mesh::Entity entity: orphanedSides)
    {
        sideKeyMap[bulkData.entity_key(entity)] = false;
    }

    communicate_side_nodes_and_check_for_attached_elements(bulkData, orphanedSides, sideKeyMap);

    std::vector<stk::mesh::Entity> badSides;
    for(std::map<stk::mesh::EntityKey, bool>::value_type data : sideKeyMap)
    {
        if(true == data.second)
            badSides.push_back(bulkData.get_entity(data.first));
    }
    return badSides;
}

std::vector<stk::mesh::Entity> get_solo_sides_without_element_on_different_proc(const stk::mesh::BulkData& bulkData)
{
    std::vector<stk::mesh::Entity> orphanedSides = get_local_solo_sides(bulkData);
    std::map<stk::mesh::EntityKey, bool> sideKeyMap;

    for(stk::mesh::Entity entity: orphanedSides)
    {
        sideKeyMap[bulkData.entity_key(entity)] = false;
    }

    communicate_local_solo_side_nodes_and_check_for_elements(bulkData, orphanedSides, sideKeyMap);

    std::vector<stk::mesh::Entity> badSides;
    for(std::map<stk::mesh::EntityKey, bool>::value_type data : sideKeyMap)
    {
        if(false == data.second)
            badSides.push_back(bulkData.get_entity(data.first));
    }
    return badSides;
}

std::vector<std::string> get_messages_for_orphaned_owned_sides(const stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& entities)
{
    std::vector<std::string> errorList;
    std::ostringstream os;
    for(const stk::mesh::Entity& entity : entities)
    {
        os.str(std::string());
        os << "ERROR: [" << bulkData.parallel_rank() << "] Side " << bulkData.entity_key(entity) << " (" << bulkData.bucket(entity).topology()
           << ") does not have upwards relations to a locally owned element. Nodes of side are {";
        unsigned num_nodes = bulkData.num_nodes(entity);
        const stk::mesh::Entity* nodes = bulkData.begin_nodes(entity);
        for(unsigned i=0;i<num_nodes;i++)
        {
            os << bulkData.entity_key(nodes[i]);
            if(i != num_nodes-1)
                os << ", ";
        }
        os << "}.\n";
        errorList.push_back(os.str());
    }
    return errorList;
}

std::vector<std::string> get_messages_for_solo_sides(const stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& entities)
{
    std::vector<std::string> errorList;
    std::ostringstream os;
    for(const stk::mesh::Entity& entity : entities)
    {
        os.str(std::string());
        os << "ERROR: [" << bulkData.parallel_rank() << "] solo side " << bulkData.entity_key(entity) << " (" << bulkData.bucket(entity).topology()
                << ") nodes = {";
        unsigned num_nodes = bulkData.num_nodes(entity);
        const stk::mesh::Entity* nodes = bulkData.begin_nodes(entity);
        for(unsigned i=0;i<num_nodes;i++)
        {
            os << bulkData.entity_key(nodes[i]);
            if(i != num_nodes-1)
                os << ", ";
        }
        os << "}.\n";
        errorList.push_back(os.str());
    }
    if(!entities.empty())
    {
      std::ofstream out("solo_faces." + std::to_string(bulkData.parallel_size()) + "." + std::to_string(bulkData.parallel_rank()),std::ios_base::app);
      for (const std::string& s : errorList) { out << s; }
      out.close();
    }
    return errorList;
}

} }
