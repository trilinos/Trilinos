/*
 * NoGhostGameofLife.cpp
 *
 *  Created on: Aug 11, 2015
 *      Author: jonchu
 */
#include "NoGhostGameofLife.hpp"
#include <stddef.h>                     // for size_t
#include <utility>                      // for pair
#include "stk_mesh/base/HashEntityAndEntityKey.hpp"  // for hash
#include "GameOfLife/GameofLifeMesh.hpp"  // for GameofLifeMesh, etc
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/GetEntities.hpp"  // for get_entities
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire
#include "stk_util/parallel/CommSparse.hpp"  // for CommSparse
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer, unpack

//public
NoGhostGameofLife::NoGhostGameofLife(GameofLifeMesh& Mesh, std::string name)
: m_bulkData(Mesh.bulk_data()),
  m_lifeField(Mesh.life_field()),
  m_neighborField(Mesh.neighbor_field()),
  m_name(name), m_stkIo(Mesh.comm())
{
    finish_construction();
}

void NoGhostGameofLife::activate_these_ids(stk::mesh::EntityIdVector& elemIds)
{
   for (stk::mesh::EntityId elemId : elemIds)
       activate_element_id(elemId);
}

void NoGhostGameofLife::run_game_of_life(int numSteps)
{
    if (0 == m_time)
        write_output_step();

    for (int timeStep = 0; timeStep < numSteps; timeStep++)
        run_game_of_life_step();
}

//test functions
stk::mesh::Entity NoGhostGameofLife::element_with_id(stk::mesh::EntityId elemId)
{
   return m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

}

bool NoGhostGameofLife::is_valid_entity(stk::mesh::Entity entity)
{
    return m_bulkData.is_valid(entity);
}

unsigned NoGhostGameofLife::num_neighbors(stk::mesh::Entity elem)
{
    return m_localElementToLocalNeighborElements[elem].size() +
            m_localElementToRemoteElementKeys[elem].size();
}

unsigned NoGhostGameofLife::num_active_elems()
{
    unsigned count = 0;
    for (stk::mesh::Entity localElem : m_elements)
        if (is_element_active(localElem))
            count++;
    return count;
}

unsigned NoGhostGameofLife::num_active_neighbors(stk::mesh::Entity elem)
{
    unsigned numActiveNeighbors = 0;
    //local elements
    for (stk::mesh::Entity localElem : m_localElementToLocalNeighborElements[elem])
        if (is_element_active(localElem))
            numActiveNeighbors++;

    //send the remote entity keys to their processors
    stk::CommSparse send(m_bulkData.parallel());
    for (int phase = 0; phase < 2; phase++)
    {
       for (stk::mesh::EntityKey remoteElemKey : m_localElementToRemoteElementKeys[elem])
       {
           int owningProc = m_remoteElementKeyToOwningProcessor[remoteElemKey];
           send.send_buffer(owningProc).pack<stk::mesh::EntityKey>(remoteElemKey);
       }
       if (0 == phase)
           send.allocate_buffers();
       else if (1 == phase)
           send.communicate();
    }

    //find if those remote (now local) elements are active and send that back
    unsigned numActiveNeighborsOnProc = 0;
    stk::CommSparse recieve(m_bulkData.parallel());

    std::unordered_map<int, unsigned> procToValue;

    for (int procNum = 0; procNum < m_numProcs; procNum++)
    {
        stk::CommBuffer& sendBuf = send.recv_buffer(procNum);
        while (sendBuf.remaining())
        {
            stk::mesh::EntityKey localElemKey;
            sendBuf.unpack<stk::mesh::EntityKey>(localElemKey);
            if (is_element_active(m_bulkData.get_entity(localElemKey)))
                numActiveNeighborsOnProc++;
        }
        procToValue[procNum] = numActiveNeighborsOnProc;
        numActiveNeighborsOnProc = 0;
    }

    for (int phase = 0; phase < 2; phase++)
    {
        for (std::pair<int, unsigned> pair : procToValue)
            recieve.send_buffer(pair.first).pack<unsigned>(pair.second);
        if (0 == phase)
            recieve.allocate_buffers();
        else if (1 == phase)
            recieve.communicate();
    }
    //get the number
    for (int procNum = 0; procNum < m_numProcs; procNum++)
    {
        stk::CommBuffer& recieveBuf= recieve.recv_buffer(procNum);
        while (recieveBuf.remaining())
        {
            unsigned remoteActiveNeighbors;
            recieveBuf.unpack<unsigned>(remoteActiveNeighbors);
            numActiveNeighbors += remoteActiveNeighbors;
        }
    }
    return numActiveNeighbors;
}

//private
bool NoGhostGameofLife::is_element_active(stk::mesh::Entity elem)
{
    return *stk::mesh::field_data(m_lifeField, elem); // field
}
void NoGhostGameofLife::activate_element(stk::mesh::Entity elem)
{
    *stk::mesh::field_data(m_lifeField, elem) = 1;
}
void NoGhostGameofLife::deactivate_element(stk::mesh::Entity elem)
{
    *stk::mesh::field_data(m_lifeField, elem) = 0;
}
void NoGhostGameofLife::finish_construction()
{
    m_numProcs = m_bulkData.parallel_size();
    m_time = 0;
    get_elements();
    create_element_connectivity_maps();
    write_output_mesh();
}
void NoGhostGameofLife::get_elements()
{
    stk::mesh::get_entities(m_bulkData, stk::topology::ELEM_RANK, m_elements);
}
void NoGhostGameofLife::create_element_connectivity_maps()
{
   create_local_element_to_local_element_map();
   create_remote_element_key_maps();
   create_local_element_to_remote_element_key_map();
}
void NoGhostGameofLife::activate_element_id(stk::mesh::EntityId elemId)
{
   stk::mesh::Entity elem = m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
   if (m_bulkData.is_valid(elem))
   {
       activate_element(elem);
       m_localActiveElements.push_back(elem);
   }
}
void NoGhostGameofLife::create_local_element_to_local_element_map()
{
   for (stk::mesh::Entity localElem : m_elements)
       create_map_for_this_element(localElem);
}
void NoGhostGameofLife::create_map_for_this_element(stk::mesh::Entity elem)
{
    const stk::mesh::Entity* localElemNodes = m_bulkData.begin_nodes(elem);
    for (unsigned nodeIndex = 0, numNodes = m_bulkData.num_nodes(elem);
            nodeIndex < numNodes; nodeIndex++)
        add_this_nodes_elements_to_this_elements_map(elem, localElemNodes[nodeIndex]);
}
void NoGhostGameofLife::add_this_nodes_elements_to_this_elements_map(stk::mesh::Entity elem,
                                                                     stk::mesh::Entity node)
{
    const stk::mesh::Entity* localElemNodeElems = m_bulkData.begin_elements(node);
    for (unsigned elemIndex = 0, numElems = m_bulkData.num_elements(node);
            elemIndex < numElems; elemIndex++)
        if (localElemNodeElems[elemIndex] != elem)
            m_localElementToLocalNeighborElements[elem].
            insert(localElemNodeElems[elemIndex]);
}
void NoGhostGameofLife::create_remote_element_key_maps()
{
    stk::CommSparse buffer(m_bulkData.parallel());
    fill_buffer_with_local_element_keys_and_remote_node_keys(buffer);
    unpack_remote_elem_key_info_from_buffer(buffer);
    create_map_of_remote_elem_keys_to_local_elements();
}
void NoGhostGameofLife::fill_buffer_with_local_element_keys_and_remote_node_keys(stk::CommSparse&
                                                                                 buffer)
{
    for (int phase = 0; phase < 2; phase++)
    {
        for (stk::mesh::Entity elem : m_elements)
            fill_buffer_with_this_elements_info(elem, buffer);
        if (0 == phase)
            buffer.allocate_buffers();
        else
            buffer.communicate();
    }
}
void NoGhostGameofLife::fill_buffer_with_this_elements_info(stk::mesh::Entity elem,
                                                         stk::CommSparse& buffer)
{
    std::unordered_map<int,std::unordered_set<stk::mesh::EntityKey,
        std::hash<stk::mesh::EntityKey>>> remoteProcessorNumberToSharedNodes;
    fill_map_with_this_elements_nodes(elem, remoteProcessorNumberToSharedNodes);
    fill_buffer_with_map_info(elem, buffer, remoteProcessorNumberToSharedNodes);
}
void NoGhostGameofLife::fill_map_with_this_elements_nodes(stk::mesh::Entity elem,
                                                          std::unordered_map<int,
                                                          std::unordered_set<stk::mesh::EntityKey,
                                                          std::hash<stk::mesh::EntityKey>>>& map)
{
    unsigned numNodes = m_bulkData.num_nodes(elem);
    const stk::mesh::Entity* nodeBegin = m_bulkData.begin_nodes(elem);
    for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
    {
        std::vector<int> sharingProcs;
        stk::mesh::EntityKey nodeKey = m_bulkData.entity_key(nodeBegin[nodeIndex]);
        m_bulkData.comm_shared_procs(nodeKey, sharingProcs);
        for (int procNum : sharingProcs)
            map[procNum].insert(nodeKey);
    }
}
void NoGhostGameofLife::fill_buffer_with_map_info(stk::mesh::Entity elem, stk::CommSparse& buffer,
                                                  std::unordered_map<int,std::unordered_set
                                                  <stk::mesh::EntityKey, std::hash
                                                  <stk::mesh::EntityKey>>>& map)
{
    for (std::pair< const int,std::unordered_set<stk::mesh::EntityKey,
            std::hash<stk::mesh::EntityKey>>>& pair : map)
    {
        int remoteProc = pair.first;
        buffer.send_buffer(remoteProc).pack<stk::mesh::EntityKey>(m_bulkData.entity_key(elem));
        buffer.send_buffer(remoteProc).pack<size_t>(pair.second.size());
        for (stk::mesh::EntityKey nodeKey : pair.second)
            buffer.send_buffer(remoteProc).pack<stk::mesh::EntityKey>(nodeKey);
    }
}
void NoGhostGameofLife::unpack_remote_elem_key_info_from_buffer(stk::CommSparse& buffer)
{
    for (int proc = 0; proc < m_numProcs; proc++)
    {
        stk::CommBuffer& buf = buffer.recv_buffer(proc);
        while (buf.remaining())
            unpack_remote_info_from_this_processor(proc, buf);
    }
}
void NoGhostGameofLife::unpack_remote_info_from_this_processor(int proc, stk::CommBuffer& buf)
{
    stk::mesh::EntityKey elemKey = stk::unpack<stk::mesh::EntityKey>(buf);
    size_t numNodes = stk::unpack<size_t>(buf);
    m_remoteElementKeys.insert(elemKey);
    m_remoteElementKeyToOwningProcessor[elemKey] = proc;
    for (unsigned nodeNum = 0; nodeNum < numNodes; nodeNum++)
        m_remoteElementKeyToLocalNodeKeys[elemKey].
        insert(stk::unpack<stk::mesh::EntityKey>(buf));
}
void NoGhostGameofLife::create_map_of_remote_elem_keys_to_local_elements()
{
    for (stk::mesh::EntityKey remoteElemKey : m_remoteElementKeys)
        for (stk::mesh::EntityKey localNodeKey : m_remoteElementKeyToLocalNodeKeys[remoteElemKey])
            map_this_remote_element_key_with_this_nodes_elements(remoteElemKey, localNodeKey);
    m_remoteElementKeyToLocalNodeKeys.clear(); //take out the trash
}
void NoGhostGameofLife::map_this_remote_element_key_with_this_nodes_elements(stk::mesh::EntityKey
                                                                             remoteKey,
                                                                             stk::mesh::EntityKey
                                                                             nodeKey)
{
    stk::mesh::Entity node = m_bulkData.get_entity(nodeKey);
    unsigned numElems = m_bulkData.num_elements(node);
    const stk::mesh::Entity* nodeElem = m_bulkData.begin_elements(node);
    for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
        m_remoteElementKeyToLocalNeighborElements[remoteKey].insert(nodeElem[elemIndex]);
}
void NoGhostGameofLife::create_local_element_to_remote_element_key_map()
{
    stk::CommSparse buffer(m_bulkData.parallel());
    fill_buffer_with_local_neighbors_of_remote_keys(buffer);
    unpack_local_and_remote_key_info_from_each_processor(buffer);
}
void NoGhostGameofLife::fill_buffer_with_local_neighbors_of_remote_keys(stk::CommSparse& buffer)
{
    for (int phase = 0; phase < 2; phase++)
    {
        for (stk::mesh::EntityKey remoteElemKey : m_remoteElementKeys)
            fill_buffer_with_local_neighbors_of_remote_element_key(remoteElemKey, buffer);
        if (0 == phase)
            buffer.allocate_buffers();
        else
            buffer.communicate();
    }
}
void NoGhostGameofLife::fill_buffer_with_local_neighbors_of_remote_element_key(stk::mesh::EntityKey
                                                                               remoteKey,
                                                                               stk::CommSparse&
                                                                               buffer)
{
    int procNum = m_remoteElementKeyToOwningProcessor[remoteKey];
    size_t numNeighbors = m_remoteElementKeyToLocalNeighborElements[remoteKey].size();
    buffer.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteKey);
    buffer.send_buffer(procNum).pack<size_t>(numNeighbors);
    for (stk::mesh::Entity localElem : m_remoteElementKeyToLocalNeighborElements[remoteKey])
        buffer.send_buffer(procNum).pack<stk::mesh::EntityKey>(m_bulkData.entity_key(localElem));
}
void NoGhostGameofLife::unpack_local_and_remote_key_info_from_each_processor(stk::CommSparse&
                                                                             buffer)
{
    for (int procRank = 0; procRank < m_numProcs; procRank++)
    {
        stk::CommBuffer& buf = buffer.recv_buffer(procRank);
        while (buf.remaining())
            unpack_local_and_remote_keys_from_buffer(buf);
    }
}
void NoGhostGameofLife::unpack_local_and_remote_keys_from_buffer(stk::CommBuffer& buf)
{
    stk::mesh::EntityKey localElemKey = stk::unpack<stk::mesh::EntityKey>(buf);
    stk::mesh::Entity localElem = m_bulkData.get_entity(localElemKey);
    size_t numRemoteNeighbors = stk::unpack<size_t>(buf);
    for (unsigned neighborNum = 0; neighborNum < numRemoteNeighbors; neighborNum++)
        m_localElementToRemoteElementKeys[localElem].insert(stk::unpack<stk::mesh::EntityKey>(buf));
}
void NoGhostGameofLife::write_output_mesh()
{
    m_name += ".e";
    m_stkIo.set_bulk_data(m_bulkData);
    m_fileHandler = m_stkIo.create_output_mesh(m_name, stk::io::WRITE_RESULTS);
    m_stkIo.add_field(m_fileHandler, m_lifeField);
    m_stkIo.write_output_mesh(m_fileHandler);
}
void NoGhostGameofLife::write_output_step()
{
    m_stkIo.begin_output_step(m_fileHandler, m_time);
    m_stkIo.write_defined_output_fields(m_fileHandler);
    m_stkIo.end_output_step(m_fileHandler);
    m_time++;
}
void NoGhostGameofLife::run_game_of_life_step()
{
    determine_elements_to_check();
    update_neighbor_values_with_local_elements();
    update_neighbor_values_with_remote_elements();
    update_element_membership();
    write_output_step();
}
void NoGhostGameofLife::determine_elements_to_check()
{
    refresh_element_maps();
    stk::CommSparse buffer(m_bulkData.parallel());
    communicate_remote_element_keys_to_check(buffer);
    recieve_local_element_keys_to_check(buffer);
}
void NoGhostGameofLife::refresh_element_maps()
{
    m_localElementsToVisit.clear();
    m_remoteElementKeysToVisit.clear();
    get_elements_to_visit();
}
void NoGhostGameofLife::get_elements_to_visit()
{
    for (stk::mesh::Entity localElem : m_localActiveElements)
    {
        m_localElementsToVisit.insert(localElem);
        for (stk::mesh::Entity localElemElem : m_localElementToLocalNeighborElements[localElem])
            m_localElementsToVisit.insert(localElemElem);
        for (stk::mesh::EntityKey remoteElemKey : m_localElementToRemoteElementKeys[localElem])
            m_remoteElementKeysToVisit.insert(remoteElemKey);
    }
}
void NoGhostGameofLife::communicate_remote_element_keys_to_check(stk::CommSparse& buffer)
{
    for (int phase = 0; phase < 2; phase++)
    {
        for (stk::mesh::EntityKey remoteElemKey : m_remoteElementKeysToVisit)
            buffer.send_buffer(m_remoteElementKeyToOwningProcessor[remoteElemKey]).
            pack<stk::mesh::EntityKey>(remoteElemKey);
        if (0 == phase)
            buffer.allocate_buffers();
        else
            buffer.communicate();
    }
}
void NoGhostGameofLife::recieve_local_element_keys_to_check(stk::CommSparse& buffer)
{
    for (int procNum = 0; procNum < m_numProcs; procNum++)
    {
        stk::CommBuffer& buf = buffer.recv_buffer(procNum);
        while (buf.remaining())
            m_localElementsToVisit.insert(m_bulkData.
            get_entity(stk::unpack<stk::mesh::EntityKey>(buf)));
    }
}
void NoGhostGameofLife::update_neighbor_values_with_local_elements()
{
    for (stk::mesh::Entity localElem : m_localElementsToVisit)
    {
        int* neighborVal = stk::mesh::field_data(m_neighborField, localElem);
        *neighborVal = 0;
        for (stk::mesh::Entity localElemElem : m_localElementToLocalNeighborElements[localElem])
            if (is_element_active(localElemElem))
                (*neighborVal)++;
    }
}
void NoGhostGameofLife::update_neighbor_values_with_remote_elements()
{
    stk::CommSparse buffer(m_bulkData.parallel());
    send_num_active_neighbors_of_remote_elem_keys(buffer);
    recieve_num_active_neighbors_of_local_elements(buffer);
}
void NoGhostGameofLife::send_num_active_neighbors_of_remote_elem_keys(stk::CommSparse& buffer)
{
    for (int phase = 0 ; phase < 2; phase++)
    {
        for (stk::mesh::EntityKey remoteElemKey : m_remoteElementKeysToVisit)
            pack_number_of_local_neighbors_of_remote_element_into_buffer(buffer, remoteElemKey);
        if (0 == phase)
            buffer.allocate_buffers();
        else if (1 == phase)
            buffer.communicate();
    }
}
void NoGhostGameofLife::pack_number_of_local_neighbors_of_remote_element_into_buffer(stk::CommSparse&
                                                                                     buffer,
                                                                                     stk::mesh::
                                                                                     EntityKey
                                                                                     remoteKey)
{
    int numActive = count_local_active_neighbors_for_remote_element_key(remoteKey);
    pack_num_active_neighbors_into_buffer(buffer, numActive, remoteKey);
}
int NoGhostGameofLife::count_local_active_neighbors_for_remote_element_key(stk::mesh::EntityKey
                                                                           remoteKey)
{
    int numActive = 0;
    for (stk::mesh::Entity localElem :
            m_remoteElementKeyToLocalNeighborElements[remoteKey])
        if (is_element_active(localElem))
            numActive++;
    return numActive;
}
void NoGhostGameofLife::pack_num_active_neighbors_into_buffer(stk::CommSparse& buffer,
                                                              int numActive, stk::mesh::EntityKey
                                                              remoteElemKey)
{
    int procNum = m_remoteElementKeyToOwningProcessor[remoteElemKey];
    buffer.send_buffer(procNum).pack<stk::mesh::EntityKey>(remoteElemKey);
    buffer.send_buffer(procNum).pack<int>(numActive);
}
void NoGhostGameofLife::recieve_num_active_neighbors_of_local_elements(stk::CommSparse& buffer)
{
    for (int procNum = 0; procNum < m_numProcs; procNum++)
    {
        stk::CommBuffer& buf = buffer.recv_buffer(procNum);
        while(buf.remaining())
            update_local_element_with_remote_neighbor_data(buf);
    }
}
void NoGhostGameofLife::update_local_element_with_remote_neighbor_data(stk::CommBuffer& buf)
{
    stk::mesh::EntityKey localElemKey = stk::unpack<stk::mesh::EntityKey>(buf);
    stk::mesh::Entity localElem = m_bulkData.get_entity(localElemKey);
    int neighborVal = stk::unpack<int>(buf);
    *stk::mesh::field_data(m_neighborField, localElem) += neighborVal;
}
void NoGhostGameofLife::update_element_membership()
{
    m_localActiveElements.clear(); // field
    for (stk::mesh::Entity localElem : m_localElementsToVisit)
    {
        switch (m_bulkData.bucket(localElem).topology())
        {
            case stk::topology::TRI_3_2D:
                update_tri_membership(localElem);
                break;
            case stk::topology::QUAD_4_2D:
                update_quad_membership(localElem);
                break;
            case stk::topology::HEX_8:
                update_hex_membership(localElem);
                break;
            default:
                ThrowRequire(false);
        }
        if (is_element_active(localElem)) // field
            m_localActiveElements.push_back(localElem);
    }
}
void NoGhostGameofLife::update_tri_membership(stk::mesh::Entity elem)
{
    switch (*stk::mesh::field_data(m_neighborField, elem))
    {
        case 2:
        case 7:
            break;
        case 3:
            activate_element(elem);
            break;
        default:
            deactivate_element(elem);
            break;
    }
}
void NoGhostGameofLife::update_quad_membership(stk::mesh::Entity elem)
{
//    maze rules
//    if (*stk::mesh::field_data(m_lifeField, elem))
//    {
//    switch (*stk::mesh::field_data(m_neighborField, elem))
//    {
//        case 3:
//            activate_element(elem);
//            break;
//        default:
//            deactivate_element(elem);
//    }
//    }
//    else
//    {
//    switch (*stk::mesh::field_data(m_neighborField, elem))
//    {
//        case 1:
//        case 2:
//        case 3:
//        case 4:
//        case 5:
//            activate_element(elem);
//            break;
//        default:
//            deactivate_element(elem);
//    }
//    }
    switch (*stk::mesh::field_data(m_neighborField, elem))
    {
        case 2:
            break;
        case 3:
            activate_element(elem);
            break;
        default:
            deactivate_element(elem);
            break;
    }
}
void NoGhostGameofLife::update_hex_membership(stk::mesh::Entity elem)
{
    switch (*stk::mesh::field_data(m_neighborField, elem))
    {
        case 4:
            break;
        case 5:
            activate_element(elem);
            break;
        default:
            deactivate_element(elem);
            break;
    }
}
