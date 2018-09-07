/*
 * GameofLife.cpp
 *
 *  Created on: Jul 25, 2015
 *      Author: Jonathan Chu
 */

#include "GameofLife.hpp"
#include "GameOfLife/GameofLifeMesh.hpp"  // for GameofLifeMesh
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/FieldParallel.hpp"  // for communicate_field_data
#include "stk_mesh/base/GetEntities.hpp"  // for get_selected_entities
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire
#include "stk_unit_test_utils/ioUtils.hpp"
/*
 * GameofLife
 */
//public
GameofLife::GameofLife(GameofLifeMesh& Mesh, std::string meshName)
:m_metaData(Mesh.meta_data()), m_bulkData(Mesh.bulk_data()),
  m_lifeField(Mesh.life_field()), m_stkIo(Mesh.comm()), m_meshName(meshName), m_time(0),
 m_elemType(Mesh.element_type()),m_activeNeighborField(Mesh.neighbor_field())
{
    get_elements();
    write_output_mesh();
    fill_neighbor_sets_of_elements();
}

void GameofLife::activate_these_ids(stk::mesh::EntityIdVector& elemIdsToActivate)
{
    activate_each_element_id_in_vector(elemIdsToActivate);
}

void GameofLife::run_game_of_life(int numSteps)
{
    if (0 == m_time)
        write_output_step();

    for (int time = 1; time <= numSteps; time++)
    {
        //std::cerr << time << std::endl;
        run_step_of_game_of_life();
    }
}

void GameofLife::put_all_nodes_in_nodeset()
{
    m_bulkData.modification_begin();

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(m_bulkData, stk::topology::NODE_RANK, nodes);
    stk::mesh::Part& nodeset1Part = m_metaData.declare_part("nodelist_1", stk::topology::NODE_RANK);
    for(stk::mesh::Entity node : nodes)
        m_bulkData.change_entity_parts(node, stk::mesh::ConstPartVector {&nodeset1Part});

    m_bulkData.modification_end();
}

void destroy_disabled_elements(stk::mesh::BulkData &bulk, const ScalarIntField &lifeField)
{
    bulk.modification_begin();

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elements);
    for(stk::mesh::Entity element : elements)
    {
        if ( *stk::mesh::field_data(lifeField, element) != 1)
        {
            stk::mesh::EntityVector nodes(bulk.begin_nodes(element), bulk.end_nodes(element));
            bulk.destroy_entity(element);
            for(unsigned j=0; j<nodes.size(); j++)
                bulk.destroy_entity(nodes[j]);
        }
    }

    bulk.modification_end();
}

void GameofLife::write_mesh()
{
    destroy_disabled_elements(m_bulkData, m_lifeField);

    put_all_nodes_in_nodeset();

    stk::io::StkMeshIoBroker stkIo(m_bulkData.parallel());
    stkIo.set_bulk_data(m_bulkData);
    size_t fh = stkIo.create_output_mesh("pic.g", stk::io::WRITE_RESULTS);
    stkIo.begin_output_step(fh, 0);
    stkIo.write_defined_output_fields(fh);
    stkIo.end_output_step(fh);
}

//protected
void GameofLife::update_this_element(stk::mesh::Entity elem)
{
    if (stk::topology::QUAD_4 == m_elemType)
        update_quad(elem);
    else if (stk::topology::HEX_8 == m_elemType)
        update_hex(elem);
    else if (stk::topology::TRIANGLE_3 == m_elemType)
        update_tri(elem);
    else
        ThrowRequire(true);
}

bool GameofLife::element_is_local(stk::mesh::Entity elem) const
{
    bool owned = false;

    if (m_bulkData.is_valid(elem) && m_bulkData.bucket(elem).owned())
        owned = true;

    return owned;
}

//test functions
unsigned GameofLife::get_num_active_elements() const
{
    unsigned count = 0;

    for (stk::mesh::Entity elem : m_elements)
        if (element_is_active(elem))
            count++;

    return count;
}

bool GameofLife::are_these_ids_active(const stk::mesh::EntityIdVector& elemIds) const
{
    bool accountedFor = true;
    for (stk::mesh::EntityId elemId : elemIds)
    {
        stk::mesh::Entity elem = m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
        if (element_is_local(elem) && !element_is_active(elem))
            accountedFor = false;
    }
    return accountedFor;
}

//private
void GameofLife::get_elements()
{
    stk::mesh::get_selected_entities(m_metaData.locally_owned_part(),
                                     m_bulkData.buckets(stk::topology::ELEM_RANK),
                                     m_elements);
    m_elemsOnProc = m_elements.size();
}
void GameofLife::write_output_mesh()
{
    m_meshName += ".e";
    m_stkIo.set_bulk_data(m_bulkData);
    m_fileHandler = m_stkIo.create_output_mesh(m_meshName, stk::io::WRITE_RESULTS);
    m_stkIo.add_field(m_fileHandler, m_lifeField);
    m_stkIo.write_output_mesh(m_fileHandler);
}

void GameofLife::fill_neighbor_sets_of_elements()
{
    for (stk::mesh::Entity elem : m_elements)
        fill_neighbor_set_of_this_element(elem);
}

void GameofLife::fill_neighbor_set_of_this_element(stk::mesh::Entity elem)
{
    unsigned numNodes = m_bulkData.num_nodes(elem);
    const stk::mesh::Entity* elemNodes = m_bulkData.begin_nodes(elem);

    for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
        add_this_nodes_elements_to_set(elemNodes[nodeIndex], elem);
}

void GameofLife::add_this_nodes_elements_to_set(stk::mesh::Entity node, stk::mesh::Entity elem)
{
    unsigned numElems = m_bulkData.num_elements(node);
    const stk::mesh::Entity* nodeElems = m_bulkData.begin_elements(node);

    for (unsigned nodeElemIndex = 0; nodeElemIndex < numElems; nodeElemIndex++)
        if (elem != nodeElems[nodeElemIndex])
            m_neighborSets[elem].insert(nodeElems[nodeElemIndex]);
}

void GameofLife::run_step_of_game_of_life()
{
    communicate_data();
    update_neighbor_fields();
    update_each_element();
    write_output_step();
}

void GameofLife::update_neighbor_fields()
{
    for (stk::mesh::Entity elem : m_elements)
        update_neighbor_field_of_this_element(elem);
}

void GameofLife::update_neighbor_field_of_this_element(stk::mesh::Entity elem)
{
    int* neighborVal = stk::mesh::field_data(m_activeNeighborField, elem);
    *neighborVal = get_num_active_neighbors_of_this_element(elem);
}

int GameofLife::get_num_active_neighbors_of_this_element(stk::mesh::Entity elem)
{
    int numActiveNeighbors = 0;

    for (stk::mesh::Entity neighborElem : m_neighborSets[elem])
        if (element_is_active(neighborElem))
            numActiveNeighbors++;

    return numActiveNeighbors;
}

void GameofLife::update_tri(stk::mesh::Entity elem)
{
    switch (*stk::mesh::field_data(m_activeNeighborField, elem))
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

void GameofLife::update_quad(stk::mesh::Entity elem)
{
    switch (*stk::mesh::field_data(m_activeNeighborField, elem))
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

void GameofLife::update_hex(stk::mesh::Entity elem)
{
    switch (*stk::mesh::field_data(m_activeNeighborField, elem))
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

void GameofLife::write_output_step()
{
    m_stkIo.begin_output_step(m_fileHandler, m_time);
    m_stkIo.write_defined_output_fields(m_fileHandler);
    m_stkIo.end_output_step(m_fileHandler);
    m_time++;
}

/*
 * PartGameofLife
 */
//public
PartGameofLife::PartGameofLife(GameofLifeMesh& Mesh, std::string meshName)
:GameofLife(Mesh, meshName), m_activePart(Mesh.active_part())
{
    m_active.push_back(&m_activePart);
}

//private
void PartGameofLife::activate_each_element_id_in_vector(stk::mesh::EntityIdVector& elemIdsToActivate)
{
    m_bulkData.modification_begin();

    for (stk::mesh::EntityId elemId : elemIdsToActivate)
        activate_this_element_id(elemId);

    m_bulkData.modification_end();
}

void PartGameofLife::activate_this_element_id(stk::mesh::EntityId elemId)
{
    stk::mesh::Entity elem = m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

    if (element_is_local(elem))
        activate_element(elem);
}

void PartGameofLife::communicate_data()
{
    // nothing here
}

bool PartGameofLife::element_is_active(stk::mesh::Entity elem) const
{
    return m_bulkData.bucket(elem).member(m_activePart);
}

void PartGameofLife::update_each_element()
{
    m_bulkData.modification_begin();

    for (stk::mesh::Entity elem : m_elements)
        update_this_element(elem);

    m_bulkData.modification_end();
}

void PartGameofLife::activate_element(stk::mesh::Entity elem)
{
    *stk::mesh::field_data(m_lifeField, elem) = 1;
    m_bulkData.change_entity_parts(elem, m_active, m_empty);
}

void PartGameofLife::deactivate_element(stk::mesh::Entity elem)
{
    *stk::mesh::field_data(m_lifeField, elem) = 0;
    m_bulkData.change_entity_parts(elem, m_empty, m_active);
}

/*
 * FieldGameofLife
 */
//public
FieldGameofLife::FieldGameofLife(GameofLifeMesh& Mesh, std::string meshName)
:GameofLife(Mesh, meshName)
{
}

//private
void FieldGameofLife::activate_each_element_id_in_vector(stk::mesh::EntityIdVector& elemIdsToActivate)
{
    for (stk::mesh::EntityId elemId : elemIdsToActivate)
        activate_this_element_id(elemId);
}

void FieldGameofLife::activate_this_element_id(stk::mesh::EntityId elemId)
{
    stk::mesh::Entity elem = m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

    if (element_is_local(elem))
        activate_element(elem);
}

void FieldGameofLife::communicate_data()
{
    communicate_field_data(m_bulkData, {&m_lifeField});
}

bool FieldGameofLife::element_is_active(stk::mesh::Entity elem) const
{
    return *stk::mesh::field_data(m_lifeField, elem);
}

void FieldGameofLife::update_each_element()
{
    for (stk::mesh::Entity elem : m_elements)
        update_this_element(elem);
}

void FieldGameofLife::activate_element(stk::mesh::Entity elem)
{
    *stk::mesh::field_data(m_lifeField, elem) = 1;
}

void FieldGameofLife::deactivate_element(stk::mesh::Entity elem)
{
    *stk::mesh::field_data(m_lifeField, elem) = 0;
}

void TodstwdGameOfLife::write_mesh()
{
    destroy_disabled_elements(m_bulkData, m_lifeField);

    stk::io::StkMeshIoBroker stkIo(m_bulkData.parallel());
    stkIo.set_bulk_data(m_bulkData);
    size_t fh = stkIo.create_output_mesh("pic.g", stk::io::WRITE_RESULTS);
    stkIo.begin_output_step(fh, 0);
    stkIo.write_defined_output_fields(fh);
    stkIo.end_output_step(fh);
}

