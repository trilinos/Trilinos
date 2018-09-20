/*
 * MeshBuilder.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: jonchu
 */
#include "MeshBuilder.hpp"
#include <stddef.h>                     // for size_t
#include <vector>                       // for allocator, vector
#include "stk_learning/MeshBuilder/CoordinateSets.hpp"
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/IossBridge.hpp"        // for put_io_part_attribute
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian2d, etc
#include "stk_mesh/base/FEMHelpers.hpp"  // for declare_element
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/GetEntities.hpp"  // for get_entities
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire

// change for meshing big things
const double zOffset = 1;

MeshBuilder::MeshBuilder(stk::ParallelMachine comm, std::string name, stk::topology::topology_t
                         elemType, int spacialDim)
:m_metaData(spacialDim), m_bulkData(m_metaData, comm, stk::mesh::BulkData::NO_AUTO_AURA),
 m_spacialDim(spacialDim), m_numProcs(stk::parallel_machine_size(comm)),
 m_procRank(stk::parallel_machine_rank(comm)), m_name(name), m_time(0)
{
    m_elemPart = &m_metaData.declare_part_with_topology("Elem Part", elemType);
    stk::io::put_io_part_attribute(m_metaData.universal_part());
    stk::io::put_io_part_attribute(*m_elemPart);
}
void MeshBuilder::commit_meta()
{
    m_metaData.commit();
}
void MeshBuilder::begin_modification()
{
    m_bulkData.modification_begin();
}
void MeshBuilder::end_modification()
{
    m_bulkData.modification_end();
}
void MeshBuilder::write_mesh()
{
    std::string timeStep = std::to_string(m_time);
    while (timeStep.size() < 8)
        timeStep = "0" + timeStep;

    stk::io::StkMeshIoBroker stkIo(m_bulkData.parallel());
    stkIo.set_bulk_data(m_bulkData);
    size_t fh = stkIo.create_output_mesh(m_name + timeStep + ".e", stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(fh);

    stkIo.begin_output_step(fh, m_time++);
    stkIo.end_output_step(fh);
}

//GoL!!!
void MeshBuilder::create_life_and_neighbor_fields(ScalarIntField*& lifeField,
                                                  ScalarIntField*& neighborField)
{
    ThrowRequire(!m_metaData.is_commit());

    lifeField = &m_metaData.declare_field<ScalarIntField>(
            stk::topology::ELEM_RANK, "lifeField");
    neighborField = &m_metaData.declare_field<ScalarIntField>(
            stk::topology::ELEM_RANK, "neighborField");
    int val = 0;
    stk::mesh::put_field_on_mesh(*lifeField, m_metaData.universal_part(), &val);
    stk::mesh::put_field_on_mesh(*neighborField, m_metaData.universal_part(), &val);
}

//test
unsigned MeshBuilder::num_elems() const
{
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(m_bulkData, stk::topology::ELEM_RANK, elements);
    return elements.size();
}

//protected
stk::mesh::Entity MeshBuilder::create_element(stk::mesh::EntityId elemId,
                                              const stk::mesh::EntityIdVector& nodeIds,
                                              int chosenProc)
{
    ThrowRequire(chosenProc < num_procs());
    stk::mesh::Entity elem = stk::mesh::Entity();

    if (elem_id_not_used(elemId))
    {
        if (m_procRank == chosenProc)
            elem = generate_element(elemId, nodeIds);
        else
            share_shared_nodes(nodeIds, chosenProc);
    }

    m_usedElemIds.insert(elemId);
    return elem;
}
void MeshBuilder::remove_element(stk::mesh::EntityId elemId)
{
    stk::mesh::Entity elem = m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
    if (m_bulkData.is_valid(elem))
        m_bulkData.destroy_entity(elem);
}

//private
void MeshBuilder::update_node_to_processor_map(const stk::mesh::EntityIdVector& nodeIds, int proc)
{
    for (const stk::mesh::EntityId nodeId : nodeIds)
        m_nodeIdToSharingProcs[nodeId].insert(proc);
}
stk::mesh::Entity MeshBuilder::generate_element(stk::mesh::EntityId elemId,
                                                const stk::mesh::EntityIdVector& nodeIds)
{
    stk::mesh::Entity elem = stk::mesh::Entity();
    elem = declare_element(m_bulkData, *m_elemPart, elemId, nodeIds);
    share_nodes(nodeIds);
    return elem;
}
bool MeshBuilder::elem_id_not_used(stk::mesh::EntityId elemId)
{
    std::unordered_set<stk::mesh::EntityId>::iterator idIter =
            m_usedElemIds.find(elemId);
    return (m_usedElemIds.end() == idIter);
}
void MeshBuilder::share_nodes(const stk::mesh::EntityIdVector& nodeIds)
{
    for (stk::mesh::EntityId nodeId : nodeIds)
    {
       stk::mesh::Entity node = m_bulkData.get_entity(stk::topology::NODE_RANK, nodeId);
       for (int procNum : m_nodeIdToSharingProcs[nodeId])
           m_bulkData.add_node_sharing(node, procNum);
    }
}
void MeshBuilder::share_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int chosenProc)
{
    update_node_to_processor_map(nodeIds, chosenProc);
    for (stk::mesh::EntityId nodeId : nodeIds)
    {
        stk::mesh::Entity node = m_bulkData.get_entity(stk::topology::NODE_RANK, nodeId);
        if (m_bulkData.is_valid(node))
            m_bulkData.add_node_sharing(node, chosenProc);
    }
}

//Quad mesh Builder
//public
QuadMeshBuilder::QuadMeshBuilder(stk::ParallelMachine comm, std::string name)
:MeshBuilder(comm, name, stk::topology::QUAD_4, 2)
{
    declare_coordinates();
}
void QuadMeshBuilder::create_element(unsigned xCoord, unsigned yCoord, int chosenProc)
{
    ElemCoordPair elemCoord(xCoord, yCoord);
    stk::mesh::Entity elem = MeshBuilder::create_element(elemCoord.elemId,
                                                         elemCoord.nodeIds, chosenProc);
    if (stk::mesh::Entity() != elem)
        label_node_coordinates(elemCoord);
}
void QuadMeshBuilder::fill_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper)
{
    std::vector<unsigned> layersPerProc;
    int layersCreated = yUpper-yLower + 1;
    unsigned basePerRow = (layersCreated)/num_procs();
    for (int procIndex = 0; procIndex < num_procs(); procIndex++)
        layersPerProc.push_back(basePerRow+(layersCreated%num_procs() > procIndex));

    unsigned bottomLayer = yLower;
    for (int procIndex = 0; procIndex < num_procs(); procIndex++)
    {
        unsigned topLayer = bottomLayer + layersPerProc[procIndex] - 1;
        fill_area_on_proc(xLower, xUpper, bottomLayer, topLayer, procIndex);
        bottomLayer = topLayer + 1;
    }
}
void QuadMeshBuilder::fill_area_randomly(unsigned xLower, unsigned xUpper,
                                         unsigned yLower, unsigned yUpper)
{
    int procCounter = 0;
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
           create_element(x, y, procCounter++%num_procs());
}
void QuadMeshBuilder::fill_area_on_proc(unsigned xLower, unsigned xUpper, unsigned yLower,
                                        unsigned yUpper, int chosenProc)
{
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
           create_element(x, y, chosenProc);
}
void QuadMeshBuilder::fill_area_with_layers(unsigned xLower, unsigned xUpper, unsigned yLower,
                                            unsigned yUpper)
{
    int rowCounter = 0;
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
            fill_area_on_proc(xLower, xUpper, y, y, rowCounter++%num_procs());
}
void QuadMeshBuilder::remove_element(unsigned xCoord, unsigned yCoord)
{
    stk::mesh::EntityId elemId = generate_two_dim_elem_id(xCoord, yCoord);
    MeshBuilder::remove_element(elemId);
}
void QuadMeshBuilder::remove_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper)
{
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
            remove_element(x, y);
}

//test
double QuadMeshBuilder::node_x_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[0];
}
double QuadMeshBuilder::node_y_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[1];
}

//private
void QuadMeshBuilder::declare_coordinates()
{
    m_coordinates = &meta_data().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian2d>>(
            stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(*m_coordinates, meta_data().universal_part(), 2, nullptr);
}
void QuadMeshBuilder::label_node_coordinates(const ElemCoordPair& elemCoords)
{
    for (unsigned nodeIndex = 0; nodeIndex < 4; nodeIndex++)
    {
       stk::mesh::Entity node = bulk_data().get_entity(stk::topology::NODE_RANK,
                                                      elemCoords.nodeIds[nodeIndex]);
       if (0 == nodeIndex)
          label_coordinates(node, elemCoords.x-1, elemCoords.y-1);
       else if (1 == nodeIndex)
          label_coordinates(node, elemCoords.x, elemCoords.y-1);
       else if (2 == nodeIndex)
          label_coordinates(node, elemCoords.x, elemCoords.y);
       else if (3 == nodeIndex)
          label_coordinates(node, elemCoords.x-1, elemCoords.y);
    }
}
void QuadMeshBuilder::label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY)
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    coord[0] = nodeX;
    coord[1] = nodeY;
}

//HexMeshBuilder
//public
HexMeshBuilder::HexMeshBuilder(stk::ParallelMachine comm, std::string name)
:MeshBuilder(comm, name, stk::topology::HEX_8, 3)
{
    declare_coordinates();
}
void HexMeshBuilder::create_element(unsigned xCoord, unsigned yCoord, unsigned zCoord, int chosenProc)
{

    ElemCoordTriple elemCoord(xCoord, yCoord, zCoord);
    stk::mesh::Entity elem = MeshBuilder::create_element(elemCoord.elemId,
                                                         elemCoord.nodeIds, chosenProc);
    if (stk::mesh::Entity() != elem)
        label_node_coordinates(elemCoord);
}
void HexMeshBuilder::fill_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                               unsigned zLower, unsigned zUpper)
{
    std::vector<unsigned> layersPerProc;
    int layersCreated = zUpper-zLower + 1;
    unsigned basePerRow = (layersCreated)/num_procs();
    for (int procIndex = 0; procIndex < num_procs(); procIndex++)
        layersPerProc.push_back(basePerRow+(layersCreated%num_procs() > procIndex));

    unsigned bottomLayer = yLower;
    for (int procIndex = 0; procIndex < num_procs(); procIndex++)
    {
        unsigned topLayer = bottomLayer + layersPerProc[procIndex] - 1;
        fill_area_on_proc(xLower, xUpper, yLower, yUpper, bottomLayer, topLayer, procIndex);
        bottomLayer = topLayer + 1;
    }
}
void HexMeshBuilder::fill_area_randomly(unsigned xLower, unsigned xUpper, unsigned yLower,
                                        unsigned yUpper, unsigned zLower, unsigned zUpper)
{
    int currentProc = 0;
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
            for (unsigned z = zLower; z <= zUpper; z++)
                create_element(x, y, z, currentProc++%num_procs());
}
void HexMeshBuilder::fill_area_on_proc(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned
                                       yUpper, unsigned zLower, unsigned zUpper, int chosenProc)
{
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
            for (unsigned z = zLower; z <= zUpper; z++)
                create_element(x, y, z, chosenProc);
}
void HexMeshBuilder::fill_area_with_layers(unsigned xLower, unsigned xUpper, unsigned yLower,
                                           unsigned yUpper, unsigned zLower, unsigned zUpper)
{
    int layerCounter = 0;
    for (unsigned z = zLower; z <= zUpper; z++)
        fill_area_on_proc(xLower, xUpper, yLower, yUpper, z, z, layerCounter++%num_procs());

}
void HexMeshBuilder::remove_element(unsigned xCoord, unsigned yCoord, unsigned zCoord)
{
    stk::mesh::EntityId elemId = generate_three_dim_elem_id(xCoord, yCoord, zCoord);
    MeshBuilder::remove_element(elemId);
}
void HexMeshBuilder::remove_area(unsigned xLower, unsigned xUpper, unsigned yLower, unsigned yUpper,
                                 unsigned zLower, unsigned zUpper)
{
    for (unsigned x = xLower; x <= xUpper; x++)
        for (unsigned y = yLower; y <= yUpper; y++)
            for (unsigned z = zLower; z <= zUpper; z++)
                remove_element(x, y, z);
}

//test function
double HexMeshBuilder::node_x_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[0];
}
double HexMeshBuilder::node_y_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[1];
}
double HexMeshBuilder::node_z_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[2];
}

//private
void HexMeshBuilder::declare_coordinates()
{
   m_coordinates = &meta_data().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d>>(
           stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field_on_mesh(*m_coordinates, meta_data().universal_part(), 3, nullptr);
}
void HexMeshBuilder::label_node_coordinates(const ElemCoordTriple& elemCoords)
{
    for (unsigned nodeIndex = 0; nodeIndex < 8; nodeIndex++)
    {
       stk::mesh::Entity node = bulk_data().get_entity(stk::topology::NODE_RANK,
                                                      elemCoords.nodeIds[nodeIndex]);
       switch (nodeIndex)
       {
           case 0:
              label_coordinates(node, elemCoords.x-1, elemCoords.y-1, elemCoords.z-1);
              break;
           case 1:
              label_coordinates(node, elemCoords.x, elemCoords.y-1, elemCoords.z-1);
              break;
           case 2:
              label_coordinates(node, elemCoords.x, elemCoords.y, elemCoords.z-1);
              break;
           case 3:
              label_coordinates(node, elemCoords.x-1, elemCoords.y, elemCoords.z-1);
              break;
           case 4:
              label_coordinates(node, elemCoords.x-1, elemCoords.y-1, elemCoords.z);
              break;
           case 5:
              label_coordinates(node, elemCoords.x, elemCoords.y-1, elemCoords.z);
              break;
           case 6:
              label_coordinates(node, elemCoords.x, elemCoords.y, elemCoords.z);
              break;
           case 7:
              label_coordinates(node, elemCoords.x-1, elemCoords.y, elemCoords.z);
              break;
       }
    }
}

void HexMeshBuilder::label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY,
                                       unsigned nodeZ)
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    coord[0] = nodeX;
    coord[1] = nodeY;
    coord[2] = nodeZ*zOffset;
}
