/*
 * MeshBuilder.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: jonchu
 */
#include "MeshPrinter.h"

MeshPrinter::MeshPrinter(stk::ParallelMachine comm, std::string name, stk::topology::topology_t
                         elemType, int spacialDim)
:m_metaData(spacialDim), m_bulkData(m_metaData, comm, stk::mesh::BulkData::NO_AUTO_AURA),
 m_numProcs(stk::parallel_machine_size(comm)), m_procRank(stk::parallel_machine_rank(comm)),
 m_name(name), m_time(0)
{
    m_elemPart = &m_metaData.declare_part_with_topology("Elem Part", elemType);
    stk::io::put_io_part_attribute(m_metaData.universal_part());
    stk::io::put_io_part_attribute(*m_elemPart);
}
void MeshPrinter::commit_meta()
{
    m_metaData.commit();
}
void MeshPrinter::begin_modification()
{
    m_bulkData.modification_begin();
}
void MeshPrinter::end_modification()
{
    m_bulkData.modification_end();
}
stk::mesh::Entity MeshPrinter::create_element(stk::mesh::EntityId elemId,
                                              const stk::mesh::EntityIdVector& nodeIds,
                                              int chosenProc)
{
    ThrowRequire(chosenProc < num_procs());
    stk::mesh::Entity elem = stk::mesh::Entity();

    if (m_procRank == chosenProc)
        elem = generate_element(elemId, nodeIds);
    else
        share_shared_nodes(nodeIds, chosenProc);

    m_usedElemIds.insert(elemId);
    return elem;
}

void MeshPrinter::write_mesh()
{
    std::string timeStep = std::to_string(m_time);
    while (timeStep.size() < 4)
        timeStep = "0" + timeStep;

    stk::io::StkMeshIoBroker stkIo(m_bulkData.parallel());
    stkIo.set_bulk_data(m_bulkData);
    size_t fh = stkIo.create_output_mesh(m_name + timeStep + ".e", stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(fh);

    stkIo.begin_output_step(fh, m_time++);
    stkIo.end_output_step(fh);
}

//private
void MeshPrinter::update_node_to_processor_map(const stk::mesh::EntityIdVector& nodeIds, int proc)
{
    for (const stk::mesh::EntityId nodeId : nodeIds)
        m_nodeIdToSharingProcs[nodeId].insert(proc);
}
stk::mesh::Entity MeshPrinter::generate_element(stk::mesh::EntityId elemId,
                                                const stk::mesh::EntityIdVector& nodeIds)
{
    stk::mesh::Entity elem = stk::mesh::Entity();
    if (elem_id_not_used(elemId))
    {
        elem = declare_element(m_bulkData, *m_elemPart, elemId, nodeIds);
        share_nodes(nodeIds);
    }
    return elem;
}
bool MeshPrinter::elem_id_not_used(stk::mesh::EntityId elemId)
{
    std::unordered_set<stk::mesh::EntityId>::iterator idIter =
            m_usedElemIds.find(elemId);
    return (m_usedElemIds.end() == idIter);
}
void MeshPrinter::share_nodes(const stk::mesh::EntityIdVector& nodeIds)
{
    for (stk::mesh::EntityId nodeId : nodeIds)
    {
       stk::mesh::Entity node = m_bulkData.get_entity(stk::topology::NODE_RANK, nodeId);
       for (int procNum : m_nodeIdToSharingProcs[nodeId])
           m_bulkData.add_node_sharing(node, procNum);
    }
}
void MeshPrinter::share_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int chosenProc)
{
    update_node_to_processor_map(nodeIds, chosenProc);
    for (stk::mesh::EntityId nodeId : nodeIds)
    {
        stk::mesh::Entity node = m_bulkData.get_entity(stk::topology::NODE_RANK, nodeId);
        if (m_bulkData.is_valid(node))
            m_bulkData.add_node_sharing(node, chosenProc);
    }
}

//Quad mesh printer
//public
QuadMeshPrinter::QuadMeshPrinter(stk::ParallelMachine comm, std::string name)
:MeshPrinter(comm, name, stk::topology::QUAD_4, 2)
{
    m_coordinates = &meta_data().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian2d>>(
            stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field(*m_coordinates, meta_data().universal_part(), 2);
}
void QuadMeshPrinter::create_element(unsigned xCoord, unsigned yCoord, int chosenProc)
{
    ElemCoordPair elemCoord(xCoord, yCoord);
    stk::mesh::Entity elem = MeshPrinter::create_element(elemCoord.elemId,
                                                         elemCoord.nodeIds, chosenProc);

    if (stk::mesh::Entity() != elem)
        label_node_coordinates(elemCoord);
}

//test
double QuadMeshPrinter::node_x_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[0];
}
double QuadMeshPrinter::node_y_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[1];
}

//private
void QuadMeshPrinter::label_node_coordinates(const ElemCoordPair& elemCoords)
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
void QuadMeshPrinter::label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY)
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    coord[0] = nodeX;
    coord[1] = nodeY;
}

//HexMeshPrinter
//public
HexMeshPrinter::HexMeshPrinter(stk::ParallelMachine comm, std::string name)
:MeshPrinter(comm, name, stk::topology::HEX_8, 3)
{
   m_coordinates = &meta_data().declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d>>(
           stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field(*m_coordinates, meta_data().universal_part(), 3);
}
void HexMeshPrinter::create_element(unsigned xCoord, unsigned yCoord, unsigned zCoord, int chosenProc)
{

    ElemCoordTriple elemCoord(xCoord, yCoord, zCoord);
    stk::mesh::Entity elem = MeshPrinter::create_element(elemCoord.elemId,
                                                         elemCoord.nodeIds, chosenProc);

    if (stk::mesh::Entity() != elem)
        label_node_coordinates(elemCoord);
}
//test function
double HexMeshPrinter::node_x_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[0];
}
double HexMeshPrinter::node_y_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[1];
}
double HexMeshPrinter::node_z_coord(stk::mesh::Entity node) const
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    return coord[2];
}
//private
void HexMeshPrinter::label_node_coordinates(const ElemCoordTriple& elemCoords)
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

void HexMeshPrinter::label_coordinates(stk::mesh::Entity node, unsigned nodeX, unsigned nodeY,
                                       unsigned nodeZ)
{
    double* const coord = stk::mesh::field_data(*m_coordinates, node);
    coord[0] = nodeX;
    coord[1] = nodeY;
    coord[2] = nodeZ;
}
