/*
 * MeshBuilder.cpp
 *
 *  Created on: Jul 10, 2015
 *      Author: Jonathan Chu
 */
#include "GameofLifeMesh.hpp"
#include <stddef.h>                     // for size_t
#include "stk_io/IossBridge.hpp"        // for put_io_part_attribute
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FEMHelpers.hpp"  // for declare_element
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/GetEntities.hpp"  // for get_entities
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_topology/topology.hpp"    // for topology, etc

// Game of Life Mesh
GameofLifeMesh::GameofLifeMesh(stk::ParallelMachine comm, stk::topology elemType, unsigned spacialDim,
                         stk::mesh::BulkData::AutomaticAuraOption auraOption)
:m_numProcs(stk::parallel_machine_size(comm)), m_procRank(stk::parallel_machine_rank(comm)),
 m_metaData(spacialDim), m_bulkData(m_metaData, comm, auraOption),
  m_comm(comm), m_elemType(elemType)
{
    declare_parts();
    declare_fields();
    put_fields_on_parts();
    put_parts_in_io();
}

// public
void GameofLifeMesh::fill_mesh()
{
    declare_element_nodes_ids();
    declare_element_ids();
    create_entities();
    number_coordinate_field();
}
//private
void GameofLifeMesh::declare_fields()
{
    m_lifeField =
            &m_metaData.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK, "Life Field");
    m_activeNeighborField =
            &m_metaData.declare_field<ScalarIntField>(stk::topology::ELEMENT_RANK,
                                                      "Neighbor Field");

    m_distFactors = &m_metaData.declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "distribution_factors");
}
void GameofLifeMesh::put_fields_on_parts()
{
    int val = 0;
    stk::mesh::put_field_on_mesh(*m_lifeField, m_metaData.universal_part(), &val);
    stk::mesh::put_field_on_mesh(*m_activeNeighborField, m_metaData.universal_part(), &val);

    double factor = 1.0;
    stk::mesh::put_field_on_mesh(*m_distFactors, *m_nodeset1, &factor);
    stk::mesh::put_field_on_mesh(*m_distFactors, *m_nodeset2, &factor);
    stk::mesh::put_field_on_mesh(*m_distFactors, *m_nodeset3, &factor);
    stk::mesh::put_field_on_mesh(*m_distFactors, *m_nodeset4, &factor);
}
void GameofLifeMesh::declare_parts()
{
    m_elemPart = &m_metaData.declare_part_with_topology("Elem_Part", m_elemType);
    stk::mesh::Part &block2 = m_metaData.declare_part_with_topology("block_2", m_elemType);
    stk::mesh::Part &block3 = m_metaData.declare_part_with_topology("block_3", m_elemType);
    stk::mesh::Part &block4 = m_metaData.declare_part_with_topology("block_4", m_elemType);

    stk::io::put_io_part_attribute(block2);
    stk::io::put_io_part_attribute(block3);
    stk::io::put_io_part_attribute(block4);

    m_activePart = &m_metaData.declare_part_with_topology("Active Part", m_elemType);

    m_nodeset1 = &m_metaData.declare_part("nodelist_1", stk::topology::NODE_RANK);
    m_nodeset2 = &m_metaData.declare_part("nodelist_2", stk::topology::NODE_RANK);
    m_nodeset3 = &m_metaData.declare_part("nodelist_3", stk::topology::NODE_RANK);
    m_nodeset4 = &m_metaData.declare_part("nodelist_4", stk::topology::NODE_RANK);
}
void GameofLifeMesh::put_parts_in_io()
{
    stk::io::put_io_part_attribute(*m_elemPart);

    stk::io::put_io_part_attribute(*m_nodeset1);
    stk::io::put_io_part_attribute(*m_nodeset2);
    stk::io::put_io_part_attribute(*m_nodeset3);
    stk::io::put_io_part_attribute(*m_nodeset4);
}
void GameofLifeMesh::declare_element_ids()
{
    m_elemIds.resize(m_elemsOnProc);
    for (unsigned index = 0; index < m_elemsOnProc; index++)
        m_elemIds[index] = index + m_elemProcOffset + 1;
}
void GameofLifeMesh::create_entities()
{
    m_bulkData.modification_begin();
    declare_entities();
    if (!only_one_active_proc())
        share_nodes_between_processors();
    m_bulkData.modification_end();
}
void GameofLifeMesh::declare_entities()
{
    for (unsigned index = 0; index < m_elemsOnProc; index++)
        declare_element(m_bulkData, *m_elemPart, m_elemIds[index], m_elemNodeIds[index]);
    stk::mesh::get_entities(m_bulkData, stk::topology::NODE_RANK, m_nodes);
}

//*
//*
//*
//*
// Two Dimensional Game ofLife Mesh
TwoDimGameofLifeMesh::TwoDimGameofLifeMesh(stk::ParallelMachine comm,
                                                     stk::topology elemType,
                                                     unsigned width, unsigned height,
                                                     stk::mesh::BulkData::AutomaticAuraOption
                                                     auraOption)
:GameofLifeMesh(comm, elemType, 2, auraOption), m_width(width), m_height(height),
        m_rowsPerProc(height/m_numProcs), m_nodesPerRow(m_width+1)
{
    declare_coordinate_field();
    meta_data().commit();
}

// private
void TwoDimGameofLifeMesh::declare_coordinate_field()
{
    m_nodeCoords = &meta_data().declare_field<stk::mesh::Field<double,stk::mesh::Cartesian2d>>(
            stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(*m_nodeCoords, meta_data().universal_part(), 2, nullptr);
}
bool TwoDimGameofLifeMesh::only_one_active_proc()
{
    return m_height < m_numProcs;
}
void TwoDimGameofLifeMesh::share_nodes_between_processors()
{
    if (m_numProcs-1!= m_procRank)
        share_top_nodes();
    if (0 != m_procRank)
        share_bottom_nodes();
}
void TwoDimGameofLifeMesh::number_coordinate_field()
{
    for (size_t index = 0, numNodes = m_nodes.size(); index<numNodes; index++)
        number_coordinate_field_of_node(index);
}
void TwoDimGameofLifeMesh::number_coordinate_field_of_node(unsigned nodeIndex)
{
    double* const coord = stk::mesh::field_data(*m_nodeCoords, m_nodes[nodeIndex]);
    coord[0] = nodeIndex%m_nodesPerRow;
    coord[1] = nodeIndex/m_nodesPerRow + m_rowsPerProc*m_procRank;
}
void TwoDimGameofLifeMesh::share_top_nodes()
{
    unsigned topNodeOffset = m_rowsPerProc*m_nodesPerRow*(m_procRank+1);
    for (unsigned
            // needed to commit after fields and to begin/end modification calls
 index = 1; index <= m_nodesPerRow; index++)
        share_node_with_this_id_to_this_processor(index+topNodeOffset, m_procRank+1);
}
void TwoDimGameofLifeMesh::share_bottom_nodes()
{
    unsigned bottomNodeOffset = m_rowsPerProc*m_nodesPerRow*(m_procRank);
    for (unsigned index = 1; index <= m_nodesPerRow; index++)
        share_node_with_this_id_to_this_processor(index+bottomNodeOffset, m_procRank-1);
}
void TwoDimGameofLifeMesh::share_node_with_this_id_to_this_processor(unsigned nodeId,
                                                                          unsigned procNum)
{
    stk::mesh::Entity node = bulk_data().get_entity(stk::topology::NODE_RANK, nodeId);
    bulk_data().add_node_sharing(node, procNum);
}
//*
//*
//*
//*
// Triangle Game of Life Mesh
TriGameofLifeMesh::TriGameofLifeMesh(stk::ParallelMachine comm, unsigned width,
                                         unsigned height,
                                         stk::mesh::BulkData::AutomaticAuraOption auraOption)
:TwoDimGameofLifeMesh(comm, stk::topology::TRIANGLE_3, width, height, auraOption)
{
    m_elemsPerRow = 2*m_width;
    m_elemProcOffset = m_procRank*m_elemsPerRow*m_rowsPerProc;
    if (m_numProcs - 1 == m_procRank)
        m_elemsOnProc = m_elemsPerRow*(m_rowsPerProc+height%m_numProcs);
    else
        m_elemsOnProc = m_elemsPerRow*m_rowsPerProc;
    fill_mesh();
}

//private
void TriGameofLifeMesh::declare_element_nodes_ids()
{
    m_elemNodeIds.resize(m_elemsOnProc);
    unsigned initialOffset = m_nodesPerRow*m_rowsPerProc*m_procRank;
    for (unsigned index = 0; index < m_elemsOnProc; index+=2)
        declare_node_ids_of_two_elements(index, initialOffset);
}
void TriGameofLifeMesh::declare_node_ids_of_two_elements(unsigned index,
                                                           unsigned initialOffset)
{
    if (index < m_elemsPerRow)
        declare_node_ids_of_two_first_row_elements(index, initialOffset);
    else
        declare_node_ids_of_two_sucessive_row_elements(index);
}
void TriGameofLifeMesh::declare_node_ids_of_two_first_row_elements(unsigned index,
                                                                     unsigned initialOffset)
{
    unsigned offset = initialOffset-index/2;
    m_elemNodeIds[index].resize(3);
    m_elemNodeIds[index+1].resize(3);
    m_elemNodeIds[index] = {index+1+offset, index+2+offset, m_nodesPerRow+1+index+offset};
    m_elemNodeIds[index+1] = {m_nodesPerRow+index+2+offset, m_nodesPerRow+index+1+offset,
                              index+2+offset};
}
void TriGameofLifeMesh::declare_node_ids_of_two_sucessive_row_elements(unsigned index)
{
    declare_node_ids_of_this_element(index);
    declare_node_ids_of_this_element(index+1);
}
void TriGameofLifeMesh::declare_node_ids_of_this_element(unsigned index)
{
    m_elemNodeIds[index].resize(3);
    for (size_t nodeIndex = 0; nodeIndex<3; nodeIndex++)
        m_elemNodeIds[index][nodeIndex] =
                m_elemNodeIds[index%m_elemsPerRow][nodeIndex]+index/m_elemsPerRow*m_nodesPerRow;
}
//*
//*
//*
//*
//Quad Game of Life Mesh
QuadGameofLifeMesh::QuadGameofLifeMesh(stk::ParallelMachine comm, unsigned width, unsigned height,
                                 stk::mesh::BulkData::AutomaticAuraOption auraOption)
:TwoDimGameofLifeMesh(comm, stk::topology::QUAD_4, width, height, auraOption)
{
    m_elemsPerRow = m_width;
    if (m_numProcs - 1 == m_procRank)
        m_elemsOnProc = m_elemsPerRow*(m_rowsPerProc+m_height%m_numProcs);
    else
        m_elemsOnProc = m_elemsPerRow*m_rowsPerProc;
    m_elemProcOffset = m_procRank*m_elemsPerRow*m_rowsPerProc;
    fill_mesh();
}

//private
void QuadGameofLifeMesh::declare_element_nodes_ids()
{
    unsigned initialOffset = m_nodesPerRow*m_procRank*m_rowsPerProc;
    m_elemNodeIds.resize(m_elemsOnProc);
    for (unsigned index = 0; index < m_elemsOnProc; index++)
        declare_node_ids_of_element(index, initialOffset);
}
void QuadGameofLifeMesh::declare_node_ids_of_element(unsigned index, unsigned initialOffset)
{
    if (index < m_elemsPerRow)
        declare_first_row_element_nodes(index, initialOffset+index);
    else
        declare_remaining_element_nodes(index);

}
void QuadGameofLifeMesh::declare_first_row_element_nodes(unsigned index, unsigned offset)
{
    m_elemNodeIds[index].resize(4);
    m_elemNodeIds[index] = {offset+1, offset+2,offset+m_nodesPerRow+2,offset+m_nodesPerRow+1};
}
void QuadGameofLifeMesh::declare_remaining_element_nodes(unsigned index)
{
    m_elemNodeIds[index].resize(4);
    for (unsigned nodeIndex = 0; nodeIndex < 4; nodeIndex++)
        m_elemNodeIds[index][nodeIndex] =
                m_elemNodeIds[index%m_elemsPerRow][nodeIndex]+index/m_elemsPerRow*m_nodesPerRow;
}
//*
//*
//*
//*
//Three Dimensional Game of Life Mesh
ThreeDimGameofLifeMesh::ThreeDimGameofLifeMesh(stk::ParallelMachine comm,
                                                         stk::topology elemType, unsigned width,
                                                         unsigned height, unsigned depth,
                                                         stk::mesh::BulkData::AutomaticAuraOption
                                                         auraOption)
:GameofLifeMesh(comm, elemType, 3, auraOption), m_width(width), m_height(height), m_depth(depth),
 m_slicesPerProc(depth/m_numProcs), m_nodeWidth(width+1), m_nodeHeight(height+1),
 m_nodesPerSlice(m_nodeWidth*m_nodeHeight)
{
    declare_coordinate_field();
    meta_data().commit();
}

//private
void ThreeDimGameofLifeMesh::declare_coordinate_field()
{
    m_nodeCoords = &meta_data().declare_field<stk::mesh::Field<double,stk::mesh::Cartesian>>(
            stk::topology::NODE_RANK, "coordinates");
    stk::mesh::put_field_on_mesh(*m_nodeCoords, meta_data().universal_part(), 3, nullptr);
}
bool ThreeDimGameofLifeMesh::only_one_active_proc()
{
    return m_depth < m_numProcs;
}
void ThreeDimGameofLifeMesh::share_nodes_between_processors()
{
    if (m_numProcs-1 != m_procRank)
        share_nodes_in_back();
    if (0 != m_procRank)
        share_nodes_in_front();
}
void ThreeDimGameofLifeMesh::share_nodes_in_back()
{
    unsigned nodeOffset = m_nodesPerSlice*m_slicesPerProc*(m_procRank+1);
    for (unsigned index = 1; index <= m_nodesPerSlice; index++)
        share_node_with_this_id_to_this_processor(index+nodeOffset, m_procRank+1);
}
void ThreeDimGameofLifeMesh::share_nodes_in_front()
{
    unsigned nodeOffset = m_nodesPerSlice*m_slicesPerProc*m_procRank;
    for (unsigned index = 1; index <= m_nodesPerSlice; index++)
        share_node_with_this_id_to_this_processor(index+nodeOffset, m_procRank-1);
}
void ThreeDimGameofLifeMesh::share_node_with_this_id_to_this_processor(unsigned nodeId,
                                                                            unsigned procNum)
{
    stk::mesh::Entity node = bulk_data().get_entity(stk::topology::NODE_RANK, nodeId);
    bulk_data().add_node_sharing(node, procNum);
}
void ThreeDimGameofLifeMesh::number_coordinate_field()
{
    for (unsigned index = 0, numNodes = m_nodes.size(); index < numNodes; index++)
        number_coordinate_field_for_node(index);
}
void ThreeDimGameofLifeMesh::number_coordinate_field_for_node(unsigned nodeIndex)
{
    double* const coord = stk::mesh::field_data(*m_nodeCoords, m_nodes[nodeIndex]);
    coord[0] = nodeIndex%m_nodesPerSlice%m_nodeWidth;
    coord[1] = nodeIndex%m_nodesPerSlice/m_nodeWidth;
    coord[2] = nodeIndex/m_nodesPerSlice + m_procRank*m_slicesPerProc;
}
//*
//*
//*
//*
//Hex Game of Life Mesh
HexGameofLifeMesh::HexGameofLifeMesh(stk::ParallelMachine comm, unsigned width, unsigned height,
                               unsigned depth, stk::mesh::BulkData::AutomaticAuraOption
                               auraOption)
:ThreeDimGameofLifeMesh(comm, stk::topology::HEX_8, width, height, depth, auraOption)
{
    m_elemsPerSlice = m_width*m_height;
    if (m_numProcs-1 == m_procRank)
        m_elemsOnProc = m_elemsPerSlice*(m_slicesPerProc+depth%m_numProcs);
    else
        m_elemsOnProc = m_elemsPerSlice*m_slicesPerProc;
    m_elemProcOffset = m_elemsPerSlice*m_slicesPerProc*m_procRank;
    fill_mesh();
}

//private
void HexGameofLifeMesh::declare_element_nodes_ids()
{
    m_elemNodeIds.resize(m_elemsOnProc);
    unsigned offset = m_nodesPerSlice*m_slicesPerProc*m_procRank;
    for (unsigned index = 0; index < m_elemsOnProc; index++)
        declare_node_ids_of_element(index, offset);
}
void HexGameofLifeMesh::declare_node_ids_of_element(unsigned index, unsigned offset)
{
    if (index < m_elemsPerSlice)
        declare_first_slice_element_node_ids(m_elemNodeIds[index], index, offset);
    else
        declare_remaining_element_node_ids(m_elemNodeIds[index],
                                           m_elemNodeIds[index%m_elemsPerSlice], index);
}
void HexGameofLifeMesh::declare_first_slice_element_node_ids(stk::mesh::EntityIdVector& V,
                                                          unsigned index,unsigned offset)
{
    V.resize(8);
    unsigned rowOffset = index/m_width;
    unsigned totalOffset = offset+rowOffset+index;
    V = {totalOffset+1+m_nodeWidth,
         totalOffset+2+m_nodeWidth,
         totalOffset+m_nodesPerSlice+2+m_nodeWidth,
         totalOffset+m_nodesPerSlice+1+m_nodeWidth,
         totalOffset+1,
         totalOffset+2,
         totalOffset+m_nodesPerSlice+2,
         totalOffset+m_nodesPerSlice+1};
}
void HexGameofLifeMesh::declare_remaining_element_node_ids(stk::mesh::EntityIdVector& newer,
                                                        stk::mesh::EntityIdVector& older,
                                                        unsigned index)
{
    newer.resize(8);
    for (unsigned nodeIndex = 0; nodeIndex < 8; nodeIndex++)
        newer[nodeIndex]= older[nodeIndex] + index/m_elemsPerSlice*m_nodesPerSlice;
}
