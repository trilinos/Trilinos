// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for sort, unique
#include <array>
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityIdVector
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine

namespace stk {
namespace mesh {
namespace fixtures {
using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;
using stk::mesh::MeshBuilder;

QuadFixture::QuadFixture( MetaData& meta, BulkData& bulk, size_t nx, size_t ny, size_t nz, size_t nid_start, size_t eid_start)
  : m_spatial_dimension(2),
    m_meta(meta),
    m_bulk_data(bulk),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    m_nx( nx ),
    m_ny( ny ),
    m_node_id_start(nid_start),
    m_elem_id_start(eid_start)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::string& coordsName,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordsName);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          bool auraOn )
  : m_spatial_dimension(2),
    m_bulk_p( build_mesh(2, pm, (auraOn ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA)) ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

void QuadFixture::node_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= m_node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id;
}

void QuadFixture::elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= m_elem_id_start;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}


void QuadFixture::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = m_nx * m_ny;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor, coordMap);
}

void QuadFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor,
    const CoordinateMapping & coordMap)
{
  {
    //sort and unique the input elements
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_ids_on_this_processor.erase(ib, ie);
  }

  m_bulk_data.modification_begin();

  {
    // Declare the elements that belong on this process

    std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(4) ;

    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0;
      elem_x_y(entity_id, ix, iy);

      elem_nodes[0] = node_id( ix   , iy );
      elem_nodes[1] = node_id( ix+1 , iy );
      elem_nodes[2] = node_id( ix+1 , iy+1 );
      elem_nodes[3] = node_id( ix   , iy+1 );

      stk::mesh::declare_element( m_bulk_data, m_elem_parts, elem_id( ix , iy ) , elem_nodes);

      for (unsigned i = 0; i<4; ++i) {
        stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , elem_nodes[i] );
        m_bulk_data.change_entity_parts(node, m_node_parts);

        DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, elem_nodes[i], node);

        STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

        // Compute and assign coordinates to the node
        unsigned nx = 0, ny = 0;
        node_x_y(elem_nodes[i], nx, ny);

        Scalar * data = stk::mesh::field_data( *m_coord_field , node );

        // The CoordinateMappings are used for 2D and 3D so make sure we give it enough space to write to.
        std::array<double, 3> temp;
        coordMap.getNodeCoordinates(temp.data(), nx, ny, 0);

        data[0] = temp[0];
        data[1] = temp[1] ;
      }
    }
  }

  m_bulk_data.modification_end();

  m_bulk_data.initialize_face_adjacent_element_graph();
}

void QuadFixture::fill_node_map(int p_rank)
{

  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny;

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  {

    std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0;
      elem_x_y(entity_id, ix, iy);

      stk::mesh::EntityId elem_nodes[4] ;

      elem_nodes[0] = node_id( ix   , iy );
      elem_nodes[1] = node_id( ix+1 , iy );
      elem_nodes[2] = node_id( ix+1 , iy+1 );
      elem_nodes[3] = node_id( ix   , iy+1 );

      for (unsigned i = 0; i<4; ++i) {
        AddToNodeProcsMMap(m_nodes_to_procs, elem_nodes[i] , p_rank);
      }
    }
  }
}

Quad9Fixture::Quad9Fixture( MetaData& meta, BulkData& bulk, size_t nx, size_t ny, size_t nz, size_t nid_start, size_t eid_start)
  : m_spatial_dimension(2),
    m_meta(meta),
    m_bulk_data(bulk),
    m_quad_part( m_meta.declare_part_with_topology("quad9_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    m_nx( nx ),
    m_ny( ny ),
    m_node_id_start(nid_start),
    m_elem_id_start(eid_start)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

Quad9Fixture::Quad9Fixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

Quad9Fixture::Quad9Fixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::string& coordsName,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordsName);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

Quad9Fixture::Quad9Fixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          bool auraOn )
  : m_spatial_dimension(2),
    m_bulk_p( build_mesh(2, pm, (auraOn ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA)) ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

void Quad9Fixture::elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= m_elem_id_start;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}

void Quad9Fixture::generate_mesh()
{
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = m_nx * m_ny;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor);
}

void Quad9Fixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor)
{
  {
    //sort and unique the input elements
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_ids_on_this_processor.erase(ib, ie);
  }

  m_bulk_data.modification_begin();

  // Declare the elements that belong on this process

  std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
  const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
  stk::mesh::EntityIdVector elemNodes ;

  for (; ib != ie; ++ib) {
    EntityId entity_id = *ib;
    elemNodes = node_ids( entity_id );
    std::vector<double> nodeCoords = node_coords(entity_id);

    stk::mesh::declare_element( m_bulk_data, m_elem_parts, entity_id , elemNodes);

    for (unsigned i = 0; i<elemNodes.size(); ++i) {
      stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , elemNodes[i] );
      m_bulk_data.change_entity_parts(node, m_node_parts);

      DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, elemNodes[i], node);

      STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

      Scalar * data = stk::mesh::field_data( *m_coord_field , node );

      data[0] = nodeCoords[2*i+0];
      data[1] = nodeCoords[2*i+1] ;
    }
  }

  m_bulk_data.modification_end();

  m_bulk_data.initialize_face_adjacent_element_graph();
}

void Quad9Fixture::fill_node_map(int p_rank)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny;

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
  const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
  for (; ib != ie; ++ib) {
    EntityId entity_id = *ib;
    stk::mesh::EntityIdVector elem_nodes = node_ids( entity_id );

    for (unsigned i = 0; i<elem_nodes.size(); ++i) {
      AddToNodeProcsMMap(m_nodes_to_procs, elem_nodes[i] , p_rank);
    }
  }
}

std::vector<double> Quad9Fixture::node_coords(EntityId elemId) const
{
  unsigned x = 0;
  unsigned y = 0;

  elem_x_y(elemId, x, y);

  std::vector<double> coordVec(18u);

  coordVec[ 0] = x;
  coordVec[ 1] = y;
  coordVec[ 2] = x+1;
  coordVec[ 3] = y;
  coordVec[ 4] = x+1;
  coordVec[ 5] = y+1;
  coordVec[ 6] = x;
  coordVec[ 7] = y+1;
  coordVec[ 8] = x+0.5;
  coordVec[ 9] = y;
  coordVec[10] = x+1;
  coordVec[11] = y+0.5;
  coordVec[12] = x+0.5;
  coordVec[13] = y+1;
  coordVec[14] = x;
  coordVec[15] = y+0.5;
  coordVec[16] = x+0.5;
  coordVec[17] = y+0.5;

  return coordVec;
}

EntityIdVector Quad9Fixture::node_ids( EntityId elemId ) const
{
  unsigned x = 0;
  unsigned y = 0;

  elem_x_y(elemId, x, y);

  unsigned nodeSouth = m_node_id_start + (m_nx+1)*(m_ny+1) + y*(3*m_nx+1) + x;
  unsigned nodeWest  = nodeSouth + m_nx + x;
  unsigned nodeEast  = nodeWest + 2;
  unsigned nodeNorth = nodeSouth + 3*m_nx + 1;
  unsigned nodeMid   = nodeWest + 1;

  unsigned nodeLowerLeft  = m_node_id_start + x + ( m_nx + 1 ) * y;
  unsigned nodeLowerRight = nodeLowerLeft + 1;
  unsigned nodeUpperLeft  = nodeLowerLeft + 1 + m_nx;
  unsigned nodeUpperRight = nodeUpperLeft + 1;

  EntityIdVector nodeVec(9u);

  nodeVec[0] = nodeLowerLeft;
  nodeVec[1] = nodeLowerRight;
  nodeVec[2] = nodeUpperRight;
  nodeVec[3] = nodeUpperLeft;
  nodeVec[4] = nodeSouth;
  nodeVec[5] = nodeEast;
  nodeVec[6] = nodeNorth;
  nodeVec[7] = nodeWest;
  nodeVec[8] = nodeMid;

  return nodeVec ;
}

namespace simple_fields {

QuadFixture::QuadFixture( MetaData& meta, BulkData& bulk, size_t nx, size_t ny, size_t nz, size_t nid_start, size_t eid_start)
  : m_spatial_dimension(2),
    m_meta(meta),
    m_bulk_data(bulk),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    m_nx( nx ),
    m_ny( ny ),
    m_node_id_start(nid_start),
    m_elem_id_start(eid_start)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::string& coordsName,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordsName);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          bool auraOn )
  : m_spatial_dimension(2),
    m_bulk_p( build_mesh(2, pm, (auraOn ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA)) ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

void QuadFixture::node_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= m_node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id;
}

void QuadFixture::elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= m_elem_id_start;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}


void QuadFixture::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = m_nx * m_ny;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor, coordMap);
}

void QuadFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor,
    const CoordinateMapping & coordMap)
{
  {
    //sort and unique the input elements
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_ids_on_this_processor.erase(ib, ie);
  }

  m_bulk_data.modification_begin();

  {
    // Declare the elements that belong on this process

    std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(4) ;

    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0;
      elem_x_y(entity_id, ix, iy);

      elem_nodes[0] = node_id( ix   , iy );
      elem_nodes[1] = node_id( ix+1 , iy );
      elem_nodes[2] = node_id( ix+1 , iy+1 );
      elem_nodes[3] = node_id( ix   , iy+1 );

      stk::mesh::declare_element( m_bulk_data, m_elem_parts, elem_id( ix , iy ) , elem_nodes);

      for (unsigned i = 0; i<4; ++i) {
        stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , elem_nodes[i] );
        m_bulk_data.change_entity_parts(node, m_node_parts);

        stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, elem_nodes[i], node);

        STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

        // Compute and assign coordinates to the node
        unsigned nx = 0, ny = 0;
        node_x_y(elem_nodes[i], nx, ny);

        Scalar * data = stk::mesh::field_data( *m_coord_field , node );

        // The CoordinateMappings are used for 2D and 3D so make sure we give it enough space to write to.
        std::array<double, 3> temp;
        coordMap.getNodeCoordinates(temp.data(), nx, ny, 0);

        data[0] = temp[0];
        data[1] = temp[1] ;
      }
    }
  }

  m_bulk_data.modification_end();

  m_bulk_data.initialize_face_adjacent_element_graph();
}

void QuadFixture::fill_node_map(int p_rank)
{

  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny;

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  {

    std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0;
      elem_x_y(entity_id, ix, iy);

      stk::mesh::EntityId elem_nodes[4] ;

      elem_nodes[0] = node_id( ix   , iy );
      elem_nodes[1] = node_id( ix+1 , iy );
      elem_nodes[2] = node_id( ix+1 , iy+1 );
      elem_nodes[3] = node_id( ix   , iy+1 );

      for (unsigned i = 0; i<4; ++i) {
        stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, elem_nodes[i] , p_rank);
      }
    }
  }
}

Quad9Fixture::Quad9Fixture( MetaData& meta, BulkData& bulk, size_t nx, size_t ny, size_t nz, size_t nid_start, size_t eid_start)
  : m_spatial_dimension(2),
    m_meta(meta),
    m_bulk_data(bulk),
    m_quad_part( m_meta.declare_part_with_topology("quad9_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    m_nx( nx ),
    m_ny( ny ),
    m_node_id_start(nid_start),
    m_elem_id_start(eid_start)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

Quad9Fixture::Quad9Fixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

Quad9Fixture::Quad9Fixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::string& coordsName,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(2).set_entity_rank_names(rank_names).create() ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordsName);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

Quad9Fixture::Quad9Fixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          bool auraOn )
  : m_spatial_dimension(2),
    m_bulk_p( build_mesh(2, pm, (auraOn ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA)) ),
    m_meta( m_bulk_p->mesh_meta_data() ),
    m_bulk_data(*m_bulk_p),
    m_quad_part( m_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_9_2D ) ),
    m_elem_parts(1, &m_quad_part),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_nx( nx ),
    m_ny( ny )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

void Quad9Fixture::elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= m_elem_id_start;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}

void Quad9Fixture::generate_mesh()
{
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = m_nx * m_ny;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor);
}

void Quad9Fixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor)
{
  {
    //sort and unique the input elements
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_ids_on_this_processor.erase(ib, ie);
  }

  m_bulk_data.modification_begin();

  // Declare the elements that belong on this process

  std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
  const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
  stk::mesh::EntityIdVector elemNodes ;

  for (; ib != ie; ++ib) {
    EntityId entity_id = *ib;
    elemNodes = node_ids( entity_id );
    std::vector<double> nodeCoords = node_coords(entity_id);

    stk::mesh::declare_element( m_bulk_data, m_elem_parts, entity_id , elemNodes);

    for (unsigned i = 0; i<elemNodes.size(); ++i) {
      stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , elemNodes[i] );
      m_bulk_data.change_entity_parts(node, m_node_parts);

      stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, elemNodes[i], node);

      STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

      Scalar * data = stk::mesh::field_data( *m_coord_field , node );

      data[0] = nodeCoords[2*i+0];
      data[1] = nodeCoords[2*i+1] ;
    }
  }

  m_bulk_data.modification_end();

  m_bulk_data.initialize_face_adjacent_element_graph();
}

void Quad9Fixture::fill_node_map(int p_rank)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny;

  const EntityId beg_elem = m_elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = m_elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
  const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
  for (; ib != ie; ++ib) {
    EntityId entity_id = *ib;
    stk::mesh::EntityIdVector elem_nodes = node_ids( entity_id );

    for (unsigned i = 0; i<elem_nodes.size(); ++i) {
      stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, elem_nodes[i] , p_rank);
    }
  }
}

std::vector<double> Quad9Fixture::node_coords(EntityId elemId) const
{
  unsigned x = 0;
  unsigned y = 0;

  elem_x_y(elemId, x, y);

  std::vector<double> coordVec(18u);

  coordVec[ 0] = x;
  coordVec[ 1] = y;
  coordVec[ 2] = x+1;
  coordVec[ 3] = y;
  coordVec[ 4] = x+1;
  coordVec[ 5] = y+1;
  coordVec[ 6] = x;
  coordVec[ 7] = y+1;
  coordVec[ 8] = x+0.5;
  coordVec[ 9] = y;
  coordVec[10] = x+1;
  coordVec[11] = y+0.5;
  coordVec[12] = x+0.5;
  coordVec[13] = y+1;
  coordVec[14] = x;
  coordVec[15] = y+0.5;
  coordVec[16] = x+0.5;
  coordVec[17] = y+0.5;

  return coordVec;
}

EntityIdVector Quad9Fixture::node_ids( EntityId elemId ) const
{
  unsigned x = 0;
  unsigned y = 0;

  elem_x_y(elemId, x, y);

  unsigned nodeSouth = m_node_id_start + (m_nx+1)*(m_ny+1) + y*(3*m_nx+1) + x;
  unsigned nodeWest  = nodeSouth + m_nx + x;
  unsigned nodeEast  = nodeWest + 2;
  unsigned nodeNorth = nodeSouth + 3*m_nx + 1;
  unsigned nodeMid   = nodeWest + 1;

  unsigned nodeLowerLeft  = m_node_id_start + x + ( m_nx + 1 ) * y;
  unsigned nodeLowerRight = nodeLowerLeft + 1;
  unsigned nodeUpperLeft  = nodeLowerLeft + 1 + m_nx;
  unsigned nodeUpperRight = nodeUpperLeft + 1;

  EntityIdVector nodeVec(9u);

  nodeVec[0] = nodeLowerLeft;
  nodeVec[1] = nodeLowerRight;
  nodeVec[2] = nodeUpperRight;
  nodeVec[3] = nodeUpperLeft;
  nodeVec[4] = nodeSouth;
  nodeVec[5] = nodeEast;
  nodeVec[6] = nodeNorth;
  nodeVec[7] = nodeWest;
  nodeVec[8] = nodeMid;

  return nodeVec ;
}
} // namespace simple_fields

} // fixtures
} // mesh
} // stk
