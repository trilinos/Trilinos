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

#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityIdVector
#include <stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/PyramidFixture.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_io/IossBridge.hpp"

namespace stk {
namespace mesh {
namespace fixtures {

PyramidFixture::PyramidFixture(MetaData& meta,
                               BulkData& bulk,
                               size_t nx,
                               size_t ny,
                               size_t nz,
                               size_t nid_start,
                               size_t eid_start)
  : m_spatial_dimension(3),
    m_nx(2*nx),
    m_ny(2*ny),
    m_nz(2*nz),
    node_id_start(nid_start),
    elem_id_start(eid_start),
    m_meta( meta),
    m_bulk_data( bulk ),
    m_elem_parts( ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    owns_mesh(false)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);
}

PyramidFixture::PyramidFixture(stk::ParallelMachine pm,
                               size_t nx,
                               size_t ny,
                               size_t nz,
                               stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : m_spatial_dimension(3),
    m_nx(2*nx),
    m_ny(2*ny),
    m_nz(2*nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(autoAuraOption)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("pyramid_part", stk::topology::PYRAMID_5) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);

}

PyramidFixture::PyramidFixture(stk::ParallelMachine pm,
                               size_t nx,
                               size_t ny,
                               size_t nz,
                               std::string coordinate_name,
                               stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : m_spatial_dimension(3),
    m_nx(2*nx),
    m_ny(2*ny),
    m_nz(2*nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(autoAuraOption)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("pyramid_part", stk::topology::PYRAMID_5) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordinate_name);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);

}

PyramidFixture::~PyramidFixture()
{

}

void PyramidFixture::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<size_t> hex_range_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_elems = m_nx * m_ny * m_nz / 8;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const size_t beg_elem = ( num_elems * p_rank ) / p_size ;
  const size_t end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( size_t i = beg_elem; i != end_elem; ++i) {
    hex_range_on_this_processor.push_back(i);
  }

  generate_mesh(hex_range_on_this_processor, coordMap);
}

void PyramidFixture::node_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  entity_id -= node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id % (m_ny+1);
  entity_id /= (m_ny+1);

  z = entity_id;
}

void PyramidFixture::hex_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  x = 2*(entity_id % (m_nx/2));
  entity_id /= (m_nx/2);

  y = 2*(entity_id % (m_ny/2));
  entity_id /= (m_ny/2);

  z = 2*entity_id;
}

void PyramidFixture::generate_mesh(std::vector<size_t> & hex_range_on_this_processor, const CoordinateMapping & coordMap)
{
  m_bulk_data.modification_begin();

  {
    int pyramid_vert[][5] = {
        {0,1,2,3,8},
        {4,7,6,5,8},
        {2,6,7,3,8},
        {4,5,1,0,8},
        {5,6,2,1,8},
        {7,4,0,3,8}};

    // Declare the elements that belong on this process
    std::vector<size_t>::iterator ib = hex_range_on_this_processor.begin();
    const std::vector<size_t>::iterator ie = hex_range_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(9);
    stk::mesh::EntityIdVector pyramid_nodes(5);

    for (; ib != ie; ++ib) {
      size_t hex_id = *ib;
      size_t ix = 0, iy = 0, iz = 0;
      hex_x_y_z(hex_id, ix, iy, iz);

      elem_nodes[0] = node_id( ix   , iy   , iz   );
      elem_nodes[1] = node_id( ix+2 , iy   , iz   );
      elem_nodes[2] = node_id( ix+2 , iy+2 , iz   );
      elem_nodes[3] = node_id( ix   , iy+2 , iz   );
      elem_nodes[4] = node_id( ix   , iy   , iz+2 );
      elem_nodes[5] = node_id( ix+2 , iy   , iz+2 );
      elem_nodes[6] = node_id( ix+2 , iy+2 , iz+2 );
      elem_nodes[7] = node_id( ix   , iy+2 , iz+2 );
      elem_nodes[8] = node_id( ix+1 , iy+1 , iz+1 );

      for (size_t pyr = 0; pyr < 6; ++pyr) {
        pyramid_nodes[0] = elem_nodes[pyramid_vert[pyr][0]];
        pyramid_nodes[1] = elem_nodes[pyramid_vert[pyr][1]];
        pyramid_nodes[2] = elem_nodes[pyramid_vert[pyr][2]];
        pyramid_nodes[3] = elem_nodes[pyramid_vert[pyr][3]];
        pyramid_nodes[4] = elem_nodes[pyramid_vert[pyr][4]];

        EntityId pyramid_id = 6*hex_id + pyr + elem_id_start;
        stk::mesh::declare_element( m_bulk_data, m_elem_parts, pyramid_id, pyramid_nodes);

        for (size_t i = 0; i<5; ++i) {
           stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , pyramid_nodes[i] );
           m_bulk_data.change_entity_parts(node, m_node_parts);

           STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
                              "This process should know about the nodes that make up its element");

           DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, pyramid_nodes[i], node);

           // Compute and assign coordinates to the node
           size_t nx = 0, ny = 0, nz = 0;
           node_x_y_z(pyramid_nodes[i], nx, ny, nz);

           Scalar * data = stk::mesh::field_data( *m_coord_field , node );

           coordMap.getNodeCoordinates(data, nx, ny, nz);

           // Since the mesh was expanded by 2x during creation, scale the coordinates accordingly
           data[0] *= 0.5;
           data[1] *= 0.5;
           data[2] *= 0.5;
         }
      }
    }
  }
  m_bulk_data.modification_end();
}

void PyramidFixture::fill_node_map(int p_rank)
{

  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny * m_nz / 8;

  const EntityId beg_elem = ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  {
    int pyramid_vert[][5] = {
        {0,1,2,3,8},
        {4,7,6,5,8},
        {2,6,7,3,8},
        {4,5,1,0,8},
        {5,6,2,1,8},
        {7,4,0,3,8}};

    // Declare the elements that belong on this process
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      size_t hex_id = *ib;
      size_t ix = 0, iy = 0, iz = 0;
      hex_x_y_z(hex_id, ix, iy, iz);

      stk::mesh::EntityId elem_nodes[9] ;
      stk::mesh::EntityId pyramid_nodes[5];

      elem_nodes[0] = node_id( ix   , iy   , iz   );
      elem_nodes[1] = node_id( ix+2 , iy   , iz   );
      elem_nodes[2] = node_id( ix+2 , iy+2 , iz   );
      elem_nodes[3] = node_id( ix   , iy+2 , iz   );
      elem_nodes[4] = node_id( ix   , iy   , iz+2 );
      elem_nodes[5] = node_id( ix+2 , iy   , iz+2 );
      elem_nodes[6] = node_id( ix+2 , iy+2 , iz+2 );
      elem_nodes[7] = node_id( ix   , iy+2 , iz+2 );
      elem_nodes[8] = node_id( ix+1 , iy+1 , iz+1 );

      for (size_t pyr = 0; pyr < 6; ++pyr) {
        pyramid_nodes[0] = elem_nodes[pyramid_vert[pyr][0]];
        pyramid_nodes[1] = elem_nodes[pyramid_vert[pyr][1]];
        pyramid_nodes[2] = elem_nodes[pyramid_vert[pyr][2]];
        pyramid_nodes[3] = elem_nodes[pyramid_vert[pyr][3]];
        pyramid_nodes[4] = elem_nodes[pyramid_vert[pyr][4]];

        for (size_t i = 0; i<5; ++i) {
          AddToNodeProcsMMap(m_nodes_to_procs, pyramid_nodes[i] , p_rank);
        }
      }
    }
  }
}

namespace simple_fields {

PyramidFixture::PyramidFixture(MetaData& meta,
                               BulkData& bulk,
                               size_t nx,
                               size_t ny,
                               size_t nz,
                               size_t nid_start,
                               size_t eid_start)
  : m_spatial_dimension(3),
    m_nx(2*nx),
    m_ny(2*ny),
    m_nz(2*nz),
    node_id_start(nid_start),
    elem_id_start(eid_start),
    m_meta( meta),
    m_bulk_data( bulk ),
    m_elem_parts( ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    owns_mesh(false)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);
}

PyramidFixture::PyramidFixture(stk::ParallelMachine pm,
                               size_t nx,
                               size_t ny,
                               size_t nz,
                               stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : m_spatial_dimension(3),
    m_nx(2*nx),
    m_ny(2*ny),
    m_nz(2*nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(autoAuraOption)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("pyramid_part", stk::topology::PYRAMID_5) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);

}

PyramidFixture::PyramidFixture(stk::ParallelMachine pm,
                               size_t nx,
                               size_t ny,
                               size_t nz,
                               std::string coordinate_name,
                               stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : m_spatial_dimension(3),
    m_nx(2*nx),
    m_ny(2*ny),
    m_nz(2*nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(autoAuraOption)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("pyramid_part", stk::topology::PYRAMID_5) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordinate_name);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);

}

PyramidFixture::~PyramidFixture()
{

}

void PyramidFixture::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<size_t> hex_range_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_elems = m_nx * m_ny * m_nz / 8;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const size_t beg_elem = ( num_elems * p_rank ) / p_size ;
  const size_t end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( size_t i = beg_elem; i != end_elem; ++i) {
    hex_range_on_this_processor.push_back(i);
  }

  generate_mesh(hex_range_on_this_processor, coordMap);
}

void PyramidFixture::node_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  entity_id -= node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id % (m_ny+1);
  entity_id /= (m_ny+1);

  z = entity_id;
}

void PyramidFixture::hex_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  x = 2*(entity_id % (m_nx/2));
  entity_id /= (m_nx/2);

  y = 2*(entity_id % (m_ny/2));
  entity_id /= (m_ny/2);

  z = 2*entity_id;
}

void PyramidFixture::generate_mesh(std::vector<size_t> & hex_range_on_this_processor, const CoordinateMapping & coordMap)
{
  m_bulk_data.modification_begin();

  {
    int pyramid_vert[][5] = {
        {0,1,2,3,8},
        {4,7,6,5,8},
        {2,6,7,3,8},
        {4,5,1,0,8},
        {5,6,2,1,8},
        {7,4,0,3,8}};

    // Declare the elements that belong on this process
    std::vector<size_t>::iterator ib = hex_range_on_this_processor.begin();
    const std::vector<size_t>::iterator ie = hex_range_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(9);
    stk::mesh::EntityIdVector pyramid_nodes(5);

    for (; ib != ie; ++ib) {
      size_t hex_id = *ib;
      size_t ix = 0, iy = 0, iz = 0;
      hex_x_y_z(hex_id, ix, iy, iz);

      elem_nodes[0] = node_id( ix   , iy   , iz   );
      elem_nodes[1] = node_id( ix+2 , iy   , iz   );
      elem_nodes[2] = node_id( ix+2 , iy+2 , iz   );
      elem_nodes[3] = node_id( ix   , iy+2 , iz   );
      elem_nodes[4] = node_id( ix   , iy   , iz+2 );
      elem_nodes[5] = node_id( ix+2 , iy   , iz+2 );
      elem_nodes[6] = node_id( ix+2 , iy+2 , iz+2 );
      elem_nodes[7] = node_id( ix   , iy+2 , iz+2 );
      elem_nodes[8] = node_id( ix+1 , iy+1 , iz+1 );

      for (size_t pyr = 0; pyr < 6; ++pyr) {
        pyramid_nodes[0] = elem_nodes[pyramid_vert[pyr][0]];
        pyramid_nodes[1] = elem_nodes[pyramid_vert[pyr][1]];
        pyramid_nodes[2] = elem_nodes[pyramid_vert[pyr][2]];
        pyramid_nodes[3] = elem_nodes[pyramid_vert[pyr][3]];
        pyramid_nodes[4] = elem_nodes[pyramid_vert[pyr][4]];

        EntityId pyramid_id = 6*hex_id + pyr + elem_id_start;
        stk::mesh::declare_element( m_bulk_data, m_elem_parts, pyramid_id, pyramid_nodes);

        for (size_t i = 0; i<5; ++i) {
           stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , pyramid_nodes[i] );
           m_bulk_data.change_entity_parts(node, m_node_parts);

           STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
                              "This process should know about the nodes that make up its element");

           stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, pyramid_nodes[i], node);

           // Compute and assign coordinates to the node
           size_t nx = 0, ny = 0, nz = 0;
           node_x_y_z(pyramid_nodes[i], nx, ny, nz);

           Scalar * data = stk::mesh::field_data( *m_coord_field , node );

           coordMap.getNodeCoordinates(data, nx, ny, nz);

           // Since the mesh was expanded by 2x during creation, scale the coordinates accordingly
           data[0] *= 0.5;
           data[1] *= 0.5;
           data[2] *= 0.5;
         }
      }
    }
  }
  m_bulk_data.modification_end();
}

void PyramidFixture::fill_node_map(int p_rank)
{

  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny * m_nz / 8;

  const EntityId beg_elem = ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  {
    int pyramid_vert[][5] = {
        {0,1,2,3,8},
        {4,7,6,5,8},
        {2,6,7,3,8},
        {4,5,1,0,8},
        {5,6,2,1,8},
        {7,4,0,3,8}};

    // Declare the elements that belong on this process
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      size_t hex_id = *ib;
      size_t ix = 0, iy = 0, iz = 0;
      hex_x_y_z(hex_id, ix, iy, iz);

      stk::mesh::EntityId elem_nodes[9] ;
      stk::mesh::EntityId pyramid_nodes[5];

      elem_nodes[0] = node_id( ix   , iy   , iz   );
      elem_nodes[1] = node_id( ix+2 , iy   , iz   );
      elem_nodes[2] = node_id( ix+2 , iy+2 , iz   );
      elem_nodes[3] = node_id( ix   , iy+2 , iz   );
      elem_nodes[4] = node_id( ix   , iy   , iz+2 );
      elem_nodes[5] = node_id( ix+2 , iy   , iz+2 );
      elem_nodes[6] = node_id( ix+2 , iy+2 , iz+2 );
      elem_nodes[7] = node_id( ix   , iy+2 , iz+2 );
      elem_nodes[8] = node_id( ix+1 , iy+1 , iz+1 );

      for (size_t pyr = 0; pyr < 6; ++pyr) {
        pyramid_nodes[0] = elem_nodes[pyramid_vert[pyr][0]];
        pyramid_nodes[1] = elem_nodes[pyramid_vert[pyr][1]];
        pyramid_nodes[2] = elem_nodes[pyramid_vert[pyr][2]];
        pyramid_nodes[3] = elem_nodes[pyramid_vert[pyr][3]];
        pyramid_nodes[4] = elem_nodes[pyramid_vert[pyr][4]];

        for (size_t i = 0; i<5; ++i) {
          stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, pyramid_nodes[i] , p_rank);
        }
      }
    }
  }
}

}

} // fixtures
} // mesh
} // stk
