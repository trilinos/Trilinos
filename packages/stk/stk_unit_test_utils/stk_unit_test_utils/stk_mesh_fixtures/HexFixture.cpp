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

#include <algorithm>                    // for sort, unique
#include <set>                          // for set
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityIdVector
#include <stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg
#include <utility>                      // for pair
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine

namespace stk {
namespace mesh {
namespace fixtures {

HexFixture::HexFixture(MetaData& meta,
                       BulkData& bulk,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       size_t nid_start,
                       size_t eid_start)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    node_id_start(nid_start),
    elem_id_start(eid_start),
    m_meta( meta),
    m_bulk_data( bulk ),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    owns_mesh(false)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::HexFixture(stk::ParallelMachine pm,
                       size_t nx,
                       size_t ny,
                       size_t nz)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(stk::mesh::BulkData::AUTO_AURA)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::HexFixture(stk::ParallelMachine pm,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       const std::string& coordinate_name)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(stk::mesh::BulkData::AUTO_AURA)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordinate_name);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::HexFixture(stk::ParallelMachine pm,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       bool auraOn)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(auraOn ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::~HexFixture()
{

}

void HexFixture::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_elems = m_nx * m_ny * m_nz ;

  m_nodes_to_procs.clear();
  for (int rank = 0; rank < static_cast<int>(p_size); ++rank) {
    fill_node_map(rank);
  }

  const EntityId beg_elem = elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor, coordMap);
}

void HexFixture::node_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  entity_id -= node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id % (m_ny+1);
  entity_id /= (m_ny+1);

  z = entity_id;
}

void HexFixture::elem_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  entity_id -= elem_id_start;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id % m_ny;
  entity_id /= m_ny;

  z = entity_id;
}

void HexFixture::fill_node_map( int p_rank)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny * m_nz ;

  const EntityId beg_elem = elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  //sort and unique the input elements
  std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
  std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

  std::sort( ib, ie);
  ib = std::unique( ib, ie);
  element_ids_on_this_processor.erase(ib, ie);

  std::set<EntityId> nodes_on_proc;

  ib = element_ids_on_this_processor.begin();
  ie = element_ids_on_this_processor.end();
  for (; ib != ie; ++ib) {
    EntityId entity_id = *ib;
    size_t ix = 0, iy = 0, iz = 0;
    elem_x_y_z(entity_id, ix, iy, iz);

    stk::mesh::EntityId elem_node[8] ;

    elem_node[0] = node_id( ix   , iy   , iz   );
    elem_node[1] = node_id( ix+1 , iy   , iz   );
    elem_node[2] = node_id( ix+1 , iy+1 , iz   );
    elem_node[3] = node_id( ix   , iy+1 , iz   );
    elem_node[4] = node_id( ix   , iy   , iz+1 );
    elem_node[5] = node_id( ix+1 , iy   , iz+1 );
    elem_node[6] = node_id( ix+1 , iy+1 , iz+1 );
    elem_node[7] = node_id( ix   , iy+1 , iz+1 );

    for (int ien = 0; ien < 8; ++ien) {
      stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, elem_node[ien], p_rank);
    }
  }
}

void HexFixture::fill_node_map(const std::map<int,std::vector<EntityId> > &parallel_distribution)
{
  m_nodes_to_procs.clear();

  std::map<int,std::vector<EntityId> >::const_iterator pd_i, pd_e;
  pd_i = parallel_distribution.begin();
  pd_e = parallel_distribution.end();
  for (; pd_i != pd_e; ++pd_i) {
    const int proc = pd_i->first;
    const std::vector<EntityId> &elements = pd_i->second;

    for (size_t e_j = 0; e_j < elements.size(); ++ e_j) {

      EntityId element_id = elements[e_j];

      size_t ix = 0, iy = 0, iz = 0;
      elem_x_y_z(element_id, ix, iy, iz);

      stk::mesh::EntityId elem_node[8] ;

      elem_node[0] = node_id( ix   , iy   , iz   );
      elem_node[1] = node_id( ix+1 , iy   , iz   );
      elem_node[2] = node_id( ix+1 , iy+1 , iz   );
      elem_node[3] = node_id( ix   , iy+1 , iz   );
      elem_node[4] = node_id( ix   , iy   , iz+1 );
      elem_node[5] = node_id( ix+1 , iy   , iz+1 );
      elem_node[6] = node_id( ix+1 , iy+1 , iz+1 );
      elem_node[7] = node_id( ix   , iy+1 , iz+1 );

      for (int ien = 0; ien < 8; ++ien) {
        stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, elem_node[ien], proc);
      }
    }
  }
}

void HexFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap)
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

    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(8);
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      size_t ix = 0, iy = 0, iz = 0;
      elem_x_y_z(entity_id, ix, iy, iz);

      elem_nodes[0] = node_id( ix   , iy   , iz   );
      elem_nodes[1] = node_id( ix+1 , iy   , iz   );
      elem_nodes[2] = node_id( ix+1 , iy+1 , iz   );
      elem_nodes[3] = node_id( ix   , iy+1 , iz   );
      elem_nodes[4] = node_id( ix   , iy   , iz+1 );
      elem_nodes[5] = node_id( ix+1 , iy   , iz+1 );
      elem_nodes[6] = node_id( ix+1 , iy+1 , iz+1 );
      elem_nodes[7] = node_id( ix   , iy+1 , iz+1 );

      stk::mesh::declare_element( m_bulk_data, m_elem_parts, elem_id( ix , iy , iz ) , elem_nodes);

      for (size_t i = 0; i<8; ++i) {
        EntityId node_id = elem_nodes[i];
        stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id );
        m_bulk_data.change_entity_parts(node, m_node_parts);

        stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, node_id, node);

        STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

        // Compute and assign coordinates to the node
        size_t nx = 0, ny = 0, nz = 0;
        node_x_y_z(elem_nodes[i], nx, ny, nz);

        Scalar * data = stk::mesh::field_data( *m_coord_field , node );

        coordMap.getNodeCoordinates(data, nx, ny, nz);
      }
    }
  }

  m_bulk_data.modification_end();
}

namespace simple_fields {

HexFixture::HexFixture(MetaData& meta,
                       BulkData& bulk,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       size_t nid_start,
                       size_t eid_start)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    node_id_start(nid_start),
    elem_id_start(eid_start),
    m_meta( meta),
    m_bulk_data( bulk ),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    owns_mesh(false)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::HexFixture(stk::ParallelMachine pm,
                       size_t nx,
                       size_t ny,
                       size_t nz)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(stk::mesh::BulkData::AUTO_AURA)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::HexFixture(stk::ParallelMachine pm,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       const std::string& coordinate_name)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(stk::mesh::BulkData::AUTO_AURA)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordinate_name);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::HexFixture(stk::ParallelMachine pm,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       bool auraOn)
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_bulk_p( MeshBuilder(pm)
                .set_aura_option(auraOn ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA)
                .set_add_fmwk_data(false)
                .set_spatial_dimension(3)
                .create() ),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("hex_part", stk::topology::HEX_8) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) )
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates");

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
}

HexFixture::~HexFixture()
{

}

void HexFixture::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_elems = m_nx * m_ny * m_nz ;

  m_nodes_to_procs.clear();
  for (int rank = 0; rank < static_cast<int>(p_size); ++rank) {
    fill_node_map(rank);
  }

  const EntityId beg_elem = elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor, coordMap);
}

void HexFixture::node_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  entity_id -= node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id % (m_ny+1);
  entity_id /= (m_ny+1);

  z = entity_id;
}

void HexFixture::elem_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const
{
  entity_id -= elem_id_start;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id % m_ny;
  entity_id /= m_ny;

  z = entity_id;
}

void HexFixture::fill_node_map( int p_rank)
{
  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny * m_nz ;

  const EntityId beg_elem = elem_id_start + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = elem_id_start + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  //sort and unique the input elements
  std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
  std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

  std::sort( ib, ie);
  ib = std::unique( ib, ie);
  element_ids_on_this_processor.erase(ib, ie);

  std::set<EntityId> nodes_on_proc;

  ib = element_ids_on_this_processor.begin();
  ie = element_ids_on_this_processor.end();
  for (; ib != ie; ++ib) {
    EntityId entity_id = *ib;
    size_t ix = 0, iy = 0, iz = 0;
    elem_x_y_z(entity_id, ix, iy, iz);

    stk::mesh::EntityId elem_node[8] ;

    elem_node[0] = node_id( ix   , iy   , iz   );
    elem_node[1] = node_id( ix+1 , iy   , iz   );
    elem_node[2] = node_id( ix+1 , iy+1 , iz   );
    elem_node[3] = node_id( ix   , iy+1 , iz   );
    elem_node[4] = node_id( ix   , iy   , iz+1 );
    elem_node[5] = node_id( ix+1 , iy   , iz+1 );
    elem_node[6] = node_id( ix+1 , iy+1 , iz+1 );
    elem_node[7] = node_id( ix   , iy+1 , iz+1 );

    for (int ien = 0; ien < 8; ++ien) {
      stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, elem_node[ien], p_rank);
    }
  }
}

void HexFixture::fill_node_map(const std::map<int,std::vector<EntityId> > &parallel_distribution)
{
  m_nodes_to_procs.clear();

  std::map<int,std::vector<EntityId> >::const_iterator pd_i, pd_e;
  pd_i = parallel_distribution.begin();
  pd_e = parallel_distribution.end();
  for (; pd_i != pd_e; ++pd_i) {
    const int proc = pd_i->first;
    const std::vector<EntityId> &elements = pd_i->second;

    for (size_t e_j = 0; e_j < elements.size(); ++ e_j) {

      EntityId element_id = elements[e_j];

      size_t ix = 0, iy = 0, iz = 0;
      elem_x_y_z(element_id, ix, iy, iz);

      stk::mesh::EntityId elem_node[8] ;

      elem_node[0] = node_id( ix   , iy   , iz   );
      elem_node[1] = node_id( ix+1 , iy   , iz   );
      elem_node[2] = node_id( ix+1 , iy+1 , iz   );
      elem_node[3] = node_id( ix   , iy+1 , iz   );
      elem_node[4] = node_id( ix   , iy   , iz+1 );
      elem_node[5] = node_id( ix+1 , iy   , iz+1 );
      elem_node[6] = node_id( ix+1 , iy+1 , iz+1 );
      elem_node[7] = node_id( ix   , iy+1 , iz+1 );

      for (int ien = 0; ien < 8; ++ien) {
        stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, elem_node[ien], proc);
      }
    }
  }
}

void HexFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap)
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

    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(8);
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      size_t ix = 0, iy = 0, iz = 0;
      elem_x_y_z(entity_id, ix, iy, iz);

      elem_nodes[0] = node_id( ix   , iy   , iz   );
      elem_nodes[1] = node_id( ix+1 , iy   , iz   );
      elem_nodes[2] = node_id( ix+1 , iy+1 , iz   );
      elem_nodes[3] = node_id( ix   , iy+1 , iz   );
      elem_nodes[4] = node_id( ix   , iy   , iz+1 );
      elem_nodes[5] = node_id( ix+1 , iy   , iz+1 );
      elem_nodes[6] = node_id( ix+1 , iy+1 , iz+1 );
      elem_nodes[7] = node_id( ix   , iy+1 , iz+1 );

      stk::mesh::declare_element( m_bulk_data, m_elem_parts, elem_id( ix , iy , iz ) , elem_nodes);

      for (size_t i = 0; i<8; ++i) {
        EntityId node_id = elem_nodes[i];
        stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id );
        m_bulk_data.change_entity_parts(node, m_node_parts);

        stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, node_id, node);

        STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

        // Compute and assign coordinates to the node
        size_t nx = 0, ny = 0, nz = 0;
        node_x_y_z(elem_nodes[i], nx, ny, nz);

        Scalar * data = stk::mesh::field_data( *m_coord_field , node );

        coordMap.getNodeCoordinates(data, nx, ny, nz);
      }
    }
  }

  m_bulk_data.modification_end();
}

} // namespace simple_fields

} // fixtures
} // mesh
} // stk
