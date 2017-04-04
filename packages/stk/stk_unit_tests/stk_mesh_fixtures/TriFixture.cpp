// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_mesh/fixtures/TriFixture.hpp>
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityIdVector
#include <stk_mesh/fixtures/FixtureNodeSharing.hpp>
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequireMsg
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_mesh/fixtures/CoordinateMapping.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { struct ConnectivityMap; } }





namespace stk {
namespace mesh {
namespace fixtures {

TriFixture::TriFixture(   stk::ParallelMachine pm
              , size_t nx
              , size_t ny
              , stk::mesh::BulkData::AutomaticAuraOption autoAuraOption
              , ConnectivityMap const* connectivity_map
            )
  : m_spatial_dimension(2),
    m_nx(nx),
    m_ny(ny),
    m_meta( m_spatial_dimension ),
    m_bulk_data(  m_meta
                , pm
                , autoAuraOption
#ifdef SIERRA_MIGRATION
                , false
#endif
                , connectivity_map
               ),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("tri_part", stk::topology::TRI_3) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( m_meta.declare_field<CoordFieldType>(stk::topology::NODE_RANK, "Coordinates") )
{

  //put coord-field on all nodes:
  put_field(
    m_coord_field,
    m_meta.universal_part(),
    m_spatial_dimension);

}

void TriFixture::generate_mesh()
{
  std::vector<size_t> quad_range_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_elems = m_nx * m_ny ;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const size_t beg_elem = ( num_elems * p_rank ) / p_size ;
  const size_t end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( size_t i = beg_elem; i != end_elem; ++i) {
    quad_range_on_this_processor.push_back(i);
  }

  generate_mesh(quad_range_on_this_processor);
}

void TriFixture::node_x_y( EntityId entity_id, size_t &x , size_t &y ) const
{
  entity_id -= 1;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id;
}

void TriFixture::quad_x_y( EntityId entity_id, size_t &x , size_t &y ) const
{
  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}

void TriFixture::generate_mesh(std::vector<size_t> & quad_range_on_this_processor)
{
  {
    //sort and unique the input elements
    std::vector<size_t>::iterator ib = quad_range_on_this_processor.begin();
    std::vector<size_t>::iterator ie = quad_range_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    quad_range_on_this_processor.erase(ib, ie);
  }

  m_bulk_data.modification_begin();

  {
    int tri_vert[][3] = { {0,1,2}, {0,2,3} };

    // Declare the elements that belong on this process
    std::vector<size_t>::iterator ib = quad_range_on_this_processor.begin();
    const std::vector<size_t>::iterator ie = quad_range_on_this_processor.end();
    stk::mesh::EntityIdVector elem_nodes(4);
    stk::mesh::EntityIdVector tri_nodes(3);

    for (; ib != ie; ++ib) {
      size_t quad_id = *ib;
      size_t ix = 0, iy = 0;
      quad_x_y(quad_id, ix, iy);

      elem_nodes[0] = node_id( ix   , iy   );
      elem_nodes[1] = node_id( ix+1 , iy   );
      elem_nodes[2] = node_id( ix+1 , iy+1 );
      elem_nodes[3] = node_id( ix   , iy+1 );

      for (size_t tri = 0; tri < 2; ++tri) {
        tri_nodes[0] = elem_nodes[tri_vert[tri][0]];
        tri_nodes[1] = elem_nodes[tri_vert[tri][1]];
        tri_nodes[2] = elem_nodes[tri_vert[tri][2]];

        EntityId tri_id = 2*quad_id + tri + 1;
        stk::mesh::declare_element( m_bulk_data, m_elem_parts, tri_id, tri_nodes);

        for (size_t i = 0; i<3; ++i) {
          stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , tri_nodes[i] );
          m_bulk_data.change_entity_parts(node, m_node_parts);

          ThrowRequireMsg( m_bulk_data.is_valid(node),
               "This process should know about the nodes that make up its element");

          DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, tri_nodes[i], node);

          // Compute and assign coordinates to the node
          size_t nx = 0, ny = 0;
          node_x_y(tri_nodes[i], nx, ny);

          Scalar * data = stk::mesh::field_data( m_coord_field , node );

          data[0] = (Scalar)nx ;
          data[1] = (Scalar)ny ;
        }
      }
    }
  }
  m_bulk_data.modification_end();
}

void TriFixture::fill_node_map(int p_rank)
{

  std::vector<EntityId> element_ids_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t num_elems = m_nx * m_ny ;

  const EntityId beg_elem = ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  {
    int tri_vert[][3] = { {0,1,2}, {0,2,3} };

    // Declare the elements that belong on this process
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      size_t quad_id = *ib;
      size_t ix = 0, iy = 0;
      quad_x_y(quad_id, ix, iy);

      stk::mesh::EntityId elem_node[4] ;
      stk::mesh::EntityId tri_node[3];

      elem_node[0] = node_id( ix   , iy    );
      elem_node[1] = node_id( ix+1 , iy    );
      elem_node[2] = node_id( ix+1 , iy+1  );
      elem_node[3] = node_id( ix   , iy+1  );

      for (size_t tri = 0; tri < 2; ++tri) {
        tri_node[0] = elem_node[tri_vert[tri][0]];
        tri_node[1] = elem_node[tri_vert[tri][1]];
        tri_node[2] = elem_node[tri_vert[tri][2]];

        for (size_t i = 0; i<3; ++i) {
          AddToNodeProcsMMap(m_nodes_to_procs, tri_node[i] , p_rank);
        }
      }
    }
  }
}

} // fixtures
} // mesh
} // stk
