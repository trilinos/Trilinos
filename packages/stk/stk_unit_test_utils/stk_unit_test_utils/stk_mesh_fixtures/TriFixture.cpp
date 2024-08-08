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

#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, etc
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <array>
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Types.hpp>      // for EntityId, EntityIdVector
#include <stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/TriFixture.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg
#include "stk_io/IossBridge.hpp"

namespace stk {
namespace mesh {
namespace fixtures {

namespace impl {

template <int DIM>
TriFixtureImpl<DIM>::TriFixtureImpl(MetaData& meta,
                                    BulkData& bulk,
                                    size_t nx,
                                    size_t ny,
                                    size_t ,
                                    size_t nid_start,
                                    size_t eid_start)
  : m_spatial_dimension(DIM),
    m_nx( nx ),
    m_ny( ny ),
    m_meta(meta),
    m_bulk_data(bulk),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("tri_part",
                                                        DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    m_node_id_start(nid_start),
    m_elem_id_start(eid_start),
    m_elem_topology( DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3),
    m_face_topology( DIM == 2 ? stk::topology::LINE_2 : stk::topology::TRI_3)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);
}

template <int DIM>
TriFixtureImpl<DIM>::TriFixtureImpl(stk::ParallelMachine pm,
                                    size_t nx,
                                    size_t ny,
                                    stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : TriFixtureImpl(pm, nx, ny, "coordinates", autoAuraOption)
{
}

template <int DIM>
TriFixtureImpl<DIM>::TriFixtureImpl(stk::ParallelMachine pm,
                                    size_t nx,
                                    size_t ny,
                                    const std::string& coordsName,
                                    stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : m_spatial_dimension(DIM),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(DIM)
                             .set_aura_option(autoAuraOption)
                             .set_add_fmwk_data(false)
                             .create() ),
    m_nx(nx),
    m_ny(ny),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("tri_part",
                                                        DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_elem_topology( DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3),
    m_face_topology( DIM == 2 ? stk::topology::LINE_2 : stk::topology::TRI_3)
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordsName);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);

}

template <int DIM>
void TriFixtureImpl<DIM>::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<size_t> quad_range_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_quads = m_nx * m_ny ;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const size_t beg_elem = ( num_quads * p_rank ) / p_size ;
  const size_t end_elem = ( num_quads * ( p_rank + 1 ) ) / p_size ;

  for ( size_t i = beg_elem; i != end_elem; ++i) {
    quad_range_on_this_processor.push_back(i);
  }

  generate_mesh(quad_range_on_this_processor, coordMap);
}

template <int DIM>
void TriFixtureImpl<DIM>::node_x_y( EntityId entity_id, size_t &x , size_t &y ) const
{
  entity_id -= m_node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id;
}

template <int DIM>
void TriFixtureImpl<DIM>::quad_x_y( EntityId entity_id, size_t &x , size_t &y ) const
{
  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}

template <int DIM>
void TriFixtureImpl<DIM>::generate_mesh(std::vector<size_t> & quad_range_on_this_processor,
                                        const CoordinateMapping & coordMap)
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

        EntityId tri_id = 2*quad_id + tri + m_elem_id_start;
        stk::mesh::declare_element( m_bulk_data, m_elem_parts, tri_id, tri_nodes);

        for (size_t i = 0; i<3; ++i) {
          stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , tri_nodes[i] );
          m_bulk_data.change_entity_parts(node, m_node_parts);

          STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
               "This process should know about the nodes that make up its element");

          DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, tri_nodes[i], node);

          // Compute and assign coordinates to the node
          size_t nx = 0, ny = 0;
          node_x_y(tri_nodes[i], nx, ny);

          Scalar * data = stk::mesh::field_data( *m_coord_field , node );

          // The CoordinateMappings are used for 2D and 3D so make sure we give it enough space to write to.
          std::array<double, 3> temp;
          coordMap.getNodeCoordinates(temp.data(), nx, ny, 0);

          data[0] = temp[0];
          data[1] = temp[1] ;
          if(DIM == 3) data[2] = 0.;
        }
      }
    }
  }
  m_bulk_data.modification_end();
}

template <int DIM>
void TriFixtureImpl<DIM>::fill_node_map(int p_rank)
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

template class TriFixtureImpl<2>;
template class TriFixtureImpl<3>;
} // impl

namespace simple_fields {
namespace impl {

template <int DIM>
TriFixtureImpl<DIM>::TriFixtureImpl(MetaData& meta,
                                    BulkData& bulk,
                                    size_t nx,
                                    size_t ny,
                                    size_t ,
                                    size_t nid_start,
                                    size_t eid_start)
  : m_spatial_dimension(DIM),
    m_nx( nx ),
    m_ny( ny ),
    m_meta(meta),
    m_bulk_data(bulk),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("tri_part",
                                                        DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_coord_field( &m_meta.declare_field<double>(stk::topology::NODE_RANK, "Coordinates") ),
    m_node_id_start(nid_start),
    m_elem_id_start(eid_start),
    m_elem_topology( DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3),
    m_face_topology( DIM == 2 ? stk::topology::LINE_2 : stk::topology::TRI_3)
{
  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);
}

template <int DIM>
TriFixtureImpl<DIM>::TriFixtureImpl(stk::ParallelMachine pm,
                                    size_t nx,
                                    size_t ny,
                                    stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : TriFixtureImpl(pm, nx, ny, "coordinates", autoAuraOption)
{
}

template <int DIM>
TriFixtureImpl<DIM>::TriFixtureImpl(stk::ParallelMachine pm,
                                    size_t nx,
                                    size_t ny,
                                    const std::string& coordsName,
                                    stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
  : m_spatial_dimension(DIM),
    m_bulk_p( MeshBuilder(pm).set_spatial_dimension(DIM)
                             .set_aura_option(autoAuraOption)
                             .set_add_fmwk_data(false)
                             .create() ),
    m_nx(nx),
    m_ny(ny),
    m_meta(m_bulk_p->mesh_meta_data()),
    m_bulk_data(*m_bulk_p),
    m_elem_parts( 1, &m_meta.declare_part_with_topology("tri_part",
                                                        DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3) ),
    m_node_parts( 1, &m_meta.declare_part_with_topology("node_part", stk::topology::NODE) ),
    m_elem_topology( DIM == 2 ? stk::topology::TRI_3_2D : stk::topology::SHELL_TRI_3),
    m_face_topology( DIM == 2 ? stk::topology::LINE_2 : stk::topology::TRI_3)
{
  m_coord_field = &m_meta.declare_field<double>(stk::topology::NODE_RANK, coordsName);

  //put coord-field on all nodes:
  put_field_on_mesh(*m_coord_field, m_meta.universal_part(), m_spatial_dimension, nullptr);
  stk::io::set_field_output_type(*m_coord_field, stk::io::FieldOutputType::VECTOR_3D);

}

template <int DIM>
void TriFixtureImpl<DIM>::generate_mesh(const CoordinateMapping & coordMap)
{
  std::vector<size_t> quad_range_on_this_processor;

  const size_t p_size = m_bulk_data.parallel_size();
  const size_t p_rank = m_bulk_data.parallel_rank();
  const size_t num_quads = m_nx * m_ny ;

  m_nodes_to_procs.clear();
  for (unsigned proc_rank = 0; proc_rank < p_size; ++proc_rank) {
    fill_node_map(proc_rank);
  }

  const size_t beg_elem = ( num_quads * p_rank ) / p_size ;
  const size_t end_elem = ( num_quads * ( p_rank + 1 ) ) / p_size ;

  for ( size_t i = beg_elem; i != end_elem; ++i) {
    quad_range_on_this_processor.push_back(i);
  }

  generate_mesh(quad_range_on_this_processor, coordMap);
}

template <int DIM>
void TriFixtureImpl<DIM>::node_x_y( EntityId entity_id, size_t &x , size_t &y ) const
{
  entity_id -= m_node_id_start;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id;
}

template <int DIM>
void TriFixtureImpl<DIM>::quad_x_y( EntityId entity_id, size_t &x , size_t &y ) const
{
  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}

template <int DIM>
void TriFixtureImpl<DIM>::generate_mesh(std::vector<size_t> & quad_range_on_this_processor,
                                        const CoordinateMapping & coordMap)
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

        EntityId tri_id = 2*quad_id + tri + m_elem_id_start;
        stk::mesh::declare_element( m_bulk_data, m_elem_parts, tri_id, tri_nodes);

        for (size_t i = 0; i<3; ++i) {
          stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , tri_nodes[i] );
          m_bulk_data.change_entity_parts(node, m_node_parts);

          STK_ThrowRequireMsg( m_bulk_data.is_valid(node),
               "This process should know about the nodes that make up its element");

          stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, tri_nodes[i], node);

          // Compute and assign coordinates to the node
          size_t nx = 0, ny = 0;
          node_x_y(tri_nodes[i], nx, ny);

          Scalar * data = stk::mesh::field_data( *m_coord_field , node );

          // The CoordinateMappings are used for 2D and 3D so make sure we give it enough space to write to.
          std::array<double, 3> temp;
          coordMap.getNodeCoordinates(temp.data(), nx, ny, 0);

          data[0] = temp[0];
          data[1] = temp[1] ;
          if(DIM == 3) data[2] = 0.;
        }
      }
    }
  }
  m_bulk_data.modification_end();
}

template <int DIM>
void TriFixtureImpl<DIM>::fill_node_map(int p_rank)
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
          stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, tri_node[i] , p_rank);
        }
      }
    }
  }
}

template class TriFixtureImpl<2>;
template class TriFixtureImpl<3>;
} // impl
} // namespace simple_fields

} // fixtures
} // mesh
} // stk
