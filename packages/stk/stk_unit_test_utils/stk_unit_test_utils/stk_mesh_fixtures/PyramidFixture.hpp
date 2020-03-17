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

#ifndef STK_MESH_FIXTURES_PYRAMID_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_PYRAMID_MESH_FIXTURE_HPP

#include <stddef.h>                     // for size_t, NULL
#include <map>                          // for multimap, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId, PartVector
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <vector>                       // for vector

#include "stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp"
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire


namespace stk {
namespace mesh {
namespace fixtures {

/**
 * A 3-dimensional X*Y*Z pyramid fixture.
 * Generates 6* X*Y*Z pyramids -- each "hex" is subdivided into 6 pyramids with an added center node.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class PyramidFixture
{
 public:
  typedef double                     Scalar ;
  typedef Field<Scalar, Cartesian>   CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  PyramidFixture(   MetaData& meta
              , BulkData& bulk
              , size_t nx
              , size_t ny
              , size_t nz
              , size_t nid_start
              , size_t eid_start
            );

  PyramidFixture(   stk::ParallelMachine pm
              , size_t nx
              , size_t ny
              , size_t nz
              , stk::mesh::BulkData::AutomaticAuraOption = stk::mesh::BulkData::AUTO_AURA
            );

  PyramidFixture(   stk::ParallelMachine pm
              , size_t nx
              , size_t ny
              , size_t nz
              , std::string coordinate_name
              , stk::mesh::BulkData::AutomaticAuraOption = stk::mesh::BulkData::AUTO_AURA
            );

  ~PyramidFixture();

  const int         m_spatial_dimension;
  const size_t      m_nx;
  const size_t      m_ny;
  const size_t      m_nz;

  const size_t      node_id_start = 1;
  const size_t      elem_id_start = 1;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1)*(m_nz+1);// + (m_nx)*(m_ny)*(m_nz)/8;
  }

  size_t num_elements() const {
    return 6*(m_nx)*(m_ny)*(m_nz)/8;
  }

 private:
  MetaData*         m_meta_p;
  BulkData*         m_bulk_p;
 public:
  MetaData&         m_meta;
  BulkData&         m_bulk_data;
  PartVector        m_elem_parts;
  PartVector        m_node_parts;
  CoordFieldType &  m_coord_field ;
  bool              owns_mesh = true;
  stk::topology     m_elem_topology = stk::topology::PYRAMID_5;
  stk::topology     m_face_topology = stk::topology::QUAD_4;


  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( size_t x , size_t y , size_t z ) const  {
    return node_id_start + x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity node( size_t x , size_t y , size_t z ) const {
    return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void hex_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<size_t> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  // When creating entities, you can tell TetFixture what parts to add
  // elements and nodes.

  template <typename Iterator>
  void add_elem_parts(Iterator itr, size_t num)
  {
    ThrowRequire(!m_meta.is_commit());
    m_elem_parts.insert(m_elem_parts.end(), itr, itr + num);
  }

  //template <typename Iterator>
  //void add_node_parts(Iterator itr, size_t num)
  //{
  //  ThrowRequire(!m_meta.is_commit());
  //  m_node_parts.insert(m_node_parts.end(), itr, itr + num);
  //}

 private:

  typedef std::multimap<EntityId, int> NodeToProcsMMap;

  NodeToProcsMMap m_nodes_to_procs;

  PyramidFixture();
  PyramidFixture( const PyramidFixture &);
  PyramidFixture & operator = (const PyramidFixture &);

  void fill_node_map( int proc_rank);
};


} // fixtures
} // mesh
} // stk
#endif
