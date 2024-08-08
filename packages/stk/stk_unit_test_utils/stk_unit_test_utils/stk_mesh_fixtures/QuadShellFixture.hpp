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

#ifndef STK_MESH_FIXTURES_QUAD_SHELL_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_QUAD_SHELL_MESH_FIXTURE_HPP

#include <map>                          // for multimap, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                       // for string
#include <vector>                       // for vector

#include "stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/BuildMesh.hpp"

namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {

/**
 * An 3-dimensional X*Y quad fixture.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class QuadShellFixture
{
public:
  typedef double Scalar;
  typedef Field<Scalar> CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */

  QuadShellFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny,
                    const std::vector<std::string>& rank_names = std::vector<std::string>() );

  QuadShellFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny, const std::string& coordsName,
                    const std::vector<std::string>& rank_names = std::vector<std::string>() );

  QuadShellFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny, bool auraOn );

  QuadShellFixture( MetaData& meta, BulkData& bulk, unsigned nx , unsigned ny, unsigned nid_start, unsigned eid_start );

  QuadShellFixture( MetaData& meta, BulkData& bulk, unsigned nx , unsigned ny, unsigned nz, unsigned nid_start, unsigned eid_start )
  : QuadShellFixture(meta, bulk, nx, ny, nid_start, eid_start) {}

  ~QuadShellFixture();

private:
  std::shared_ptr<BulkData>     m_bulk_p;

public:
  const unsigned                m_spatial_dimension;
  MetaData&                     m_meta ;
  BulkData&                     m_bulk_data ;
  Part &                        m_quad_part ;
  PartVector                    m_elem_parts;
  PartVector                    m_node_parts;
  CoordFieldType *              m_coord_field ;
  const unsigned                m_node_id_start = 1;
  const unsigned                m_elem_id_start = 1;
  bool                          owns_mesh = true;
  stk::topology                 m_elem_topology = stk::topology::SHELL_QUAD_4;
  stk::topology                 m_face_topology = stk::topology::QUAD_4;
  const unsigned                m_nx ;
  const unsigned                m_ny ;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1);
  }

  size_t num_elements() const {
    return (m_nx)*(m_ny);
  }

  /**
   * Thinking in terms of rows and columns of nodes, get the id of the node in
   * the (x, y) position.
   */
  EntityId node_id( unsigned x , unsigned y ) const
    { return m_node_id_start + x + ( m_nx + 1 ) * y ; }

  /**
   * Thinking in terms of rows and columns of elements, get the id of the
   * element in the (x, y) position.
   */
  EntityId elem_id( unsigned x , unsigned y ) const
    { return m_elem_id_start + x + m_nx * y ; }

  /**
   * Thinking in terms of rows and columns of nodes, get the node in
   * the (x, y) position. Return NULL if this process doesn't know about
   * this node.
   */
  Entity node( unsigned x , unsigned y ) const
  { return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y) ); }

  /**
   * Thinking in terms of rows and columns of elements, get the element in
   * the (x, y) position. Return NULL if this process doesn't know about
   * this element.
   */
  Entity elem( unsigned x , unsigned y ) const
  { return m_bulk_data.get_entity( stk::topology::ELEMENT_RANK, elem_id(x, y)); }

  /**
   * Thinking in terms of a 2D grid of nodes, compute the (x, y) position
   * of a node given it's id.
   */
  void node_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const;

  /**
   * Thinking in terms of a 2D grid of elements, compute the (x, y) position
   * of an element given it's id.
   */
  void elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping() );

private:
  typedef std::multimap<EntityId, int> NodeToProcsMMap;

  NodeToProcsMMap m_nodes_to_procs;

  QuadShellFixture();
  QuadShellFixture( const QuadShellFixture & );
  QuadShellFixture & operator=( const QuadShellFixture & );

  void fill_node_map( int proc_rank);
};

namespace simple_fields {

/**
 * An 3-dimensional X*Y quad fixture.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
QuadShellFixture
{
public:
  typedef double Scalar;
  typedef Field<Scalar> CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */

  QuadShellFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny,
                    const std::vector<std::string>& rank_names = std::vector<std::string>() );

  QuadShellFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny, const std::string& coordsName,
                    const std::vector<std::string>& rank_names = std::vector<std::string>() );

  QuadShellFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny, bool auraOn );

  QuadShellFixture( MetaData& meta, BulkData& bulk, unsigned nx , unsigned ny, unsigned nid_start, unsigned eid_start );

  QuadShellFixture( MetaData& meta, BulkData& bulk, unsigned nx , unsigned ny, unsigned nz, unsigned nid_start, unsigned eid_start )
  : QuadShellFixture(meta, bulk, nx, ny, nid_start, eid_start) {}

  ~QuadShellFixture();

private:
  std::shared_ptr<BulkData>     m_bulk_p;

public:
  const unsigned                m_spatial_dimension;
  MetaData&                     m_meta ;
  BulkData&                     m_bulk_data ;
  Part &                        m_quad_part ;
  PartVector                    m_elem_parts;
  PartVector                    m_node_parts;
  CoordFieldType *              m_coord_field ;
  const unsigned                m_node_id_start = 1;
  const unsigned                m_elem_id_start = 1;
  bool                          owns_mesh = true;
  stk::topology                 m_elem_topology = stk::topology::SHELL_QUAD_4;
  stk::topology                 m_face_topology = stk::topology::QUAD_4;
  const unsigned                m_nx ;
  const unsigned                m_ny ;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1);
  }

  size_t num_elements() const {
    return (m_nx)*(m_ny);
  }

  /**
   * Thinking in terms of rows and columns of nodes, get the id of the node in
   * the (x, y) position.
   */
  EntityId node_id( unsigned x , unsigned y ) const
    { return m_node_id_start + x + ( m_nx + 1 ) * y ; }

  /**
   * Thinking in terms of rows and columns of elements, get the id of the
   * element in the (x, y) position.
   */
  EntityId elem_id( unsigned x , unsigned y ) const
    { return m_elem_id_start + x + m_nx * y ; }

  /**
   * Thinking in terms of rows and columns of nodes, get the node in
   * the (x, y) position. Return NULL if this process doesn't know about
   * this node.
   */
  Entity node( unsigned x , unsigned y ) const
  { return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y) ); }

  /**
   * Thinking in terms of rows and columns of elements, get the element in
   * the (x, y) position. Return NULL if this process doesn't know about
   * this element.
   */
  Entity elem( unsigned x , unsigned y ) const
  { return m_bulk_data.get_entity( stk::topology::ELEMENT_RANK, elem_id(x, y)); }

  /**
   * Thinking in terms of a 2D grid of nodes, compute the (x, y) position
   * of a node given it's id.
   */
  void node_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const;

  /**
   * Thinking in terms of a 2D grid of elements, compute the (x, y) position
   * of an element given it's id.
   */
  void elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping() );

private:
  typedef std::multimap<EntityId, int> NodeToProcsMMap;

  NodeToProcsMMap m_nodes_to_procs;

  QuadShellFixture();
  QuadShellFixture( const QuadShellFixture & );
  QuadShellFixture & operator=( const QuadShellFixture & );

  void fill_node_map( int proc_rank);
};

} // namespace simple_fields

} // fixtures
} // mesh
} // stk
#endif
