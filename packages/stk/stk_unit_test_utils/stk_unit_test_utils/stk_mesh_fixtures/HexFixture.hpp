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

#ifndef STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stddef.h>
#include <map>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <string>
#include <vector>

#include "stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace mesh {
namespace fixtures {

/**
 * A 3-dimensional X*Y*Z hex fixture.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class HexFixture
{
public:
  static std::string name() { return "HexFixture"; }
  typedef double        Scalar;
  typedef Field<Scalar> CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  HexFixture(MetaData& meta,
             BulkData& bulk,
             size_t nx,
             size_t ny,
             size_t nz,
             size_t nid_start,
             size_t eid_start);

  HexFixture(stk::ParallelMachine pm,
             size_t nx,
             size_t ny,
             size_t nz);

  HexFixture(stk::ParallelMachine pm,
             size_t nx,
             size_t ny,
             size_t nz,
             const std::string& coordinate_name);

  HexFixture(stk::ParallelMachine pm,
             size_t nx,
             size_t ny,
             size_t nz,
             bool auraOn);

  const int                     m_spatial_dimension;
  const size_t                  m_nx;
  const size_t                  m_ny;
  const size_t                  m_nz;
  const size_t                  node_id_start = 1;
  const size_t                  elem_id_start = 1;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1)*(m_nz+1);
  }

  size_t num_elements() const {
    return (m_nx)*(m_ny)*(m_nz);
  }

 private:
  std::shared_ptr<BulkData> m_bulk_p;
 public:
  MetaData&                     m_meta;
  BulkData&                     m_bulk_data;
  PartVector                    m_elem_parts;
  PartVector                    m_node_parts;
  CoordFieldType *              m_coord_field ;
  bool owns_mesh = true;
  stk::topology     m_elem_topology = stk::topology::HEX_8;
  stk::topology     m_face_topology = stk::topology::QUAD_4;

  ~HexFixture();

  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( size_t x , size_t y , size_t z ) const  {
    return node_id_start + x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  EntityId elem_id( size_t x , size_t y , size_t z ) const  {
    return elem_id_start + x + m_nx * ( y + m_ny * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity node( size_t x , size_t y , size_t z ) const {
    return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the elements in the
   * (x, y, z) position. Return NULL if this process doesn't know about this
   * element.
   */
  Entity elem( size_t x , size_t y , size_t z ) const {
    return m_bulk_data.get_entity( stk::topology::ELEMENT_RANK, elem_id(x, y, z) );
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
  void elem_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void fill_node_map( int proc_rank);

  void fill_node_map(const std::map<int,std::vector<EntityId> > &parallel_distribution);

  // When creating entities, you can tell HexFixture what parts to add
  // elements and nodes.

  template <typename Iterator>
  void add_elem_parts(Iterator itr, size_t num)
  {
    STK_ThrowRequire(!m_meta.is_commit());
    m_elem_parts.insert(m_elem_parts.end(), itr, itr + num);
  }

  //template <typename Iterator>
  //void add_node_parts(Iterator itr, size_t num)
  //{
  //  STK_ThrowRequire(!m_meta.is_commit());
  //  m_node_parts.insert(m_node_parts.end(), itr, itr + num);
  //}

 private:
  NodeToProcsMMap m_nodes_to_procs;

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator=(const HexFixture &);
};

namespace simple_fields {

/**
 * A 3-dimensional X*Y*Z hex fixture.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class HexFixture
{
public:
  static std::string name() { return "HexFixture"; }
  typedef double        Scalar;
  typedef Field<Scalar> CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
  HexFixture(MetaData& meta,
             BulkData& bulk,
             size_t nx,
             size_t ny,
             size_t nz,
             size_t nid_start,
             size_t eid_start);

  STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
  HexFixture(stk::ParallelMachine pm,
             size_t nx,
             size_t ny,
             size_t nz);

  STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
  HexFixture(stk::ParallelMachine pm,
             size_t nx,
             size_t ny,
             size_t nz,
             const std::string& coordinate_name);

  STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
  HexFixture(stk::ParallelMachine pm,
             size_t nx,
             size_t ny,
             size_t nz,
             bool auraOn);

  const int                     m_spatial_dimension;
  const size_t                  m_nx;
  const size_t                  m_ny;
  const size_t                  m_nz;
  const size_t                  node_id_start = 1;
  const size_t                  elem_id_start = 1;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1)*(m_nz+1);
  }

  size_t num_elements() const {
    return (m_nx)*(m_ny)*(m_nz);
  }

 private:
  std::shared_ptr<BulkData> m_bulk_p;
 public:
  MetaData&                     m_meta;
  BulkData&                     m_bulk_data;
  PartVector                    m_elem_parts;
  PartVector                    m_node_parts;
  CoordFieldType *              m_coord_field ;
  bool owns_mesh = true;
  stk::topology     m_elem_topology = stk::topology::HEX_8;
  stk::topology     m_face_topology = stk::topology::QUAD_4;

  ~HexFixture();

  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( size_t x , size_t y , size_t z ) const  {
    return node_id_start + x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  EntityId elem_id( size_t x , size_t y , size_t z ) const  {
    return elem_id_start + x + m_nx * ( y + m_ny * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity node( size_t x , size_t y , size_t z ) const {
    return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the elements in the
   * (x, y, z) position. Return NULL if this process doesn't know about this
   * element.
   */
  Entity elem( size_t x , size_t y , size_t z ) const {
    return m_bulk_data.get_entity( stk::topology::ELEMENT_RANK, elem_id(x, y, z) );
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
  void elem_x_y_z( EntityId entity_id, size_t &x , size_t &y , size_t &z ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void fill_node_map( int proc_rank);

  void fill_node_map(const std::map<int,std::vector<EntityId> > &parallel_distribution);

  // When creating entities, you can tell HexFixture what parts to add
  // elements and nodes.

  template <typename Iterator>
  void add_elem_parts(Iterator itr, size_t num)
  {
    STK_ThrowRequire(!m_meta.is_commit());
    m_elem_parts.insert(m_elem_parts.end(), itr, itr + num);
  }

  //template <typename Iterator>
  //void add_node_parts(Iterator itr, size_t num)
  //{
  //  STK_ThrowRequire(!m_meta.is_commit());
  //  m_node_parts.insert(m_node_parts.end(), itr, itr + num);
  //}

 private:
  NodeToProcsMMap m_nodes_to_procs;

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator=(const HexFixture &);
};

} // namespace simple_fields

} // fixtures
} // mesh
} // stk
#endif
