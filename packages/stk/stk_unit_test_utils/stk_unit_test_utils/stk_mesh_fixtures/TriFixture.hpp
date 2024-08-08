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

#ifndef STK_MESH_FIXTURES_TRI_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_TRI_MESH_FIXTURE_HPP

#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire
#include <map>                          // for multimap, etc
#include <memory>
#include <stddef.h>                     // for size_t, NULL
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId, PartVector
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <vector>                       // for vector

#include "stk_unit_test_utils/stk_mesh_fixtures/CoordinateMapping.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk {
namespace mesh {
namespace fixtures {

/**
 * A 2-dimensional X*Y tri fixture.
 * Generates 2* X*Y triangles -- each "quad" is subdivided into 2 triangles.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
namespace impl {

template <int DIM>
class TriFixtureImpl
{
 public:
  typedef double        Scalar;
  typedef Field<Scalar> CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  // The nz argument is ignored, only here to allow use in templated functions that also operate on 3D fixtures
  TriFixtureImpl(MetaData& meta, BulkData& bulk, size_t nx, size_t ny, size_t nz, size_t nid_start, size_t eid_start);

  TriFixtureImpl(stk::ParallelMachine pm,
                 size_t nx,
                 size_t ny,
                 stk::mesh::BulkData::AutomaticAuraOption = stk::mesh::BulkData::AUTO_AURA);

  TriFixtureImpl(stk::ParallelMachine pm,
                 size_t nx,
                 size_t ny,
                 const std::string& coordsName,
                 stk::mesh::BulkData::AutomaticAuraOption = stk::mesh::BulkData::AUTO_AURA);

  const int         m_spatial_dimension;

 private:
  std::shared_ptr<BulkData> m_bulk_p;

 public:
  const size_t      m_nx;
  const size_t      m_ny;
  MetaData &        m_meta;
  BulkData &        m_bulk_data;
  PartVector        m_elem_parts;
  PartVector        m_node_parts;
  CoordFieldType *  m_coord_field;
  const size_t      m_node_id_start = 1;
  const size_t      m_elem_id_start = 1;
  stk::topology     m_elem_topology;
  stk::topology     m_face_topology;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1);
  }

  size_t num_elements() const {
    return 2*(m_nx)*(m_ny);
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( size_t x , size_t y ) const  {
    return m_node_id_start + x + ( m_nx + 1 ) * y;
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity node( size_t x , size_t y ) const {
    return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y) );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y( EntityId entity_id, size_t &x , size_t &y ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void quad_x_y( EntityId entity_id, size_t &x , size_t &y ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<size_t> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping());

 private:
  typedef std::multimap<EntityId, int> NodeToProcsMMap;

  NodeToProcsMMap m_nodes_to_procs;

  TriFixtureImpl();
  TriFixtureImpl( const TriFixtureImpl &);
  TriFixtureImpl & operator=(const TriFixtureImpl &);

  void fill_node_map( int proc_rank);
};

} // impl

using TriFixture = impl::TriFixtureImpl<2>;
using TriShellFixture = impl::TriFixtureImpl<3>;

namespace simple_fields {
namespace impl {

template <int DIM>
class TriFixtureImpl
{
 public:
  typedef double        Scalar;
  typedef Field<Scalar> CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  // The nz argument is ignored, only here to allow use in templated functions that also operate on 3D fixtures
  TriFixtureImpl(MetaData& meta, BulkData& bulk, size_t nx, size_t ny, size_t nz, size_t nid_start, size_t eid_start);

  TriFixtureImpl(stk::ParallelMachine pm,
                 size_t nx,
                 size_t ny,
                 stk::mesh::BulkData::AutomaticAuraOption = stk::mesh::BulkData::AUTO_AURA);

  TriFixtureImpl(stk::ParallelMachine pm,
                 size_t nx,
                 size_t ny,
                 const std::string& coordsName,
                 stk::mesh::BulkData::AutomaticAuraOption = stk::mesh::BulkData::AUTO_AURA);

  const int         m_spatial_dimension;

 private:
  std::shared_ptr<BulkData> m_bulk_p;

 public:
  const size_t      m_nx;
  const size_t      m_ny;
  MetaData &        m_meta;
  BulkData &        m_bulk_data;
  PartVector        m_elem_parts;
  PartVector        m_node_parts;
  CoordFieldType *  m_coord_field;
  const size_t      m_node_id_start = 1;
  const size_t      m_elem_id_start = 1;
  stk::topology     m_elem_topology;
  stk::topology     m_face_topology;

  size_t num_nodes() const {
    return (m_nx+1)*(m_ny+1);
  }

  size_t num_elements() const {
    return 2*(m_nx)*(m_ny);
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( size_t x , size_t y ) const  {
    return m_node_id_start + x + ( m_nx + 1 ) * y;
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity node( size_t x , size_t y ) const {
    return m_bulk_data.get_entity( stk::topology::NODE_RANK , node_id(x, y) );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y( EntityId entity_id, size_t &x , size_t &y ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void quad_x_y( EntityId entity_id, size_t &x , size_t &y ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<size_t> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping());

 private:
  typedef std::multimap<EntityId, int> NodeToProcsMMap;

  NodeToProcsMMap m_nodes_to_procs;

  TriFixtureImpl();
  TriFixtureImpl( const TriFixtureImpl &);
  TriFixtureImpl & operator=(const TriFixtureImpl &);

  void fill_node_map( int proc_rank);
};

} // impl

using TriFixture
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead") = impl::TriFixtureImpl<2>;
using TriShellFixture
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")  = impl::TriFixtureImpl<3>;

} // namespace simple_fields

} // fixtures
} // mesh
} // stk
#endif
