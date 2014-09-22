/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_QUAD_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_QUAD_MESH_FIXTURE_HPP

#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityId
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }




namespace stk {
namespace mesh {
namespace fixtures {

/**
 * An 2-dimensional X*Y quad fixture.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class QuadFixture
{
 public:
  typedef int Scalar ;
  typedef Field<Scalar, Cartesian>    CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  QuadFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny, const std::vector<std::string>& rank_names = std::vector<std::string>() );

  ~QuadFixture() {}

  const unsigned                m_spatial_dimension;
  MetaData                      m_meta ;
  BulkData                      m_bulk_data ;
  Part &                        m_quad_part ;
  CoordFieldType &              m_coord_field ;
  const unsigned                m_nx ;
  const unsigned                m_ny ;

  /**
   * Thinking in terms of rows and columns of nodes, get the id of the node in
   * the (x, y) position.
   */
  EntityId node_id( unsigned x , unsigned y ) const
    { return 1 + x + ( m_nx + 1 ) * y ; }

  /**
   * Thinking in terms of rows and columns of elements, get the id of the
   * element in the (x, y) position.
   */
  EntityId elem_id( unsigned x , unsigned y ) const
    { return 1 + x + m_nx * y ; }

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
  void generate_mesh();

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor );
 private:

  typedef std::multimap<EntityId, int> NodeToProcsMMap;

  NodeToProcsMMap m_nodes_to_procs;

  QuadFixture();
  QuadFixture( const QuadFixture & );
  QuadFixture & operator = ( const QuadFixture & );

  void fill_node_map( int proc_rank);
};

} // fixtures
} // mesh
} // stk
#endif
