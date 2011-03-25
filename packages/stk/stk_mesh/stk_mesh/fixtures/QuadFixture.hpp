/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_QUAD_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_QUAD_MESH_FIXTURE_HPP

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

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
  typedef Field<Scalar*,ElementNode>  CoordGatherFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  QuadFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny );

  ~QuadFixture() {}

  const unsigned                m_spatial_dimension;
  fem::FEMMetaData              m_fem_meta ;
  BulkData                      m_bulk_data ;
  Part &                        m_quad_part ;
  CoordFieldType &              m_coord_field ;
  CoordGatherFieldType &        m_coord_gather_field ;
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
  Entity * node( unsigned x , unsigned y ) const
  { return m_bulk_data.get_entity( fem::NODE_RANK , node_id(x, y) ); }

  /**
   * Thinking in terms of rows and columns of elements, get the element in
   * the (x, y) position. Return NULL if this process doesn't know about
   * this element.
   */
  Entity * elem( unsigned x , unsigned y ) const
  { return m_bulk_data.get_entity( m_fem_meta.element_rank(), elem_id(x, y)); }

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

 private:
  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor );

  QuadFixture();
  QuadFixture( const QuadFixture & );
  QuadFixture & operator = ( const QuadFixture & );
};

} // fixtures
} // mesh
} // stk
#endif
