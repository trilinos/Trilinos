/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

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
  typedef double                     Scalar ;
  typedef Field<Scalar, Cartesian>   CoordFieldType;
  typedef Field<Scalar*,ElementNode> CoordGatherFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  HexFixture(stk::ParallelMachine pm, unsigned nx, unsigned ny, unsigned nz);

  const int                     m_spatial_dimension;
  const unsigned                m_nx;
  const unsigned                m_ny;
  const unsigned                m_nz;
  fem::FEMMetaData              m_fem_meta;
  BulkData                      m_bulk_data;
  Part &                        m_hex_part;
  Part &                        m_node_part;
  CoordFieldType &              m_coord_field ;
  CoordGatherFieldType &        m_coord_gather_field ;

  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( unsigned x , unsigned y , unsigned z ) const  {
    return 1 + x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  EntityId elem_id( unsigned x , unsigned y , unsigned z ) const  {
    return 1 + x + m_nx * ( y + m_ny * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity * node( unsigned x , unsigned y , unsigned z ) const {
    return m_bulk_data.get_entity( fem::FEMMetaData::NODE_RANK , node_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the elements in the
   * (x, y, z) position. Return NULL if this process doesn't know about this
   * element.
   */
  Entity * elem( unsigned x , unsigned y , unsigned z ) const {
    return m_bulk_data.get_entity( m_fem_meta.element_rank(), elem_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y_z( EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void elem_x_y_z( EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh();

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor );

 private:

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator = (const HexFixture &);
};


} // fixtures
} // mesh
} // stk
#endif
