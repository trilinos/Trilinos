/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <unit_tests/UnitTestBoxFixture.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

//----------------------------------------------------------------------------

BoxFixture::~BoxFixture()
{}

BoxFixture::BoxFixture( stk::ParallelMachine pm ,
                        unsigned elem_nx , unsigned elem_ny , unsigned elem_nz )
  : m_meta_data( stk::mesh::fem_entity_rank_names() ),
    m_bulk_data( m_meta_data , pm ),
    m_elem_block( m_meta_data.declare_part("block", stk::mesh::Element) ),
    m_coord_field( m_meta_data.declare_field<CoordFieldType>("Coordinates") ),
    m_elem_nx( elem_nx ),
    m_elem_ny( elem_ny ),
    m_elem_nz( elem_nz )
{
  typedef shards::Hexahedron<8> Hex8 ;
  enum { SpatialDim = 3 };
  enum { NodesPerElem = Hex8::node_count };

  // Set topology of the element block part
  stk::mesh::set_cell_topology<shards::Hexahedron<8> >(m_elem_block);

  //put coord-field on all nodes:
  stk::mesh::put_field( m_coord_field, stk::mesh::Node, m_meta_data.universal_part(), SpatialDim );

  // Meta data is complete
  m_meta_data.commit();

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = elem_nx * elem_ny * elem_nz ;

  const unsigned beg_elem = ( num_elems * p_rank ) / p_size ;
  const unsigned end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  m_bulk_data.modification_begin();

  {
    unsigned i_elem = 0 ;

    for ( unsigned iz = 0 ; iz < elem_nz ; ++iz ) {
    for ( unsigned iy = 0 ; iy < elem_ny ; ++iy ) {
    for ( unsigned ix = 0 ; ix < elem_nx ; ++ix , ++i_elem ) {
      if ( beg_elem <= i_elem && i_elem < end_elem ) {

        stk::mesh::EntityId elem_node[8] ;

        elem_node[0] = node_id( ix   , iy   , iz );
        elem_node[1] = node_id( ix+1 , iy   , iz );
        elem_node[2] = node_id( ix+1 , iy   , iz+1 );
        elem_node[3] = node_id( ix   , iy   , iz+1 );
        elem_node[4] = node_id( ix   , iy+1 , iz );
        elem_node[5] = node_id( ix+1 , iy+1 , iz );
        elem_node[6] = node_id( ix+1 , iy+1 , iz );
        elem_node[7] = node_id( ix   , iy+1 , iz+1 );

        stk::mesh::declare_element( m_bulk_data, m_elem_block, elem_id( ix , iy , iz ) , elem_node);
      }
    }
  }
  }
  }

  m_bulk_data.modification_end();

  for ( unsigned iz = 0 ; iz <= elem_nz ; ++iz ) {
  for ( unsigned iy = 0 ; iy <= elem_ny ; ++iy ) {
  for ( unsigned ix = 0 ; ix <= elem_nx ; ++ix ) {
    stk::mesh::Entity * const entity =
      m_bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy,iz) );

    if ( entity ) {
      Scalar * data = stk::mesh::field_data( m_coord_field , *entity );

      data[0] = ix ;
      data[1] = iy ;
      data[2] = - iz ;
    }
  }
  }
  }
}

//----------------------------------------------------------------------------

QuadFixture::~QuadFixture()
{}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                        unsigned elem_nx , unsigned elem_ny )
  : m_meta_data( stk::mesh::fem_entity_rank_names() ),
    m_bulk_data( m_meta_data , pm ),
    m_elem_block( m_meta_data.declare_part("block", stk::mesh::Element) ),
    m_coord_field( m_meta_data.declare_field<CoordFieldType>("Coordinates") ),
    m_elem_nx( elem_nx ),
    m_elem_ny( elem_ny )
{
  typedef shards::Quadrilateral<4> Quad4 ;
  enum { SpatialDim = 2 };
  enum { NodesPerElem = Quad4::node_count };

  // Set topology of the element block part
  stk::mesh::set_cell_topology< Quad4 >(m_elem_block);

  //put coord-field on all nodes:
  stk::mesh::put_field( m_coord_field, stk::mesh::Node, m_meta_data.universal_part(), SpatialDim );

  // Meta data is complete
  m_meta_data.commit();

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = elem_nx * elem_ny ;

  const unsigned beg_elem = ( num_elems * p_rank ) / p_size ;
  const unsigned end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  m_bulk_data.modification_begin();

  {
    unsigned i_elem = 0 ;

    for ( unsigned iy = 0 ; iy < elem_ny ; ++iy ) {
    for ( unsigned ix = 0 ; ix < elem_nx ; ++ix , ++i_elem ) {
      if ( beg_elem <= i_elem && i_elem < end_elem ) {

        stk::mesh::EntityId elem_node[4] ;

        elem_node[0] = node_id( ix   , iy );
        elem_node[1] = node_id( ix+1 , iy );
        elem_node[2] = node_id( ix+1 , iy+1 );
        elem_node[3] = node_id( ix   , iy+1 );

        stk::mesh::declare_element( m_bulk_data, m_elem_block, elem_id( ix , iy ) , elem_node);
      }
    }
  }
  }

  m_bulk_data.modification_end();

  for ( unsigned iy = 0 ; iy <= elem_ny ; ++iy ) {
  for ( unsigned ix = 0 ; ix <= elem_nx ; ++ix ) {
    stk::mesh::Entity * const entity =
      m_bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy) );

    if ( entity ) {
      Scalar * data = stk::mesh::field_data( m_coord_field , *entity );

      data[0] = ix ;
      data[1] = iy ;
    }
  }
  }
}

