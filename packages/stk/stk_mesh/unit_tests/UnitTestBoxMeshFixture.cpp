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

#include <unit_tests/UnitTestBoxMeshFixture.hpp>
#include <unit_tests/stk_utest_macros.hpp>

BoxMeshFixture::~BoxMeshFixture()
{}

BoxMeshFixture::BoxMeshFixture( stk::ParallelMachine pm )
  : m_meta_data( stk::mesh::fem_entity_rank_names() ),
    m_bulk_data( m_meta_data , pm ),
    m_elem_block( m_meta_data.declare_part("block1", stk::mesh::Element) ),
    m_elem_block2 ( m_meta_data.declare_part("block2", stk::mesh::Element) ),
    m_elem_block3 ( m_meta_data.declare_part("block3", stk::mesh::Element) ),
    m_coord_field( m_meta_data.declare_field<CoordFieldType>("Coordinates") ),
    m_coord_gather_field(
      m_meta_data.declare_field<CoordGatherFieldType>("GatherCoordinates") ),
    m_quad_field( m_meta_data.declare_field<QuadFieldType>("Quad") ),
    m_basis_field( m_meta_data.declare_field<BasisFieldType>("Basis") )
{
  typedef shards::Hexahedron<8> Hex8 ;
  enum { SpatialDim = 3 };
  enum { NodesPerElem = Hex8::node_count };

  // Set topology of the element block part
  stk::mesh::set_cell_topology<shards::Hexahedron<8> >(m_elem_block);
  stk::mesh::set_cell_topology<shards::Particle >(m_elem_block3);

  //put coord-field on all nodes:
  stk::mesh::put_field( m_coord_field, stk::mesh::Node, m_meta_data.universal_part(), SpatialDim );

  //put coord-gather-field on all elements:
  stk::mesh::put_field( m_coord_gather_field, stk::mesh::Element, m_meta_data.universal_part(), NodesPerElem);

  // Field relation so coord-gather-field on elements points
  // to coord-field of the element's nodes
  m_meta_data.declare_field_relation( m_coord_gather_field, stk::mesh::element_node_stencil<shards::Hexahedron<8> >, m_coord_field);
  m_meta_data.declare_field_relation( m_coord_gather_field, stk::mesh::element_node_stencil<void>, m_coord_field);
  m_meta_data.declare_field_relation( m_coord_gather_field, stk::mesh::element_node_lock_stencil<void>, m_coord_field);

  // Meta data is complete
  m_meta_data.commit();

  STKUNIT_EXPECT_TRUE( NULL != stk::mesh::ElementNode::tag().name() );
  STKUNIT_EXPECT_TRUE( NULL != stk::mesh::QuadratureTag::tag().name() );
  STKUNIT_EXPECT_TRUE( NULL != stk::mesh::BasisTag::tag().name() );
}

void BoxMeshFixture::fill_mesh()
{
  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = 8 ;

  for ( int iz = 0 ; iz < 3 ; ++iz ) {
  for ( int iy = 0 ; iy < 3 ; ++iy ) {
  for ( int ix = 0 ; ix < 3 ; ++ix ) {
    m_node_id[iz][iy][ix] = 1 + ix + 3 * ( iy + 3 * iz );
    m_node_coord[iz][iy][ix][0] = ix ;
    m_node_coord[iz][iy][ix][1] = iy ;
    m_node_coord[iz][iy][ix][2] = - iz ;
  }
  }
  }

  const unsigned beg_elem = ( num_elems * p_rank ) / p_size ;
  const unsigned end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  m_bulk_data.modification_begin();

  unsigned elem = 0 ;
  stk::mesh::EntityId face_id = 0;

  for ( int iz = 0 ; iz < 2 ; ++iz ) {
  for ( int iy = 0 ; iy < 2 ; ++iy ) {
  for ( int ix = 0 ; ix < 2 ; ++ix , ++elem ) {
    if ( beg_elem <= elem && elem < end_elem ) {
      stk::mesh::EntityId elem_id = 1 + ix + 3 * ( iy + 3 * iz );
      stk::mesh::EntityId elem_node[8] ;
      Scalar * node_coord[8] ;
      face_id = elem_id;

      elem_node[0] = m_node_id[iz  ][iy  ][ix  ] ;
      elem_node[1] = m_node_id[iz  ][iy  ][ix+1] ;
      elem_node[2] = m_node_id[iz+1][iy  ][ix+1] ;
      elem_node[3] = m_node_id[iz+1][iy  ][ix  ] ;
      elem_node[4] = m_node_id[iz  ][iy+1][ix  ] ;
      elem_node[5] = m_node_id[iz  ][iy+1][ix+1] ;
      elem_node[6] = m_node_id[iz+1][iy+1][ix+1] ;
      elem_node[7] = m_node_id[iz+1][iy+1][ix  ] ;

      node_coord[0] = m_node_coord[iz  ][iy  ][ix  ] ;
      node_coord[1] = m_node_coord[iz  ][iy  ][ix+1] ;
      node_coord[2] = m_node_coord[iz+1][iy  ][ix+1] ;
      node_coord[3] = m_node_coord[iz+1][iy  ][ix  ] ;
      node_coord[4] = m_node_coord[iz  ][iy+1][ix  ] ;
      node_coord[5] = m_node_coord[iz  ][iy+1][ix+1] ;
      node_coord[6] = m_node_coord[iz+1][iy+1][ix+1] ;
      node_coord[7] = m_node_coord[iz+1][iy+1][ix  ] ;

      stk::mesh::Entity& element = stk::mesh::declare_element( m_bulk_data, m_elem_block, elem_id, elem_node);

      stk::mesh::PairIterRelation rel = element.relations(stk::mesh::Node);

      //Test0.1 on declare_element_side
      {
        int ok = 0 ;
        try {

          stk::mesh::Entity & face = stk::mesh::declare_element_side( m_bulk_data, face_id, element, 0, &m_elem_block);
          stk::mesh::PairIterRelation rel2 = face.relations(stk::mesh::Node);

        }
        catch( const std::exception & x ) {
          ok = 1 ;
          std::cout << "UnitTestBoxMeshFixture CORRECTLY caught error for : "
                    << x.what()
                    << std::endl ;
        }

        if ( ! ok ) {
          throw std::runtime_error("test 0.1 UnitTestBoxMeshFixture FAILED to catch error for stk::mesh::declare_element_side");
        }
      }

      for(int i=0; i< 8 ; ++i ) {
        stk::mesh::Entity& node = *rel[i].entity();
        Scalar* data = stk::mesh::field_data( m_coord_field, node);
        data[0] = node_coord[i][0];
        data[1] = node_coord[i][1];
        data[2] = node_coord[i][2];
      }

      //Test1 on declare_element_side for local side id no > no of sides
      {
        int ok = 0 ;
        try {
          stk::mesh::Entity& elementSide = stk::mesh::declare_element_side(m_bulk_data, 6, element, 9, &m_elem_block);
          stk::mesh::PairIterRelation rel2 = elementSide.relations(stk::mesh::Node);
        }
        catch( const std::exception & x ) {
          ok = 1 ;
          std::cout << "UnitTestBoxMeshFixture CORRECTLY caught error for : "
                    << x.what()
                    << std::endl ;
        }
        if ( ! ok ) {
          throw std::runtime_error("UnitTestBoxMeshFixture FAILED to catch error for stk::mesh::declare_element_side");
        }
      }
    }
  }
  }
  }

  m_bulk_data.modification_end();

  for ( int iz = 0 ; iz < 3 ; ++iz ) {
  for ( int iy = 0 ; iy < 3 ; ++iy ) {
  for ( int ix = 0 ; ix < 3 ; ++ix ) {
    // Find node
    stk::mesh::EntityId node_id = 1 + ix + 3 * ( iy + 3 * iz );
    m_nodes[iz][iy][ix] = m_bulk_data.get_entity( stk::mesh::Node , node_id );
  }
  }
  }

  for ( int iz = 0 ; iz < 2 ; ++iz ) {
  for ( int iy = 0 ; iy < 2 ; ++iy ) {
  for ( int ix = 0 ; ix < 2 ; ++ix , ++elem ) {
    stk::mesh::EntityId elem_id = 1 + ix + 3 * ( iy + 3 * iz );
    // Find element
    m_elems[iz][iy][ix] = m_bulk_data.get_entity( stk::mesh::Element , elem_id );
  }
  }
  }

  //Tests on declare_element and set_cell_topology
  const unsigned beg_elem2 = ( num_elems * p_rank ) / p_size ;
  const unsigned end_elem2 = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  unsigned elem2 = 0 ;

  for ( int iz2 = 0 ; iz2 < 2 ; ++iz2 ) {
  for ( int iy2 = 0 ; iy2 < 2 ; ++iy2 ) {
  for ( int ix2 = 0 ; ix2 < 2 ; ++ix2 , ++elem2 ) {
    if ( beg_elem2 <= elem2 && elem2 < end_elem2 ) {
      stk::mesh::EntityId elem_id2 = 1 + ix2 + 3 * ( iy2 + 3 * iz2 );
      stk::mesh::EntityId elem_node2[8] ;
      Scalar * node_coord2[8] ;

      elem_node2[0] = m_node_id[iz2  ][iy2  ][ix2  ] ;
      elem_node2[1] = m_node_id[iz2  ][iy2  ][ix2+1] ;
      elem_node2[2] = m_node_id[iz2+1][iy2  ][ix2+1] ;
      elem_node2[3] = m_node_id[iz2+1][iy2  ][ix2  ] ;
      elem_node2[4] = m_node_id[iz2  ][iy2+1][ix2  ] ;
      elem_node2[5] = m_node_id[iz2  ][iy2+1][ix2+1] ;
      elem_node2[6] = m_node_id[iz2+1][iy2+1][ix2+1] ;
      elem_node2[7] = m_node_id[iz2+1][iy2+1][ix2  ] ;

      node_coord2[0] = m_node_coord[iz2  ][iy2  ][ix2  ] ;
      node_coord2[1] = m_node_coord[iz2  ][iy2  ][ix2+1] ;
      node_coord2[2] = m_node_coord[iz2+1][iy2  ][ix2+1] ;
      node_coord2[3] = m_node_coord[iz2+1][iy2  ][ix2  ] ;
      node_coord2[4] = m_node_coord[iz2  ][iy2+1][ix2  ] ;
      node_coord2[5] = m_node_coord[iz2  ][iy2+1][ix2+1] ;
      node_coord2[6] = m_node_coord[iz2+1][iy2+1][ix2+1] ;
      node_coord2[7] = m_node_coord[iz2+1][iy2+1][ix2  ] ;

      //Test2 on declare_element ok
      {
        int ok = 0 ;
        try {

          stk::mesh::Entity& element2 = stk::mesh::declare_element(m_bulk_data, m_elem_block2, elem_id2, elem_node2);
          stk::mesh::PairIterRelation rel2 = element2.relations(stk::mesh::Node);

        }
        catch( const std::exception & x ) {
          ok = 1 ;
          std::cout << "UnitTestBoxMeshFixture CORRECTLY caught error for : "
                    << x.what()
                    << std::endl ;
        }
        if ( ! ok ) {
          throw std::runtime_error("UnitTestBoxMeshFixture FAILED to catch error for stk::mesh::declare_element");
        }
      }
    }
  }
  }
  }

  //Test3 on set_cell_topology ok
  {
    int ok = 0 ;
    try {
      set_cell_topology( m_elem_block2, NULL );
    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestBoxMeshFixture CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }
    if ( ! ok ) {
      throw std::runtime_error("UnitTestBoxMeshFixture FAILED to catch error for stk::mesh::set_cell_topology for a part and singleton = NULL");
    }
  }
}
