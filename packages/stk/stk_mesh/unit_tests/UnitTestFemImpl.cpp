/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/femImpl/FiniteElementMeshImpl.hpp>
#include <stk_mesh/femImpl/PartCellTopologyMap.hpp>

#include <stk_mesh/fem/FiniteElementMesh.hpp>


namespace {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

STKUNIT_UNIT_TEST( UnitTestPartCellTopologyMap , construct )
{
  for ( unsigned spatial_dimension = 1 ; spatial_dimension < 4 ; ++spatial_dimension ) {

    stk::mesh::MetaData meta_data( stk::mesh::impl::finite_element_mesh_entity_rank_names( spatial_dimension ) );

    stk::mesh::impl::PartCellTopologyMap part_cell( meta_data , spatial_dimension );

    // Ordinary topologies:

    stk::mesh::Part * nil   = 0 ;
    stk::mesh::Part * node  = part_cell.get_part( * shards::getCellTopologyData< shards::Node >() , NULL );

    stk::mesh::Part * line2 = part_cell.get_part( * shards::getCellTopologyData< shards::Line<2> >() , NULL );
    stk::mesh::Part * line3 = part_cell.get_part( * shards::getCellTopologyData< shards::Line<3> >() , NULL );

    stk::mesh::Part * tri3 = part_cell.get_part( * shards::getCellTopologyData< shards::Triangle<3> >() , NULL );
    stk::mesh::Part * tri6 = part_cell.get_part( * shards::getCellTopologyData< shards::Triangle<6> >() , NULL );
    stk::mesh::Part * tri4 = part_cell.get_part( * shards::getCellTopologyData< shards::Triangle<4> >() , NULL );

    stk::mesh::Part * quad4 = part_cell.get_part( * shards::getCellTopologyData< shards::Quadrilateral<4> >() , NULL );
    stk::mesh::Part * quad8 = part_cell.get_part( * shards::getCellTopologyData< shards::Quadrilateral<8> >() , NULL );
    stk::mesh::Part * quad9 = part_cell.get_part( * shards::getCellTopologyData< shards::Quadrilateral<9> >() , NULL );

    stk::mesh::Part * tet4  = part_cell.get_part( * shards::getCellTopologyData< shards::Tetrahedron<4> >() , NULL );
    stk::mesh::Part * tet10 = part_cell.get_part( * shards::getCellTopologyData< shards::Tetrahedron<10> >() , NULL );
    stk::mesh::Part * tet8  = part_cell.get_part( * shards::getCellTopologyData< shards::Tetrahedron<8> >() , NULL );

    stk::mesh::Part * pyr5  = part_cell.get_part( * shards::getCellTopologyData< shards::Pyramid<5> >() , NULL );
    stk::mesh::Part * pyr13 = part_cell.get_part( * shards::getCellTopologyData< shards::Pyramid<13> >() , NULL );
    stk::mesh::Part * pyr14 = part_cell.get_part( * shards::getCellTopologyData< shards::Pyramid<14> >() , NULL );

    stk::mesh::Part * wedge6  = part_cell.get_part( * shards::getCellTopologyData< shards::Wedge<6> >() , NULL );
    stk::mesh::Part * wedge15 = part_cell.get_part( * shards::getCellTopologyData< shards::Wedge<15> >() , NULL );
    stk::mesh::Part * wedge18 = part_cell.get_part( * shards::getCellTopologyData< shards::Wedge<18> >() , NULL );

    stk::mesh::Part * hex8  = part_cell.get_part( * shards::getCellTopologyData< shards::Hexahedron<8> >() , NULL );
    stk::mesh::Part * hex20 = part_cell.get_part( * shards::getCellTopologyData< shards::Hexahedron<20> >() , NULL );
    stk::mesh::Part * hex27 = part_cell.get_part( * shards::getCellTopologyData< shards::Hexahedron<27> >() , NULL );


    stk::mesh::Part * particle = part_cell.get_part( * shards::getCellTopologyData< shards::Particle >() , NULL );

    stk::mesh::Part * beam2 = part_cell.get_part( * shards::getCellTopologyData< shards::Beam<2> >() , NULL );
    stk::mesh::Part * beam3 = part_cell.get_part( * shards::getCellTopologyData< shards::Beam<3> >() , NULL );

    stk::mesh::Part * shellLine2 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellLine<2> >() , NULL );
    stk::mesh::Part * shellLine3 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellLine<3> >() , NULL );

    stk::mesh::Part * shellTri3 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellTriangle<3> >() , NULL );
    stk::mesh::Part * shellTri6 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellTriangle<6> >() , NULL );

    stk::mesh::Part * shellQuad4 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellQuadrilateral<4> >() , NULL );
    stk::mesh::Part * shellQuad8 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellQuadrilateral<8> >() , NULL );
    stk::mesh::Part * shellQuad9 = part_cell.get_part( * shards::getCellTopologyData< shards::ShellQuadrilateral<9> >() , NULL );


    STKUNIT_EXPECT_TRUE( node );
    STKUNIT_EXPECT_TRUE( line2 );
    STKUNIT_EXPECT_TRUE( line3 );
    STKUNIT_EXPECT_TRUE( particle );

    STKUNIT_EXPECT_EQUAL( 0u , node->primary_entity_rank() );
    STKUNIT_EXPECT_EQUAL( 1u , line2->primary_entity_rank() );
    STKUNIT_EXPECT_EQUAL( 1u , line3->primary_entity_rank() );
    STKUNIT_EXPECT_EQUAL( spatial_dimension , particle->primary_entity_rank() );

    if ( 1 < spatial_dimension ) {
      STKUNIT_EXPECT_TRUE( tri3 );
      STKUNIT_EXPECT_TRUE( tri6 );
      STKUNIT_EXPECT_TRUE( tri4 );
      STKUNIT_EXPECT_TRUE( quad4 );
      STKUNIT_EXPECT_TRUE( quad8 );
      STKUNIT_EXPECT_TRUE( quad9 );

      STKUNIT_EXPECT_TRUE( beam2 );
      STKUNIT_EXPECT_TRUE( beam3 );

      STKUNIT_EXPECT_EQUAL( 2u , tri3->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 2u , tri6->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 2u , tri4->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 2u , quad4->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 2u , quad8->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 2u , quad9->primary_entity_rank() );

      STKUNIT_EXPECT_EQUAL( spatial_dimension , beam2->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( spatial_dimension , beam3->primary_entity_rank() );
    }
    else {
      STKUNIT_EXPECT_EQUAL( nil , tri3 );
      STKUNIT_EXPECT_EQUAL( nil , tri6 );
      STKUNIT_EXPECT_EQUAL( nil , tri4 );
      STKUNIT_EXPECT_EQUAL( nil , quad4 );
      STKUNIT_EXPECT_EQUAL( nil , quad8 );
      STKUNIT_EXPECT_EQUAL( nil , quad9 );
      STKUNIT_EXPECT_EQUAL( nil , beam2 );
      STKUNIT_EXPECT_EQUAL( nil , beam3 );
    }

    if ( 2 == spatial_dimension ) {
      STKUNIT_EXPECT_EQUAL( 2u , shellLine2->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 2u , shellLine3->primary_entity_rank() );
    }
    else {
      STKUNIT_EXPECT_EQUAL( nil , shellLine2 );
      STKUNIT_EXPECT_EQUAL( nil , shellLine3 );
    }

    if ( 2 < spatial_dimension ) {
      STKUNIT_EXPECT_TRUE( tet4 );
      STKUNIT_EXPECT_TRUE( tet10 );
      STKUNIT_EXPECT_TRUE( tet8 );

      STKUNIT_EXPECT_TRUE( pyr5 );
      STKUNIT_EXPECT_TRUE( pyr13 );
      STKUNIT_EXPECT_TRUE( pyr14 );

      STKUNIT_EXPECT_TRUE( wedge6 );
      STKUNIT_EXPECT_TRUE( wedge15 );
      STKUNIT_EXPECT_TRUE( wedge18 );

      STKUNIT_EXPECT_TRUE( hex8 );
      STKUNIT_EXPECT_TRUE( hex20 );
      STKUNIT_EXPECT_TRUE( hex27 );

      STKUNIT_EXPECT_TRUE( shellTri3 );
      STKUNIT_EXPECT_TRUE( shellTri6 );

      STKUNIT_EXPECT_TRUE( shellQuad4 );
      STKUNIT_EXPECT_TRUE( shellQuad8 );
      STKUNIT_EXPECT_TRUE( shellQuad9 );

      STKUNIT_EXPECT_EQUAL( 3u , tet4->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , tet10->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , tet8->primary_entity_rank() );

      STKUNIT_EXPECT_EQUAL( 3u , pyr5->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , pyr13->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , pyr14->primary_entity_rank() );

      STKUNIT_EXPECT_EQUAL( 3u , wedge6->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , wedge15->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , wedge18->primary_entity_rank() );

      STKUNIT_EXPECT_EQUAL( 3u , hex8->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , hex20->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , hex27->primary_entity_rank() );

      STKUNIT_EXPECT_EQUAL( 3u , shellTri3->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , shellTri6->primary_entity_rank() );

      STKUNIT_EXPECT_EQUAL( 3u , shellQuad4->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , shellQuad8->primary_entity_rank() );
      STKUNIT_EXPECT_EQUAL( 3u , shellQuad9->primary_entity_rank() );
    }
    else {
      STKUNIT_EXPECT_EQUAL( nil , tet4 );
      STKUNIT_EXPECT_EQUAL( nil , tet10 );
      STKUNIT_EXPECT_EQUAL( nil , tet8 );

      STKUNIT_EXPECT_EQUAL( nil , pyr5 );
      STKUNIT_EXPECT_EQUAL( nil , pyr13 );
      STKUNIT_EXPECT_EQUAL( nil , pyr14 );

      STKUNIT_EXPECT_EQUAL( nil , wedge6 );
      STKUNIT_EXPECT_EQUAL( nil , wedge15 );
      STKUNIT_EXPECT_EQUAL( nil , wedge18 );

      STKUNIT_EXPECT_EQUAL( nil , hex8 );
      STKUNIT_EXPECT_EQUAL( nil , hex20 );
      STKUNIT_EXPECT_EQUAL( nil , hex27 );

      STKUNIT_EXPECT_EQUAL( nil , shellTri3 );
      STKUNIT_EXPECT_EQUAL( nil , shellTri6 );

      STKUNIT_EXPECT_EQUAL( nil , shellQuad4 );
      STKUNIT_EXPECT_EQUAL( nil , shellQuad8 );
      STKUNIT_EXPECT_EQUAL( nil , shellQuad9 );
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestPartCellTopologyMap , part_supersets )
{
  static const char method[] = "UnitTestPartCellTopologyMap::part_supersets" ;
  const unsigned spatial_dimension = 3 ;

  stk::mesh::MetaData meta_data( stk::mesh::impl::finite_element_mesh_entity_rank_names( spatial_dimension ) );

  stk::mesh::impl::PartCellTopologyMap part_cell( meta_data , spatial_dimension );

  const CellTopologyData * const top_hex8 = shards::getCellTopologyData< shards::Hexahedron<8> >();
  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  stk::mesh::Part & hex8 = * part_cell.get_part( *top_hex8 , method );
  stk::mesh::Part & tet4 = * part_cell.get_part( *top_tet4 , method );

  stk::mesh::Part & part_hex_1 = meta_data.declare_part( "block_1" , 3 );
  stk::mesh::Part & part_hex_2 = meta_data.declare_part( "block_2" , 3 );
  stk::mesh::Part & part_hex_3 = meta_data.declare_part( "block_3" , 3 );
  stk::mesh::Part & part_tet_1 = meta_data.declare_part( "block_4" , 3 );
  stk::mesh::Part & part_tet_2 = meta_data.declare_part( "block_5" , 3 );

  meta_data.declare_part_subset( hex8 , part_hex_1 );
  meta_data.declare_part_subset( hex8 , part_hex_2 );
  meta_data.declare_part_subset( hex8 , part_hex_3 );

  meta_data.declare_part_subset( tet4 , part_tet_1 );
  meta_data.declare_part_subset( tet4 , part_tet_2 );

  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( hex8 , method ) );
  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( part_hex_1 , method ) );
  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( part_hex_2 , method ) );
  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( part_hex_3 , method ) );

  STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( tet4 , method ) );
  STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( part_tet_1 , method ) );
  STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::impl::PartCellTopologyMap::get_cell_topology( part_tet_2 , method ) );

}

STKUNIT_UNIT_TEST( UnitTestPartCellTopologyMap , errors )
{
  const unsigned spatial_dimension = 2 ;

  stk::mesh::MetaData meta_data( stk::mesh::impl::finite_element_mesh_entity_rank_names( spatial_dimension ) );

  stk::mesh::impl::PartCellTopologyMap part_cell( meta_data , spatial_dimension );

  const CellTopologyData * const top_quad4 = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  const CellTopologyData * const top_tri3 = shards::getCellTopologyData< shards::Triangle<3> >();
  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  STKUNIT_ASSERT_THROW( part_cell.declare_part( *top_quad4 , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( part_cell.declare_part( *top_quad4 , 1 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( part_cell.declare_part( *top_quad4 , 3 ) , std::runtime_error );

  STKUNIT_ASSERT_THROW( part_cell.declare_part( *top_tri3 , 8 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( part_cell.declare_part( *top_tet4 , 2 ) , std::runtime_error );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

STKUNIT_UNIT_TEST( UnitTestPartCellTopologyMap , bucket )
{
  stk::ParallelMachine machine = MPI_COMM_WORLD ;

  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  const unsigned spatial_dimension = 3 ;

  stk::mesh::FiniteElementMesh<> mesh( spatial_dimension , machine );

  stk::mesh::Part & block = mesh.declare_part< shards::Tetrahedron<4> >( "block_1" );

  mesh.metaData.commit();

  mesh.bulkData.modification_begin();

  stk::mesh::Entity * element = NULL ;

  if ( mesh.bulkData.parallel_rank() == 0 ) {
    int node_ids[4] = { 1 , 2 , 3 , 4 };

    element = & declare_element( mesh.bulkData , block , 1 , node_ids );
  }

  mesh.bulkData.modification_end();

  if ( element ) {
    STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::get_cell_topology( *element ) );
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace

