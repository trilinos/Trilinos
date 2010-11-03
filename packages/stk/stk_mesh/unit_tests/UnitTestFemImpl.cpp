/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fem/TopologicalMetaData.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stdexcept>

using stk::mesh::MetaData;
using stk::mesh::TopologicalMetaData;

namespace {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * Really just a meta-data, topological-meta-data pair. They share the same
 * spatial dimension and the meta-data is associated with the
 * topological-meta-data.
 */
class MetaTopoPair {
public:
  MetaData            m_meta_data ;
  TopologicalMetaData m_topo_data ;

  MetaTopoPair( unsigned spatial_dimension )
    : m_meta_data(TopologicalMetaData::entity_rank_names(spatial_dimension))
    , m_topo_data(m_meta_data, spatial_dimension)
    {}

  int get_entity_rank( const CellTopologyData * top ) const
    { return m_topo_data.get_entity_rank( top ); }
};

STKUNIT_UNIT_TEST( UnitTestTopologicalMetaData , entity_rank )
{
  MetaTopoPair test1D( 1u );
  MetaTopoPair test2D( 2u );
  MetaTopoPair test3D( 3u );

  STKUNIT_EXPECT_EQUAL( test1D.m_topo_data.spatial_dimension , 1u );
  STKUNIT_EXPECT_EQUAL( test1D.m_topo_data.node_rank , 0u );
  STKUNIT_EXPECT_EQUAL( test1D.m_topo_data.edge_rank , 0u );
  STKUNIT_EXPECT_EQUAL( test1D.m_topo_data.side_rank , 0u );
  STKUNIT_EXPECT_EQUAL( test1D.m_topo_data.element_rank , 1u );
  STKUNIT_EXPECT_EQUAL( test1D.m_topo_data.patch_rank , 2u );

  STKUNIT_EXPECT_EQUAL( test2D.m_topo_data.spatial_dimension , 2u );
  STKUNIT_EXPECT_EQUAL( test2D.m_topo_data.node_rank , 0u );
  STKUNIT_EXPECT_EQUAL( test2D.m_topo_data.edge_rank , 1u );
  STKUNIT_EXPECT_EQUAL( test2D.m_topo_data.side_rank , 1u );
  STKUNIT_EXPECT_EQUAL( test2D.m_topo_data.element_rank , 2u );
  STKUNIT_EXPECT_EQUAL( test2D.m_topo_data.patch_rank , 3u );

  STKUNIT_EXPECT_EQUAL( test3D.m_topo_data.spatial_dimension , 3u );
  STKUNIT_EXPECT_EQUAL( test3D.m_topo_data.node_rank , 0u );
  STKUNIT_EXPECT_EQUAL( test3D.m_topo_data.edge_rank , 1u );
  STKUNIT_EXPECT_EQUAL( test3D.m_topo_data.side_rank , 2u );
  STKUNIT_EXPECT_EQUAL( test3D.m_topo_data.element_rank , 3u );
  STKUNIT_EXPECT_EQUAL( test3D.m_topo_data.patch_rank , 4u );
}


STKUNIT_UNIT_TEST( UnitTestTopologicalMetaData , cellTopology )
{
  const CellTopologyData * node  = shards::getCellTopologyData< shards::Node >();

  const CellTopologyData * line2 = shards::getCellTopologyData< shards::Line<2> >();
  const CellTopologyData * line3 = shards::getCellTopologyData< shards::Line<3> >();

  const CellTopologyData * tri3 = shards::getCellTopologyData< shards::Triangle<3> >();
  const CellTopologyData * tri6 = shards::getCellTopologyData< shards::Triangle<6> >();
  const CellTopologyData * tri4 = shards::getCellTopologyData< shards::Triangle<4> >();

  const CellTopologyData * quad4 = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  const CellTopologyData * quad8 = shards::getCellTopologyData< shards::Quadrilateral<8> >();
  const CellTopologyData * quad9 = shards::getCellTopologyData< shards::Quadrilateral<9> >();

  const CellTopologyData * tet4  = shards::getCellTopologyData< shards::Tetrahedron<4> >();
  const CellTopologyData * tet10 = shards::getCellTopologyData< shards::Tetrahedron<10> >();
  const CellTopologyData * tet8  = shards::getCellTopologyData< shards::Tetrahedron<8> >();

  const CellTopologyData * pyr5  = shards::getCellTopologyData< shards::Pyramid<5> >();
  const CellTopologyData * pyr13 = shards::getCellTopologyData< shards::Pyramid<13> >();
  const CellTopologyData * pyr14 = shards::getCellTopologyData< shards::Pyramid<14> >();

  const CellTopologyData * wedge6  = shards::getCellTopologyData< shards::Wedge<6> >();
  const CellTopologyData * wedge15 = shards::getCellTopologyData< shards::Wedge<15> >();
  const CellTopologyData * wedge18 = shards::getCellTopologyData< shards::Wedge<18> >();

  const CellTopologyData * hex8  = shards::getCellTopologyData< shards::Hexahedron<8> >();
  const CellTopologyData * hex20 = shards::getCellTopologyData< shards::Hexahedron<20> >();
  const CellTopologyData * hex27 = shards::getCellTopologyData< shards::Hexahedron<27> >();


  const CellTopologyData * particle = shards::getCellTopologyData< shards::Particle >();

  const CellTopologyData * beam2 = shards::getCellTopologyData< shards::Beam<2> >();
  const CellTopologyData * beam3 = shards::getCellTopologyData< shards::Beam<3> >();

  const CellTopologyData * shellLine2 = shards::getCellTopologyData< shards::ShellLine<2> >();
  const CellTopologyData * shellLine3 = shards::getCellTopologyData< shards::ShellLine<3> >();

  const CellTopologyData * shellTri3 = shards::getCellTopologyData< shards::ShellTriangle<3> >();
  const CellTopologyData * shellTri6 = shards::getCellTopologyData< shards::ShellTriangle<6> >();

  const CellTopologyData * shellQuad4 = shards::getCellTopologyData< shards::ShellQuadrilateral<4> >();
  const CellTopologyData * shellQuad8 = shards::getCellTopologyData< shards::ShellQuadrilateral<8> >();
  const CellTopologyData * shellQuad9 = shards::getCellTopologyData< shards::ShellQuadrilateral<9> >();


  for ( unsigned spatial_dimension = 1 ;
        spatial_dimension < 4 ; ++spatial_dimension ) {

    MetaTopoPair test( spatial_dimension );

    STKUNIT_EXPECT_EQUAL( 0 , test.get_entity_rank( node ) );
    STKUNIT_EXPECT_EQUAL( 1 , test.get_entity_rank( line2 ) );
    STKUNIT_EXPECT_EQUAL( 1 , test.get_entity_rank( line3 ) );
    STKUNIT_EXPECT_EQUAL( (int) spatial_dimension , test.get_entity_rank( particle ) );
    if ( 1 < spatial_dimension ) {
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( tri3 ) );
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( tri6 ) );
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( tri4 ) );
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( quad4 ) );
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( quad8 ) );
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( quad9 ) );

      STKUNIT_EXPECT_EQUAL( (int) spatial_dimension , test.get_entity_rank( beam2 ) );
      STKUNIT_EXPECT_EQUAL( (int) spatial_dimension , test.get_entity_rank( beam3 ) );
    }
    else {
      STKUNIT_ASSERT_THROW( test.get_entity_rank( tri3 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( tri6 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( tri4 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( quad4 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( quad8 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( quad9 ) , std::runtime_error );

      STKUNIT_ASSERT_THROW( test.get_entity_rank( beam2 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( beam3 ) , std::runtime_error );
    }

    if ( 2 == spatial_dimension ) {
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( shellLine2 ) );
      STKUNIT_EXPECT_EQUAL( 2 , test.get_entity_rank( shellLine3 ) );
    }
    else {
      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellLine2 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellLine3 ) , std::runtime_error );
    }

    if ( 2 < spatial_dimension ) {
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( tet4 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( tet10 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( tet8 ) );

      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( pyr5 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( pyr13 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( pyr14 ) );

      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( wedge6 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( wedge15 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( wedge18 ) );

      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( hex8 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( hex20 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( hex27 ) );

      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( shellTri3 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( shellTri6 ) );

      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( shellQuad4 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( shellQuad8 ) );
      STKUNIT_EXPECT_EQUAL( 3 , test.get_entity_rank( shellQuad9 ) );
    }
    else {
      STKUNIT_ASSERT_THROW( test.get_entity_rank( tet4 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( tet10 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( tet8 ) , std::runtime_error );

      STKUNIT_ASSERT_THROW( test.get_entity_rank( pyr5 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( pyr13 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( pyr14 ) , std::runtime_error );

      STKUNIT_ASSERT_THROW( test.get_entity_rank( wedge6 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( wedge15 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( wedge18 ) , std::runtime_error );

      STKUNIT_ASSERT_THROW( test.get_entity_rank( hex8 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( hex20 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( hex27 ) , std::runtime_error );

      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellTri3 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellTri6 ) , std::runtime_error );

      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellQuad4 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellQuad8 ) , std::runtime_error );
      STKUNIT_ASSERT_THROW( test.get_entity_rank( shellQuad9 ) , std::runtime_error );
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestTopologicalMetaData , part_supersets )
{
  static const char method[] = "UnitTestPartCellTopologyMap::part_supersets" ;
  const unsigned spatial_dimension = 3 ;

  MetaTopoPair test( spatial_dimension );

  const CellTopologyData * const top_hex8 = shards::getCellTopologyData< shards::Hexahedron<8> >();
  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  stk::mesh::Part & hex8 = test.m_topo_data.declare_part( top_hex8->name , top_hex8 );
  stk::mesh::Part & tet4 = test.m_topo_data.declare_part( top_tet4->name , top_tet4 );

  stk::mesh::Part & part_hex_1 = test.m_meta_data.declare_part( "block_1" , spatial_dimension );
  stk::mesh::Part & part_hex_2 = test.m_meta_data.declare_part( "block_2" , spatial_dimension );
  stk::mesh::Part & part_hex_3 = test.m_meta_data.declare_part( "block_3" , spatial_dimension );
  stk::mesh::Part & part_tet_1 = test.m_meta_data.declare_part( "block_4" , spatial_dimension );
  stk::mesh::Part & part_tet_2 = test.m_meta_data.declare_part( "block_5" , spatial_dimension );

  test.m_meta_data.declare_part_subset( hex8 , part_hex_1 );
  test.m_meta_data.declare_part_subset( hex8 , part_hex_2 );
  test.m_meta_data.declare_part_subset( hex8 , part_hex_3 );

  test.m_meta_data.declare_part_subset( tet4 , part_tet_1 );
  test.m_meta_data.declare_part_subset( tet4 , part_tet_2 );

  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::TopologicalMetaData::get_cell_topology( hex8 , method ) );
  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::TopologicalMetaData::get_cell_topology( part_hex_1 , method ) );
  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::TopologicalMetaData::get_cell_topology( part_hex_2 , method ) );
  STKUNIT_EXPECT_EQUAL( top_hex8 , stk::mesh::TopologicalMetaData::get_cell_topology( part_hex_3 , method ) );

  STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::TopologicalMetaData::get_cell_topology( tet4 , method ) );
  STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::TopologicalMetaData::get_cell_topology( part_tet_1 , method ) );
  STKUNIT_EXPECT_EQUAL( top_tet4 , stk::mesh::TopologicalMetaData::get_cell_topology( part_tet_2 , method ) );

}

STKUNIT_UNIT_TEST( UnitTestPartCellTopologyMap , errors )
{
  const unsigned spatial_dimension = 2 ;

  MetaTopoPair test( spatial_dimension );

  const CellTopologyData * const top_quad4 = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  const CellTopologyData * const top_tri3 = shards::getCellTopologyData< shards::Triangle<3> >();
  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  STKUNIT_ASSERT_THROW( test.m_topo_data.declare_cell_topology( top_quad4 , 0 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( test.m_topo_data.declare_cell_topology( top_quad4 , 1 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( test.m_topo_data.declare_cell_topology( top_quad4 , 3 ) , std::runtime_error );

  STKUNIT_ASSERT_THROW( test.m_topo_data.declare_cell_topology( top_tri3 , 3 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( test.m_topo_data.declare_cell_topology( top_tet4 , 3 ) , std::runtime_error );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*

class FiniteElementMesh {
public:

  typedef stk::mesh::MetaData            MetaData ;
  typedef stk::mesh::TopologicalMetaData TopologicalMetaData ;
  typedef stk::mesh::BulkData            BulkData ;

  MetaData m_meta_data ;
  TopologicalMetaData m_topo_data ;
  BulkData     bulkData ;

  FiniteElementMesh( unsigned spatial_dimension ,
                     stk::ParallelMachine machine )
    : m_meta_data( TopologicalMetaData::entity_rank_names( spatial_dimension ) )
    , m_topo_data( metaData , spatial_dimension )
    , bulkData( metaData , machine )
    {}
};

STKUNIT_UNIT_TEST( UnitTestPartCellTopologyMap , bucket )
{
  stk::ParallelMachine machine = MPI_COMM_WORLD ;

  const CellTopologyData * const top_tet4 =
    shards::getCellTopologyData< shards::Tetrahedron<4> >();

  const unsigned spatial_dimension = 3 ;

  FiniteElementMesh mesh( spatial_dimension , machine );

  stk::mesh::Part & block =
    mesh.m_topo_data.declare_part< shards::Tetrahedron<4> >( "block_1" );

  mesh.m_meta_data.commit();

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

*/

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace

