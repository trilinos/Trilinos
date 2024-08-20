// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#define HAVE_SHARDS_DEBUG

#include <iostream>
#include <stdexcept>
#include <Shards_CellTopology.hpp>
#include <Shards_BasicTopologies.hpp>

namespace {

#define REQUIRE( S ) \
  if ( ! ( S ) ) { throw std::runtime_error(std::string(#S)); }

#define REQUIRE_EX( S ) \
  { bool flag = true ; \
    try { S ; } catch( const std::exception & ) { flag = false ; } \
    if ( flag ) { throw std::runtime_error(std::string(#S)); } \
  }

template< class Traits , unsigned Dim , unsigned Ord >
void test_subcell( const shards::CellTopology & parent )
{
  typedef typename Traits::template subcell< Dim , Ord > Subcell ;
  typedef typename Subcell::topology SubcellTraits ;

  std::cout << "    subcell(" << Dim << "," << Ord << ") = "
            << parent.getCellTopologyData(Dim,Ord)->name << std::endl ;

  REQUIRE( SubcellTraits::key          == parent.getKey(Dim,Ord) )
  REQUIRE( SubcellTraits::side_count   == parent.getSideCount(Dim,Ord) )
  REQUIRE( SubcellTraits::node_count   == parent.getNodeCount(Dim,Ord) )
  REQUIRE( SubcellTraits::edge_count   == parent.getEdgeCount(Dim,Ord) )
  REQUIRE( SubcellTraits::vertex_count == parent.getVertexCount(Dim,Ord) )
}


template< class Traits , unsigned Dim , unsigned Count >
struct test_all_subcell
{
  test_all_subcell( const shards::CellTopology & parent )
  {
    test_all_subcell< Traits , Dim , Count - 1 > previous( parent );
    test_subcell< Traits , Dim , Count - 1 >( parent );
  }
};

template< class Traits , unsigned Dim >
struct test_all_subcell<Traits,Dim,1>
{
  test_all_subcell( const shards::CellTopology & parent )
  {
    test_subcell< Traits , Dim , 0 >( parent );
  }
};

template< class Traits , unsigned Dim >
struct test_all_subcell<Traits,Dim,0>
{
  test_all_subcell( const shards::CellTopology & )
  { }
};

template< class Traits >
void test_cell()
{
  typedef typename Traits::template subcell< 0 > Subcell_0 ;
  typedef typename Traits::template subcell< 1 > Subcell_1 ;
  typedef typename Traits::template subcell< 2 > Subcell_2 ;
  typedef typename Traits::template subcell< 3 > Subcell_3 ;
  typedef typename Traits::template subcell< Traits::dimension > SelfSubcell ;
  typedef typename SelfSubcell::topology SelfTraits ;

  const CellTopologyData * const cell_data =
    shards::getCellTopologyData< Traits >();

  std::cout << "  Testing " << cell_data->name << std::endl ;

  const shards::CellTopology top( cell_data );

  enum { same_type = shards::SameType< typename Traits::Traits , SelfTraits >::value };

  REQUIRE( same_type );
  REQUIRE( cell_data            == top.getCellTopologyData() )
  REQUIRE( Traits::key          == top.getKey() )
  REQUIRE( Traits::side_count   == top.getSideCount() )
  REQUIRE( Traits::node_count   == top.getNodeCount() )
  REQUIRE( Traits::edge_count   == top.getEdgeCount() )
  REQUIRE( Traits::vertex_count == top.getVertexCount() )
  REQUIRE( Traits::permutation_count == top.getNodePermutationCount() )
  REQUIRE( Subcell_0::count     == top.getSubcellCount(0) )
  REQUIRE( Subcell_1::count     == top.getSubcellCount(1) )
  REQUIRE( Subcell_2::count     == top.getSubcellCount(2) )
  REQUIRE( Subcell_3::count     == top.getSubcellCount(3) )

  REQUIRE_EX( top.getSubcellCount(4) )
  REQUIRE_EX( top.getSubcellCount(5) )

  test_all_subcell< Traits , 0 , Subcell_0::count > testd0( top );
  test_all_subcell< Traits , 1 , Subcell_1::count > testd1( top );
  test_all_subcell< Traits , 2 , Subcell_2::count > testd2( top );
  test_all_subcell< Traits , 3 , Subcell_3::count > testd3( top );

  if ( 3 == cell_data->dimension ) {
    for ( unsigned i = 0 ; i < cell_data->subcell_count[2] ; ++i ) {
      const CellTopologyData * const face_top =
        cell_data->subcell[2][i].topology ;
      for ( unsigned j = 0 ; j < face_top->subcell_count[1] ; ++j ) {
        REQUIRE( 0 <= mapCellFaceEdge( cell_data , i , j ) );
      }
    }
  }
}

template< class Traits >
void test_permutation( const unsigned p )
{
  const CellTopologyData * const cell_data =
    shards::getCellTopologyData< Traits >();

  int node_expected[ Traits::node_count ];
  int node_permuted[ Traits::node_count ];

std::cout << "Test " << cell_data->name << " perm " << p << " of " << Traits::permutation_count << " = {" ;
std::cout.flush();

  REQUIRE( p < Traits::permutation_count );

  for ( unsigned i = 0 ; i < Traits::node_count ; ++i ) {
    node_expected[i] = ( i + 1 ) * 100 ;
  }

  for ( unsigned i = 0 ; i < Traits::node_count ; ++i ) {
    const unsigned ip = cell_data->permutation[p].node[i];

    std::cout << " " << ip ;
    std::cout.flush();

    REQUIRE( ip < Traits::node_count );

    node_permuted[i] = node_expected[ ip ];
  }

  std::cout << " } inverse {" ;

  for ( unsigned i = 0 ; i < Traits::node_count ; ++i ) {
    const unsigned ip = cell_data->permutation_inverse[p].node[i];

    std::cout << " " << ip ;
    std::cout.flush();
  }

  std::cout << " }" << std::endl ;
  std::cout.flush();

  for ( unsigned i = 0 ; i < Traits::node_count ; ++i ) {
    const unsigned ip = cell_data->permutation_inverse[p].node[i];
    REQUIRE( ip < Traits::node_count );
    REQUIRE( node_permuted[ip] == node_expected[i] );
  }

  const int find_p =
    shards::findPermutation( * cell_data , node_expected , node_permuted );

  REQUIRE( (int) p == find_p );
}

void local_test_cell_topology()
{
  test_cell< shards::Node >();
  test_cell< shards::Particle >();

  test_cell< shards::Line<2> >();
  test_cell< shards::ShellLine<2> >();
  test_cell< shards::ShellLine<3> >();
  test_cell< shards::Beam<2> >();
  test_cell< shards::Beam<3> >();

  test_cell< shards::Triangle<3> >();
  test_cell< shards::Triangle<6> >();
  test_cell< shards::ShellTriangle<3> >();
  test_cell< shards::ShellTriangle<6> >();

  test_cell< shards::Quadrilateral<4> >();
  test_cell< shards::Quadrilateral<8> >();
  test_cell< shards::Quadrilateral<9> >();
  test_cell< shards::ShellQuadrilateral<4> >();
  test_cell< shards::ShellQuadrilateral<8> >();
  test_cell< shards::ShellQuadrilateral<9> >();

  test_cell< shards::Pentagon<5> >();
  test_cell< shards::Hexagon<6> >();

  test_cell< shards::Tetrahedron<4> >();
  test_cell< shards::Tetrahedron<10> >();
  test_cell< shards::Tetrahedron<11> >();

  test_cell< shards::Pyramid<5> >();
  test_cell< shards::Pyramid<13> >();
  test_cell< shards::Pyramid<14> >();

  test_cell< shards::Wedge<6> >();
  test_cell< shards::Wedge<15> >();
  test_cell< shards::Wedge<18> >();

  test_cell< shards::Hexahedron<8> >();
  test_cell< shards::Hexahedron<20> >();
  test_cell< shards::Hexahedron<27> >();

  for ( unsigned i = 0 ; i < 2 ; ++i ) {
    test_permutation< shards::Line<2> >( i );
  }
  for ( unsigned i = 0 ; i < 2 ; ++i ) {
    test_permutation< shards::Line<3> >( i );
  }

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    test_permutation< shards::Triangle<3> >( i );
  }
  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    test_permutation< shards::Triangle<6> >( i );
  }
  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    test_permutation< shards::Triangle<4> >( i );
  }

  for ( unsigned i = 0 ; i < 8 ; ++i ) {
    test_permutation< shards::Quadrilateral<4> >( i );
  }
  for ( unsigned i = 0 ; i < 8 ; ++i ) {
    test_permutation< shards::Quadrilateral<8> >( i );
  }
  for ( unsigned i = 0 ; i < 8 ; ++i ) {
    test_permutation< shards::Quadrilateral<9> >( i );
  }
}

}

void test_shards_cell_topology()
{
  static const char method[] = "test_shards_cell_topology" ;

  try {
    local_test_cell_topology();
    std::cout << method << "\n" << "End Result: TEST PASSED" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << method << "\n" << "End Result: TEST FAILED: " << x.what() << std::endl ;
    throw x ;
  }
  catch( ... ) {
    std::cout << method << "\n" << "End Result: TEST FAILED: <unknown>" << std::endl ;
    throw ;
  }
}

