//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <utility>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <BoxElemPart.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

void box_partition( const unsigned global_size ,
                    const unsigned global_rank ,
                    const unsigned global_box[][2] ,
                          unsigned box[][2] )
{
  box[0][0] = global_box[0][0] ; box[0][1] = global_box[0][1] ;
  box[1][0] = global_box[1][0] ; box[1][1] = global_box[1][1] ;
  box[2][0] = global_box[2][0] ; box[2][1] = global_box[2][1] ;

  unsigned ip = 0 ;
  unsigned np = global_size ;

  while ( 1 < np ) {

    // P = [ ip + j * portion , ip + ( j + 1 ) * portion )

    unsigned jip , jup ;

    {
      const unsigned part = ( 0 == ( np % 5 ) ) ? 5 : (
                            ( 0 == ( np % 3 ) ) ? 3 : 2 );

      const unsigned portion = np / part ;

      if ( 2 < part || global_rank < ip + portion ) {
        jip = portion * size_t( double( global_rank - ip ) / double(portion) );
        jup = jip + portion ;
      }
      else {
        jip = portion ;
        jup = np ;
      }
    }

    // Choose axis with largest count:

    const unsigned nb[3] = {
      box[0][1] - box[0][0] ,
      box[1][1] - box[1][0] ,
      box[2][1] - box[2][0] };

    const unsigned axis = nb[2] > nb[1] ? ( nb[2] > nb[0] ? 2 : 0 )
                                        : ( nb[1] > nb[0] ? 1 : 0 );

    box[ axis ][1] = box[ axis ][0] + unsigned( double(nb[axis]) * ( double(jup) / double(np) ));
    box[ axis ][0] = box[ axis ][0] + unsigned( double(nb[axis]) * ( double(jip) / double(np) ));

    np = jup - jip ;
    ip = ip + jip ;
  }
}

} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

void BoxElemPart::local( const unsigned  rank ,
                               unsigned  uses_elem[][2] ,
                               unsigned  owns_node[][2] ,
                               unsigned  uses_node[][2] ) const
{
  if ( BoxElemPart::DecomposeElem == m_decompose ) {

    Kokkos::Example::box_partition( m_global_size , rank , m_global_elem_box , uses_elem );

    for ( unsigned i = 0 ; i < 3 ; ++i ) {
      owns_node[i][0] = uses_elem[i][0] ;
      owns_node[i][1] = uses_elem[i][1] + ( m_global_elem_box[i][1] == uses_elem[i][1] ? 1 : 0 );
    }
  }
  else {

    const unsigned global_vert[3][2] =
      { { 0 , m_global_elem_box[0][1] + 1 },
        { 0 , m_global_elem_box[1][1] + 1 },
        { 0 , m_global_elem_box[2][1] + 1 } };

    Kokkos::Example::box_partition( m_global_size , rank , global_vert , owns_node );

    for ( unsigned i = 0 ; i < 3 ; ++i ) {
      uses_elem[i][0] = global_vert[i][0] == owns_node[i][0] ? owns_node[i][0] : owns_node[i][0] - 1 ;
      uses_elem[i][1] = global_vert[i][1] == owns_node[i][1] ? owns_node[i][1] - 1 : owns_node[i][1] ;
    }
  }

  for ( unsigned i = 0 ; i < 3 ; ++i ) {
    uses_node[i][0] = uses_elem[i][0] ;
    uses_node[i][1] = uses_elem[i][1] + 1 ;
  }

  if ( BoxElemPart::ElemQuadratic == m_elem_order ) {
    for ( unsigned i = 0 ; i < 3 ; ++i ) {
      owns_node[i][0] = 2 * owns_node[i][0] ;
      uses_node[i][0] = 2 * uses_node[i][0] ;
      owns_node[i][1] = 2 * owns_node[i][1] - 1 ;
      uses_node[i][1] = 2 * uses_node[i][1] - 1 ;
    }
  }
}

BoxElemPart::BoxElemPart(
  const BoxElemPart::ElemOrder elem_order ,
  const BoxElemPart::Decompose decompose ,
  const unsigned global_size ,
  const unsigned global_rank ,
  const unsigned elem_nx ,
  const unsigned elem_ny ,
  const unsigned elem_nz )
{
  m_global_size = global_size ;
  m_global_rank = global_rank ;
  m_elem_order  = elem_order ;
  m_decompose   = decompose ;

  m_global_elem_box[0][0] = 0 ; m_global_elem_box[0][1] = elem_nx ;
  m_global_elem_box[1][0] = 0 ; m_global_elem_box[1][1] = elem_ny ;
  m_global_elem_box[2][0] = 0 ; m_global_elem_box[2][1] = elem_nz ;

  m_global_node_box[0][0] = 0 ; m_global_node_box[0][1] = 0 ;
  m_global_node_box[1][0] = 0 ; m_global_node_box[1][1] = 0 ;
  m_global_node_box[2][0] = 0 ; m_global_node_box[2][1] = 0 ;

  if ( ElemLinear == elem_order ) {
    m_global_node_box[0][1] = elem_nx + 1 ;
    m_global_node_box[1][1] = elem_ny + 1 ;
    m_global_node_box[2][1] = elem_nz + 1 ;
  }
  else if ( ElemQuadratic == elem_order ) {
    m_global_node_box[0][1] = 2 * elem_nx + 1 ;
    m_global_node_box[1][1] = 2 * elem_ny + 1 ;
    m_global_node_box[2][1] = 2 * elem_nz + 1 ;
  }

  //----------------------------------------

  local( m_global_rank , m_uses_elem_box , m_owns_node_box[0] , m_uses_node_box );

  const size_t global_node_count = Kokkos::Example::box_count( m_global_node_box );
  const size_t global_elem_count = Kokkos::Example::box_count( m_global_elem_box );

  //----------------------------------------

  size_t elem_count = Kokkos::Example::box_count( m_uses_elem_box );
  size_t node_count = Kokkos::Example::box_count( m_owns_node_box[0] );

  m_owns_node[0][0] = global_rank ;
  m_owns_node[0][1] = node_count ;
  m_owns_node_count = 1 ;
  m_send_node_count = 0 ;

  for ( unsigned rr = 1 ; rr < m_global_size ; ++rr ) {

    const unsigned rank = ( m_global_rank + rr ) % m_global_size ;

    unsigned elem_box[3][2] , o_node_box[3][2] , u_node_box[3][2] ;

    // Boxes for process 'rank'
    local( rank , elem_box , o_node_box , u_node_box );

    // Box that this process uses but is owned by process 'rank'
    Kokkos::Example::box_intersect( m_owns_node_box[ m_owns_node_count ] , m_uses_node_box , o_node_box );

    m_owns_node[ m_owns_node_count ][1] = Kokkos::Example::box_count( m_owns_node_box[ m_owns_node_count ] );

    if ( m_owns_node[ m_owns_node_count ][1] ) {
      m_owns_node[ m_owns_node_count ][0] = rank ;
      ++m_owns_node_count ;
    }

    // Box that this process owns and is used by process 'rank'
    Kokkos::Example::box_intersect( m_send_node_box[ m_send_node_count ] , m_owns_node_box[0] , u_node_box );

    m_send_node[ m_send_node_count ][1] = Kokkos::Example::box_count( m_send_node_box[ m_send_node_count ] );

    if ( m_send_node[ m_send_node_count ][1] ) {
      m_send_node[ m_send_node_count ][0] = rank ;
      ++m_send_node_count ;
    }

    // Error checking:

    unsigned test_box[3][2] ;

    elem_count += Kokkos::Example::box_count( elem_box );
    node_count += Kokkos::Example::box_count( o_node_box );

    {
      Kokkos::Example::box_intersect( test_box , m_owns_node_box[0] , o_node_box );

      if ( Kokkos::Example::box_count( test_box ) ) {
        std::cout << "owns_node[" << m_global_rank << "]{"
                  << " [" << m_owns_node_box[0][0][0] << "," << m_owns_node_box[0][0][1] << ")"
                  << " [" << m_owns_node_box[0][1][0] << "," << m_owns_node_box[0][1][1] << ")"
                  << " [" << m_owns_node_box[0][2][0] << "," << m_owns_node_box[0][2][1] << ")"
                  << "} intersects"
                  << " owns_node[" << rank << "]{"
                  << " [" << o_node_box[0][0] << "," << o_node_box[0][1] << ")"
                  << " [" << o_node_box[1][0] << "," << o_node_box[1][1] << ")"
                  << " [" << o_node_box[2][0] << "," << o_node_box[2][1] << ")"
                  << "}" << std::endl ;
      }
    }

    if ( DecomposeElem == decompose ) {

      Kokkos::Example::box_intersect( test_box , m_uses_elem_box , elem_box );

      if ( Kokkos::Example::box_count( test_box ) ) {

        std::cout << "ElemBox[" << m_global_rank << "]{"
                  << " [" << m_uses_elem_box[0][0] << "," << m_uses_elem_box[0][1] << ")"
                  << " [" << m_uses_elem_box[1][0] << "," << m_uses_elem_box[1][1] << ")"
                  << " [" << m_uses_elem_box[2][0] << "," << m_uses_elem_box[2][1] << ")"
                  << "} intersects"
                  << " ElemBox[" << rank << "]{"
                  << " [" << elem_box[0][0] << "," << elem_box[0][1] << ")"
                  << " [" << elem_box[1][0] << "," << elem_box[1][1] << ")"
                  << " [" << elem_box[2][0] << "," << elem_box[2][1] << ")"
                  << "}" << std::endl ;
      }
    }
  }

  // Sentinal values at the end of the owns and send lists:

  m_owns_node[ m_owns_node_count ][0] = ~0u ;
  m_owns_node[ m_owns_node_count ][1] = ~0u ;
  m_owns_node_box[ m_owns_node_count ][0][0] = 0u ; m_owns_node_box[ m_owns_node_count ][0][0] = ~0u ;
  m_owns_node_box[ m_owns_node_count ][1][0] = 0u ; m_owns_node_box[ m_owns_node_count ][1][0] = ~0u ;
  m_owns_node_box[ m_owns_node_count ][2][0] = 0u ; m_owns_node_box[ m_owns_node_count ][2][0] = ~0u ;

  m_send_node[ m_send_node_count ][0] = ~0u ;
  m_send_node[ m_send_node_count ][1] = ~0u ;
  m_send_node_box[ m_send_node_count ][0][0] = 0u ; m_send_node_box[ m_send_node_count ][0][0] = ~0u ;
  m_send_node_box[ m_send_node_count ][1][0] = 0u ; m_send_node_box[ m_send_node_count ][1][0] = ~0u ;
  m_send_node_box[ m_send_node_count ][2][0] = 0u ; m_send_node_box[ m_send_node_count ][2][0] = ~0u ;

  {
    size_t count = 0 ;
    for ( unsigned i = 0 ; i < m_owns_node_count ; ++i ) {
      count += m_owns_node[i][1] ;
    }
    if ( count != Kokkos::Example::box_count( m_uses_node_box ) ) {
      std::cout << "Node uses count = " << Kokkos::Example::box_count( m_uses_node_box )
                << " error count = " << count << std::endl ;
    }
  }

  if ( global_node_count != node_count ) {
    std::cout << "Node count = " << global_node_count << " overlap error count = " << node_count << std::endl ;
  }

  if ( DecomposeElem == decompose && global_elem_count != elem_count ) {
    std::cout << "Elem count = " << global_elem_count << " overlap error count = " << elem_count << std::endl ;
  }
}

void BoxElemPart::print( std::ostream & s ) const
{
  s << "BoxElemPart P[" << m_global_rank << ":" << m_global_size << "]"
    << std::endl
    << "  elem_box {"
    << " [" << m_uses_elem_box[0][0] << "," << m_uses_elem_box[0][1] << ")"
    << " [" << m_uses_elem_box[1][0] << "," << m_uses_elem_box[1][1] << ")"
    << " [" << m_uses_elem_box[2][0] << "," << m_uses_elem_box[2][1] << ")"
    << " } / {"
    << " [" << m_global_elem_box[0][0] << "," << m_global_elem_box[0][1] << ")"
    << " [" << m_global_elem_box[1][0] << "," << m_global_elem_box[1][1] << ")"
    << " [" << m_global_elem_box[2][0] << "," << m_global_elem_box[2][1] << ")"
    << " }"
    << std::endl
    << "  node_box {"
    << " [" << m_owns_node_box[0][0][0] << "," << m_owns_node_box[0][0][1] << ")"
    << " [" << m_owns_node_box[0][1][0] << "," << m_owns_node_box[0][1][1] << ")"
    << " [" << m_owns_node_box[0][2][0] << "," << m_owns_node_box[0][2][1] << ")"
    << " } / {"
    << " [" << m_uses_node_box[0][0] << "," << m_uses_node_box[0][1] << ")"
    << " [" << m_uses_node_box[1][0] << "," << m_uses_node_box[1][1] << ")"
    << " [" << m_uses_node_box[2][0] << "," << m_uses_node_box[2][1] << ")"
    << " } / {"
    << " [" << m_global_node_box[0][0] << "," << m_global_node_box[0][1] << ")"
    << " [" << m_global_node_box[1][0] << "," << m_global_node_box[1][1] << ")"
    << " [" << m_global_node_box[2][0] << "," << m_global_node_box[2][1] << ")"
    << " }"
    << std::endl ;

  for ( unsigned i = 1 ; i < m_owns_node_count ; ++i ) {
    s << "  P[" << m_owns_node[i][0] << "]"
      << " recv node_box {"
      << " [" << m_owns_node_box[i][0][0] << "," << m_owns_node_box[i][0][1] << ")"
      << " [" << m_owns_node_box[i][1][0] << "," << m_owns_node_box[i][1][1] << ")"
      << " [" << m_owns_node_box[i][2][0] << "," << m_owns_node_box[i][2][1] << ")"
      << " }"
      << std::endl ;
  }

  for ( unsigned i = 0 ; i < m_send_node_count ; ++i ) {
    s << "  P[" << m_send_node[i][0] << "]"
      << " send node_box {"
      << " [" << m_send_node_box[i][0][0] << "," << m_send_node_box[i][0][1] << ")"
      << " [" << m_send_node_box[i][1][0] << "," << m_send_node_box[i][1][1] << ")"
      << " [" << m_send_node_box[i][2][0] << "," << m_send_node_box[i][2][1] << ")"
      << " }"
      << std::endl ;
  }
}

} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------


