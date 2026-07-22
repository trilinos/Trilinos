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

#ifndef KOKKOS_BOXELEMPART_HPP
#define KOKKOS_BOXELEMPART_HPP

#include <utility>
#include <ostream>
#include <Kokkos_Macros.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

KOKKOS_INLINE_FUNCTION
void box_intersect( unsigned box[][2] ,
                    const unsigned boxA[][2] ,
                    const unsigned boxB[][2] )
{
  for ( unsigned i = 0 ; i < 3 ; ++i ) {
    box[i][0] = boxA[i][0] > boxB[i][0] ? boxA[i][0] : boxB[i][0] ;
    box[i][1] = boxA[i][1] < boxB[i][1] ? boxA[i][1] : boxB[i][1] ;
    if ( box[i][0] > box[i][1] ) box[i][1] = box[i][0] ;
  }
}

KOKKOS_INLINE_FUNCTION
size_t box_count( const unsigned box[][2] )
{
  return size_t( box[0][1] - box[0][0] ) *
         size_t( box[1][1] - box[1][0] ) *
         size_t( box[2][1] - box[2][0] );
}

KOKKOS_INLINE_FUNCTION
void box_ghost_layer( const unsigned global_box[][2] ,
                      const unsigned local_box[][2] ,
                      const unsigned ghost_layer ,
                            unsigned ghost_box[][2] )
{
  for ( unsigned i = 0 ; i < 3 ; ++i ) {
    ghost_box[i][0] = global_box[i][0] + ghost_layer > local_box[i][0] ? global_box[i][0] : local_box[i][0] - ghost_layer ;
    ghost_box[i][1] = global_box[i][1] < local_box[i][1] + ghost_layer ? global_box[i][1] : local_box[i][1] + ghost_layer ;
  }
}

void box_partition( const unsigned global_size ,
                    const unsigned global_rank ,
                    const unsigned global_box[][2] ,
                          unsigned box[][2] );

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

/** \brief Partition a box of hexahedral elements among subdomains.
 *
 *  Nodes are ordered locally as follows:
 *    { owned_by[ this_process ] ,
 *      owned_by[ neighbor_process[0] ] ,
 *      owned_by[ neighbor_process[1] ] ,
 *      owned_by[ neighbor_process[2] ] ,
 *      ... };
 */
class BoxElemPart {
public:

  enum Decompose { DecomposeNode , DecomposeElem };
  enum ElemOrder { ElemLinear , ElemQuadratic };

  BoxElemPart( const ElemOrder elem_order ,
               const Decompose decompose ,
               const unsigned global_size ,
               const unsigned global_rank ,
               const unsigned elem_nx ,
               const unsigned elem_ny ,
               const unsigned elem_nz );

  KOKKOS_INLINE_FUNCTION
  size_t global_elem_count() const
    { return Kokkos::Example::box_count( m_global_elem_box ); }

  KOKKOS_INLINE_FUNCTION
  size_t global_node_count() const
    { return Kokkos::Example::box_count( m_global_node_box ); }

  KOKKOS_INLINE_FUNCTION
  size_t uses_elem_count() const
    { return Kokkos::Example::box_count( m_uses_elem_box ); }

  KOKKOS_INLINE_FUNCTION
  size_t owns_node_count() const
    { return Kokkos::Example::box_count( m_owns_node_box[0] ); }

  KOKKOS_INLINE_FUNCTION
  size_t uses_node_count() const
    { return Kokkos::Example::box_count( m_uses_node_box ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  size_t uses_elem_offset( const unsigned ix ,
                           const unsigned iy ,
                           const unsigned iz ) const
  {
    return size_t( ix - m_uses_elem_box[0][0] ) + size_t( m_uses_elem_box[0][1] - m_uses_elem_box[0][0] ) * (
           size_t( iy - m_uses_elem_box[1][0] ) + size_t( m_uses_elem_box[1][1] - m_uses_elem_box[1][0] ) * (
           size_t( iz - m_uses_elem_box[2][0] ) ) );
  }

  KOKKOS_INLINE_FUNCTION
  void uses_elem_coord( size_t lid , unsigned c[] ) const
  {
    const unsigned nx = m_uses_elem_box[0][1] - m_uses_elem_box[0][0] ;
    const unsigned ny = m_uses_elem_box[1][1] - m_uses_elem_box[1][0] ;

    c[0] = m_uses_elem_box[0][0] + lid % nx ; lid /= nx ;
    c[1] = m_uses_elem_box[1][0] + lid % ny ; lid /= ny ;
    c[2] = m_uses_elem_box[2][0] + lid ;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  unsigned global_coord_max( unsigned axis ) const
  { return m_global_node_box[axis][1] - 1 ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  void local_node_coord( size_t lid , unsigned coord[] ) const
  {
    // Local id within an 'owns' block (has sentinal)
    unsigned j = 0 ;
    while ( m_owns_node[j][1] <= lid ) { lid -= m_owns_node[j][1] ; ++j ; }

    // Map to global coordinates:
    const unsigned nx = m_owns_node_box[j][0][1] - m_owns_node_box[j][0][0] ;
    const unsigned ny = m_owns_node_box[j][1][1] - m_owns_node_box[j][1][0] ;

    coord[0] = m_owns_node_box[j][0][0] + lid % nx ; lid /= nx ;
    coord[1] = m_owns_node_box[j][1][0] + lid % ny ; lid /= ny ;
    coord[2] = m_owns_node_box[j][2][0] + lid ;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned local_node_id( const unsigned c[] ) const
  {
    // Find which 'owns' block and accumulate the offset of this block:
    size_t lid = 0 ;
    unsigned j = 0 ;
    while ( ! ( m_owns_node_box[j][0][0] <= c[0] && c[0] < m_owns_node_box[j][0][1] &&
                m_owns_node_box[j][1][0] <= c[1] && c[1] < m_owns_node_box[j][1][1] &&
                m_owns_node_box[j][2][0] <= c[2] && c[2] < m_owns_node_box[j][2][1] ) ) {
      
      lid += m_owns_node[j][1] ;
      ++j ;
    }

    // Map offset to the block plus offset within the block:
    return lid +
           size_t( c[0] - m_owns_node_box[j][0][0] ) + size_t( m_owns_node_box[j][0][1] - m_owns_node_box[j][0][0] ) * (
           size_t( c[1] - m_owns_node_box[j][1][0] ) + size_t( m_owns_node_box[j][1][1] - m_owns_node_box[j][1][0] ) * (
           size_t( c[2] - m_owns_node_box[j][2][0] ) ) );
  }

  KOKKOS_INLINE_FUNCTION
  size_t global_node_id( const unsigned c[] ) const
  {
    return size_t( c[0] - m_global_node_box[0][0] ) + size_t( m_global_node_box[0][1] - m_global_node_box[0][0] ) * (
           size_t( c[1] - m_global_node_box[1][0] ) + size_t( m_global_node_box[1][1] - m_global_node_box[1][0] ) * (
           size_t( c[2] - m_global_node_box[2][0] ) ) );
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  unsigned recv_node_msg_count() const { return m_owns_node_count - 1 ; }

  KOKKOS_INLINE_FUNCTION
  unsigned recv_node_rank(  unsigned msg ) const { return m_owns_node[msg+1][0] ; }

  KOKKOS_INLINE_FUNCTION
  unsigned recv_node_count( unsigned msg ) const { return m_owns_node[msg+1][1] ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  unsigned send_node_msg_count() const { return m_send_node_count ; }

  KOKKOS_INLINE_FUNCTION
  unsigned send_node_rank(  unsigned msg ) const { return m_send_node[msg][0] ; }

  KOKKOS_INLINE_FUNCTION
  unsigned send_node_count( unsigned msg ) const { return m_send_node[msg][1] ; }

  KOKKOS_INLINE_FUNCTION
  unsigned send_node_id_count() const
  {
    unsigned count = 0 ;
    for ( unsigned i = 0 ; i < m_send_node_count ; ++i ) {
      count += m_send_node[i][1] ;
    }
    return count ;
  }

  KOKKOS_INLINE_FUNCTION
  unsigned send_node_id( unsigned item ) const
  {
    // Find which send list this send item is in:
    unsigned j = 0 ;
    while ( m_send_node[j][1] <= item ) { item -= m_send_node[j][1] ; ++j ; }

    // Map to global coordinate:
    const unsigned nx = m_send_node_box[j][0][1] - m_send_node_box[j][0][0] ;
    const unsigned ny = m_send_node_box[j][1][1] - m_send_node_box[j][1][0] ;

    unsigned c[3] ;

    c[0] = m_send_node_box[j][0][0] + item % nx ; item /= nx ;
    c[1] = m_send_node_box[j][1][0] + item % ny ; item /= ny ;
    c[2] = m_send_node_box[j][2][0] + item ;

    // Map to local id:
    return size_t( c[0] - m_owns_node_box[0][0][0] ) + size_t( m_owns_node_box[0][0][1] - m_owns_node_box[0][0][0] ) * (
           size_t( c[1] - m_owns_node_box[0][1][0] ) + size_t( m_owns_node_box[0][1][1] - m_owns_node_box[0][1][0] ) * (
           size_t( c[2] - m_owns_node_box[0][2][0] ) ) );
  }

  //----------------------------------------

  void print( std::ostream & s ) const ;

private:

  // Maximum number of processes in a neighborhood, including this process
  enum { PROC_NEIGH_MAX = 32 };

  void local( const unsigned  rank ,
                    unsigned  uses_elem[][2] ,
                    unsigned  owns_node[][2] ,
                    unsigned  uses_node[][2] ) const ;

  unsigned  m_global_size ;
  unsigned  m_global_rank ;
  Decompose m_decompose ;
  ElemOrder m_elem_order ;

  unsigned m_global_elem_box[3][2] ;
  unsigned m_global_node_box[3][2] ;
  unsigned m_uses_elem_box[3][2] ;
  unsigned m_uses_node_box[3][2] ;

  // [ processor rank , count ]
  unsigned m_owns_node_box[ PROC_NEIGH_MAX ][3][2] ;
  unsigned m_owns_node[     PROC_NEIGH_MAX ][2] ;
  unsigned m_owns_node_count ;

  unsigned m_send_node_box[ PROC_NEIGH_MAX ][3][2] ;
  unsigned m_send_node[     PROC_NEIGH_MAX ][2] ;
  unsigned m_send_node_count ;
};

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_BOXELEMPART_HPP */

