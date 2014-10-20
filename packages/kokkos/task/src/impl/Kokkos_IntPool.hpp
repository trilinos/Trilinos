/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_INTPOOL_HPP
#define KOKKOS_IMPL_INTPOOL_HPP

#include <iostream>

#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class Device >
class IntPool {
private:

  enum { SHIFT  = 8 };
  enum { BLOCK  = 1 << SHIFT /* = 256 */ };
  enum { MASK   = BLOCK - 1 };
  enum { STRIDE = BLOCK + 1 }; // 0 == m_pool.dimension_0() % STRIDE

  // { { entries[BLOCK] , offset_to_free_list }* , noop }
  //
  View<int*,Device> m_pool ;
  unsigned          m_size ;

  // Offset to entry at begining of block's free list.
  // Precondition: 0 <= i && i < m_size
  KOKKOS_FORCEINLINE_FUNCTION
  int block_list_begin( const int i ) const
    { return ( i & ~int(MASK) ) + ( i >> SHIFT ) + BLOCK ; }

  KOKKOS_FORCEINLINE_FUNCTION
  int entry( const int i ) const
    { return i + ( i >> SHIFT ); }

public:

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  bool is_claimed( const int i ) const
    {
      bool result = 0 <= i && i < int(m_size) ;

      if ( result ) {
        // Verify that the block is OK and the index is not in the free list

        const int j   = entry(i);
        const int end = block_list_begin( i );
        const int beg = end - BLOCK ;

        int k = m_pool[ end ];
        for ( int count = 0 ; count < BLOCK && 0 <= k ; ++count , k = m_pool[ k ] ) {
          if ( k < beg || end <= k || j == k ) { result = false ; }
        }
        if ( -1 != k ) { result = false ; }
      }
      return result ;
    }

  template < typename iType >
  KOKKOS_FORCEINLINE_FUNCTION
  int & operator[]( const iType & i ) const
    { return m_pool[ entry(i) ]; }

  //----------------------------------------

  enum ClaimStatus { SUCCESS /* Claimed an entry */
                   , FAIL    /* Failed, try again with new hint */
                   , FULL    /* Full, no entries to claim */
                   };

  // Claim the next available entry.
  // Input 'value' is a hint where to look for an available entry.
  // If SUCCESS then 'value' is the claimed entry.
  // If FAIL    then 'value' is a recommended hint for the next attempt.
  // If FULL    then no entries are available.
  KOKKOS_INLINE_FUNCTION
  ClaimStatus claim( int & value ) const
    {
      int volatile * const pool = m_pool.ptr_on_device();
      ClaimStatus result = FULL ;

      // Free list begin of the block for this hint.

      for ( unsigned k = 0 , jfree = block_list_begin( unsigned(value) % m_size ) ;
            k < m_pool.dimension_0() ; k += STRIDE ) {
        const int j = pool[jfree] ;
        if ( 0 <= j ) {
          if ( atomic_compare_exchange_strong( pool + jfree , j , pool[j] ) ) {
            // pool[jfree] = pool[j]
            result = SUCCESS ;
            pool[j] = 0 ;
          }
          else {
            result = FAIL ;
          }
          value = j - j / STRIDE ;
          break ;
        }
        if ( m_pool.dimension_0() < ( jfree += STRIDE ) ) jfree = BLOCK ;
      }

      return result ;
    }

  KOKKOS_INLINE_FUNCTION
  void release( const int i ) const
    {
      int volatile * const pool = m_pool.ptr_on_device();

      const int j     = entry(i);
      const int jfree = block_list_begin(i);

      while ( ! atomic_compare_exchange_strong( pool + jfree , ( pool[j] = pool[jfree] ) , j ) );
    }
 
  //-----------------------------------

  void print( std::ostream & s ) const
    {
      for ( unsigned i = 0 ; i < m_pool.dimension_0() ; ++i ) {
        if ( 0 == i % STRIDE ) s << std::endl ;
        s << " " << m_pool[i] ;
      }
      s << std::endl ;
    }

  //-----------------------------------

  struct Init {
    typedef Device device_type ;
    View<int*,Device> m_pool ;
    unsigned          m_size ;

    Init( const View<int*,Device> v , const unsigned n )
      : m_pool(v) , m_size(n) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const int j ) const
      {
        const int knext = 1 + ( j % STRIDE );
        const int jnext = ( knext <  BLOCK ) ? j + 1 : (
                          ( knext == BLOCK ) ? -1 : j - BLOCK );
        m_pool[j] = ( jnext - jnext / STRIDE ) < int(m_size) ? jnext : -1 ;
      }
  };

  IntPool( int capacity )
    : m_pool( Kokkos::ViewAllocateWithoutInitializing("IntPool"), 1 + ( ( capacity + MASK ) >> SHIFT ) * STRIDE )
    , m_size( capacity )
    { Kokkos::parallel_for( m_pool.dimension_0() , Init(m_pool,m_size) ); }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_INTPOOL_HPP */

