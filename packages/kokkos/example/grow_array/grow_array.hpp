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

#ifndef EXAMPLE_GROW_ARRAY
#define EXAMPLE_GROW_ARRAY

#include <stdlib.h>

#include <KokkosCore_config.h>
#include <Kokkos_Atomic.hpp>

//----------------------------------------------------------------------------

namespace Example {

template< class Device >
struct GrowArrayFunctor {

  typedef Device device_type ;

  enum { SHIFT = sizeof(int) == 8 ? 6 : 5 };
  enum { MASK  = ( 1 << SHIFT ) - 1 };

  const Kokkos::View<int*,Device>  m_search_flags ; // bit flags for values to append
  const Kokkos::View<int*,Device>  m_search_array ; // array to append values
  const Kokkos::View<int,Device>   m_search_count ; // offset
  const int m_search_team_chunk ;

  GrowArrayFunctor( int array_length , int search_length , int print = 1 )
    : m_search_flags( "flags" , ( search_length + MASK ) & ~unsigned(MASK) )
    , m_search_array( "array" , array_length )
    , m_search_count( "count" )
    , m_search_team_chunk( 2048 )
    {
      typename Kokkos::View<int,Device>::HostMirror  count = Kokkos::create_mirror_view( m_search_count );
      typename Kokkos::View<int*,Device>::HostMirror flags = Kokkos::create_mirror_view( m_search_flags );

      for ( int i = 0 ; i < array_length ; ++i ) {
        const long int index = ( lrand48() * search_length ) >> 31 ;
        flags( index >> SHIFT ) |= ( 1 << ( index & MASK ) );
      }

      Kokkos::deep_copy( m_search_flags , flags );

      Kokkos::ParallelWorkRequest work ;

      // Each team works on 'm_search_team_chunk' span of the search_length
      work.team_size   = Device::team_recommended();
      work.league_size = ( search_length + m_search_team_chunk - 1 ) / m_search_team_chunk ;

      Kokkos::parallel_for( work , *this );

      typename Kokkos::View<int*,Device>::HostMirror results = Kokkos::create_mirror_view( m_search_array );

      Kokkos::deep_copy( results , m_search_array );
      Kokkos::deep_copy( count ,   m_search_count );

      int error_count = 0 ;
      for ( int i = 0 ; i < *count ; ++i ) {
        const int index = results(i);
        const bool flag = flags( index >> SHIFT ) & ( 1 << ( index & MASK ) );
        if ( ! flag ) {
          if ( print ) std::cerr << "( " << i << " , " << index << " )";
          ++error_count ;
        }
      }

      if ( error_count ) {
        std::cerr << std::endl << "Example::GrowArrayFunctor( " << array_length
                  << " , " << search_length
                  << " ) error_count( " << error_count << " )"
                  << std::endl ;
      }
    }

  KOKKOS_INLINE_FUNCTION
  bool flag_is_set( const int index ) const
    {
      // 64 or 32 bit integer:

      const int j = index >> SHIFT ;
      const int k = 1 << ( index & MASK );
      const int s = ( j < int(m_search_flags.dimension_0()) ) && ( 0 != ( m_search_flags(j) & k ) );

      return s ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( Device dev ) const
    {
      enum { LOCAL_BUFFER_LENGTH = 16 };

      int local_buffer[ LOCAL_BUFFER_LENGTH ] ;
      int local_count = 0 ;

      // Each team searches 'm_search_team_chunk' indices

            int search_team_begin = dev.league_rank() * m_search_team_chunk ;
      const int search_team_end   = search_team_begin + m_search_team_chunk ;

      int k = 0 ;

      while ( search_team_begin < search_team_end ) {

        const int search_index = search_team_begin + dev.team_rank();

        // If this thread should include the search index.
        if ( flag_is_set(search_index) ) {
          local_buffer[ local_count ] = search_index ;
          ++local_count ;
        }

        search_team_begin += dev.team_size(); // Striding team by team size
        ++k ;

        // Write buffer if end of search or a thread might have filled its buffer.
        if ( k == LOCAL_BUFFER_LENGTH || ! ( search_team_begin < search_team_end ) ) {

          // Team's exclusive scan of threads' contributions, with global offset.
          const int team_offset = dev.team_scan( local_count , & *m_search_count );

          // Copy locally buffered entries into global array:
          for ( int i = 0 ; i < local_count ; ++i ) {
            m_search_array( team_offset + i ) = local_buffer[i] ;
          }

          k = 0 ;
          local_count = 0 ;
        }
      }
    }
};






} // namespace Example

//----------------------------------------------------------------------------

#endif /* #ifndef EXAMPLE_GROW_ARRAY */

