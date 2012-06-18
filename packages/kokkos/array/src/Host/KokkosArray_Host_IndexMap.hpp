/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_HOST_INDEXMAP_HPP
#define KOKKOS_HOST_INDEXMAP_HPP

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_IndexMapRight_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class ValueType ,
          class IndexMapOutput ,
          class IndexMapInput >
struct HostIndexMapDeepCopy : public HostThreadWorker<void>
{
  typedef Host::size_type                      size_type ;
  typedef Impl::MemoryView< ValueType , Host > data_type ;

private:

  data_type       output_data ;
  data_type       input_data ;
  IndexMapOutput  output_map ;
  IndexMapInput   input_map ;

  const size_type N1 ;
  const size_type N2 ;
  const size_type N3 ;
  const size_type N4 ;
  const size_type N5 ;
  const size_type N6 ;
  const size_type N7 ;

  void copy8( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < N7 ; ++i7 ) {
      const size_type i_out = output_map(i0,i1,i2,i3,i4,i5,i6,i7);
      const size_type i_in  = input_map( i0,i1,i2,i3,i4,i5,i6,i7);
      output_data[i_out] = input_data[i_in];
    }}}}}}}}
  }

  void copy7( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
      const size_type i_out = output_map(i0,i1,i2,i3,i4,i5,i6);
      const size_type i_in  = input_map( i0,i1,i2,i3,i4,i5,i6);
      output_data[i_out] = input_data[i_in];
    }}}}}}}
  }

  void copy6( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
      const size_type i_out = output_map(i0,i1,i2,i3,i4,i5);
      const size_type i_in  = input_map( i0,i1,i2,i3,i4,i5);
      output_data[i_out] = input_data[i_in];
    }}}}}}
  }

  void copy5( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
      const size_type i_out = output_map(i0,i1,i2,i3,i4);
      const size_type i_in  = input_map( i0,i1,i2,i3,i4);
      output_data[i_out] = input_data[i_in];
    }}}}}
  }

  void copy4( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
      const size_type i_out = output_map(i0,i1,i2,i3);
      const size_type i_in  = input_map( i0,i1,i2,i3);
      output_data[i_out] = input_data[i_in];
    }}}}
  }

  void copy3( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
      const size_type i_out = output_map(i0,i1,i2);
      const size_type i_in  = input_map( i0,i1,i2);
      output_data[i_out] = input_data[i_in];
    }}}
  }

  void copy2( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
      const size_type i_out = output_map(i0,i1);
      const size_type i_in  = input_map( i0,i1);
      output_data[i_out] = input_data[i_in];
    }}
  }

  void copy1( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
      output_data[i0] = input_data[i0];
    }
  }

  void execute_on_thread( HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( output_map.dimension(0) );

    switch( output_map.rank() ) {
    case 1: copy1( range ); break ;
    case 2: copy2( range ); break ;
    case 3: copy3( range ); break ;
    case 4: copy4( range ); break ;
    case 5: copy5( range ); break ;
    case 6: copy6( range ); break ;
    case 7: copy7( range ); break ;
    case 8: copy8( range ); break ;
    }
  }

  HostIndexMapDeepCopy(
    const data_type & oData , const IndexMapOutput & oMap ,
    const data_type & iData , const IndexMapInput & iMap )
    : output_data( oData ), input_data( iData )
    , output_map(  oMap ),  input_map(  iMap )
    , N1( 1 < output_data.rank() ? output_data.dimension(1) : 0 )
    , N2( 2 < output_data.rank() ? output_data.dimension(2) : 0 )
    , N3( 3 < output_data.rank() ? output_data.dimension(3) : 0 )
    , N4( 4 < output_data.rank() ? output_data.dimension(4) : 0 )
    , N5( 5 < output_data.rank() ? output_data.dimension(5) : 0 )
    , N6( 6 < output_data.rank() ? output_data.dimension(6) : 0 )
    , N7( 7 < output_data.rank() ? output_data.dimension(7) : 0 )
    {}

public:

  static inline
  void deep_copy( const data_type & oData , const IndexMapOutput & oMap ,
                  const data_type & iData , const IndexMapInput & iMap )
  {
    typedef MemoryManager< Host > memory_manager ;

    memory_manager::disable_memory_view_tracking();

    const HostIndexMapDeepCopy driver( oData , oMap , iData , iMap );

    memory_manager::enable_memory_view_tracking();

    HostThreadWorker<void>::execute( driver );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_HOST_INDEXMAP_HPP */


