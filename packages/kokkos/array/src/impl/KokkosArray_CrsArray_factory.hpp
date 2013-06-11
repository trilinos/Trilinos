/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOSARRAY_IMPL_CRSARRAY_FACTORY_HPP
#define KOKKOSARRAY_IMPL_CRSARRAY_FACTORY_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class Arg1Type , class Arg2Type , typename SizeType >
inline
typename CrsArray< DataType , Arg1Type , Arg2Type , SizeType >::HostMirror
create_mirror_view( const CrsArray<DataType,Arg1Type,Arg2Type,SizeType > & view ,
                    typename Impl::enable_if< ViewTraits<DataType,Arg1Type,Arg2Type,void>::is_hostspace >::type * = 0 )
{
  return view ;
}

template< class DataType , class Arg1Type , class Arg2Type , typename SizeType >
inline
typename CrsArray< DataType , Arg1Type , Arg2Type , SizeType >::HostMirror
create_mirror_view( const CrsArray<DataType,Arg1Type,Arg2Type,SizeType > & view ,
                    typename Impl::enable_if< ! ViewTraits<DataType,Arg1Type,Arg2Type,void>::is_hostspace >::type * = 0 )
{
  // Force copy:
  typedef Impl::ViewAssignment< Impl::LayoutDefault > alloc ;
  typedef CrsArray< DataType , Arg1Type , Arg2Type , SizeType >  crsarray_type ;

  typename crsarray_type::HostMirror               tmp ;
  typename crsarray_type::row_map_type::HostMirror tmp_row_map ;

  // Allocation to match:
  alloc( tmp_row_map , view.row_map );
  alloc( tmp.entries , view.entries );

  // Assignment of 'const' from 'non-const'
  tmp.row_map = tmp_row_map ;

  // Deep copy:
  deep_copy( tmp_row_map , view.row_map );
  deep_copy( tmp.entries , view.entries );

  return tmp ;
}

template< class DataType , class Arg1Type , class Arg2Type , typename SizeType >
inline
typename CrsArray< DataType , Arg1Type , Arg2Type , SizeType >::HostMirror
create_mirror( const CrsArray<DataType,Arg1Type,Arg2Type,SizeType > & view )
{
#if KOKKOSARRAY_MIRROR_VIEW_OPTIMIZE
  // Allow choice via type:
  return create_mirror_view( view );
#else
  // Force copy:
  typedef Impl::ViewAssignment< Impl::LayoutDefault > alloc ;
  typedef CrsArray< DataType , Arg1Type , Arg2Type , SizeType >  crsarray_type ;

  typename crsarray_type::HostMirror               tmp ;
  typename crsarray_type::row_map_type::HostMirror tmp_row_map ;

  // Allocation to match:
  alloc( tmp_row_map , view.row_map );
  alloc( tmp.entries , view.entries );

  // Assignment of 'const' from 'non-const'
  tmp.row_map = tmp_row_map ;

  // Deep copy:
  deep_copy( tmp_row_map , view.row_map );
  deep_copy( tmp.entries , view.entries );

  return tmp ;
#endif
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class CrsArrayType , class InputSizeType >
inline
typename CrsArrayType::crsarray_type
create_crsarray( const std::string & label ,
                 const std::vector< InputSizeType > & input )
{
  typedef CrsArrayType                  output_type ;
  typedef std::vector< InputSizeType >  input_type ;

  typedef typename output_type::entries_type   entries_type ;

  typedef View< typename output_type::size_type [] ,
                typename output_type::array_layout ,
                typename output_type::device_type > work_type ;

  output_type output ;

  // Create the row map:

  const size_t length = input.size();

  {
    work_type row_work( "tmp" , length + 1 );

    typename work_type::HostMirror row_work_host =
      create_mirror_view( row_work );

    size_t sum = 0 ;
    row_work_host[0] = 0 ;
    for ( size_t i = 0 ; i < length ; ++i ) {
      row_work_host[i+1] = sum += input[i];
    }

    deep_copy( row_work , row_work_host );

    output.entries   = entries_type( label , sum );
    output.row_map   = row_work ;
  }

  return output ;
}

//----------------------------------------------------------------------------

template< class CrsArrayType , class InputSizeType >
inline
typename CrsArrayType::crsarray_type
create_crsarray( const std::string & label ,
                 const std::vector< std::vector< InputSizeType > > & input )
{
  typedef CrsArrayType                                output_type ;
  typedef std::vector< std::vector< InputSizeType > > input_type ;
  typedef typename output_type::entries_type          entries_type ;
  typedef typename output_type::size_type             size_type ;

  typedef typename
    Impl::assert_shape_is_rank_one< typename entries_type::shape_type >::type
      ok_rank ;

  typedef View< typename output_type::size_type [] ,
                typename output_type::array_layout ,
                typename output_type::device_type > work_type ;

  output_type output ;

    // Create the row map:

  const size_t length = input.size();

  {
    work_type row_work( "tmp" , length + 1 );

    typename work_type::HostMirror row_work_host =
      create_mirror_view( row_work );

    size_t sum = 0 ;
    row_work_host[0] = 0 ;
    for ( size_t i = 0 ; i < length ; ++i ) {
      row_work_host[i+1] = sum += input[i].size();
    }

    deep_copy( row_work , row_work_host );

    output.entries   = entries_type( label , sum );
    output.row_map   = row_work ;
  }

  // Fill in the entries:
  {
    typename entries_type::HostMirror host_entries =
      create_mirror_view( output.entries );

    size_t sum = 0 ;
    for ( size_t i = 0 ; i < length ; ++i ) {
      for ( size_t j = 0 ; j < input[i].size() ; ++j , ++sum ) {
        host_entries( sum ) = input[i][j] ;
      }
    }

    deep_copy( output.entries , host_entries );
  }

  return output ;
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_IMPL_CRSARRAY_FACTORY_HPP */

