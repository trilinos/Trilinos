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

#ifndef KOKKOS_HOST_CRSMAP_HPP
#define KOKKOS_HOST_CRSMAP_HPP

namespace Kokkos {
namespace Impl {

template< typename SizeType >
class CreateCrsMap< CrsMap< Host , SizeType > , CrsMap< Host , SizeType > > {
public:

  typedef CrsMap< Host , SizeType > type ;

  template< typename IteratorType >
  static
  type create( const std::string & label ,
               const IteratorType second_count_begin ,
               const IteratorType second_count_end )
  {
    type m ;

    SizeType first_count = 0 ;
    SizeType total_count = 0 ;

    for ( IteratorType i = second_count_begin ; i != second_count_end ; ++i ) {
      ++first_count ;
      total_count += *i ;
    }

    m.m_first_count = first_count ;
    m.m_total_count = total_count ;
    m.m_memory.allocate( first_count + 1 , label );

    SizeType * const offset = m.m_memory.ptr_on_device();

    first_count = 0 ;
    total_count = 0 ;
    offset[0] = 0 ;

    for ( IteratorType i = second_count_begin ; i != second_count_end ; ++i ) {
      offset[ ++first_count ] = ( total_count += *i );
    }

    return m ;
  }

  static
  type create( const type & rhs ) { return rhs ; }
};

template< typename DstSizeType , typename SrcSizeType >
class CreateCrsMap< CrsMap< Host , DstSizeType > ,
                    CrsMap< Host , SrcSizeType > > {
public:

  typedef CrsMap< Host , DstSizeType > type ;

  static
  type create( const CrsMap<Host,SrcSizeType> & rhs )
  {
    type m ;

    m.m_first_count = (DstSizeType) rhs.m_first_count ;
    m.m_total_count = (DstSizeType) rhs.m_total_count ;
    m.m_memory.allocate( m.m_first_count + 1 , std::string() );

    DstSizeType * const dst = m  .m_memory.ptr_on_device();
    SrcSizeType * const src = rhs.m_memory.ptr_on_device();

    HostParallelCopy<DstSizeType,SrcSizeType>(
      m  .m_memory.ptr_on_device() ,
      rhs.m_memory.ptr_on_device() ,
      m  .m_first_count + 1 );

    return m ;
  }
};

//----------------------------------------------------------------------------

template<>
class CreateCrsMap< Host , Host > {
public:
  typedef CrsMap< Host , Host::size_type > type ;

  template< typename IteratorType >
  static
  type create( const std::string & label ,
               const IteratorType second_count_begin ,
               const IteratorType second_count_end )
  {
    return CreateCrsMap< type >
      ::create( label , second_count_begin , second_count_end );
  }
};

template< typename SrcSizeType >
class CreateCrsMap< Host , CrsMap< Host , SrcSizeType > > {
public:
  typedef CrsMap< Host , Host::size_type > type ;

  static
  type create( const CrsMap< Host , SrcSizeType > & rhs )
  {
    return CreateCrsMap< type , CrsMap< Host , SrcSizeType > >::create( rhs );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_HOST_CRSMAP_HPP */

