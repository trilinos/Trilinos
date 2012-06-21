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

#ifndef KOKKOS_IMPL_ARRAY_FACTORY_HPP
#define KOKKOS_IMPL_ARRAY_FACTORY_HPP

#include <vector>
#include <impl/KokkosArray_MemoryView.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class ViewType >
inline
ViewType create( const std::string & label )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , void >::create( label );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<1> >::create( label , n0 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<2> >
    ::create( label , n0 , n1 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 ,
                                             const size_t n2 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<3> >
    ::create( label , n0 , n1 , n2 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 ,
                                             const size_t n2 ,
                                             const size_t n3 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<4> >
    ::create( label , n0 , n1 , n2 , n3 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 ,
                                             const size_t n2 ,
                                             const size_t n3 ,
                                             const size_t n4 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<5> >
    ::create( label , n0 , n1 , n2 , n3 , n4 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 ,
                                             const size_t n2 ,
                                             const size_t n3 ,
                                             const size_t n4 ,
                                             const size_t n5 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<6> >
    ::create( label , n0 , n1 , n2 , n3 , n4 , n5 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 ,
                                             const size_t n2 ,
                                             const size_t n3 ,
                                             const size_t n4 ,
                                             const size_t n5 ,
                                             const size_t n6 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<7> >
    ::create( label , n0 , n1 , n2 , n3 , n4 , n5 , n6 );
}

template< class ViewType >
inline
ViewType create( const std::string & label , const size_t n0 , 
                                             const size_t n1 ,
                                             const size_t n2 ,
                                             const size_t n3 ,
                                             const size_t n4 ,
                                             const size_t n5 ,
                                             const size_t n6 ,
                                             const size_t n7 )
{
  typename StaticAssertSame<
             ViewType ,
             View< typename ViewType::data_type ,
                   typename ViewType::layout_spec ,
                   typename ViewType::device_type > >::type view_type ;

  return Impl::Factory< view_type , unsigned<8> >
    ::create( label , n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< class DataType ,
          class LayoutDst , class DeviceDst ,
          class LayoutSrc , DeviceSrc >
inline
void deep_copy( const View<DataType,LayoutDst,DeviceDst> & dst ,
                const View<DataType,LayoutSrc,DeviceSrc> & src )
{
  typedef View<DataType,LayoutDst,DeviceDst> dst_type ;
  typedef View<DataType,LayoutSrc,DeviceSrc> src_type ;

  if ( dst.operator!=(src) ) {

    Impl::array_require_equal_dimension( dst.dimension(0) , src.dimension(0) );

    Impl::Factory< dst_type , src_type >::deep_copy( dst , src );
  }
}

} /* namespace KokkosArray */

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class ArrayType , class DeviceOutput >
struct Factory< Array< ArrayType , DeviceOutput > , MirrorUseView >
{
  typedef Array< ArrayType , DeviceOutput > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create( const Array< ArrayType , DeviceInput > & input )
  {
    typedef Array< ArrayType , DeviceInput > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_ARRAY_FACTORY_HPP */

