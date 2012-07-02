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

#ifndef KOKKOS_IMPL_VIEW_FACTORY_HPP
#define KOKKOS_IMPL_VIEW_FACTORY_HPP

#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , void >::create( label );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<1> >::create( label , n0 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<2> >
    ::create( label , n0 , n1 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<3> >
    ::create( label , n0 , n1 , n2 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<4> >
    ::create( label , n0 , n1 , n2 , n3 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<5> >
    ::create( label , n0 , n1 , n2 , n3 , n4 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 , const size_t n5 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<6> >
    ::create( label , n0 , n1 , n2 , n3 , n4 , n5 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 , const size_t n5 ,
                                    const size_t n6 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<7> >
    ::create( label , n0 , n1 , n2 , n3 , n4 , n5 , n6 );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 , const size_t n5 ,
                                    const size_t n6 , const size_t n7 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  return Impl::Factory< view_type , Impl::unsigned_<8> >
    ::create( label , n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< class DataTypeDst , class LayoutDst , class DeviceDst ,
          class DataTypeSrc , class LayoutSrc , class DeviceSrc >
inline
void deep_copy( const View<DataTypeDst,LayoutDst,DeviceDst> & dst ,
                const View<DataTypeSrc,LayoutSrc,DeviceSrc> & src )
{
  typedef View<DataTypeDst,LayoutDst,DeviceDst> dst_type ;
  typedef View<DataTypeSrc,LayoutSrc,DeviceSrc> src_type ;

  Impl::Factory< dst_type , src_type >::deep_copy( dst , src );
}

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DataType , class LayoutType , class DeviceOutput >
struct Factory< View< DataType , LayoutType , DeviceOutput > , MirrorUseView >
{
  typedef View< DataType , LayoutType , DeviceOutput > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create( const View< DataType, LayoutType, DeviceInput > & input )
  {
    typedef View< DataType , LayoutType , DeviceInput > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_VIEW_FACTORY_HPP */

