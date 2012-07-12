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

#include <iostream>

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

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<0> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create();

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<1> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<2> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 , const size_t n2 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<3> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1,n2);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 , const size_t n2 , const size_t n3 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<4> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1,n2,n3);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 , const size_t n2 , const size_t n3 ,
        const size_t n4 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<5> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1,n2,n3,n4);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 , const size_t n2 , const size_t n3 ,
        const size_t n4 , const size_t n5 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<6> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1,n2,n3,n4,n5);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 , const size_t n2 , const size_t n3 ,
        const size_t n4 , const size_t n5 , const size_t n6 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<7> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1,n2,n3,n4,n5,n6);

  return Impl::Factory< view_type , void >::create( label , shape );
}

template< class ViewType >
inline
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label ,
        const size_t n0 , const size_t n1 , const size_t n2 , const size_t n3 ,
        const size_t n4 , const size_t n5 , const size_t n6 , const size_t n7 )
{
  typedef View< typename ViewType::data_type ,
                typename ViewType::layout_type ,
                typename ViewType::device_type > view_type ;

  typedef typename view_type::shape_type    shape_type ;
  typedef typename view_type::memory_space  memory_space ;
  typedef  Impl::Factory< shape_type , memory_space > shape_factory ;

  typedef typename
    Impl::StaticAssertSame< Impl::unsigned_<8> ,
                            Impl::unsigned_< shape_type::rank_dynamic >
                          >::type ok_rank ;

  const shape_type shape = shape_factory::create(n0,n1,n2,n3,n4,n5,n6,n7);

  return Impl::Factory< view_type , void >::create( label , shape );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Vector of a multivector */
template< typename ValueType , class DeviceType >
inline
View< ValueType[] , LayoutLeft , DeviceType >
view( const View< ValueType[][0] , LayoutLeft , DeviceType > & input ,
      const size_t J )
{
  typedef typename
    Impl::StaticAssert< 0 == Impl::rank<ValueType>::value >::type ok_rank ;

  typedef View< ValueType[][0] , LayoutLeft , DeviceType > input_type ;

  return Impl::Factory< void , input_type >::view( input , J );
}

/** \brief  Value within a multivector */
template< typename ValueType , class DeviceType >
inline
View< ValueType , LayoutLeft , DeviceType >
view( const View< ValueType[][0] , LayoutLeft , DeviceType > & input ,
      const size_t I , const size_t J )
{
  typedef typename
    Impl::StaticAssert< 0 == Impl::rank<ValueType>::value >::type ok_rank ;

  typedef View< ValueType[][0] , LayoutLeft , DeviceType > input_type ;

  return Impl::Factory< void , input_type >::view( input , I , J );
}

/** \brief  Value within a vector */
template< typename ValueType , class DeviceType >
inline
View< ValueType , LayoutLeft , DeviceType >
view( const View< ValueType[] , LayoutLeft , DeviceType > & input ,
      const size_t I )
{
  typedef typename
    Impl::StaticAssert< 0 == Impl::rank<ValueType>::value >::type ok_rank ;

  typedef View< ValueType[] , LayoutLeft , DeviceType > input_type ;

  return Impl::Factory< void , input_type >::view( input , I );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Deep copy compatible arrays */

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

//----------------------------------------------------------------------------
/** \brief  Deep copy a single value */

template< class DataType , class LayoutType , class DeviceType >
inline
void deep_copy( const View<DataType,LayoutType,DeviceType> & dst ,
                const DataType & src )
{
  typedef View< DataType , LayoutType , DeviceType >  dst_type ;
  typedef DataType                                    src_type ;
  typedef typename dst_type::shape_type               shape_type ;
  typedef typename dst_type::value_type               value_type ;

  typedef typename
    Impl::assert_shape_is_rank< shape_type , 0 >::type ok_rank ;

  typedef typename
    Impl::StaticAssertAssignable< value_type , DataType >::type ok_assign ;

  Impl::Factory< dst_type , src_type >::deep_copy( dst , src );
}

/** \brief  Deep copy a single value */
template< class DataType , class LayoutType , class DeviceType >
inline
void deep_copy( DataType & dst , 
                const View<DataType,LayoutType,DeviceType> & src )
{
  typedef DataType                                    dst_type ;
  typedef View< DataType , LayoutType , DeviceType >  src_type ;
  typedef typename src_type::shape_type               shape_type ;
  typedef typename src_type::value_type               value_type ;

  typedef typename
    Impl::assert_shape_is_rank< shape_type , 0 >::type ok_rank ;

  typedef typename
    Impl::StaticAssertAssignable< DataType , value_type >::type ok_assign ;

  Impl::Factory< dst_type , src_type >::deep_copy( dst , src );
}

//----------------------------------------------------------------------------

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


/** \brief  Create subviews of multivectors */
template< typename ValueType , class Device >
struct Factory< void , View< ValueType[][0] , LayoutLeft , Device > >
{
  typedef typename StaticAssert< rank<ValueType>::value == 0 >::type ok_rank ;

  // The multivector type:
  typedef View< ValueType[][0] , LayoutLeft , Device > input_type ;

  typedef typename Device::memory_space memory_space ;
  typedef typename input_type::shape_map shape_map ;

  inline static
  View< ValueType[] , LayoutLeft , Device >
  view( const input_type & input , const size_t J )
  {
    View< ValueType[] , LayoutLeft , Device > output ;

    if ( input ) {
      output.m_shape.N0 = input.m_shape.N0 ;
      output.m_ptr_on_device = input.m_ptr_on_device +
                               shape_map::offset( input.m_shape , 0 , J );
      memory_space::increment( output.m_ptr_on_device );
    }

    return output ;
  }

  inline static
  View< ValueType , LayoutLeft , Device >
  view( const input_type & input , const size_t I , size_t J )
  {

    View< ValueType , LayoutLeft , Device > output ;

    if ( input ) {
      output.m_ptr_on_device = input.m_ptr_on_device +
                               shape_map::offset( input.m_shape , I , J );

      memory_space::increment( output.m_ptr_on_device );
    }

    return output ;
  }
};

/** \brief  Create subview of vectors */
template< typename ValueType , class Device >
struct Factory< void , View< ValueType[] , LayoutLeft , Device > >
{
  typedef typename StaticAssert< rank<ValueType>::value == 0 >::type ok_rank ;

  // The multivector type:
  typedef View< ValueType[] , LayoutLeft , Device > input_type ;

  typedef typename Device::memory_space memory_space ;
  typedef typename input_type::shape_map shape_map ;

  inline static
  View< ValueType , LayoutLeft , Device >
  view( const input_type & input , const size_t I )
  {
    View< ValueType , LayoutLeft , Device > output ;

    if ( input ) {
      output.m_ptr_on_device = input.m_ptr_on_device +
                               shape_map::offset( input.m_shape , I );

      memory_space::increment( output.m_ptr_on_device );
    }

    return output ;
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_VIEW_FACTORY_HPP */

