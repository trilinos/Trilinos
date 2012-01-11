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

#ifndef KOKKOS_VALUE_HPP
#define KOKKOS_VALUE_HPP

#include <cstddef>
#include <string>
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Plain-old-data value allocated on a compute device.
 */
template< typename ValueType , class DeviceType >
class Value {
public:
  typedef ValueType  value_type ;
  typedef DeviceType device_type ;

  typedef Value< value_type , void /* Host */ > HostMirror ;

  /*------------------------------------------------------------------*/
  /** \brief  Access value */
  value_type & operator* () const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  Value();

  /** \brief  Construct a view of the array */
  Value( const Value & rhs );

  /** \brief  Assign to a view of the rhs.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  Value & operator = ( const Value & rhs );

  /**  \brief  Destroy this view of the value.
   *           If the last view then allocated memory is deallocated.
   */
  ~Value();

  /*------------------------------------------------------------------*/
  /** \brief  Allow the Value to be a parallel reduce
   *          'finalize functor' that assigns the reduced value
   *          on the device.
   */
  void operator()( const value_type & rhs ) const ;
};

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
Value< ValueType , DeviceType >
create_labeled_value( const std::string & label );

template< typename ViewType >
Value< typename ViewType::value_type , typename ViewType::device_type >
create_labeled_value( const std::string & label );

template< typename ValueType , class DeviceType >
Value< ValueType , DeviceType >
create_value();

template< typename ViewType >
Value< typename ViewType::value_type , typename ViewType::device_type >
create_value();

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
typename Value< ValueType , DeviceType >::HostMirror
create_mirror( const Value< ValueType , DeviceType > & v );

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
void deep_copy( const Value< ValueType , DeviceType > & dst ,
                const ValueType & src );

template< typename ValueType , class DeviceType >
void deep_copy( ValueType & dst ,
                const Value< ValueType , DeviceType > & src );

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const Value< ValueType , DeviceDst > & dst ,
                const Value< ValueType , DeviceSrc > & src );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
Value< ValueType , DeviceType >
create_labeled_value( const std::string & label )
{ return Value< ValueType , DeviceType >( label ); }

template< typename ValueType , class DeviceType >
inline
Value< ValueType , DeviceType >
create_value()
{ return create_labeled_value< ValueType , DeviceType >( std::string() ); }

template< typename ViewType >
inline
Value< typename ViewType::value_type , typename ViewType::device_type >
create_labeled_value( const std::string & label )
{
  return create_labeled_value< typename ViewType::value_type ,
                               typename ViewType::device_type >( label );
}

template< typename ViewType >
inline
Value< typename ViewType::value_type , typename ViewType::device_type >
create_value()
{
  return create_labeled_value< typename ViewType::value_type ,
                               typename ViewType::device_type >(std::string());
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Impl {

template< typename ValueType , class Device >
class CreateMirror< Value< ValueType , Device > , true /* view */ >
{
public:
  typedef  Value< ValueType , Device >            View ;
  typedef  typename Value< ValueType , Device >::HostMirror  HostMirror ;

  static
  HostMirror create( const View & v ) { return HostMirror( v ); }
};

template< typename ValueType , class Device >
class CreateMirror< Value< ValueType , Device > , false /* copy */ >
{
public:
  typedef  Value< ValueType , Device >            View ;
  typedef  typename Value< ValueType , Device >::HostMirror  HostMirror ;

  static
  HostMirror create( const View & )
    { return create_labeled_value< HostMirror >( std::string() ); }
};

} // namespace Impl

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
typename Value< ValueType , DeviceType >::HostMirror
create_mirror( const Value< ValueType , DeviceType > & v )
{
  typedef Value< ValueType , DeviceType >     device_view ;
  typedef typename device_view::HostMirror      host_view ;
  typedef typename host_view::device_type     host_device ;
  typedef typename host_device::memory_space  host_memory ;
  typedef typename DeviceType::memory_space   device_memory ;

  enum { optimize = Impl::SameType< device_memory , host_memory >::value
#if defined( KOKKOS_MIRROR_VIEW_OPTIMIZE )
                    && KOKKOS_MIRROR_VIEW_OPTIMIZE
#endif
       };

  return Impl::CreateMirror< device_view , optimize >::create( v );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst >
inline
void deep_copy( const Value<ValueType,DeviceDst> & dst ,
                const ValueType & src )
{
  typedef Value< ValueType , DeviceDst > dst_type ;
  typedef ValueType src_type ;
  Impl::DeepCopy<dst_type,src_type>::run( dst , src );
}

template< typename ValueType , class DeviceSrc >
inline
void deep_copy( ValueType & dst ,
                const Value<ValueType,DeviceSrc> & src )
{
  typedef ValueType dst_type ;
  typedef Value< ValueType , DeviceSrc > src_type ;
  Impl::DeepCopy<dst_type,src_type>::run( dst , src );
}

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const Value< ValueType , DeviceDst > & dst ,
                const Value< ValueType , DeviceSrc > & src )
{
  typedef Value< ValueType , DeviceDst > dst_type ;
  typedef Value< ValueType , DeviceSrc > src_type ;
  if ( dst.operator != ( src ) ) {
    Impl::DeepCopy<dst_type,src_type>::run( dst , src );
  }
}

//----------------------------------------------------------------------------

} //Kokkos

#endif /* KOKKOS_VALUE_HPP */


