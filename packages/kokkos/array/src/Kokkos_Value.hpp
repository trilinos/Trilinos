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

template< class Type >
inline
Value< typename Type::value_type , typename Type::device_type >
create_value()
{
  typedef Value< typename Type::value_type , typename Type::device_type > type ;
  return Impl::Factory<type,void>::create();
}

template< class Type >
inline
Value< typename Type::value_type , typename Type::device_type >
create_value( const std::string & label )
{
  typedef Value< typename Type::value_type , typename Type::device_type > type ;
  return Impl::Factory<type,void>::create( label );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst >
inline
void deep_copy( const Value<ValueType,DeviceDst> & dst ,
                const ValueType & src )
{
  typedef Value< ValueType , DeviceDst > dst_type ;
  typedef ValueType src_type ;
  Impl::Factory<dst_type,src_type>::deep_copy( dst , src );
}

template< typename ValueType , class DeviceSrc >
inline
void deep_copy( ValueType & dst ,
                const Value<ValueType,DeviceSrc> & src )
{
  typedef ValueType dst_type ;
  typedef Value< ValueType , DeviceSrc > src_type ;
  Impl::Factory<dst_type,src_type>::deep_copy( dst , src );
}

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const Value< ValueType , DeviceDst > & dst ,
                const Value< ValueType , DeviceSrc > & src )
{
  typedef Value< ValueType , DeviceDst > dst_type ;
  typedef Value< ValueType , DeviceSrc > src_type ;
  if ( dst.operator != ( src ) ) {
    Impl::Factory<dst_type,src_type>::deep_copy( dst , src );
  }
}

//----------------------------------------------------------------------------

} //Kokkos

#include <impl/Kokkos_Value_factory.hpp>

#endif /* KOKKOS_VALUE_HPP */


