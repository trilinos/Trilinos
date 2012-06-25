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
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

namespace KokkosArray {

/// \class Value
/// \brief View of a single value allocated on a compute device.
///
/// \tparam ValueType The type of the value.  Must be a "plain old
///   data" type.
/// \tparam DeviceType Type of the compute device.
///
/// This class wraps a single datum of type ValueType.  That datum's
/// memory lives in the compute device's memory.  Depending on the
/// type of compute device, that may mean that the host processor
/// can't access the value directly (without first copying to a host
/// mirror Value).
template< typename ValueType , class DeviceType >
class Value {
public:
  typedef ValueType  value_type ;
  typedef DeviceType device_type ;

  /// \brief The type of a Value that resides in host memory.
  ///
  /// If you want to read on the host a Value stored on a compute
  /// device, you'll need to do a deep_copy from the compute device's
  /// Value to a Value of this type.  A Value of this type is stored
  /// on the host, and its raw value can be accessed directly via its
  /// operator* method (see below).
  typedef Value< value_type , void /* Host */ > HostMirror ;

  /// \brief Access the value directly.
  ///
  /// \warning If this is not a host Value (i.e., if the DeviceType is
  ///   not a host device), then this method may only be called in a
  ///   compute kernel.  If this is a host Value, then you may call
  ///   this method either inside or outside of a compute kernel.
  value_type & operator* () const ;

  //! Construct a NULL view (a view of no value).
  Value();

  //! Construct a view of the given value.
  Value( const Value & rhs );

  /// \brief Replace this Value's view with a view of rhs.
  ///
  /// If this Value is the last view, then allocated memory is
  /// deallocated.
  Value & operator = ( const Value & rhs );

  /// \brief Destroy this view of the value.
  ///
  /// If this is the last view, then allocated memory is deallocated.
  ~Value();

  /// \brief In a kernel, assign rhs to the value stored in *this.
  ///
  /// This method allows the Value to work as a parallel reduce
  /// 'finalize functor' that assigns the reduced value on the device.
  ///
  /// \warning If the DeviceType is not a host device, then this
  ///   method may only be called in a compute kernel, unless the
  ///   device type's memory is accessible from the host.
  void operator()( const value_type & rhs ) const ;
};

//----------------------------------------------------------------------------

/// \fn create_value()
/// \brief Nonmember constructor of a Value.
/// \tparam Type Specialization of Value (see above).
template< class Type >
inline
Value< typename Type::value_type , typename Type::device_type >
create_value()
{
  typedef Value< typename Type::value_type , typename Type::device_type > type ;
  return Impl::Factory<type,void>::create();
}

/// \fn create_value()
/// \brief Nonmember constructor of a Value, with a label.
/// \tparam Type Specialization of Value (see above).
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

/// \brief Copy a raw ValueType (a plain old datum) to a Value.
/// \relatesalso Value
template< typename ValueType , class DeviceDst >
inline
void deep_copy( const Value<ValueType,DeviceDst> & dst ,
                const ValueType & src )
{
  typedef Value< ValueType , DeviceDst > dst_type ;
  typedef ValueType src_type ;
  Impl::Factory<dst_type,src_type>::deep_copy( dst , src );
}

/// \brief Copy a Value to a raw ValueType (a plain old datum) in host memory.
/// \relatesalso Value
template< typename ValueType , class DeviceSrc >
inline
void deep_copy( ValueType & dst ,
                const Value<ValueType,DeviceSrc> & src )
{
  typedef ValueType dst_type ;
  typedef Value< ValueType , DeviceSrc > src_type ;
  Impl::Factory<dst_type,src_type>::deep_copy( dst , src );
}

/// \brief Copy a Value to another Value.
/// \relatesalso Value
///
/// The two Values may have different device types.
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

} //KokkosArray

#include <impl/KokkosArray_Value_factory.hpp>

#endif /* KOKKOS_VALUE_HPP */


