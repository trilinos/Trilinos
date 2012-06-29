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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef KOKKOS_HOST_VALUE_HPP
#define KOKKOS_HOST_VALUE_HPP

#include <string.h>

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_Value_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< Value< ValueType , Host > , void >
{
  typedef Value< ValueType , Host > output_type ;

  static inline
  output_type create( const std::string & label )
  {
    output_type output ;
    output.m_memory = KokkosArray::create< typename output_type::view_type >( label );
    memset( output.m_memory.ptr_on_device() , 0 , sizeof(ValueType) );
    return output ;
  }

  static inline
  output_type create() { return create( std::string() ); }
};

template< typename ValueType , class Device >
struct Factory< Value< ValueType , HostMapped< Device > > , void >
{
  typedef Value< ValueType , HostMapped< Device > > output_type ;

  static inline
  output_type create( const std::string & label )
  {
    output_type output ;
    output.m_memory.allocate( 1 , label );
    memset( output.ptr_on_device() , 0 , sizeof(ValueType) );
    return output ;
  }

  static inline
  output_type create() { return create( std::string() ); }
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< Value< ValueType , Host > ,
                Value< ValueType , Host > >
{
  static inline
  void deep_copy( const Value< ValueType , Host > & output ,
                  const Value< ValueType , Host > & input )
  { *output = *input ; }
};

template< typename ValueType >
struct Factory< Value< ValueType , Host > , ValueType >
{
  static inline
  void deep_copy( const Value< ValueType , Host > & output ,
                  const ValueType & input )
  { *output = input ; }
};

template< typename ValueType >
struct Factory< ValueType , Value< ValueType , Host > >
{
  static inline
  void deep_copy( ValueType & output ,
                  const Value< ValueType , Host > & input )
  { output = *input ; }
};

//----------------------------------------------------------------------------

template< typename ValueType , class Device >
struct Factory< Value< ValueType , HostMapped< Device > > ,
                Value< ValueType , HostMapped< Device > > >
{
  static inline
  void deep_copy( const Value< ValueType , HostMapped< Device > > & output ,
                  const Value< ValueType , HostMapped< Device > > & input )
  { *output = *input ; }
};

template< typename ValueType , class Device >
struct Factory< Value< ValueType , Host > ,
                Value< ValueType , HostMapped< Device > > >
{
  static inline
  void deep_copy( const Value< ValueType , Host > & output ,
                  const Value< ValueType , HostMapped< Device > > & input )
  { *output = *input ; }
};

template< typename ValueType , class Device >
struct Factory< Value< ValueType , HostMapped< Device > > ,
                Value< ValueType , Host > >
{
  static inline
  void deep_copy( const Value< ValueType , HostMapped< Device > > & output ,
                  const Value< ValueType , Host > & input )
  { *output = *input ; }
};

template< typename ValueType , class Device >
struct Factory< Value< ValueType , HostMapped< Device > > , ValueType >
{
  static inline
  void deep_copy( const Value< ValueType , HostMapped< Device > > & output ,
                  const ValueType & input )
  { *output = input ; }
};

template< typename ValueType , class Device >
struct Factory< ValueType , Value< ValueType , HostMapped< Device > > >
{
  static inline
  void deep_copy( ValueType & output ,
                  const Value< ValueType , HostMapped< Device > > & input )
  { output = *input ; }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_HOST_VALUE_HPP */

