/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOS_REDUCEOPERATOR_HPP
#define KOKKOS_REDUCEOPERATOR_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  The reduce operator is the aggregate of
 *          ValueOper::value_type
 *          ValueOper::init( update )
 *          ValueOper::join( update , input )
 */
template < class ValueOper ,
           class ValueType = typename ValueOper::value_type >
struct ReduceOperator
{
private:

  ReduceOperator();
  ReduceOperator & operator = ( const ReduceOperator & );

public:

  typedef ValueType    value_type ;
  typedef value_type & reference_type ;
  typedef value_type * pointer_type ;

  inline static
  unsigned value_size( const ValueOper & )
    { return sizeof(value_type); }

  KOKKOS_INLINE_FUNCTION
  explicit ReduceOperator( const ValueOper & ) {}

  KOKKOS_INLINE_FUNCTION
  unsigned value_size() const { return sizeof(value_type); }

  KOKKOS_INLINE_FUNCTION
  unsigned value_count() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  reference_type reference( void * const p ) const { return *((value_type*)p); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  reference_type reference( void * const p , const iType & i ) const
    { return ((value_type*)p)[i]; }

  KOKKOS_INLINE_FUNCTION
  void init( void * update ) const
    { ValueOper::init( *((value_type*)update) ); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void init( void * update , const iType & i ) const
    { ValueOper::init( ((value_type*)update)[i] ); }

  KOKKOS_INLINE_FUNCTION
  void join( volatile void * update , const volatile void * input ) const
    {
      typedef       volatile value_type * vvp ;
      typedef const volatile value_type * cvvp ;

      ValueOper::join( *vvp(update) , *cvvp(input) );
    }

  template< unsigned N >
  KOKKOS_INLINE_FUNCTION
  void join( volatile void * update ) const
    {
      typedef       volatile value_type * vvp ;
      typedef const volatile value_type * cvvp ;

      for ( unsigned i = 1 ; i < N ; ++i ) {
        ValueOper::join( *vvp(update) , cvvp(update)[i] );
      }
    }

  KOKKOS_INLINE_FUNCTION
  void final( void * ) const {}
};

//----------------------------------------------------------------------------
/** \brief  The reduce operator is the aggregate of
 *          ValueOper::value_type
 *          ValueOper::init( update , value_count )
 *          ValueOper::join( update , input , value_count )
 */
template < class ValueOper , typename MemberType >
struct ReduceOperator< ValueOper , MemberType[] >
{
private:

  ReduceOperator();
  ReduceOperator & operator = ( const ReduceOperator & );

  const unsigned int m_value_count ;

public:

  typedef MemberType   value_type[] ;
  typedef MemberType * reference_type ;
  typedef MemberType * pointer_type ;

  inline static
  unsigned value_size( const ValueOper & f )
    { return sizeof(MemberType) * f.value_count ; }

  explicit ReduceOperator( const ValueOper & f )
    : m_value_count( f.value_count ) {}

  KOKKOS_INLINE_FUNCTION
  unsigned value_size() const
    { return sizeof(MemberType) * m_value_count ; }

  KOKKOS_INLINE_FUNCTION
  unsigned value_count() const
    { return m_value_count ; }

  KOKKOS_INLINE_FUNCTION
  reference_type reference( void * const p ) const
    { return (reference_type)p; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  reference_type reference( void * const p , const iType & i ) const
    { return ((reference_type)p) + m_value_count * i ; }

  KOKKOS_INLINE_FUNCTION
  void init( void * update ) const
    {
      ValueOper::init( (reference_type) update , m_value_count );
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void init( void * update , const iType & i ) const
    {
      ValueOper::init( ((reference_type) update) + m_value_count * i  , m_value_count );
    }

  KOKKOS_INLINE_FUNCTION
  void join( volatile void * update , const volatile void * input ) const
    {
      typedef       volatile MemberType * vvp ;
      typedef const volatile MemberType * cvvp ;

      ValueOper::join( vvp(update) , cvvp(input) , m_value_count );
    }

  template< unsigned N >
  KOKKOS_INLINE_FUNCTION
  void join( volatile void * update ) const
    {
      typedef       volatile MemberType * vvp ;
      typedef const volatile MemberType * cvvp ;

      for ( unsigned i = 1 ; i < N ; ++i ) {
        ValueOper::join( vvp(update) , cvvp(update) + m_value_count * i , m_value_count );
      }
    }

  KOKKOS_INLINE_FUNCTION
  void final( void * ) const {}
};

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_REDUCEOPERATOR_HPP */

