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

#ifndef KOKKOSARRAY_REDUCEOPERATOR_HPP
#define KOKKOSARRAY_REDUCEOPERATOR_HPP

#include <KokkosArray_Macros.hpp>
#include <KokkosArray_View.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Define a facade for a finalize functor type
 *          and verify that the reduction value type is compatible.
 *
 *  By default the finalize functor type is passed through.
 *  A compatible View type is given a finalize functor facade.
 */
template< class FinalizeType , class ValueType >
struct ReduceOperatorFinalize {

  typedef typename
    StaticAssertSame< typename FinalizeType::value_type , ValueType >
      ::type ok_type ;

  typedef FinalizeType type ;

  // Only instantiate this function if it is called
  template< class F >
  inline static
  unsigned value_count( const F & f )
    {
      typedef
        typename StaticAssertSame< FinalizeType , F >::type ok_functor_type ;

      return f.value_count ;
    }
};

template < class DataType , class LayoutType , class DeviceType ,
           class ScalarType >
struct ReduceOperatorFinalize< View< DataType , LayoutType , DeviceType > ,
                               ScalarType >
{
  typedef View< DataType , LayoutType , DeviceType > view_type ;

  typedef typename
    StaticAssertSame< typename view_type::scalar_type , ScalarType >
      ::type ok_type ;

  struct type {

    const view_type m_view ;

    typedef typename view_type::device_type  device_type ;
    typedef ScalarType                       value_type ;

    type( const view_type & v ) : m_view( v ) {}

    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void operator()( const ScalarType & input ) const
      { *m_view = input ; }
  };
};

template < class DataType , class LayoutType , class DeviceType ,
           class ScalarType >
struct ReduceOperatorFinalize< View< DataType , LayoutType , DeviceType > ,
                               ScalarType[] >
{
  typedef View< DataType , LayoutType , DeviceType > view_type ;

  typedef typename
    StaticAssertSame< typename view_type::scalar_type , ScalarType >
      ::type ok_type ;

  typedef typename
    StaticAssert< view_type::Rank == 1 >::type ok_rank ;

  inline static
  unsigned value_count( const view_type & v ) { return v.dimension_0(); }

  struct type {

    const view_type m_view ;

    typedef typename view_type::device_type  device_type ;
    typedef ScalarType                       value_type[] ;
    const unsigned                           value_count ;

    explicit
    type( const view_type & v )
      : m_view( v ), value_count( v.dimension_0() ) {}

    KOKKOSARRAY_INLINE_DEVICE_FUNCTION
    void operator()( const ScalarType input[] ) const
    { for ( unsigned i = 0 ; i < value_count ; ++i ) m_view(i) = input[i] ; }
  };
};

//----------------------------------------------------------------------------
/** \brief  The reduce operator is the aggregate of
 *          ValueOper::value_type
 *          ValueOper::init( update )
 *          ValueOper::join( update , input )
 *          FinalizeFunctor::operator()( input )
 */
template < class ValueOper ,
           class FinalizeFunctor ,
           class ValueType  = typename ValueOper::value_type >
struct ReduceOperator
{
private:

  ReduceOperator();
  ReduceOperator & operator = ( const ReduceOperator & );

  typedef ReduceOperatorFinalize< FinalizeFunctor , ValueType >
    wrapper_type ;

  const typename wrapper_type::type  m_finalize ;

public:

  typedef ValueType    value_type ;
  typedef value_type & reference_type ;

  inline static
  unsigned value_size( const FinalizeFunctor & )
    { return sizeof(value_type); }

  KOKKOSARRAY_INLINE_FUNCTION
  explicit ReduceOperator( const FinalizeFunctor & finalize )
    : m_finalize( finalize ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned value_size() const { return sizeof(value_type); }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void join( volatile void * update , const volatile void * input ) const
    {
      typedef       volatile value_type * vvp ;
      typedef const volatile value_type * cvvp ;

      ValueOper::join( *vvp(update) , *cvvp(input) );
    }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  reference_type init( void * update ) const
    {
      reference_type ref = *((value_type*) update);
      ValueOper::init( ref );
      return ref ;
    }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void finalize( const void * input ) const
    { m_finalize( *( (const value_type *) input ) ); }
};

//----------------------------------------------------------------------------
/** \brief  The reduce operator is the aggregate of
 *          ValueOper::value_type
 *          ValueOper::init( update , value_count )
 *          ValueOper::join( update , input , value_count )
 *          FinalizeFunctor::operator()( input )
 */
template < class ValueOper ,
           class FinalizeFunctor ,
           typename MemberType >
struct ReduceOperator< ValueOper , FinalizeFunctor , MemberType[] >
{
private:

  ReduceOperator();
  ReduceOperator & operator = ( const ReduceOperator & );

  typedef ReduceOperatorFinalize< FinalizeFunctor , MemberType[] >
    wrapper_type ;

  const typename wrapper_type::type  m_finalize ;

public:

  typedef MemberType   value_type[] ;
  typedef MemberType * reference_type ;

  inline static
  unsigned value_size( const FinalizeFunctor & f )
    { return sizeof(MemberType) * wrapper_type::value_count( f ); }

  explicit ReduceOperator( const FinalizeFunctor & finalize )
    : m_finalize( finalize ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned value_size() const
    { return sizeof(MemberType) * m_finalize.value_count ; }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void join( volatile void * update , const volatile void * input ) const
    {
      typedef       volatile MemberType * vvp ;
      typedef const volatile MemberType * cvvp ;

      ValueOper::join( vvp(update) , cvvp(input) , m_finalize.value_count );
    }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  reference_type init( void * update ) const
    {
      reference_type ref = (reference_type) update ;
      ValueOper::init( ref , m_finalize.value_count );
      return ref ;
    }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void finalize( const void * input ) const
    { m_finalize( (const MemberType *) input ); }
};

} // namespace Impl
} // namespace KokkosArray

#endif /* KOKKOSARRAY_REDUCEOPERATOR_HPP */

