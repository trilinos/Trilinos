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

#ifndef KOKKOSARRAY_PARALLELREDUCE_HPP
#define KOKKOSARRAY_PARALLELREDUCE_HPP

#include <KokkosArray_Macros.hpp>
#include <cstddef>

//----------------------------------------------------------------------------

namespace KokkosArray {

/** \brief  Run 'functor' in parallel and reduce FunctorType::value_type.
 *
 *  If 'value_type' is a plain-old-data then call:
 *    FunctorType::operator()( int_type , value_type & ) 
 *    FunctorType::init( value_type & );
 *    FunctorType::join( volatile value_type & update ,
 *                       const volatile value_type & input );
 *
 *  If 'value_type' is 'type[]' and 'type' is plain-old-data then call:
 *    FunctorType::operator()( int_type , value_type [] ) 
 *    FunctorType::init( value_type [] , int_type count );
 *    FunctorType::join( volatile value_type update [] ,
 *                       const volatile value_type input [] ,
 *                       int_type count );
 */
template< class FunctorType >
typename FunctorType::value_type
parallel_reduce( const size_t        work_count ,
                 const FunctorType & functor );

template< class FunctorType , class FinalizeType >
void parallel_reduce( const size_t         work_count ,
                      const FunctorType  & functor ,
                      const FinalizeType & finalize );

} // namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {

/** \brief  Multiple functor parallel reduction.
 *
 *  Create and execute a collection of reduction functors
 *  that contribute to a common final result.
 */
template< class ReduceOper ,
          class FinalizeType = typename ReduceOper::value_type ,
          class DeviceType   = typename ReduceOper::device_type >
class MultiFunctorParallelReduce {
public:
  typedef typename DeviceType::size_type size_type ;

  FinalizeType result ;

  MultiFunctorParallelReduce();

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & );

  void execute();
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Inlined implementation:
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class FunctorType  /* parallel work operator */ ,
          class ReduceOper   /* value_type, init, join */ ,
          class FinalizeType /* serial finalization of reduction value */ ,
          class DeviceType >
class ParallelReduce ;

template< typename ValueType , class DeviceType >
class FunctorAssignment ;

template< class FunctorType ,
          class ReduceOper ,
          class FinalizeType ,
          class DeviceType >
class MultiFunctorParallelReduceMember ;

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

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

  typedef typename
    StaticAssertSame< ValueType , 
                      typename FinalizeFunctor::value_type >
      ::type ok_finalize ;

  FinalizeFunctor m_finalize ;

public:

  typedef ValueType    value_type ;
  typedef value_type & reference_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  explicit ReduceOperator( const FinalizeFunctor & finalize )
  : m_finalize( finalize )
  {}

  KOKKOSARRAY_INLINE_FUNCTION
  unsigned value_size() const { return sizeof(value_type); }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void join( void * update , const void * input ) const
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
  {
    typedef const value_type * cvp ;
    m_finalize( *cvp(input) );
  }
};

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

  typedef typename
    StaticAssertSame< MemberType[] , 
                      typename FinalizeFunctor::value_type >
      ::type ok_finalize ;

  FinalizeFunctor m_finalize ;

public:

  typedef MemberType   value_type[] ;
  typedef MemberType * reference_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  explicit ReduceOperator( const FinalizeFunctor & finalize )
    : m_finalize( finalize )
    {}

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  ReduceOperator( const ReduceOperator & rhs )
    : m_finalize( rhs.m_finalize ) {}
  
  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  unsigned value_size() const
    { return sizeof(MemberType) * m_finalize.value_count ; }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void join( void * update , const void * input ) const
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
  {
    typedef const MemberType * cvp ;
    m_finalize( *cvp(input) );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {

template< class FunctorType >
typename FunctorType::value_type
parallel_reduce( const size_t work_count ,
                 const FunctorType & functor )
{
  typedef typename FunctorType::device_type device_type ;
  typedef typename FunctorType::value_type  value_type ;

  value_type result ;

  { // Destruction of 'tmp' guarantees data is assigned to 'result'
    typedef Impl::FunctorAssignment< value_type , device_type > Finalize ;

    Finalize tmp( result );

    Impl::ParallelReduce< FunctorType, FunctorType, Finalize, device_type >
      ( work_count , functor , tmp );
  }

  return result ;
}

template< class FunctorType , class FinalizeType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      const FinalizeType & finalize )
{
  typedef typename FunctorType::device_type device_type ;

  Impl::ParallelReduce< FunctorType, FunctorType, FinalizeType, device_type >
    ( work_count , functor , finalize );
}

} // namespace KokkosArray

//----------------------------------------------------------------------------

#endif /* KOKKOSARRAY_PARALLELREDUCE_HPP */

