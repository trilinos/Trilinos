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

#include <cstddef>
#include <sstream>
#include <KokkosArray_ParallelFor.hpp>
#include <impl/KokkosArray_Error.hpp>
#include <impl/KokkosArray_ReduceOperator.hpp>

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

template< class FunctorType >
typename FunctorType::value_type
vector_parallel_reduce( const size_t        work_count ,
                        const FunctorType & functor );

template< class FunctorType , class FinalizeType >
void vector_parallel_reduce( const size_t         work_count ,
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
          class ValueOper    /* value_type, init, join */ ,
          class FinalizeType /* serial finalization of reduction value */ ,
          class DeviceType ,
          class WorkSpec = void >
class ParallelReduce ;

template< class FunctorType ,
          class ValueOper ,
          class FinalizeType ,
          class DeviceType >
class MultiFunctorParallelReduceMember ;

template< typename ValueType , class DeviceType >
class ParallelReduceFunctorValue ;

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class FunctorType >
typename FunctorType::value_type
parallel_reduce( const size_t work_count ,
                 const FunctorType & functor )
{
  typedef typename FunctorType::device_type device_type ;
  typedef typename FunctorType::value_type  value_type ;

  typedef Impl::ParallelReduceFunctorValue< value_type , device_type >
    FinalizeType ; 

  const FinalizeType finalize ; 

  Impl::ParallelReduce< FunctorType, FunctorType, FinalizeType, device_type >
    ( work_count , functor , finalize );

  return finalize.result();
}

template< class FunctorType >
typename FunctorType::value_type
vector_parallel_reduce( const size_t work_count ,
                        const FunctorType & functor )
{
  typedef typename FunctorType::device_type device_type ;
  typedef typename FunctorType::value_type  value_type ;

  typedef Impl::ParallelReduceFunctorValue< value_type , device_type >
    FinalizeType ; 

  const FinalizeType finalize ; 

  Impl::ParallelReduce< FunctorType, FunctorType, FinalizeType, device_type , Impl::VectorParallel >
    ( work_count , functor , finalize );

  return finalize.result();
}

//----------------------------------------------------------------------------

template< class FunctorType , class FinalizeType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      const FinalizeType & finalize )
{
  typedef typename FunctorType::device_type device_type ;

  Impl::ParallelReduce< FunctorType, FunctorType, FinalizeType, device_type >
    ( work_count , functor , finalize );
}

//----------------------------------------------------------------------------

template< class FunctorType , typename MemberType >
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      MemberType value[] ,
                      const unsigned count )
{
  typedef typename FunctorType::device_type device_type ;
  typedef typename FunctorType::value_type  value_type ;

  typedef
    typename Impl::StaticAssertSame< value_type , MemberType[] >::type
      ok_match_type ;

  if ( functor.value_count != count ) {
    std::ostringstream msg ;
    msg << "KokkosArray::parallel_reduce( <array_type> ) ERROR "
        << "given incompatible array lengths: functor.value_count("
        << functor.value_count << ") != count(" << count << ")" ;
    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }

  typedef Impl::ParallelReduceFunctorValue< value_type , device_type >
    FinalizeType ; 

  const FinalizeType finalize( functor.value_count );

  Impl::ParallelReduce< FunctorType, FunctorType, FinalizeType, device_type >
    ( work_count , functor , finalize );

  finalize.result( value );
}

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------

#endif /* KOKKOSARRAY_PARALLELREDUCE_HPP */

