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

#ifndef KOKKOS_THREADS_PARALLEL_HPP
#define KOKKOS_THREADS_PARALLEL_HPP

#include <vector>
// #include <algorithm>

#include <KokkosArray_ParallelFor.hpp>
#include <KokkosArray_ParallelReduce.hpp>

#include <impl/KokkosArray_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , Kokkos::Threads , void > {
public:

  struct Adapter {

    const FunctorType m_functor ;

    inline explicit Adapter( const FunctorType & f ) : m_functor(f) {}

    inline void operator()( Kokkos::Threads dev ) const { m_functor(dev); }
    inline void join( volatile void * , volatile const void * ) const {}
    inline void finalize( void * ) const {}
  };

  ParallelFor( const FunctorType & functor )
    { Kokkos::Impl::ThreadsExec::template start<Adapter>( functor ); }
};

template< class FunctorType >
class ParallelFor< FunctorType , Kokkos::Threads , size_t >
{
public:

  struct Adapter {
    FunctorType functor ;
    size_t      nwork ;

    inline void operator()( Kokkos::Threads dev ) const
    {
      const std::pair< size_t , size_t > range = dev.work_range( nwork );
// #pragma ivdep
      for ( size_t iwork = range.first ; iwork < range.second ; ++iwork ) {
        functor(iwork);
      }
    }

    inline void join( volatile void * , volatile const void * ) const {}
    inline void finalize( void * ) const {}

    Adapter( const FunctorType & f , const size_t n )
      : functor(f), nwork(n) {}
  };

  ParallelFor( const size_t NWork , const FunctorType & functor )
  {
    if ( NWork ) {
      Kokkos::Impl::ThreadsExec::template start<Adapter>( functor , NWork );
    }
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename ValueType >
class ParallelReduceFunctorValue< ValueType , Kokkos::Threads >
{
public:
  typedef ValueType value_type ;

  ParallelReduceFunctorValue() {}

  inline void operator()( const value_type & ) const {}

  void result( value_type & value ) const
  {
    value_type * const ptr = (value_type*) Kokkos::Impl::ThreadsExec::root_reduce_scratch();
    value = *ptr ;
  }
};

template< typename MemberType >
class ParallelReduceFunctorValue< MemberType[] , Kokkos::Threads >
{
public:
  typedef MemberType    value_type[] ;
  const HostSpace::size_type value_count ;

  inline void operator()( const MemberType [] ) const {}

  explicit
  ParallelReduceFunctorValue( HostSpace::size_type n )
    : value_count(n)
    {}

  void result( value_type result ) const
  {
    MemberType * const ptr = (MemberType *) Kokkos::Impl::ThreadsExec::root_reduce_scratch();

    for ( HostSpace::size_type i = 0 ; i < value_count ; ++i ) result[i] = ptr[i] ;
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class FinalizeType >
class ParallelReduce< Kokkos::Threads , FunctorType , FinalizeType , void >
{
public:

  ParallelReduce( const FunctorType & functor , const FinalizeType & finalize )
    {
      typedef ReduceOperator< FunctorType , FinalizeType > ReduceOp ;
      Kokkos::Impl::ThreadsExec::template start<ReduceOp>( functor , finalize );
    }
};

template< class FunctorType , class FinalizeType >
class ParallelReduce< Kokkos::Threads , FunctorType , FinalizeType , size_t >
{
public:

  ParallelReduce( const size_t NWork ,
                  const FunctorType & functor ,
                  const FinalizeType & finalize )
    {
      typedef ReduceOperator< FunctorType , FinalizeType , size_t > ReduceOp ;
      Kokkos::Impl::ThreadsExec::template start<ReduceOp>( functor , finalize , NWork );
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADS_PARALLEL_HPP */

