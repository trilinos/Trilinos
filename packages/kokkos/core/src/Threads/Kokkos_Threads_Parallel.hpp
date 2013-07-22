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

#ifndef KOKKOS_THREADS_PARALLEL_HPP
#define KOKKOS_THREADS_PARALLEL_HPP

// #include <vector>
// #include <algorithm>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelFor< FunctorType , ParallelWorkRequest , Kokkos::Threads >
{
private:

  const FunctorType  m_functor ;
        bool         m_active ;

public:

  inline ParallelFor( const FunctorType & f , const ParallelWorkRequest & )
    : m_functor(f), m_active(true)
    { Kokkos::Impl::ThreadsExec::start( *this ); }

  void wait() { if ( m_active ) { Kokkos::Impl::ThreadsExec::fence(); m_active = false ; } }

  ~ParallelFor() { wait(); }

  inline void apply( Kokkos::Threads dev , void * ) const { m_functor(dev); }

  inline void join( volatile void * , volatile const void * ) const {}
  inline void init( void * ) const {}
  inline void final( void * ) const {}
};


template< class FunctorType >
class ParallelFor< FunctorType , size_t , Kokkos::Threads >
{
private:

  const FunctorType m_functor ;
  const size_t      m_nwork ;
        bool        m_active ;

public:

  ParallelFor( const FunctorType & f , const size_t n )
    : m_functor(f), m_nwork(n), m_active(true)
    { Kokkos::Impl::ThreadsExec::start( *this ); }

  void wait() { if ( m_active ) { Kokkos::Impl::ThreadsExec::fence(); m_active = false ; } }

  ~ParallelFor() { wait(); }

  inline void apply( Kokkos::Threads dev , void * ) const
  {
    const std::pair< size_t , size_t > range = dev.work_range( m_nwork );
// #pragma ivdep
    for ( size_t iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_functor(iwork);
    }
  }

  inline void join( volatile void * , volatile const void * ) const {}
  inline void init( void * ) const {}
  inline void final( void * ) const {}
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
struct ThreadsReduceAdapter ;

template< class FunctorType >
struct ThreadsReduceAdapter< FunctorType , Kokkos::ParallelWorkRequest >
  : public Kokkos::Impl::ReduceAdapter< FunctorType >
{
  typedef Kokkos::Impl::ReduceAdapter< FunctorType > base_type ;

  inline ThreadsReduceAdapter( const FunctorType & f , const Kokkos::ParallelWorkRequest & )
    : base_type(f) {}

  inline void apply( Kokkos::Threads dev , void * ptr ) const
  {
    typename base_type::reference_type update = base_type::reference( ptr );

    base_type::m_functor(dev,update);
  }
};


template< class FunctorType >
struct ThreadsReduceAdapter< FunctorType , size_t > : public Kokkos::Impl::ReduceAdapter< FunctorType >
{
  typedef Kokkos::Impl::ReduceAdapter< FunctorType > base_type ;

  const size_t m_nwork ;

  ThreadsReduceAdapter( const FunctorType & f , const size_t n )
    : base_type(f), m_nwork(n) {}

  inline void apply( Kokkos::Threads dev , void * ptr ) const
  {
    typename base_type::reference_type update = base_type::reference( ptr );

    const std::pair< size_t , size_t > range = dev.work_range( m_nwork );
// #pragma ivdep
    for ( size_t iwork = range.first ; iwork < range.second ; ++iwork ) {
      base_type::m_functor(iwork,update);
    }
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelReduce< FunctorType , WorkSpec , Kokkos::Threads >
{
private:

  typedef ThreadsReduceAdapter< FunctorType , WorkSpec > ReduceAdapter ;
  typedef typename ReduceAdapter::pointer_type  pointer_type ;

  ReduceAdapter  m_reduce ;
  pointer_type   m_result_ptr ;

public:

  // Create new functor in the asynchronous functor memory space
  // and then launch it.
  ParallelReduce( const FunctorType & functor , const WorkSpec & work ,
                  const pointer_type result_ptr = 0 )
    : m_reduce( functor , work )
    , m_result_ptr( result_ptr )
    {
      if ( result_ptr ) Kokkos::Impl::ThreadsExec::acquire( this );

      Kokkos::Impl::ThreadsExec::start( *this );
    }

  void wait()
  {
    if ( m_result_ptr ) {

      Kokkos::Impl::ThreadsExec::fence();

      const pointer_type * const data = (const pointer_type) Kokkos::Impl::ThreadsExec::root_reduce_scratch();

      for ( unsigned i = 0 ; i < m_reduce.value_count() ; ++i ) { m_result_ptr[i] = data[i]; }

      m_result_ptr = 0 ;

      Kokkos::Impl::ThreadsExec::release( this );
    }
  }

  ~ParallelReduce() { wait(); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADS_PARALLEL_HPP */

