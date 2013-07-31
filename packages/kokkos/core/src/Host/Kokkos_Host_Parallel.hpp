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

#ifndef KOKKOS_HOST_PARALLEL_HPP
#define KOKKOS_HOST_PARALLEL_HPP

#include <vector>
#include <algorithm>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <impl/Kokkos_StaticAssert.hpp>

#include <Host/Kokkos_Host_Thread.hpp>

//only enable for intel 13 or better
#define KOKKOS_ENABLE_SIMD     defined( __INTEL_COMPILER ) && ! defined( __MIC__ ) \
                                 && (__INTEL_COMPILER > 1299)

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Base class for a parallel driver executing on a thread pool. */

struct HostThreadWorker {
public:

  virtual ~HostThreadWorker() {}

  void execute() const ;
  void execute_serial() const ;

  /** \brief  Virtual method called on threads.
   *
   *  The function must call either 'HostThread::end_barrier'
   *  or 'HostThread::end_reduce' as the ending statement.
   */
  virtual void execute_on_thread( HostThread & ) const = 0 ;

  //-----------------------------------

protected:

  HostThreadWorker() {}

private:

  HostThreadWorker( const HostThreadWorker & );
  HostThreadWorker & operator = ( const HostThreadWorker & );
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelFor< FunctorType , WorkSpec , Host > : public HostThreadWorker
{
public:

  typedef Host::size_type  size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  void execute_on_thread( HostThread & thread ) const
  {
#if KOKKOS_ENABLE_SIMD
    enum { vectorize = is_same<WorkSpec,VectorParallel>::value && 1 < HostSpace::WORK_ALIGNMENT };
#else
    enum { vectorize = 0 };
#endif
    const std::pair< size_type , size_type > range =
      thread.work_range( m_work_count );

    if ( ! vectorize ) {
      for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
        m_work_functor( iwork );
      }
    }
#if KOKKOS_ENABLE_SIMD
    else {
#pragma simd
#pragma ivdep
      for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
        m_work_functor( iwork );
      }
    }
#endif

    // Required end-of-function barrier
    thread.end_barrier();
  }

  ParallelFor( const FunctorType & functor , const size_t work_count )
    : m_work_functor( functor )
    , m_work_count( work_count )
    { HostThreadWorker::execute(); }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class HostMultiFunctorParallelForMember ;

template<>
class HostMultiFunctorParallelForMember<void> {
public:

  virtual ~HostMultiFunctorParallelForMember() {}

  virtual void apply( HostThread & ) const = 0 ;
};

template< class FunctorType >
class HostMultiFunctorParallelForMember
  : public HostMultiFunctorParallelForMember<void> {
public:
  typedef Host::size_type size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ~HostMultiFunctorParallelForMember() {}

  HostMultiFunctorParallelForMember(
    const FunctorType & work_functor ,
    const size_type work_count )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}

  void apply( HostThread & this_thread ) const
  {
    const std::pair< size_type , size_type > range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork );
    }
  }
};

} // namespace Impl

template<>
class MultiFunctorParallelFor< Host >
  : public Impl::HostThreadWorker
{
public:
  typedef Host::size_type               size_type ;
  typedef Impl::HostMultiFunctorParallelForMember<void>  worker_type ;

  typedef std::vector< worker_type * > MemberContainer ;

  typedef MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->apply( this_thread );
    }

    this_thread.end_barrier();
  }

public:

  MultiFunctorParallelFor() : m_member_functors() {}

  ~MultiFunctorParallelFor()
  {
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  {
    typedef Impl::HostMultiFunctorParallelForMember<FunctorType> member_work_type ;

    worker_type * const m = new member_work_type( functor , work_count );

    m_member_functors.push_back( m );
  }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class WorkSpec >
class ParallelReduce< FunctorType , WorkSpec , Host > : public HostThreadWorker
{
public:

  typedef          Host::size_type      size_type ;
  typedef ReduceAdapter< FunctorType >  reduce_oper ;

  typedef typename reduce_oper::pointer_type pointer_type ;

  const reduce_oper   m_reduce ;
  const size_type     m_work_count ;

  void execute_on_thread( HostThread & this_thread ) const
  {
#if KOKKOS_ENABLE_SIMD
    enum { work_align = is_same<WorkSpec,VectorParallel>::value &&
                        power_of_two<HostSpace::WORK_ALIGNMENT>::value
                      ? HostSpace::WORK_ALIGNMENT : 1 };
    enum { work_mask  = work_align - 1 };
#else
    enum { work_align = 1 };
    enum { work_mask  = 0 };
#endif

    // Iterate this thread's work

    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_work_count );

    if ( ! work_mask ) {
      // This thread's reduction value, initialized
      m_reduce.init( this_thread.reduce_data() );

      for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
        m_reduce.m_functor( iwork , m_reduce.reference( this_thread.reduce_data() ) );
      }
    }
#if KOKKOS_ENABLE_SIMD
    else {

#pragma simd
#pragma ivdep
      for ( size_type j = 0 ; j < HostSpace::WORK_ALIGNMENT ; ++j ) {
        m_reduce.init( this_thread.reduce_data() , j );
      }

#pragma simd vectorlength(work_align)
#pragma ivdep
      for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
        m_reduce.m_functor( iwork , m_reduce.reference( this_thread.reduce_data() , iwork & work_mask ) );
      }

      m_reduce.template join< HostSpace::WORK_ALIGNMENT >( this_thread.reduce_data() );
    }
#endif

    // End the routine with a reduction.
    this_thread.end_reduce( m_reduce );
  }

  ParallelReduce( const FunctorType  & functor ,
                  const size_type      work_count ,
                  pointer_type         result = 0 )
    : m_reduce( functor )
    , m_work_count( work_count )
    {
      Host::resize_reduce_scratch( m_reduce.value_size() );
      HostThreadWorker::execute();
      m_reduce.final( Host::root_reduce_scratch() );

      if ( result ) {
        pointer_type ptr = (pointer_type) Host::root_reduce_scratch();
        for ( unsigned i = 0 ; i < m_reduce.value_count() ; ++i ) result[i] = ptr[i] ;
      }
    }

  void wait() {}
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ReduceOper >
class HostMultiFunctorParallelReduceMember ;

template< class ReduceOper >
class HostMultiFunctorParallelReduceMember<void,ReduceOper> {
public:
  virtual ~HostMultiFunctorParallelReduceMember() {}

  virtual void apply( HostThread & , const ReduceOper & ) const = 0 ;
};


template< class FunctorType , class ReduceOper >
class HostMultiFunctorParallelReduceMember
  : public HostMultiFunctorParallelReduceMember<void,ReduceOper> {
public:
  typedef Host::size_type size_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ~HostMultiFunctorParallelReduceMember() {}

  HostMultiFunctorParallelReduceMember(
    const FunctorType & work_functor ,
    const size_type work_count )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}

  // virtual method
  void apply( HostThread & this_thread ,
              const ReduceOper & reduce ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork , reduce.reference( this_thread.reduce_data() ) );
    }
  }
};

} // namespace Impl

template< class ReduceOper >
class MultiFunctorParallelReduce< ReduceOper , Host > : public Impl::HostThreadWorker
{
public:

  typedef Impl::ReduceAdapter< ReduceOper > reduce_oper ;
  typedef          Host::size_type         size_type ;
  typedef typename ReduceOper::value_type  value_type ;
  typedef Impl::HostMultiFunctorParallelReduceMember<void,reduce_oper> worker_type ;

  typedef std::vector< worker_type * > MemberContainer ;

  typedef typename MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;
  reduce_oper     m_reduce ;

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    // This thread's reduction value, initialized
    m_reduce.init( this_thread.reduce_data() );

    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->apply( this_thread , m_reduce );
    }

    // End the routine with a reduction
    this_thread.end_reduce( m_reduce );
  }

public:

  explicit MultiFunctorParallelReduce( const ReduceOper & oper )
    : m_member_functors()
    , m_reduce( oper )
    { }

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  {
    typedef Impl::HostMultiFunctorParallelReduceMember<FunctorType,reduce_oper> member_work_type ;

    worker_type * const m = new member_work_type( functor , work_count );

    m_member_functors.push_back( m );
  }

  void execute() const
  {
    Host::resize_reduce_scratch( m_reduce.value_size() );
    Impl::HostThreadWorker::execute();
    m_reduce.final( Host::root_reduce_scratch() );
  }

  void output( typename reduce_oper::pointer_type result )
  {
    typename reduce_oper::pointer_type ptr =
      (typename reduce_oper::pointer_type) Host::root_reduce_scratch();
    for ( unsigned i = 0 ; i < m_reduce.value_count() ; ++i ) result[i] = ptr[i] ;
  }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename DstType , typename SrcType  >
class HostParallelCopy : public HostThreadWorker
{
public:

        DstType * const m_dst ;
  const SrcType * const m_src ;
  const HostSpace::size_type m_count ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    const HostThread::work_range_type range =
      this_thread.work_range( m_count );
    DstType * const x_end = m_dst + range.second ;
    DstType *       x     = m_dst + range.first ;
    const SrcType * y     = m_src + range.first ;

    for ( ; x_end != x ; ++x , ++y ) { *x = (DstType) *y ; }

    this_thread.end_barrier();
  }

  HostParallelCopy( DstType * dst , const SrcType * src ,
                    HostSpace::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    { HostThreadWorker::execute(); }
};

template< typename DstType >
class HostParallelFill : public HostThreadWorker
{
public:

  DstType * const m_dst ;
  const DstType   m_src ;
  const HostSpace::size_type m_count ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    const HostThread::work_range_type range =
      this_thread.work_range( m_count );
    DstType * const x_end = m_dst + range.second ;
    DstType *       x     = m_dst + range.first ;

    for ( ; x_end != x ; ++x ) { *x = m_src ; }

    this_thread.end_barrier();
  }

  template< typename SrcType >
  HostParallelFill( DstType * dst , const SrcType & src ,
                    HostSpace::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    { HostThreadWorker::execute(); }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_HOST_PARALLEL_HPP */

