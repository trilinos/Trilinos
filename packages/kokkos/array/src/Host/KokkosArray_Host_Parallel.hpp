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

#ifndef KOKKOSARRAY_HOST_PARALLEL_HPP
#define KOKKOSARRAY_HOST_PARALLEL_HPP

namespace KokkosArray {
namespace Impl {

class HostInternal ;

//----------------------------------------------------------------------------
/** \brief  A thread within the pool. */

class HostThread {
public:

  typedef Host::size_type size_type ;

  inline size_type rank() const { return m_thread_rank ; }

  inline size_type gang_rank() const { return m_gang_rank ; }

  inline size_type worker_rank() const { return m_worker_rank ; }

  //----------------------------------------------------------------------
  /** \brief  Compute a range of work for this thread's rank */

  inline
  std::pair< size_type , size_type >
    work_range( const size_type work_count ) const
  {
    const size_type reverse_rank    = m_thread_count - ( m_thread_rank + 1 );
    const size_type work_per_thread = ( work_count + m_thread_count - 1 )
                                      / m_thread_count ;
    const size_type work_previous   = work_per_thread * reverse_rank ;
    const size_type work_end        = work_count > work_previous
                                    ? work_count - work_previous : 0 ;

    return std::pair<size_type,size_type>(
      ( work_end > work_per_thread ?
        work_end - work_per_thread : 0 ) , work_end );
  }

  //----------------------------------------------------------------------
  /** \brief  This thread waits for each fan-in thread in the barrier.
   *
   *  A parallel work function must call either barrier or reduce
   *  on all threads at the end of the work function.
   *  Entry condition: in the Active   state
   *  Exit  condition: in the Inactive state
   */
  void barrier()
  {
    // The 'wait' function repeatedly polls the 'thread' state
    // which may reside in a different NUMA region.
    // Thus the fan is intra-node followed by inter-node
    // to minimize inter-node memory access.

    for ( unsigned i = 0 ; i < m_fan_count ; ++i ) {
      m_fan[i]->wait( HostThread::ThreadActive );
    }

    if ( m_thread_rank ) {
      set( HostThread::ThreadInactive );
    }
  }

  //----------------------------------------------------------------------
  /** \brief  This thread participates in the fan-in reduction.
   *
   *  A parallel work function must call either barrier or reduce
   *  on all threads at the end of the work function.
   *  Entry condition: in the Active   state
   *  Exit  condition: in the Inactive state
   */
  template< class ReduceTraits >
  inline
  void reduce( typename ReduceTraits::value_type & update )
  {
    typedef typename ReduceTraits::value_type value_type ;

    // Fan-in reduction of other threads' reduction data.
    // 1) Wait for source thread to complete its work and
    //    set its own state to 'Reducing'.
    // 2) Join source thread reduce data.
    // 3) Release source thread's reduction data and
    //    set the source thread's state to 'Inactive' state.

    for ( unsigned i = 0 ; i < m_fan_count ; ++i ) {
      // Wait until the source thread is finished with its work
      // and enters the reducing state.
      // Join the source thread's reduction data into this thread.
      // Release the source thread.

      m_fan[i]->wait( HostThread::ThreadActive );

      ReduceTraits::join( update, *((const value_type *) m_fan[i]->m_reduce) );

      m_fan[i]->set( HostThread::ThreadInactive );
    }

    if ( m_thread_rank ) {
      // If this is not the root thread then it will give its
      // reduction data to another thread.
      // Set the reduction data and then set the 'Reducing' state.
      // Wait for the other thread to claim reduction data and
      // deactivate this thread.

      m_reduce = & update ;
      set(  HostThread::ThreadReducing );
      wait( HostThread::ThreadReducing );
      m_reduce = NULL ;
    }
  }

  //----------------------------------------------------------------------

private:

  ~HostThread();
  HostThread();

  /** \brief States of a worker thread */
  enum State { ThreadTerminating ///<  Exists, termination in progress
             , ThreadInactive    ///<  Exists, waiting for work
             , ThreadActive      ///<  Exists, performing work
             , ThreadReducing    ///<  Exists, waiting for reduction
             };

  void set(  const State flag ) { m_state = flag ; }
  void wait( const State flag );

  static const unsigned max_fan_count = 16 ;
  static const unsigned max_thread_count = 1 << max_fan_count ;

  HostThread         *  m_fan[ max_fan_count ] ;
  unsigned              m_fan_count ;
  unsigned              m_thread_rank ;
  unsigned              m_thread_count ;
  unsigned              m_gang_rank ;
  unsigned              m_gang_count ;
  unsigned              m_worker_rank ;
  unsigned              m_worker_count ;
  const void * volatile m_reduce ;    ///< Reduction memory
  long         volatile m_state ;     ///< Thread control flag

  friend class HostInternal ;
};

//----------------------------------------------------------------------------
/** \brief  Base class for a parallel driver executing on a thread pool. */

template< class ValueType = void > class HostThreadWorker ;

template<>
class HostThreadWorker<void> {
public:

  /** \brief  Virtual method called on threads */
  virtual void execute_on_thread( HostThread & ) const = 0 ;

  virtual ~HostThreadWorker() {}

protected:

  HostThreadWorker() {}

  static void execute( const HostThreadWorker & );

private:

  HostThreadWorker( const HostThreadWorker & );
  HostThreadWorker & operator = ( const HostThreadWorker & );
};

template< class ValueType >
class HostThreadWorker {
public:

  /** \brief  Virtual method called on threads */
  virtual void execute_on_thread( HostThread & , ValueType & ) const = 0 ;

  virtual ~HostThreadWorker() {}

protected:

  HostThreadWorker() {}

private:

  HostThreadWorker( const HostThreadWorker & );
  HostThreadWorker & operator = ( const HostThreadWorker & );
};

//----------------------------------------------------------------------------

template< typename DstType , typename SrcType  >
class HostParallelCopy : public HostThreadWorker<void> {
public:

        DstType * const m_dst ;
  const SrcType * const m_src ;
  const Host::size_type m_count ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( m_count );
    DstType * const x_end = m_dst + range.second ;
    DstType *       x     = m_dst + range.first ;
    const SrcType * y     = m_src + range.first ;

    for ( ; x_end != x ; ++x , ++y ) { *x = (DstType) *y ; }

    this_thread.barrier();
  }

  HostParallelCopy( DstType * dst , const SrcType * src ,
                    Host::size_type count )
    : HostThreadWorker<void>()
    , m_dst( dst ), m_src( src ), m_count( count )
    { HostThreadWorker<void>::execute( *this ); }
};

template< typename ValueType >
struct DeepCopy<ValueType,Host::memory_space,Host::memory_space> {
  DeepCopy( ValueType * dst , const ValueType * src , size_t count )
  {
    HostParallelCopy< ValueType , ValueType >( dst , src , count );
  }
};


template< typename DstType >
class HostParallelFill : public HostThreadWorker<void> {
public:

  DstType * const m_dst ;
  const DstType   m_src ;
  const Host::size_type m_count ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( m_count );
    DstType * const x_end = m_dst + range.second ;
    DstType *       x     = m_dst + range.first ;

    for ( ; x_end != x ; ++x ) { *x = m_src ; }

    this_thread.barrier();
  }

  template< typename SrcType >
  HostParallelFill( DstType * dst , const SrcType & src ,
                    Host::size_type count )
    : HostThreadWorker<void>()
    , m_dst( dst ), m_src( src ), m_count( count )
    { HostThreadWorker<void>::execute( *this ); }
};

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HOST_PARALLEL_HPP */

