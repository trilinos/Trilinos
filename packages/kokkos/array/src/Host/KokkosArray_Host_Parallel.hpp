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

  std::pair< size_type , size_type >
    work_range( const size_type work_count ) const ;

  //----------------------------------------------------------------------
  /** \brief  This thread waits for each fan-in thread in the barrier.
   *          Once all threads have fanned in then fan-out reactivation.
   *
   *  A parallel work function must call barrier on all threads.
   */
  void barrier();

  //----------------------------------------------------------------------

  inline
  void * reduce_data() const { return m_reduce ; }

  //----------------------------------------------------------------------
  /** \brief  This thread participates in the fan-in reduction.
   *
   */
  template< class ReduceOper >
  inline
  void reduce( const ReduceOper & reduce )
    {
      // Fan-in reduction of other threads' reduction data.

      for ( unsigned i = 0 ; i < m_fan_count ; ++i ) {

        // Wait for source thread to complete its work and
        // set its own state to 'Rendezvous'.
        m_fan[i]->wait( HostThread::ThreadActive );

        // Join source thread reduce data.
        reduce.join( m_reduce , m_fan[i]->m_reduce );

        // Reset the source thread to 'Active' state.
        m_fan[i]->set( HostThread::ThreadActive );
      }

      if ( m_thread_rank ) {
        // If this is not the root thread then it will give its
        // reduction data to another thread.
        // Set the 'Rendezvous' state.
        // Wait for the other thread to process reduction data
        // and then reactivate this thread.

        set(  HostThread::ThreadRendezvous );
        wait( HostThread::ThreadRendezvous );
      }
      else {
        reduce.finalize( m_reduce );
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
             , ThreadRendezvous  ///<  Exists, waiting in a barrier or reduce
             };

  void set(  const State flag ) { m_state = flag ; }
  void wait( const State flag );

  static const unsigned max_fan_count = 16 ;
  static const unsigned max_thread_count = 1 << max_fan_count ;

  HostThread   *  m_fan[ max_fan_count ] ;
  size_type       m_fan_count ;
  size_type       m_thread_rank ;
  size_type       m_thread_count ;
  size_type       m_gang_rank ;
  size_type       m_gang_count ;
  size_type       m_worker_rank ;
  size_type       m_worker_count ;
  void         *  m_reduce ;    ///< Reduction memory
  long   volatile m_state ;     ///< Thread control flag

  friend class HostInternal ;
  friend class HostThreadWorker ;
};

//----------------------------------------------------------------------------
/** \brief  Base class for a parallel driver executing on a thread pool. */

struct HostThreadWorker {
public:

  virtual ~HostThreadWorker() {}

  void execute() const ;

  /** \brief  Virtual method called on threads */
  virtual void execute_on_thread( HostThread & ) const = 0 ;

  /** \brief Wait for fanin/fanout threads to deactivate themselves. */
  void fanin_deactivation( HostThread & thread ) const ;

protected:

  HostThreadWorker() {}

private:

  HostThreadWorker( const HostThreadWorker & );
  HostThreadWorker & operator = ( const HostThreadWorker & );
};

//----------------------------------------------------------------------------

template< class WorkerType >
struct HostParallelLaunch {
private:

  struct ThreadWorker : public HostThreadWorker {
    const WorkerType & m_worker ;

    void execute_on_thread( HostThread & thread ) const
      { m_worker( thread ); }

    ThreadWorker( const WorkerType & worker )
      : m_worker( worker ) {}
  };

public:

  HostParallelLaunch( const WorkerType & worker )
    { ThreadWorker( worker ).execute(); }
};

//----------------------------------------------------------------------------

template< typename DstType , typename SrcType  >
class HostParallelCopy {
public:

        DstType * const m_dst ;
  const SrcType * const m_src ;
  const Host::size_type m_count ;

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( m_count );
    DstType * const x_end = m_dst + range.second ;
    DstType *       x     = m_dst + range.first ;
    const SrcType * y     = m_src + range.first ;

    for ( ; x_end != x ; ++x , ++y ) { *x = (DstType) *y ; }
  }

  HostParallelCopy( DstType * dst , const SrcType * src ,
                    Host::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    { HostParallelLaunch< HostParallelCopy >( *this ); }
};

template< typename DstType >
class HostParallelFill {
public:

  DstType * const m_dst ;
  const DstType   m_src ;
  const Host::size_type m_count ;

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( m_count );
    DstType * const x_end = m_dst + range.second ;
    DstType *       x     = m_dst + range.first ;

    for ( ; x_end != x ; ++x ) { *x = m_src ; }
  }

  template< typename SrcType >
  HostParallelFill( DstType * dst , const SrcType & src ,
                    Host::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    { HostParallelLaunch< HostParallelFill >( *this ); }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HOST_PARALLEL_HPP */

