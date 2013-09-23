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

#ifndef KOKKOS_HOST_THREADDATA_HPP
#define KOKKOS_HOST_THREADDATA_HPP

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Traits.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \name  Internal interface to threading runtime. */
/**@{ */

/** \brief  Wait for *flag != value */
void host_thread_wait( volatile int * const flag , const int value );

/** \brief  Wait for *flag != value and yield resources while waiting. */
void host_thread_wait_yield( volatile int * const flag , const int value );

/**@} */

//----------------------------------------------------------------------------
/** \brief  A thread within the pool. */

class HostThread {
public:

  typedef HostSpace::size_type size_type ;

  inline size_type rank() const { return m_thread_rank ; }
  inline size_type count() const { return m_thread_count ; }

  inline size_type gang_rank()    const { return m_gang_rank ; }
  inline size_type gang_count()   const { return m_gang_count ; }
  inline size_type worker_rank()  const { return m_worker_rank ; }
  inline size_type worker_count() const { return m_worker_count ; }

  //----------------------------------------------------------------------
  /** \brief  Compute a range of work for this thread's rank */

  typedef std::pair< size_type , size_type > work_range_type ;

  inline
  work_range_type work_range( const size_type work_count ) const
    {
      enum { work_align = HostSpace::WORK_ALIGNMENT };
      enum { work_shift = power_of_two< work_align >::value };
      enum { work_mask  = work_align - 1 };

      // unit_of_work_count = ( work_count + work_mask ) >> work_shift
      // unit_of_work_per_thread = ( unit_of_work_count + thread_count - 1 ) / thread_count
      // work_per_thread = unit_of_work_per_thread * work_align

      const size_type work_per_thread =
        ( ( ( ( work_count + work_mask ) >> work_shift ) + m_thread_count - 1 ) / m_thread_count ) << work_shift ;

      const size_type work_begin = std::min( m_thread_rank * work_per_thread , work_count );
      const size_type work_end   = std::min( work_begin + work_per_thread , work_count );

      return work_range_type( work_begin , work_end );
    }

  //----------------------------------------------------------------------

  inline
  void * reduce_data() const
    {
#if defined( __INTEL_COMPILER )
__assume_aligned(m_reduce,MEMORY_ALIGNMENT);
#endif
      return m_reduce ;
    }

  //----------------------------------------------------------------------

  inline
  size_type fan_count() { return m_fan_count ; }

  inline
  HostThread & fan( unsigned i ) const { return *m_fan[i]; }

  //----------------------------------------------------------------------
  /** \brief End-of-function barrier */
  inline
  void end_barrier() const
  {
    for ( size_type i = 0 ; i < m_fan_count ; ++i ) {
      host_thread_wait( & m_fan[i]->m_state , HostThread::ThreadActive );
    }
  }

  /** \brief  End-of-function reduction */
  template< class FunctorType >
  void end_reduce( const FunctorType & functor ) const
  {
    typedef ReduceAdapter< FunctorType > Reduce ;

    for ( size_type i = 0 ; i < m_fan_count ; ++i ) {
      host_thread_wait( & m_fan[i]->m_state , HostThread::ThreadActive );
      Reduce::join( functor , m_reduce , m_fan[i]->m_reduce );
    }
  }

  //-----------------------------------

  inline
  void barrier() const
  {
    // Wait until fan-in thread enters the 'Rendezvous' state
    for ( size_type i = 0 ; i < m_fan_count ; ++i ) {
      host_thread_wait( & m_fan[i]->m_state , HostThread::ThreadActive );
    }

    // If not the root thread then enter 'Rendezvous' state
    if ( m_thread_rank ) {
      volatile int & state = const_cast<volatile int &>( m_state );
      state = HostThread::ThreadRendezvous ;
      host_thread_wait( & state , HostThread::ThreadRendezvous );
    }

    // Reset threads to the active state via fan-out.
    for ( size_type i = m_fan_count ; 0 < i ; ) {
      m_fan[--i]->m_state = HostThread::ThreadActive ;
    }
  }

  template< class FunctorType >
  inline
  void reduce( const FunctorType & functor ) const
  {
    typedef ReduceAdapter< FunctorType > Reduce ;

    // Fan-in reduction of other threads' reduction data.

    for ( size_type i = 0 ; i < m_fan_count ; ++i ) {
      // Wait for source thread to complete its work and
      // set its own state to 'Rendezvous'.
      host_thread_wait( & m_fan[i]->m_state , HostThread::ThreadActive );

      // Join source thread reduce data.
      Reduce::join( functor , m_reduce , m_fan[i]->m_reduce );
    }

    if ( m_thread_rank ) {
      // If this is not the root thread then it will give its
      // reduction data to another thread.
      // Set the 'Rendezvous' state.
      // Wait for the other thread to process reduction data
      // and then reactivate this thread.

      volatile int & state = const_cast<volatile int &>( m_state );
      state = HostThread::ThreadRendezvous ;
      host_thread_wait( & state , HostThread::ThreadRendezvous );
    }

    // Reset threads to the active state via fan-out.
    for ( size_type i = m_fan_count ; 0 < i ; ) {
      m_fan[--i]->m_state = HostThread::ThreadActive ;
    }
  }

  //----------------------------------------------------------------------

  inline ~HostThread();
  inline HostThread();

  //----------------------------------------------------------------------

  inline static
  HostThread * get_thread( const unsigned entry )
    { return m_thread[ entry ]; }

  static
  void set_thread( const unsigned rank , HostThread * );

  /**  */
  void set_topology( const unsigned thread_rank , const unsigned thread_count ,
                     const unsigned gang_rank ,   const unsigned gang_count ,
                     const unsigned worker_rank , const unsigned worker_count );


  /** \brief  Setup relationships between threads.
   *
   *  1) Order threads by group, worker, and
   *     pointer as a last resort.
   *  2) Assign thread_rank, gang_rank, and worker_rank
   *     according to ordering.
   *  3) Assign hierarchical fan in/out relationships:
   *     intra-gang followed by inter-gang
   */
  static
  void set_thread_relationships();

  static
  void clear_thread_relationships();

  static
  HostThread * clear_thread( const unsigned entry );

  /** \brief  Resize the reduction scratch space.
   *
   *  1) Must be called by the proper thread for proper first-touch.
   *  2) Must be called while thread is exclusively executing.
   */
  void resize_reduce( unsigned size );

  //----------------------------------------------------------------------

  static const unsigned max_fan_count = 16 ;
  static const unsigned max_thread_count = 1 << max_fan_count ;

  /** \brief States of a worker thread */
  enum State { ThreadTerminating ///<  Exists, termination in progress
             , ThreadInactive    ///<  Exists, waiting for work
             , ThreadActive      ///<  Exists, performing work
             , ThreadRendezvous  ///<  Exists, waiting in a barrier or reduce
             };

  int  volatile m_state ;     ///< Thread control flag

private:

  size_type     m_fan_count ;
  size_type     m_thread_rank ;
  size_type     m_thread_count ;
  size_type     m_gang_rank ;
  size_type     m_gang_count ;
  size_type     m_worker_rank ;
  size_type     m_worker_count ;
  void        * m_reduce ;    ///< Reduction memory
  HostThread  * m_fan[ max_fan_count ] ;

  static HostThread * m_thread[ max_thread_count ];
  static int          m_relations ;

  HostThread( const HostThread & );
  HostThread & operator = ( const HostThread & );

  static void warn_destroy_with_reduce();
};

inline
HostThread::HostThread()
  : m_state( HostThread::ThreadInactive )
  , m_fan_count(0)
  , m_thread_rank(0)
  , m_thread_count(1)
  , m_gang_rank(0)
  , m_gang_count(1)
  , m_worker_rank(0)
  , m_worker_count(1)
  , m_reduce(0)
{}

inline
HostThread::~HostThread()
{
  if ( m_reduce ) warn_destroy_with_reduce();
}

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_HOST_THREADDATA_HPP */

