/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef KOKKOS_DEVICENUMA_PARALLEL_HPP
#define KOKKOS_DEVICENUMA_PARALLEL_HPP

namespace Kokkos {
namespace Impl {

class DeviceNUMAInternal ;

//----------------------------------------------------------------------------
/** \brief  A thread within the pool. */

class DeviceNUMAThread {
public:

  typedef DeviceNUMA::size_type size_type ;

  inline size_type rank() const { return m_thread_rank ; }

  /** \brief  This thread waits for each fan-in thread in the barrier.
   *
   *  All threads must call this function.
   *  Entry condition: in the Active   state
   *  Exit  condition: in the Inactive state
   */
  void barrier()
  {
    DeviceNUMAThread * const thread_beg = m_fan_begin ;
    DeviceNUMAThread *       thread     = m_fan_end ;

    while ( thread_beg < thread ) {
      (--thread)->wait( DeviceNUMAThread::ThreadActive );
    }

    set( DeviceNUMAThread::ThreadInactive );
  }

  /** \brief  This thread participates in the fan-in reduction.
   *
   *  All threads must call this function.
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

    DeviceNUMAThread * const thread_beg    = m_fan_begin ;
    DeviceNUMAThread *       thread_source = m_fan_end ;

    while ( thread_beg < thread_source ) {
      --thread_source ;

      // Wait until the source thread is finished with its work.
      thread_source->wait( DeviceNUMAThread::ThreadActive );

      // Join the source thread's reduction
      ReduceTraits::join( update ,
                          *((const value_type *) thread_source->m_reduce ) );

      thread_source->m_reduce = NULL ;
      thread_source->set( DeviceNUMAThread::ThreadInactive );
    }

    if ( m_thread_rank ) {
      // If this is not the root thread then it will give its
      // reduction data to another thread.
      // Set the reduction data and then set the 'Reducing' state.
      // Wait for the other thread to claim reduction data and
      // deactivate this thread.

      m_reduce = & update ;
      set(  DeviceNUMAThread::ThreadReducing );
      wait( DeviceNUMAThread::ThreadReducing );
    }
  }

  inline
  std::pair< size_type , size_type >
    work_range( const size_type work_count ) const
  {
    const size_type work_per_thread = ( work_count + m_thread_count - 1 ) / m_thread_count ;
    const size_type work_end        = work_count - work_per_thread * m_thread_reverse_rank ;

    return std::pair<size_type,size_type>( work_end > work_per_thread ?
                                           work_end - work_per_thread : 0 , work_end );
  }

  void driver();

private:

  ~DeviceNUMAThread() {}

  DeviceNUMAThread()
    : m_fan_begin( NULL )
    , m_fan_end(   NULL )
    , m_thread_count( 0 )
    , m_thread_rank( 0 )
    , m_thread_reverse_rank( 0 )
    , m_reduce( NULL )
    , m_state( 0 )
    {}

  /** \brief States of a worker thread */
  enum State { ThreadNull = 0    ///<  Does not exist
             , ThreadTerminating ///<  Exists, termination in progress
             , ThreadInactive    ///<  Exists, waiting for work
             , ThreadActive      ///<  Exists, performing work
             , ThreadReducing    ///<  Exists, waiting for reduction
             };

  void set(  const State flag ) { m_state = flag ; }
  void wait( const State flag );

  DeviceNUMAThread    * m_fan_begin ; ///< Begin of thread fan in
  DeviceNUMAThread    * m_fan_end ;   ///< End of thread fan in
  unsigned              m_thread_count ;
  unsigned              m_thread_rank ;
  unsigned              m_thread_reverse_rank ;
  const void * volatile m_reduce ;    ///< Reduction memory
  long         volatile m_state ;     ///< Thread control flag

  friend class DeviceNUMAInternal ;
};

//----------------------------------------------------------------------------
/** \brief  Base class for a parallel driver executing on a thread pool. */

class DeviceNUMAWorker {
public:

  /** \brief  Virtual method called on threads */
  virtual void execute_on_thread( DeviceNUMAThread & ) const = 0 ;

  virtual ~DeviceNUMAWorker() {}

protected:

  DeviceNUMAWorker() {}

  static void execute( const DeviceNUMAWorker & );

private:

  DeviceNUMAWorker( const DeviceNUMAWorker & );
  DeviceNUMAWorker & operator = ( const DeviceNUMAWorker & );
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_DEVICENUMA_PARALLEL_HPP */

