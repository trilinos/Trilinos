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

#ifndef KOKKOS_DEVICEPTHREAD_IMPL_HPP
#define KOKKOS_DEVICEPTHREAD_IMPL_HPP

namespace Kokkos {
namespace Impl {

class DevicePthreadPool ;

//----------------------------------------------------------------------------
/** \brief  Base class for a parallel driver executing on pthread pool. */

class DevicePthreadWorker {
protected:

  /** \brief  Virtual method called by the pthread */
  virtual void execute_on_thread( DevicePthreadController & ) const ;

  virtual ~DevicePthreadWorker() {}

  DevicePthreadWorker() {}

  static DevicePthread::size_type work_per_thread( DevicePthread::size_type );

private:

  DevicePthreadWorker( const DevicePthreadWorker & );
  DevicePthreadWorker & operator = ( const DevicePthreadWorker & );

  friend class DevicePthreadPool ;
};

//----------------------------------------------------------------------------
/** \brief  Control information for a pthread */

class DevicePthreadController {
public:

  DevicePthreadController() {}
  ~DevicePthreadController() {}

  typedef DevicePthread::size_type size_type ;

  inline size_type rank() const { return m_rank ; }

  /** \brief  This thread participate in the thread barrier.
   *
   *  All threads must call this function.
   *  Entry condition: in the Active   state
   *  Exit  condition: in the Inactive state
   */
  void barrier()
  {
    DevicePthreadController * const thread_beg = this[0].m_thread_fan ;
    DevicePthreadController *       thread     = this[1].m_thread_fan ;

    while ( thread_beg < thread ) {
      (--thread)->wait( DevicePthreadController::Active );
    }

    set( DevicePthreadController::Inactive );
  }

  /** \brief  This thread participates in the thread reduction.
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
    //    set its own state to 'OnHold'.
    // 2) Join source thread reduce data.
    // 3) Release source thread's reduction data and
    //    set the source thread's state to 'Inactive' state.

    DevicePthreadController * const thread_begin  = this[0].m_thread_fan ;
    DevicePthreadController *       thread_source = this[1].m_thread_fan ;

    while ( thread_begin < thread_source ) {
      --thread_source ;

      thread_source->wait( DevicePthreadController::Active );

      ReduceTraits::join( update , *((const value_type *) thread_source->m_reduce ) );

      thread_source->m_reduce = NULL ;
      thread_source->set( DevicePthreadController::Inactive );
    }

    if ( m_rank ) {
      // If this is not the root thread then it will give its
      // reduction data to another thread.
      // Set the reduction data and then set the 'OnHold' state.
      // Wait for the other thread to claim reduction data and
      // deactivate this thread.

      m_reduce = & update ;
      set(  DevicePthreadController::OnHold );
      wait( DevicePthreadController::OnHold );
    }
  }

private:

  enum Control { Inactive = 0 , Active = 1 , OnHold = 2 };

  void set(  const Control flag ) { m_control = flag ; }
  void wait( const Control flag );

  DevicePthreadController * m_thread_fan ; ///< Begining of fan in threads
  long                      m_rank ;       ///< Rank for this thread's work
  const void *     volatile m_reduce ;     ///< Reduction memory
  long             volatile m_control ;    ///< Thread control flag

  friend class DevicePthreadPool ;
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_DEVICEPTHREAD_IMPL_HPP */

