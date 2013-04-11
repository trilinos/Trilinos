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

#include <Host/KokkosArray_Host_Thread.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>

namespace KokkosArray {
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
   *  The function must call either 'end_barrier'
   *  or 'end_reduce' as the ending statement.
   */
  virtual void execute_on_thread( HostThread & ) const = 0 ;

  /** \brief End-of-function barrier */
  inline
  void end_barrier( HostThread & thread ) const
  {
    const unsigned n = thread.fan_count();

    for ( unsigned i = 0 ; i < n ; ++i ) {
      host_thread_wait( & thread.fan(i).m_state , HostThread::ThreadActive );
    }
  }

  /** \brief  End-of-function reduction */
  template< class ReduceOper >
  void end_reduce( HostThread & thread , const ReduceOper & reduce ) const
  {
    const unsigned n = thread.fan_count();

    for ( unsigned i = 0 ; i < n ; ++i ) {
      HostThread & th = thread.fan(i);
      host_thread_wait( & th.m_state , HostThread::ThreadActive );
      reduce.join( thread.reduce_data() , th.reduce_data() );
    }
  }

  //-----------------------------------

  inline
  void barrier( HostThread & thread ) const
  {
    const unsigned n = thread.fan_count();

    // Wait until fan-in thread enters the 'Rendezvous' state
    for ( unsigned i = 0 ; i < n ; ++i ) {
      host_thread_wait( & thread.fan(i).m_state , HostThread::ThreadActive );
    }

    // If not the root thread then enter 'Rendezvous' state
    if ( thread.rank() ) {
      thread.m_state = HostThread::ThreadRendezvous ;
      host_thread_wait( & thread.m_state , HostThread::ThreadRendezvous );
    }

    // Reset threads to the active state via fan-out.
    for ( unsigned i = n ; 0 < i ; ) {
      thread.fan(--i).m_state = HostThread::ThreadActive ;
    }
  }

  //-----------------------------------

  template< class ReduceOper >
  inline
  void reduce( HostThread & thread , const ReduceOper & reduce_op ) const
  {
    const unsigned n = thread.fan_count();

    // Fan-in reduction of other threads' reduction data.

    for ( unsigned i = 0 ; i < n ; ++i ) {
      HostThread & th = thread.fan(i);

      // Wait for source thread to complete its work and
      // set its own state to 'Rendezvous'.
      host_thread_wait( & th.m_state , HostThread::ThreadActive );

      // Join source thread reduce data.
      reduce_op.join( thread.reduce_data() , th.reduce_data() );
    }

    if ( thread.rank() ) {
      // If this is not the root thread then it will give its
      // reduction data to another thread.
      // Set the 'Rendezvous' state.
      // Wait for the other thread to process reduction data
      // and then reactivate this thread.

      thread.m_state = HostThread::ThreadRendezvous ;
      host_thread_wait( & thread.m_state , HostThread::ThreadRendezvous );
    }

    // Reset threads to the active state via fan-out.
    for ( unsigned i = n ; 0 < i ; ) {
      thread.fan(--i).m_state = HostThread::ThreadActive ;
    }
  }

protected:

  HostThreadWorker() {}

private:

  HostThreadWorker( const HostThreadWorker & );
  HostThreadWorker & operator = ( const HostThreadWorker & );
};

//----------------------------------------------------------------------------

template< typename DstType , typename SrcType  >
class HostParallelCopy
  : public HostThreadWorker
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

    end_barrier( this_thread );
  }

  HostParallelCopy( DstType * dst , const SrcType * src ,
                    HostSpace::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    { HostThreadWorker::execute(); }
};

template< typename DstType >
class HostParallelFill
  : public HostThreadWorker
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

    end_barrier( this_thread );
  }

  template< typename SrcType >
  HostParallelFill( DstType * dst , const SrcType & src ,
                    HostSpace::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    { HostThreadWorker::execute(); }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HOST_PARALLEL_HPP */

