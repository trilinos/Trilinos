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

namespace KokkosArray {
namespace Impl {

class HostInternal ;

void host_wait( volatile int * const state , const int value );

//----------------------------------------------------------------------------
/** \brief  This thread waits for each fan-in thread in the barrier.
 *          Once all threads have fanned in then fan-out reactivation.
 *
 *  A parallel work function must call barrier on all threads.
 */
void host_barrier( HostThread & );

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

