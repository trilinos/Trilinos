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

#ifndef KOKKOS_THREADSEXEC_HPP
#define KOKKOS_THREADSEXEC_HPP

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class ThreadsExec {
public:

  enum { MAX_FAN_COUNT    = 16 };
  enum { MAX_THREAD_COUNT = 1 << MAX_FAN_COUNT };

  /** \brief States of a worker thread */
  enum { Terminating ///<  Termination in progress
       , Inactive    ///<  Exists, waiting for work
       , Active      ///<  Exists, performing work
       , Rendezvous  ///<  Exists, waiting in a barrier or reduce
       };

private:

  friend class Kokkos::Threads ;

  void        * m_reduce ;    ///< Reduction memory
  void        * m_shared ;    ///< Shared memory
  int           m_state ;
  int           m_fan_size ;
  int           m_team_rank ;
  int           m_team_size ;
  int           m_league_rank ;
  int           m_league_size ;
  int           m_thread_rank ; // m_team_rank + m_team_size * m_league_rank
  int           m_thread_size ; // m_team_size * m_league_size
  ThreadsExec * m_fan[ MAX_FAN_COUNT ] ;

  static void activate_threads();
  static void global_lock();
  static void global_unlock();
  static void wait( volatile int & , const int );
  static void wait_yield( volatile int & , const int );
  static bool spawn();

  static void execute_sleep( Threads , const void * );
  static void execute_reduce_resize( Threads , const void * );
  static void execute_shared_resize( Threads , const void * );

  ThreadsExec( const ThreadsExec & );
  ThreadsExec & operator = ( const ThreadsExec & );

  static void execute_serial( void (*)( Threads , const void * ) );

  /** \brief  Function executed by each thread */
  template< class FunctorAdapter >
  static void execute( Threads dev , const void * arg )
  {
    const FunctorAdapter & f = * ((const FunctorAdapter *) arg );

    // A reduction must initialize reduction value.
    f.init( dev.m_exec.m_reduce );

    f.apply( dev , dev.m_exec.m_reduce );

    const int n = dev.m_exec.m_fan_size ;

    for ( int i = 0 ; i < n ; ++i ) {
      wait( dev.m_exec.m_fan[i]->m_state , ThreadsExec::Active );

      // A reduction must join reduction values.
      f.join( dev.m_exec.m_reduce , dev.m_exec.m_fan[i]->m_reduce );
    }

    if ( ! dev.m_exec.m_thread_rank ) {
      // A reduction may perform serial finalization of reduction value.
      f.final( dev.m_exec.m_reduce );
    }

    dev.m_exec.m_state = ThreadsExec::Inactive ;
  }

  static void start( void (*)( Threads dev , const void * ) , const void * );

public:

  /** \brief  Wait for previous asynchronous functor to
   *          complete and release the Threads device.
   *          Acquire the Threads device and start this functor.
   */
  template< class FunctorAdapter >
  static void start( const FunctorAdapter & f )
    { start( & ThreadsExec::template execute<FunctorAdapter> , & f ); }

  static void acquire( const void * ); ///<  Acquire the Threads device
  static void release( const void * ); ///<  Release the Threads device

  static void driver(void);

  ~ThreadsExec();
  ThreadsExec();

  static void set_threads_relationships( const std::pair<unsigned,unsigned> team_topo ,
                                         ThreadsExec * threads[] );

  static void resize_reduce_scratch( size_t );
  static void resize_shared_scratch( size_t );

  static void * root_reduce_scratch();

  static bool is_process();

  static void verify_is_process( const std::string & );

  static void initialize( const std::pair<unsigned,unsigned> team_topo ,
                                std::pair<unsigned,unsigned> core_topo );

  static void finalize();

  static void print_configuration( std::ostream & , const bool detail = false );

  //------------------------------------

  static bool sleep();
  static bool wake();
  static void fence();
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {

inline void Threads::initialize( 
  const std::pair<unsigned,unsigned> team_league_size ,
  const std::pair<unsigned,unsigned> hardware_topology )
{
  Impl::ThreadsExec::initialize( team_league_size , hardware_topology );
}

inline void Threads::finalize()
{
  Impl::ThreadsExec::finalize();
}

inline void Threads::print_configuration( std::ostream & s )
{
  Impl::ThreadsExec::print_configuration( s );
}

inline int Threads::league_rank() const
{ return m_exec.m_league_rank ; }

inline int Threads::league_size() const
{ return m_exec.m_league_size ; }

inline int Threads::team_rank() const
{ return m_exec.m_team_rank ; }

inline int Threads::team_size() const
{ return m_exec.m_team_size ; }

inline
std::pair< size_t , size_t >
Threads::work_range( const size_t work_count ) const
{
  enum { work_align = Kokkos::HostSpace::WORK_ALIGNMENT };
  enum { work_shift = Kokkos::Impl::power_of_two< work_align >::value };
  enum { work_mask  = work_align - 1 };

  // unit_of_work_count = ( work_count + work_mask ) >> work_shift
  // unit_of_work_per_thread = ( unit_of_work_count + thread_count - 1 ) / thread_count
  // work_per_thread = unit_of_work_per_thread * work_align

  const size_t work_per_thread =
    ( ( ( ( work_count + work_mask ) >> work_shift ) + m_exec.m_thread_size - 1 ) / m_exec.m_thread_size ) << work_shift ;

  const size_t work_begin = std::min( m_exec.m_thread_rank * work_per_thread , work_count );
  const size_t work_end   = std::min( work_begin + work_per_thread , work_count );

  return std::pair< size_t , size_t >( work_begin , work_end );
}

} /* namespace Kokkos */

#endif /* #define KOKKOS_THREADSEXEC_HPP */

