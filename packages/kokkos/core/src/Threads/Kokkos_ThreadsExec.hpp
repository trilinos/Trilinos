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

//----------------------------------------------------------------------------

template< class > struct ThreadsExecAdapter ;

//----------------------------------------------------------------------------

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

  template< class > friend struct ThreadsExecAdapter ;

  friend class Kokkos::Threads ;

  void        * m_reduce ;    ///< Reduction memory
  void        * m_shared ;    ///< Shared memory
  int volatile  m_state ;
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

  static void execute_sleep( ThreadsExec & , const void * );
  static void execute_reduce_resize( ThreadsExec & , const void * );
  static void execute_shared_resize( ThreadsExec & , const void * );

  ThreadsExec( const ThreadsExec & );
  ThreadsExec & operator = ( const ThreadsExec & );

  static void execute_serial( void (*)( ThreadsExec & , const void * ) );

public:

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

  Threads threads() { return Threads(*this); }

  inline
  std::pair< size_t , size_t >
  work_range( const size_t work_count ) const
  {
    enum { work_align = Kokkos::HostSpace::WORK_ALIGNMENT };
    enum { work_shift = Kokkos::Impl::power_of_two< work_align >::value };
    enum { work_mask  = work_align - 1 };

    // unit_of_work_count = ( work_count + work_mask ) >> work_shift
    // unit_of_work_per_thread = ( unit_of_work_count + thread_count - 1 ) / thread_count
    // work_per_thread = unit_of_work_per_thread * work_align

    const size_t work_per_thread =
      ( ( ( ( work_count + work_mask ) >> work_shift ) + m_thread_size - 1 ) / m_thread_size ) << work_shift ;

    const size_t work_begin = std::min( m_thread_rank * work_per_thread , work_count );
    const size_t work_end   = std::min( work_begin + work_per_thread , work_count );

    return std::pair< size_t , size_t >( work_begin , work_end );
  }

  static void * execute( void (*)( ThreadsExec & , const void * ) , const void * );

  /** \brief  Wait for previous asynchronous functor to
   *          complete and release the Threads device.
   *          Acquire the Threads device and start this functor.
   */
  static void start( void (*)( ThreadsExec & , const void * ) , const void * );

  static void fence();

  static bool sleep();
  static bool wake();
};

//----------------------------------------------------------------------------

template< class FunctorType >
struct ThreadsExecAdapter< ParallelFor< FunctorType , size_t , Threads > >
{
  const FunctorType m_func ;
  const int         m_work ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ThreadsExecAdapter & self = * ((const ThreadsExecAdapter *) arg );

    const std::pair<size_t,size_t> work = exec.work_range( self.m_work );

    for ( size_t iwork = work.first ; iwork < work.second ; ++iwork ) {
      self.m_func( iwork );
    }

    {
      const int n = exec.m_fan_size ;

      for ( int i = 0 ; i < n ; ++i ) {
        ThreadsExec::wait( exec.m_fan[i]->m_state , ThreadsExec::Active );
      }
    }
  }

  ThreadsExecAdapter( const FunctorType & functor , const size_t work )
    : m_func( functor ), m_work( work ) {}
};

//----------------------------------------------------------------------------

template< class FunctorType >
struct ThreadsExecAdapter< ParallelFor< FunctorType , ParallelWorkRequest , Threads > >
{
  const FunctorType m_func ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ThreadsExecAdapter & self = * ((const ThreadsExecAdapter *) arg );

    self.m_func( exec.threads() );

    {
      const int n = exec.m_fan_size ;

      for ( int i = 0 ; i < n ; ++i ) {
        ThreadsExec::wait( exec.m_fan[i]->m_state , ThreadsExec::Active );
      }
    }
  }

  ThreadsExecAdapter( const FunctorType & functor , const ParallelWorkRequest & work )
    : m_func( functor ) {}
};

//----------------------------------------------------------------------------

template< class FunctorType >
struct ThreadsExecAdapter< ParallelReduce< FunctorType , size_t , Threads > >
{
  typedef ReduceAdapter< FunctorType > Reduce ;

  const FunctorType m_func ;
  const int         m_work ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ThreadsExecAdapter & self = * ((const ThreadsExecAdapter *) arg );

    typename Reduce::reference_type update = Reduce::reference( exec.m_reduce );

    self.m_func.init( update ); // Initialize thread-local value

    const std::pair<size_t,size_t> work = exec.work_range( self.m_work );

    for ( size_t iwork = work.first ; iwork < work.second ; ++iwork ) {
      self.m_func( iwork , update );
    }

    {
      const int n = exec.m_fan_size ;

      for ( int i = 0 ; i < n ; ++i ) {

        ThreadsExec & fan = * exec.m_fan[i] ;

        ThreadsExec::wait( fan.m_state , ThreadsExec::Active );

        self.m_func.join( update , Reduce::reference( fan.m_reduce ) );
      }
    }
  }

  ThreadsExecAdapter( const FunctorType & functor , const size_t work )
    : m_func( functor ), m_work( work ) {}
};

//----------------------------------------------------------------------------

template< class FunctorType >
struct ThreadsExecAdapter< ParallelReduce< FunctorType , ParallelWorkRequest , Threads > >
{
  typedef ReduceAdapter< FunctorType > Reduce ;

  const FunctorType  m_func ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ThreadsExecAdapter & self = * ((const ThreadsExecAdapter *) arg );

    typename Reduce::reference_type update = Reduce::reference( exec.m_reduce );

    self.m_func.init( update ); // Initialize thread-local value

    self.m_func( exec.threads() , update );

    {
      const int n = exec.m_fan_size ;

      for ( int i = 0 ; i < n ; ++i ) {

        ThreadsExec & fan = * exec.m_fan[i] ;

        ThreadsExec::wait( fan.m_state , ThreadsExec::Active );

        self.m_func.join( update , Reduce::reference( fan.m_reduce ) );
      }
    }
  }

  ThreadsExecAdapter( const FunctorType & functor , const ParallelWorkRequest & )
    : m_func( functor ) {}
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

inline void Threads::print_configuration( std::ostream & s , bool detail )
{
  Impl::ThreadsExec::print_configuration( s , detail );
}

inline void Threads::fence()
{ Impl::ThreadsExec::fence() ; }

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
{ return m_exec.work_range( work_count ); }

} /* namespace Kokkos */

#endif /* #define KOKKOS_THREADSEXEC_HPP */

