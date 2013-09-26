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

#include <string>
#include <iostream>
#include <stdexcept>
#include <Kokkos_OpenMP.hpp>
#include <Kokkos_hwloc.hpp>

namespace Kokkos {

namespace {

std::pair<unsigned,unsigned> s_coordinates[ Impl::HostThread::max_thread_count ] ;

//----------------------------------------------------------------------------

unsigned bind_host_thread()
{
  const std::pair<unsigned,unsigned> current = hwloc::get_this_thread_coordinate();
  const unsigned thread_count = (unsigned) omp_get_num_threads();

  unsigned i = 0 ;

  // Match one of the requests:
  for ( i = 0 ; i < thread_count && current != s_coordinates[i] ; ++i );

  if ( thread_count == i ) {
    // Match the NUMA request:
    for ( i = 0 ; i < thread_count && current.first != s_coordinates[i].first ; ++i );
  }

  if ( thread_count == i ) {
    // Match any unclaimed request:
    for ( i = 0 ; i < thread_count && ~0u == s_coordinates[i].first  ; ++i );
  }

  if ( i < thread_count ) {
    if ( ! hwloc::bind_this_thread( s_coordinates[i] ) ) i = thread_count ;
  }

  if ( i < thread_count ) {

#if 0
    if ( current != s_coordinates[i] ) {
      std::cout << "  Kokkos::OpenMP rebinding omp_thread["
                << omp_get_thread_num()
                << "] from ("
                << current.first << "," << current.second
                << ") to ("
                << s_coordinates[i].first << ","
                << s_coordinates[i].second
                << ")" << std::endl ;
    }
#endif

    s_coordinates[i].first  = ~0u ;
    s_coordinates[i].second = ~0u ;
  }

  return i ;
}

//----------------------------------------------------------------------------

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

Impl::HostThread * OpenMP::m_host_threads[ Impl::HostThread::max_thread_count ];

//----------------------------------------------------------------------------

void OpenMP::assert_ready( const char * const function )
{
  const bool error_not_initialized = 0 == m_host_threads[0] ;
  const bool error_in_parallel     = 0 != omp_in_parallel();

  if ( error_not_initialized || error_in_parallel ) {
    std::string msg(function);
    msg.append(" ERROR" );
    if ( error_not_initialized ) {
      msg.append(" : Not initialized");
    }
    if ( error_in_parallel ) {
      msg.append(" : Already within an OMP parallel region");
    }
    throw std::runtime_error(msg);
  }
}

//----------------------------------------------------------------------------

void OpenMP::initialize( const unsigned team_count ,
                         const unsigned threads_per_team ,
                         const unsigned use_numa_count ,
                         const unsigned use_cores_per_numa )
{
  const bool ok_inactive = 0 == m_host_threads[0] ;
  const bool ok_serial   = 0 == omp_in_parallel();

  if ( ok_inactive && ok_serial ) {

    // If user specifies threads_per_team then kill existing threads
    // allocate new threads.

    if ( threads_per_team ) {
      omp_set_num_threads( team_count * threads_per_team ); // Spawn threads
    }

    const unsigned thread_count = (unsigned) omp_get_max_threads();

    const unsigned numa_count       = std::min( use_numa_count ,     hwloc::get_available_numa_count() );
    const unsigned cores_per_numa   = std::min( use_cores_per_numa , hwloc::get_available_cores_per_numa() );
    const unsigned threads_per_core = hwloc::get_available_threads_per_core();

    // If there are more threads than "allowed" cores
    // then omp threads have already been bound (or overallocated)
    // and there is no opportunity to improve locality.

    const bool bind_threads = thread_count <= numa_count * cores_per_numa * threads_per_core ;

    //------------------------------------

    if ( bind_threads ) {

      {
        const std::pair<unsigned,unsigned> team_topo( team_count , threads_per_team );
        const std::pair<unsigned,unsigned> core_topo( numa_count , cores_per_numa );
        const std::pair<unsigned,unsigned> master_coord = hwloc::get_this_thread_coordinate();

        Impl::host_thread_mapping( team_topo , core_topo , core_topo , master_coord , s_coordinates );
      }

      const std::pair<unsigned,unsigned> master_core = s_coordinates[0] ;

      s_coordinates[0] = std::pair<unsigned,unsigned>(~0u,~0u);

#pragma omp parallel
      {
#pragma omp critical
        {
          if ( thread_count != (unsigned) omp_get_num_threads() ) {
            Kokkos::Impl::throw_runtime_exception( "omp_get_max_threads() != omp_get_num_threads()" );
          }
          if ( 0 != omp_get_thread_num() ) {
            const unsigned bind_rank = bind_host_thread();
            Impl::HostThread * const th = new Impl::HostThread();
            Impl::HostThread::set_thread( bind_rank , th );
            m_host_threads[ omp_get_thread_num() ] = th ;
          }
        }
// END #pragma omp critical
      }
// END #pragma omp parallel

      // Bind master thread last
      hwloc::bind_this_thread( master_core );

      {
        Impl::HostThread * const th = new Impl::HostThread();
        Impl::HostThread::set_thread( 0 , th );
        m_host_threads[ 0 ] = th ;
      }
    }
    //------------------------------------
    else {

#pragma omp parallel
      {
#pragma omp critical
        {
          if ( thread_count != (unsigned) omp_get_num_threads() ) {
            Kokkos::Impl::throw_runtime_exception( "omp_get_max_threads() != omp_get_num_threads()" );
          }

          const unsigned rank = (unsigned) omp_get_thread_num() ;

          Impl::HostThread * const th = new Impl::HostThread();
          Impl::HostThread::set_thread( rank , th );
          m_host_threads[ rank ] = th ;
        }
// END #pragma omp critical
      }
// END #pragma omp parallel

    }
    //------------------------------------

    // Set the thread's ranks and counts
    for ( unsigned thread_rank = 0 ; thread_rank < thread_count ; ++thread_rank ) {

      unsigned team_rank    = 0 ;
      unsigned worker_count = 0 ;
      unsigned worker_rank  = 0 ;

      // Distribute threads among gangs:

      // thread_count = k * bin + ( team_count - k ) * ( bin + 1 )
      const unsigned bin  = thread_count / team_count ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = team_count * bin1 - thread_count ;
      const unsigned part = k * bin ;

      if ( thread_rank < part ) {
        team_rank    = thread_rank / bin ;
        worker_rank  = thread_rank % bin ;
        worker_count = bin ;
      }
      else {
        team_rank    = k + ( thread_rank - part ) / bin1 ;
        worker_rank  = ( thread_rank - part ) % bin1 ;
        worker_count = bin1 ;
      }

      Impl::HostThread::get_thread( thread_rank )->
        set_topology( thread_rank , thread_count ,
                      team_rank ,   team_count ,
                      worker_rank , worker_count );
    }

    Impl::HostThread::set_thread_relationships();
  }
  else {
    std::ostringstream msg ;

    msg << "Kokkos::OpenMP::initialize() FAILED" ;

    if ( ! ok_inactive ) {
      msg << " : Device is already active" ;
    }
    if ( ! ok_serial ) {
      msg << " : Called within an OMP parallel region" ;
    }

    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }
}

void OpenMP::finalize()
{
  assert_ready("Kokkos::OpenMP::finalize");

  resize_reduce_scratch(0);

  for ( unsigned i = 0 ; i < Impl::HostThread::max_thread_count ; ++i ) {

    m_host_threads[i] = 0 ;

    Impl::HostThread * const th = Impl::HostThread::clear_thread( i );

    if ( th ) { delete th ; }
  }

  hwloc::unbind_this_thread();
}

//----------------------------------------------------------------------------

void OpenMP::resize_reduce_scratch( unsigned size )
{
  static unsigned m_reduce_size = 0 ;

  const unsigned rem = size % Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size ) || ( m_reduce_size < size ) ) {

#pragma omp parallel
    {
#pragma omp critical
      {
        get_host_thread()->resize_reduce(size);
      }
    }

    m_reduce_size = size ;
  }
}

void * OpenMP::root_reduce_scratch()
{
  Impl::HostThread * const th = m_host_threads[0];
  return th ? th->reduce_data() : 0 ;
}

} // namespace Kokkos


