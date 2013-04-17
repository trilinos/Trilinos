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

#include <string>
#include <iostream>
#include <stdexcept>
#include <KokkosArray_OpenMP.hpp>
#include <KokkosArray_hwloc.hpp>

namespace KokkosArray {

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

#if 1
    if ( current != s_coordinates[i] ) {
      std::cout << "  KokkosArray::OpenMP rebinding omp_thread["
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

void thread_mapping( const unsigned gang_count )
{
  const std::pair<unsigned,unsigned> core_topo   = hwloc::get_core_topology();
  const std::pair<unsigned,unsigned> master_core = hwloc::get_this_thread_coordinate();

  const unsigned thread_count = (unsigned) omp_get_max_threads();

  // Map threads to cores, gangs, and workers
  for ( unsigned thread_rank = 0 ; thread_rank < thread_count ; ++thread_rank ) {

    unsigned gang_rank    = 0 ;
    unsigned worker_count = 0 ;
    unsigned worker_rank  = 0 ;

    // Distribute threads among gangs:
    {
      // thread_count = k * bin + ( gang_count - k ) * ( bin + 1 )
      const unsigned bin  = thread_count / gang_count ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = gang_count * bin1 - thread_count ;
      const unsigned part = k * bin ;

      if ( thread_rank < part ) {
        gang_rank    = thread_rank / bin ;
        worker_rank  = thread_rank % bin ;
        worker_count = bin ;
      }
      else {
        gang_rank    = k + ( thread_rank - part ) / bin1 ;
        worker_rank  = ( thread_rank - part ) % bin1 ;
        worker_count = bin1 ;
      }
    }

    // Distribute gangs amont NUMA regions:
    {
      // gang_count = k * bin + ( #NUMA - k ) * ( bin + 1 )
      const unsigned bin  = gang_count / core_topo.first ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = core_topo.first * bin1 - gang_count ;
      const unsigned part = k * bin ;

      s_coordinates[ thread_rank ].first =
        ( gang_rank < part )
        ? ( gang_rank / bin )
        : ( k + ( gang_rank - part ) / bin1 );
    }

    // Distribute threads of the gang
    {
      // worker_count = k * bin + ( (#CORE/NUMA) - k ) * ( bin + 1 )
      const unsigned bin  = worker_count / core_topo.second ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = core_topo.second * bin1 - worker_count ;
      const unsigned part = k * bin ;

      s_coordinates[thread_rank].second =
        ( worker_rank < part )
        ? ( worker_rank / bin )
        : ( k + ( worker_rank - part ) / bin1 );
    }
  }

  //------------------------------------
  // Don't move the master thread:
  // Reserve entry #0 for the master thread by trading out entry #0
  // with the best-fit entry for the master threads current location.
  {
    unsigned i = 0 ;

    // First try for an exact match:
    for ( ; i < thread_count && master_core != s_coordinates[i] ; ++i );

    if ( i == thread_count ) {
      // Exact match failed: take the first entry in the NUMA region.
      for ( i = 0 ; i < thread_count &&
                    master_core.first != s_coordinates[i].first ; ++i );
    }

    if ( i == thread_count ) i = 0 ;

    s_coordinates[i] = s_coordinates[0] ;
    s_coordinates[0] = master_core ;
  }

  //------------------------------------

#if 0
  for ( unsigned thread_rank = 0 ; thread_rank < thread_count ; ++thread_rank ) {
    std::cout << "KokkosArray::OpenMP::thread_mapping"
              << " rank[" << thread_rank << "]"
              << " core(" << s_coordinates[thread_rank].first
              << "," << s_coordinates[thread_rank].second
              << ")" << std::endl ;
  }
#endif

}

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

void OpenMP::initialize( const unsigned gang_count )
{
  const bool ok_inactive = 0 == m_host_threads[0] ;
  const bool ok_serial   = 0 == omp_in_parallel();

  if ( ok_inactive && ok_serial ) {

    const unsigned thread_count = (unsigned) omp_get_max_threads();

    thread_mapping( gang_count );

    const std::pair<unsigned,unsigned> master_core = s_coordinates[0] ;

    s_coordinates[0] = std::pair<unsigned,unsigned>(~0u,~0u);

#pragma omp parallel
    {
#pragma omp critical
      {
        if ( thread_count != (unsigned) omp_get_num_threads() ) {
          KokkosArray::Impl::throw_runtime_exception( "omp_get_max_threads() != omp_get_num_threads()" );
        }
        if ( 0 != omp_get_thread_num() ) {
          const unsigned bind_rank = bind_host_thread();
          Impl::HostThread * const th = new Impl::HostThread();
          Impl::HostThread::set_thread( bind_rank , th );
          m_host_threads[ omp_get_thread_num() ] = th ;
        }
      }
    }
// END #pragma omp parallel

    // Bind master thread last
    hwloc::bind_this_thread( master_core );

    {
      Impl::HostThread * const th = new Impl::HostThread();
      Impl::HostThread::set_thread( 0 , th );
      m_host_threads[ 0 ] = th ;
    }

    // Set the thread's ranks and counts
    for ( unsigned thread_rank = 0 ; thread_rank < thread_count ; ++thread_rank ) {

      unsigned gang_rank    = 0 ;
      unsigned worker_count = 0 ;
      unsigned worker_rank  = 0 ;

      // Distribute threads among gangs:

      // thread_count = k * bin + ( gang_count - k ) * ( bin + 1 )
      const unsigned bin  = thread_count / gang_count ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = gang_count * bin1 - thread_count ;
      const unsigned part = k * bin ;

      if ( thread_rank < part ) {
        gang_rank    = thread_rank / bin ;
        worker_rank  = thread_rank % bin ;
        worker_count = bin ;
      }
      else {
        gang_rank    = k + ( thread_rank - part ) / bin1 ;
        worker_rank  = ( thread_rank - part ) % bin1 ;
        worker_count = bin1 ;
      }

      Impl::HostThread::get_thread( thread_rank )->
        set_topology( thread_rank , thread_count ,
                      gang_rank ,   gang_count ,
                      worker_rank , worker_count );
    }

    Impl::HostThread::set_thread_relationships();
  }
  else {
    std::ostringstream msg ;

    msg << "KokkosArray::OpenMP::initialize() FAILED" ;

    if ( ! ok_inactive ) {
      msg << " : Device is already active" ;
    }
    if ( ! ok_serial ) {
      msg << " : Called within an OMP parallel region" ;
    }

    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }
}

void OpenMP::finalize()
{
  assert_ready("KokkosArray::OpenMP::finalize");

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

  const unsigned rem = size % HostSpace::MEMORY_ALIGNMENT ;

  if ( rem ) size += HostSpace::MEMORY_ALIGNMENT - rem ;

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

} // namespace KokkosArray


