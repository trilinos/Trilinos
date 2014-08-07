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

#include <stdio.h>
#include <limits>
#include <iostream>
#include <Kokkos_OpenMP.hpp>
#include <Kokkos_hwloc.hpp>
#include <impl/Kokkos_Error.hpp>
#include <iostream>

#ifdef KOKKOS_HAVE_OPENMP

namespace Kokkos {
namespace Impl {
namespace {

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel();

int kokkos_omp_in_critical_region = ( Kokkos::HostSpace::register_in_parallel( kokkos_omp_in_parallel ) , 0 );

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel()
{
#ifndef __CUDA_ARCH__
  return omp_in_parallel() && ! kokkos_omp_in_critical_region ;
#else
  return 0;
#endif
}

unsigned s_threads_per_core = 0 ;
unsigned s_threads_per_numa = 0 ;
bool s_using_hwloc = false;

} // namespace
} // namespace Impl
} // namespace Kokkos


namespace Kokkos {
namespace Impl {

int OpenMPexec::m_map_rank[ OpenMPexec::MAX_THREAD_COUNT ] = { 0 };

OpenMPexec * OpenMPexec::m_pool[ OpenMPexec::MAX_THREAD_COUNT ] = { 0 };

#if 0
void OpenMPexecTeamMember::init( const int league_size , const int team_size )
{
  // Execution is using device-team interface:
  const unsigned pool_size = omp_get_num_threads();

  // Round up team size to be a multiple of threads per core:
  const unsigned team_alloc_core = s_threads_per_core * ( ( team_size + s_threads_per_core - 1 ) / s_threads_per_core );

  // Number of teams which can be allocated:
  const unsigned team_count = pool_size / team_alloc_core ;

  // Number of threads to allocate per team:
  const unsigned team_alloc = pool_size / team_count ;

  const unsigned pool_rank_rev = pool_size - ( m_exec.m_pool_rank + 1 );
  const unsigned team_rank_rev = pool_rank_rev % team_alloc ;

  // May be using fewer threads per team than a multiple of threads per core,
  // some threads will idle.

  if ( int(team_rank_rev) < team_size ) {
    const size_t pool_league_size     = pool_size     / team_alloc ;
    const size_t pool_league_rank_rev = pool_rank_rev / team_alloc ;
    const size_t pool_league_rank     = pool_league_size - ( pool_league_rank_rev + 1 );

    m_team_base   = m_exec.m_pool + team_alloc * pool_league_rank_rev ;

    m_team_shared = execution_space::
     scratch_memory_space( ( (char*) (*m_team_base)->scratch_thread() ) + TEAM_REDUCE_SIZE , m_team_shmem );

    m_team_size        = team_size ;
    m_team_rank        = team_size - ( team_rank_rev + 1 );
    m_league_size      = league_size ;
    m_league_rank      = ( league_size *  pool_league_rank    ) / pool_league_size ;
    m_league_end       = ( league_size * (pool_league_rank+1) ) / pool_league_size ;
  }
}
#endif

void OpenMPexec::verify_is_process( const char * const label )
{
  if ( omp_in_parallel() ) {
    std::string msg( label );
    msg.append( " ERROR: in parallel" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

void OpenMPexec::verify_initialized( const char * const label )
{
  if ( 0 == m_pool[0] ) {
    std::string msg( label );
    msg.append( " ERROR: not initialized" );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

void OpenMPexec::clear_scratch()
{
#pragma omp parallel
  {
    const int rank_rev = m_map_rank[ omp_get_thread_num() ];

#pragma omp critical
    {
      kokkos_omp_in_critical_region = 1 ;

      m_pool[ rank_rev ]->~OpenMPexec();
      HostSpace::decrement( m_pool[ rank_rev ] );
      m_pool[ rank_rev ] = 0 ;

      kokkos_omp_in_critical_region = 0 ;
    }
/* END #pragma omp critical */
  }
/* END #pragma omp parallel */
}

void OpenMPexec::resize_scratch( size_t reduce_size , size_t thread_size )
{
  enum { ALIGN_MASK = Kokkos::Impl::MEMORY_ALIGNMENT - 1 };
  enum { ALLOC_EXEC = ( sizeof(OpenMPexec) + ALIGN_MASK ) & ~ALIGN_MASK };

  const size_t old_reduce_size = m_pool[0] ? m_pool[0]->m_scratch_reduce_end : 0 ;
  const size_t old_thread_size = m_pool[0] ? m_pool[0]->m_scratch_thread_end - m_pool[0]->m_scratch_reduce_end : 0 ;

  reduce_size = ( reduce_size + ALIGN_MASK ) & ~ALIGN_MASK ;
  thread_size = ( thread_size + ALIGN_MASK ) & ~ALIGN_MASK ;

  // Requesting allocation and old allocation is too small:

  const bool allocate = ( old_reduce_size < reduce_size ) ||
                        ( old_thread_size < thread_size );

  if ( allocate ) {
    if ( reduce_size < old_reduce_size ) { reduce_size = old_reduce_size ; }
    if ( thread_size < old_thread_size ) { thread_size = old_thread_size ; }
  }

  const size_t alloc_size = allocate ? ALLOC_EXEC + reduce_size + thread_size : 0 ;
  const int    pool_size  = omp_get_max_threads();

  if ( allocate ) {

    clear_scratch();

#pragma omp parallel
    {
      const int rank_rev = m_map_rank[ omp_get_thread_num() ];
      const int rank     = pool_size - ( rank_rev + 1 );

#pragma omp critical
      {
        kokkos_omp_in_critical_region = 1 ;

        m_pool[ rank_rev ] =
          (OpenMPexec *) HostSpace::allocate( "openmp_scratch" , typeid(unsigned char) , 1 , alloc_size );
        new( m_pool[ rank_rev ] ) OpenMPexec( rank , ALLOC_EXEC , reduce_size , thread_size );

        kokkos_omp_in_critical_region = 0 ;
      }
/* END #pragma omp critical */
    }
/* END #pragma omp parallel */
  }
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

void OpenMP::scratch_memory_space::get_shmem_error()
{
  Kokkos::Impl::throw_runtime_exception( std::string("OpenMPexec::get_shmem FAILED : exceeded shared memory size" ) );
}

KOKKOS_FUNCTION
unsigned OpenMP::team_max()
{
#ifndef __CUDA_ARCH__
  Impl::OpenMPexec::verify_initialized("Kokkos::OpenMP::team_max" );
  Impl::OpenMPexec::verify_is_process("Kokkos::OpenMP::team_max" );

  return Impl::s_threads_per_numa ;
#else
  return 0;
#endif
}

KOKKOS_FUNCTION
unsigned OpenMP::team_recommended()
{
#ifndef __CUDA_ARCH__
  Impl::OpenMPexec::verify_initialized("Kokkos::OpenMP::team_recommended" );
  Impl::OpenMPexec::verify_is_process("Kokkos::OpenMP::team_recommended" );

  return Impl::s_threads_per_core ;
#else
  return 0;
#endif
}

//----------------------------------------------------------------------------

int OpenMP::is_initialized()
{ return 0 != Impl::OpenMPexec::m_pool[0]; }

void OpenMP::initialize( unsigned thread_count ,
                         unsigned use_numa_count ,
                         unsigned use_cores_per_numa )
{
  if(thread_count==0) thread_count = omp_get_max_threads();
  const bool is_initialized = 0 != Impl::OpenMPexec::m_pool[0] ;

  bool thread_spawn_failed = false ;

  if ( ! is_initialized ) {

    // Use hwloc thread pinning if concerned with locality.
    // If spreading threads across multiple NUMA regions.
    // If hyperthreading is enabled.
    Impl::s_using_hwloc = hwloc::available() && (
                            ( 1 < Kokkos::hwloc::get_available_numa_count() ) ||
                            ( 1 < Kokkos::hwloc::get_available_threads_per_core() ) );

    std::pair<unsigned,unsigned> threads_coord[ Impl::OpenMPexec::MAX_THREAD_COUNT ];

    if(Impl::s_using_hwloc)
      hwloc::thread_mapping( "Kokkos::OpenMP::initialize" ,
                           false /* do not allow asynchronous */ ,
                           thread_count ,
                           use_numa_count ,
                           use_cores_per_numa ,
                           threads_coord );

    // Spawn threads:

    omp_set_num_threads( thread_count );

    // Verify OMP interaction:
    if ( int(thread_count) != omp_get_max_threads() ) {
      thread_spawn_failed = true ;
    }

    // Verify spawning and bind threads:
#pragma omp parallel
    {
#pragma omp critical
      {
        if ( int(thread_count) != omp_get_num_threads() ) {
          thread_spawn_failed = true ;
        }

        // Call to 'bind_this_thread' is not thread safe so place this whole block in a critical region.
        // Call to 'new' may not be thread safe as well.

        // Reverse the rank for threads so that the scan operation reduces to the highest rank thread.

        const unsigned omp_rank    = omp_get_thread_num();
        const unsigned thread_r    = Impl::s_using_hwloc ? Kokkos::hwloc::bind_this_thread( thread_count , threads_coord ) : omp_rank ;

        Impl::OpenMPexec::m_map_rank[ omp_rank ] = thread_r ;
      }
/* END #pragma omp critical */
    }
/* END #pragma omp parallel */

    if ( ! thread_spawn_failed ) {
      Impl::s_threads_per_numa = Impl::s_using_hwloc ? thread_count / use_numa_count : thread_count;
      Impl::s_threads_per_core = Impl::s_using_hwloc ? thread_count / ( use_numa_count * use_cores_per_numa ) : 1;

      Impl::OpenMPexec::resize_scratch( 1024 , 1024 );
    }
  }

  if ( is_initialized || thread_spawn_failed ) {
    std::string msg("Kokkos::OpenMP::initialize ERROR");

    if ( is_initialized ) { msg.append(" : already initialized"); }
    if ( thread_spawn_failed ) { msg.append(" : failed spawning threads"); }

    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

//----------------------------------------------------------------------------

void OpenMP::finalize()
{
  Impl::OpenMPexec::verify_initialized( "OpenMP::finalize" );
  Impl::OpenMPexec::verify_is_process( "OpenMP::finalize" );

  Impl::OpenMPexec::clear_scratch();

  omp_set_num_threads(0);

  if(Impl::s_using_hwloc)
    hwloc::unbind_this_thread();
}

//----------------------------------------------------------------------------

void OpenMP::print_configuration( std::ostream & s , const bool detail )
{
  Impl::OpenMPexec::verify_is_process( "OpenMP::print_configuration" );

  s << "Kokkos::OpenMP" ;

#if defined( KOKKOS_HAVE_OPENMP )
  s << " KOKKOS_HAVE_OPENMP" ;
#endif
#if defined( KOKKOS_HAVE_HWLOC )

  const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

  s << " hwloc[" << numa_count << "x" << cores_per_numa << "x" << threads_per_core << "]"
    << " hwloc_binding_" << ( Impl::s_using_hwloc ? "enabled" : "disabled" )
    ;
#endif

  const bool is_initialized = 0 != Impl::OpenMPexec::m_pool[0] ;

  if ( is_initialized ) {
    s << " threads[" << omp_get_max_threads() << "]"
      << " threads_per_numa[" << Impl::s_threads_per_numa << "]"
      << " threads_per_core[" << Impl::s_threads_per_core << "]"
      << std::endl ;

    if ( detail ) {
      std::vector< std::pair<unsigned,unsigned> > coord( omp_get_max_threads() );

#pragma omp parallel
      {
#pragma omp critical
        {
          coord[ omp_get_thread_num() ] = hwloc::get_this_thread_coordinate();
        }
/* END #pragma omp critical */
      }
/* END #pragma omp parallel */

      for ( unsigned i = 0 ; i < coord.size() ; ++i ) {
        s << "  thread omp_rank[" << i << "]"
          << " kokkos_rank[" << Impl::OpenMPexec::m_map_rank[ i ] << "]"
          << " hwloc_coord[" << coord[i].first << "." << coord[i].second << "]"
          << std::endl ;
      }
    }
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

} // namespace Kokkos

#endif //KOKKOS_HAVE_OPENMP
