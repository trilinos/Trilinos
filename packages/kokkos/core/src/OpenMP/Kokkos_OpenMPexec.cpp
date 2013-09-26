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

#include <Kokkos_OpenMP.hpp>
#include <Kokkos_hwloc.hpp>

namespace Kokkos {
namespace Impl {
namespace {



} // namespace
} // namespace Impl
} // namespace Kokkos


namespace Kokkos {
namespace Impl {

OpenMPexec * OpenMPexec::m_thread[ OpenMPexec::MAX_THREAD_COUNT ] = { 0 };

OpenMPexec::OpenMPexec( const unsigned league_rank ,
                        const unsigned league_size ,
                        const unsigned team_rank ,
                        const unsigned team_size )
  : m_reduce(0)
  , m_shared(0)
  , m_shared_end(0)
  , m_shared_iter(0)
  , m_state( OpenMPexec::Active )
  , m_fan_team_size(0)
  , m_team_rank( team_rank )
  , m_team_size( team_size )
  , m_init_league_rank( league_rank )
  , m_init_league_size( league_size )
  , m_work_league_rank( league_rank )
  , m_work_league_end(  league_rank + 1 )
  , m_work_league_size( league_size )
{
  for ( int i = 0 ; i < MAX_FAN_COUNT ; ++i ) { m_fan_team[i] = 0 ; }
}


void OpenMPexec::verify_is_process( const char * const label )
{
  if ( omp_in_parallel() ) {
    std::string msg( label );
    msg.append( " ERROR: in parallel" );
    Kokkos::Impl::throw_runtime_exception( label );
  }
}

void OpenMPexec::verify_initialized( const char * const label )
{
  if ( 1 < omp_get_max_threads() && 0 == m_thread[0] ) {
    std::string msg( label );
    msg.append( " ERROR: not initialized" );
    Kokkos::Impl::throw_runtime_exception( label );
  }
}

void OpenMPexec::initialize( const unsigned team_count ,
                             const unsigned threads_per_team ,
                             const unsigned numa_count ,
                             const unsigned cores_per_numa )
{
  if ( omp_in_parallel() || 1 < omp_get_max_threads() ) {
    Kokkos::Impl::throw_runtime_exception("Kokkos::OpenMP ERROR: initializing but threads are already spawned" );
  }

  const unsigned thread_count = team_count * threads_per_team ;

  if ( thread_count <= 1 ) return ;

  const bool hwloc_avail = Kokkos::hwloc::available();

  const std::pair<unsigned,unsigned>
    hwloc_core_topo( Kokkos::hwloc::get_available_numa_count() ,
                     Kokkos::hwloc::get_available_cores_per_numa() );

  std::pair<unsigned,unsigned> team_topology( team_count , threads_per_team );
  std::pair<unsigned,unsigned> use_core_topology( numa_count , cores_per_numa );
  std::pair<unsigned,unsigned> master_coord = Kokkos::hwloc::get_this_thread_coordinate();

  std::pair<unsigned,unsigned> threads_coord[ OpenMPexec::MAX_THREAD_COUNT ];

  if ( hwloc_avail ) {

    if ( 0 == use_core_topology.first && 0 == use_core_topology.second ) {
      use_core_topology = Kokkos::hwloc::use_core_topology( thread_count );
    }

    Kokkos::hwloc::thread_mapping( team_topology , use_core_topology , hwloc_core_topo , master_coord , threads_coord );
  }

  // Spawn threads:

  omp_set_num_threads( thread_count );

  // Bind threads and allocate thread data:

#pragma omp parallel
  {
    if ( thread_count != (unsigned) omp_get_num_threads() ) {
      Kokkos::Impl::throw_runtime_exception( "Kokkos::OpenMP ERROR: thread_count != omp_get_num_threads()" );
    }
    const unsigned omp_rank    = omp_get_thread_num();
    const unsigned thread_rank = hwloc_avail ? Kokkos::hwloc::bind_this_thread( thread_count , threads_coord ) : omp_rank ;
    const unsigned league_rank = thread_rank / threads_per_team ;
    const unsigned team_rank   = thread_rank % threads_per_team ;
    m_thread[ omp_rank ] = new Impl::OpenMPexec( league_rank , team_count , team_rank , threads_per_team );
  }
/* END #pragma omp parallel */

  // Set threads' fan_team relationships:

#pragma omp parallel
  {
    OpenMPexec & th = * m_thread[ omp_get_thread_num() ];

    for ( int n = 1 ; ( th.m_team_rank + n < th.m_team_size ) &&
                      ( 0 == ( n & th.m_team_rank ) ) ; n <<= 1 ) {

      // Find thread with equal league_rank and team_rank == th.m_team_rank + n

      for ( unsigned i = 0 ; i < thread_count ; ++i ) {
        if ( m_thread[i]->m_init_league_rank == th.m_init_league_rank &&
             m_thread[i]->m_team_rank        == th.m_team_rank + n ) {
          th.m_fan_team[ th.m_fan_team_size++ ] = m_thread[i] ;
          break ;
        }
      }
    }
  }
/* END #pragma omp parallel */

  resize_reduce_scratch( 4096 );
  resize_shared_scratch( 4096 );
}

void OpenMPexec::finalize()
{
  verify_is_process( "OpenMP::finalize" );

  resize_reduce_scratch(0);
  resize_shared_scratch(0);

  for ( int i = 0 ; i < OpenMPexec::MAX_THREAD_COUNT ; ++i ) {
    if ( m_thread[i] ) { delete m_thread[i] ; }
    m_thread[i] = 0 ;
  }

  omp_set_num_threads(0);

  hwloc::unbind_this_thread();
}

void OpenMPexec::resize_reduce_scratch( size_t size )
{
  static size_t s_size = 0 ;

  verify_is_process( "OpenMP::resize_reduce_scratch" );

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size && 0 != s_size ) || s_size < size ) {

#pragma omp parallel
    {
      OpenMPexec & th = * m_thread[ omp_get_thread_num() ];

#pragma omp critical
      {
        if ( th.m_reduce ) {
          HostSpace::decrement( th.m_reduce );
          th.m_reduce = 0 ;
        }

        if ( size ) {
          th.m_reduce = HostSpace::allocate( "reduce_scratch_space" , typeid(unsigned char) , 1 , size );
        }
      }
/* END #pragma omp critical */
    }
/* END #pragma omp parallel */
  }

  s_size = size ;
}

void OpenMPexec::resize_shared_scratch( size_t size )
{
  static size_t s_size = 0 ;

  verify_is_process( "OpenMP::resize_shared_scratch" );

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size && 0 != s_size ) || s_size < size ) {

#pragma omp parallel
    {
      OpenMPexec & th = * m_thread[ omp_get_thread_num() ];

      if ( 0 == th.m_team_rank ) {
#pragma omp critical
        {
          if ( th.m_reduce ) {
            HostSpace::decrement( th.m_shared );
            th.m_shared = 0 ;
          }

          if ( size ) {
            th.m_shared = HostSpace::allocate( "reduce_shared_space" , typeid(unsigned char) , 1 , size );
          }
        }
/* END #pragma omp critical */
        // Push to threads in the same team

        for ( int i = 0 ; i < omp_get_num_threads() ; ++i ) {
          if ( th.m_init_league_rank == m_thread[i]->m_init_league_rank ) {
            m_thread[i]->m_shared = th.m_shared ;
          }
        }
      }
    }
/* END #pragma omp parallel */
  }

  s_size = size ;
}

} // namespace Impl
} // namespace Kokkos

