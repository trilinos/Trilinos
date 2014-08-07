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

#ifndef KOKKOS_OPENMPEXEC_HPP
#define KOKKOS_OPENMPEXEC_HPP

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_spinwait.hpp>

#include <Kokkos_Atomic.hpp>

namespace Kokkos {
namespace Impl {

class OpenMPexecTeamMember ;

//----------------------------------------------------------------------------
/** \brief  Data for OpenMP thread execution */

class OpenMPexec {
public:

  enum { MAX_THREAD_COUNT = 4096 };

private:

  static int          m_map_rank[ MAX_THREAD_COUNT ];
  static OpenMPexec * m_pool[ MAX_THREAD_COUNT ]; // Indexed by: m_pool_rank_rev

  friend class Kokkos::Impl::OpenMPexecTeamMember ;
  friend class Kokkos::OpenMP ;

  int const  m_pool_rank ;
  int const  m_pool_rank_rev ;
  int const  m_scratch_exec_end ;
  int const  m_scratch_reduce_end ;
  int const  m_scratch_thread_end ;

  int volatile  m_barrier_state ;

  OpenMPexec();
  OpenMPexec( const OpenMPexec & );
  OpenMPexec & operator = ( const OpenMPexec & );

  static void clear_scratch();

public:

  static inline int pool_size() { return omp_get_max_threads(); }
  inline int pool_rank() const { return m_pool_rank ; }
  inline int pool_rank_rev() const { return m_pool_rank_rev ; }

  void * scratch_reduce() const { return ((char *) this) + m_scratch_exec_end ; }
  void * scratch_thread() const { return ((char *) this) + m_scratch_reduce_end ; }

  static void resize_scratch( size_t reduce_size , size_t thread_size );

  static inline
  void * scratch_reduce( const int pool_rank_rev )
    { return m_pool[ pool_rank_rev ]->scratch_reduce(); }

  ~OpenMPexec() {}

  OpenMPexec( const int poolRank 
            , const int scratch_exec_size
            , const int scratch_reduce_size
            , const int scratch_thread_size )
    : m_pool_rank( poolRank )
    , m_pool_rank_rev( pool_size() - ( poolRank + 1 ) )
    , m_scratch_exec_end( scratch_exec_size )
    , m_scratch_reduce_end( m_scratch_exec_end   + scratch_reduce_size )
    , m_scratch_thread_end( m_scratch_reduce_end + scratch_thread_size )
    , m_barrier_state(0)
    {}

  static void finalize();

  static void initialize( const unsigned  team_count ,
                          const unsigned threads_per_team ,
                          const unsigned numa_count ,
                          const unsigned cores_per_numa );

  static void verify_is_process( const char * const );
  static void verify_initialized( const char * const );

  inline static
  OpenMPexec * get_thread_omp() { return m_pool[ m_map_rank[ omp_get_thread_num() ] ]; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class OpenMPexecTeamMember {
private:

  enum { TEAM_REDUCE_SIZE = 16 };

  /** \brief  Thread states for team synchronization */
  enum { Active = 0 , Rendezvous = 1 };

  typedef Kokkos::OpenMP                         execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  Impl::OpenMPexec    & m_exec ;
  OpenMPexec * const *  m_team_base ;
  scratch_memory_space  m_team_shared ;
  int                   m_team_shmem ;
  int                   m_team_rank_rev ;
  int                   m_team_rank ;
  int                   m_team_size ;
  int                   m_league_rank ;
  int                   m_league_end ;
  int                   m_league_size ;

  // Fan-in team threads, root of the fan-in which does not block returns true
  KOKKOS_INLINE_FUNCTION
  bool team_fan_in() const
    {
      const int rank_rev = m_team_size - ( m_team_rank + 1 );

      for ( int n = 1 , j ; ( ( j = rank_rev + n ) < m_team_size ) && ! ( rank_rev & n ) ; n <<= 1 ) {
        Impl::spinwait( m_team_base[j]->m_barrier_state , Active );
      }

      if ( rank_rev ) {
        m_exec.m_barrier_state = Rendezvous ;
        Impl::spinwait( m_exec.m_barrier_state , Rendezvous );
      }

      return 0 == rank_rev ;
    }

  KOKKOS_INLINE_FUNCTION
  void team_fan_out() const
    {
      const int rank_rev = m_team_size - ( m_team_rank + 1 );

      for ( int n = 1 , j ; ( ( j = rank_rev + n ) < m_team_size ) && ! ( rank_rev & n ) ; n <<= 1 ) {
        m_team_base[j]->m_barrier_state = Active ;
      }
    }

public:

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
    {
      if ( 1 < m_team_size ) {
        team_fan_in();
        team_fan_out();
      }
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename ArgType >
  KOKKOS_INLINE_FUNCTION ArgType team_scan( const ArgType & value , ArgType * const global_accum ) const
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(ArgType) < TEAM_REDUCE_SIZE , ArgType , void >::type type ;

      volatile type * const work_value  = ((type*) m_exec.scratch_thread());

      *work_value = value ;

      memory_fence();

      if ( team_fan_in() ) {
        // The last thread to synchronize returns true, all other threads wait for team_fan_out()
        // m_team_base[0]                 == highest ranking team member
        // m_team_base[ m_team_size - 1 ] == lowest ranking team member
        //
        // 1) copy from lower to higher rank, initialize lowest rank to zero
        // 2) prefix sum from lowest to highest rank, skipping lowest rank

        type accum = 0 ;

        if ( global_accum ) {
          for ( int i = m_team_size ; i-- ; ) {
            type & val = *((type*) m_team_base[i]->scratch_thread());
            accum += val ;
          }
          accum = atomic_fetch_add( global_accum , accum );
        }

        for ( int i = m_team_size ; i-- ; ) {
          type & val = *((type*) m_team_base[i]->scratch_thread());
          const type offset = accum ;  
          accum += val ;
          val = offset ;
        }

        memory_fence();
      }

      team_fan_out();

      return *work_value ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const
    { return this-> template team_scan<Type>( value , 0 ); }

  //----------------------------------------
  // Private for the driver

private:

  void init( const int league_size , const int team_size );

  typedef execution_space::scratch_memory_space space ;

public:

  template< class WorkArgTag >
  inline
  OpenMPexecTeamMember( Impl::OpenMPexec & exec
                      , const TeamPolicy< execution_space , WorkArgTag > & team
                      , const int shmem_size
                      )
    : m_exec( exec )
    , m_team_base(0)
    , m_team_shared(0,0)
    , m_team_shmem( shmem_size )
    , m_team_rank_rev(0)
    , m_team_rank(0)
    , m_team_size( team.team_size() )
    , m_league_rank(0)
    , m_league_end(0)
    , m_league_size( team.league_size() )
    {
      const int pool_rank_rev        = exec.pool_rank_rev();
      const int pool_team_rank_rev   = pool_rank_rev % team.team_alloc();
      const int pool_league_rank_rev = pool_rank_rev / team.team_alloc();
      const int league_iter_end      = team.league_size() - pool_league_rank_rev * team.team_iter();

      if ( pool_team_rank_rev < m_team_size && 0 < league_iter_end ) {
        m_team_base     = exec.m_pool + team.team_alloc() * pool_league_rank_rev ;
        m_team_rank_rev = pool_team_rank_rev ;
        m_team_rank     = m_team_size - ( m_team_rank_rev + 1 );
        m_league_end    = league_iter_end ;
        m_league_rank   = league_iter_end > team.team_iter() ? league_iter_end - team.team_iter() : 0 ;
        m_team_shared   = space( ( (char*) (*m_team_base)->scratch_thread() ) + TEAM_REDUCE_SIZE , m_team_shmem );
      }
    }

  bool valid() const
    { return m_league_rank < m_league_end ; }

  void next()
    {
      if ( ++m_league_rank < m_league_end ) {
        team_barrier();
        m_team_shared = space( ( (char*) (*m_team_base)->scratch_thread() ) + TEAM_REDUCE_SIZE , m_team_shmem );
      }
    }

  static inline int team_reduce_size() { return TEAM_REDUCE_SIZE ; }
};

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template < class WorkArgTag >
class TeamPolicy< Kokkos::OpenMP , WorkArgTag > {
private:

  int m_league_size ;
  int m_team_size ;
  int m_team_alloc ;
  int m_team_iter ;

  inline void init( const int league_size_request
                  , const int team_size_request )
    {
      const int team_max   = Kokkos::OpenMP::team_max();
      const int team_grain = Kokkos::OpenMP::team_recommended();
      const int pool_size  = omp_get_max_threads();

      m_league_size = league_size_request ;

      m_team_size = team_size_request < team_max ?
                    team_size_request : team_max ;

      // Round team size up to a multiple of 'team_gain'
      const int team_size_grain = team_grain * ( ( m_team_size + team_grain - 1 ) / team_grain );
      const int team_count      = pool_size / team_size_grain ;

      // Constraint : pool_size = m_team_alloc * team_count
      m_team_alloc = pool_size / team_count ;

      // Maxumum number of iterations each team will take:
      m_team_iter  = ( m_league_size + team_count - 1 ) / team_count ;
    }

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::OpenMP             execution_space ; ///< Execution space

  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size_request , int team_size_request )
    { init( league_size_request , team_size_request ); }

  TeamPolicy( int league_size_request , int team_size_request )
    { init( league_size_request , team_size_request ); }

  inline int team_alloc() const { return m_team_alloc ; }
  inline int team_iter()  const { return m_team_iter ; }

  typedef Impl::OpenMPexecTeamMember member_type ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

KOKKOS_INLINE_FUNCTION
unsigned OpenMP::hardware_thread_id() {
#ifndef __CUDA_ARCH__
  return omp_get_thread_num();
#else
  return 0;
#endif
}

KOKKOS_INLINE_FUNCTION
unsigned OpenMP::max_hardware_threads() {
#ifndef __CUDA_ARCH__
  return omp_get_max_threads();
#else
  return 1;
#endif
}

} // namespace Kokkos

#endif /* #ifndef KOKKOS_OPENMPEXEC_HPP */

