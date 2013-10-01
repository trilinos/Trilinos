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

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Data for OpenMP thread execution */

class OpenMPexec {
public:

  enum { MAX_FAN_COUNT    = 16 };
  enum { MAX_THREAD_COUNT = 1 << MAX_FAN_COUNT };
  enum { VECTOR_LENGTH    = 8 };

  /** \brief  Thread states for team barrier */
  enum { Active , Rendezvous };

private:

  friend class Kokkos::OpenMP ;

  void        * m_reduce ;    ///< Reduction memory
  void        * m_shared ;    ///< Shared memory
  int           m_shared_end ;
  int           m_shared_iter ;
  int volatile  m_state ;
  int           m_fan_team_size ;
  int           m_team_rank ;
  int           m_team_size ;
  int           m_init_league_rank ;
  int           m_init_league_size ;

  int           m_work_league_rank ;
  int           m_work_league_end ;
  int           m_work_league_size ;

  OpenMPexec  * m_fan_team[ MAX_FAN_COUNT ];

  static OpenMPexec * m_thread[ MAX_THREAD_COUNT ];

  OpenMPexec();
  OpenMPexec( const OpenMPexec & );
  OpenMPexec & operator = ( const OpenMPexec & );

public:

  ~OpenMPexec();

  OpenMPexec( const unsigned league_rank ,
              const unsigned league_size ,
              const unsigned team_rank ,
              const unsigned team_size );

  static void finalize();

  static void initialize( const unsigned  team_count ,
                          const unsigned threads_per_team ,
                          const unsigned numa_count ,
                          const unsigned cores_per_numa );

  static void verify_is_process( const char * const );
  static void verify_initialized( const char * const );

  static void resize_reduce_scratch( size_t );
  static void resize_shared_scratch( size_t );

  inline static
  OpenMPexec * get_thread( const unsigned entry )
    { return m_thread[ entry ]; }

  template< class FunctorType >
  inline
  typename ReduceAdapter< FunctorType >::reference_type
  reduce_reference( const FunctorType & ) const
  {
    return ReduceAdapter< FunctorType >::reference( m_reduce );
  }

  template< class FunctorType >
  inline
  typename ReduceAdapter< FunctorType >::reference_type
  reduce_reference( const FunctorType & f , int i ) const
  {
    typedef ReduceAdapter< FunctorType > Reduce ;
    typedef typename Reduce::pointer_type Pointer ;

    return Reduce::reference( Pointer(m_reduce) + i * Reduce::value_count( f ) );
  }

  template< class FunctorType >
  inline
  typename ReduceAdapter< FunctorType >::pointer_type
  reduce_pointer( const FunctorType & ) const
  {
    return typename ReduceAdapter< FunctorType >::pointer_type(m_reduce);
  }

  template< class FunctorType >
  inline
  typename ReduceAdapter< FunctorType >::pointer_type
  reduce_pointer( const FunctorType & f , int i ) const
  {
    typedef ReduceAdapter< FunctorType > Reduce ;
    typedef typename Reduce::pointer_type Pointer ;

    return Pointer(m_reduce) + i * Reduce::value_count( f );
  }

  //----------------------------------------------------------------------
  /** \brief  Compute a range of work for this thread's rank */

  inline
  std::pair< size_t , size_t >
  work_range( const size_t work_count ) const
  {
    typedef integral_constant< size_t , VECTOR_LENGTH - 1 > work_mask ;

    const size_t thread_size = m_team_size * m_work_league_size ;

    // work per thread rounded up and aligned to vector length:

    const size_t work_per_thread =
      ( ( ( work_count + thread_size - 1 ) / thread_size ) + work_mask::value ) & ~(work_mask::value);

    const size_t work_begin = std::min( work_count , work_per_thread * ( m_team_rank + m_team_size * m_work_league_rank ) );
    const size_t work_end   = std::min( work_count , work_per_thread + work_begin );

    return std::pair< size_t , size_t >( work_begin , work_end );
  }

  //----------------------------------------------------------------------

  void * get_shmem( const int );

  void team_barrier()
    {
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        spinwait( m_fan_team[i]->m_state , OpenMPexec::Active );
      }
      if ( m_team_rank ) {
        m_state = Rendezvous ;
        spinwait( m_state , OpenMPexec::Rendezvous );
      }
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        m_fan_team[i]->m_state = OpenMPexec::Active ;
      }
    }

  inline
  void team_work_init( int work_league_size )
    {
      const int work_per_team = ( work_league_size + m_init_league_size - 1 ) / m_init_league_size ;
      m_work_league_rank = std::min( work_league_size , work_per_team * m_init_league_rank );
      m_work_league_end  = std::min( work_league_size , work_per_team + m_work_league_rank );
      m_work_league_size = work_league_size ;
    }

  inline
  bool team_work_avail()
    {
      m_shared_iter = 0 ;
      const bool avail = m_work_league_rank < m_work_league_end ;
      if ( ! avail ) {
        m_work_league_rank = m_init_league_rank ;
        m_work_league_end  = m_init_league_rank + 1 ;
        m_work_league_size = m_init_league_size ;
      }
      return avail ;
    }

  inline
  void team_work_next()
    { if ( ++m_work_league_rank < m_work_league_end ) team_barrier(); }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline OpenMP::OpenMP( Impl::OpenMPexec & e ) : m_exec(e) {}

inline int OpenMP::league_rank() const { return m_exec.m_work_league_rank ; }
inline int OpenMP::league_size() const { return m_exec.m_work_league_size ; }
inline int OpenMP::team_rank() const { return m_exec.m_team_rank ; }
inline int OpenMP::team_size() const { return m_exec.m_team_size ; }

inline void OpenMP::team_barrier() { m_exec.team_barrier() ; }

inline void * OpenMP::get_shmem( const int size ) { return m_exec.get_shmem(size) ; }

} // namespace Kokkos

#endif /* #ifndef KOKKOS_OPENMPEXEC_HPP */

