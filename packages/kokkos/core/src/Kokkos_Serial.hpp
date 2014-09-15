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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Serial.hpp
/// \brief Declaration and definition of Kokkos::Serial device.

#ifndef KOKKOS_SERIAL_HPP
#define KOKKOS_SERIAL_HPP

#include <cstddef>
#include <iosfwd>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/// \class Serial
/// \brief Kokkos device for non-parallel execution
///
/// A "device" represents a parallel execution model.  It tells Kokkos
/// how to parallelize the execution of kernels in a parallel_for or
/// parallel_reduce.  For example, the Threads device uses Pthreads or
/// C++11 threads on a CPU, the OpenMP device uses the OpenMP language
/// extensions, and the Cuda device uses NVIDIA's CUDA programming
/// model.  The Serial device executes "parallel" kernels
/// sequentially.  This is useful if you really do not want to use
/// threads, or if you want to explore different combinations of MPI
/// and shared-memory parallel programming models.
class Serial {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! The tag (what type of kokkos_object is this).
  typedef Impl::ExecutionSpaceTag  kokkos_tag ;

  //! The device type (same as this class).
  typedef Serial                device_type ;
  typedef Serial                execution_space ;
  //! The size_type typedef best suited for this device.
  typedef HostSpace::size_type  size_type ;
  //! This device's preferred memory space.
  typedef HostSpace             memory_space ;
  //! This device's preferred array layout.
  typedef LayoutRight           array_layout ;
  /// \brief This device's host mirror type.
  ///
  /// Serial is a host device, so the host mirror type is the same as
  /// the device type itself.
  typedef Serial                host_mirror_device_type ;

  /// \brief  Scratch memory space
  typedef ScratchMemorySpace< Kokkos::Serial >  scratch_memory_space ;

  //@}

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.
  ///
  /// For the Serial device, this method <i>always</i> returns false,
  /// because parallel_for or parallel_reduce with the Serial device
  /// always execute sequentially.
  inline static int in_parallel() { return false ; }

  /** \brief  Set the device in a "sleep" state.
   *
   * This function sets the device in a "sleep" state in which it is
   * not ready for work.  This may consume less resources than if the
   * device were in an "awake" state, but it may also take time to
   * bring the device from a sleep state to be ready for work.
   *
   * \return True if the device is in the "sleep" state, else false if
   *   the device is actively working and could not enter the "sleep"
   *   state.
   */
  static bool sleep();

  /// \brief Wake the device from the 'sleep' state so it is ready for work.
  ///
  /// \return True if the device is in the "ready" state, else "false"
  ///  if the device is actively working (which also means that it's
  ///  awake).
  static bool wake();

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void fence() {}

  static void initialize( unsigned threads_count = 1 ,
                          unsigned use_numa_count = 0 ,
                          unsigned use_cores_per_numa = 0 ,
                          bool allow_asynchronous_threadpool = false) {
    (void) threads_count;
    (void) use_numa_count;
    (void) use_cores_per_numa;
    (void) allow_asynchronous_threadpool;
  }

  static int is_initialized() { return 1 ; }

  //! Free any resources being consumed by the device.
  static void finalize() {}

  //! Print configuration information to the given output stream.
  static void print_configuration( std::ostream & , const bool detail = false );

  //--------------------------------------------------------------------------

  inline static int thread_pool_size( int = 0 ) { return 1 ; }
  KOKKOS_INLINE_FUNCTION static int thread_pool_rank() { return 0 ; }

  //--------------------------------------------------------------------------

  inline static unsigned hardware_thread_id() { return thread_pool_rank(); }
  inline static unsigned max_hardware_threads() { return thread_pool_size(0); }

  static inline int team_max()         { return thread_pool_size(1) ; }
  static inline int team_recommended() { return thread_pool_size(2); }

  //--------------------------------------------------------------------------

  static void * scratch_memory_resize( unsigned reduce_size , unsigned shared_size );

  //--------------------------------------------------------------------------
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template<>
struct VerifyExecutionCanAccessMemorySpace
  < Kokkos::Serial::memory_space
  , Kokkos::Serial::scratch_memory_space
  >
{
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

namespace SerialImpl {

struct Sentinel {

  void *   m_scratch ;
  unsigned m_reduce_end ;
  unsigned m_shared_end ;

  Sentinel();
  ~Sentinel();
  static Sentinel & singleton();
};

inline
unsigned align( unsigned n );
}
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

class SerialTeamMember {
private:
  typedef Kokkos::ScratchMemorySpace< Kokkos::Serial > scratch_memory_space ;
  const scratch_memory_space  m_space ;
  const int                   m_league_rank ;
  const int                   m_league_size ;

  SerialTeamMember & operator = ( const SerialTeamMember & );

public:

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space & team_shmem() const { return m_space ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {}

  template< class JoinOp >
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & ) const
    { return value ; }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const
    {
      const Type tmp = global_accum ? *global_accum : Type(0) ;
      if ( global_accum ) { *global_accum += value ; }
      return tmp ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & ) const
    { return Type(0); }

  //----------------------------------------
  // Execution space specific:

  SerialTeamMember( int arg_league_rank
                  , int arg_league_size 
                  , int arg_shared_size
                  );
};

template<unsigned int VectorLength>
class SerialTeamVectorMember {
private:
  typedef Kokkos::ScratchMemorySpace< Kokkos::Serial > scratch_memory_space ;
  const scratch_memory_space  m_space ;
  const int                   m_league_rank ;
  const int                   m_league_size ;

  SerialTeamVectorMember & operator = ( const SerialTeamVectorMember & );

public:

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space & team_shmem() const { return m_space ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {}

  template< class JoinOp >
  typename JoinOp::value_type team_reduce( const typename JoinOp::value_type & value
                                         , const JoinOp & ) const
    { return value ; }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const
    {
      const Type tmp = global_accum ? *global_accum : Type(0) ;
      if ( global_accum ) { *global_accum += value ; }
      return tmp ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & ) const
    { return Type(0); }

#ifdef KOKKOS_HAVE_CXX11

  /** \brief  Guarantees execution of op() with only a single vector lane of this thread. */
  template< class Operation >
  KOKKOS_INLINE_FUNCTION void vector_single(const Operation & op) const {
    op();
  }

  /** \brief  Intra-thread vector parallel for. Executes op(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
   * This functionality requires C++11 support.*/
  template< typename iType, class Operation, typename ValueType >
  KOKKOS_INLINE_FUNCTION void vector_par_for(const iType n, const Operation & op) const {
    #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
    #pragma ivdep
    #endif
    for(int i=0; i<n ; i++) {
      op(i);
    }
  }

  /** \brief  Intra-thread vector parallel reduce. Executes op(iType i, ValueType & val) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
   * val is performed and put into result. This functionality requires C++11 support.*/
  template< typename iType, class Operation, typename ValueType >
  KOKKOS_INLINE_FUNCTION void vector_par_reduce(const iType n, const Operation & op, ValueType& result) const {

    result = ValueType();

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for(int i=0; i<n ; i++) {
      ValueType tmp = ValueType();
      op(i,tmp);
      result+=tmp;
    }
  }

  /** \brief  Intra-thread vector parallel reduce. Executes op(iType i, ValueType & val) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
   * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
   * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
   * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
   * '1 for *'). This functionality requires C++11 support.*/
  template< typename iType, class Operation, typename ValueType, class JoinType >
  KOKKOS_INLINE_FUNCTION void vector_par_reduce(const iType n, const Operation & op, ValueType& init_result, const JoinType & join) const {

    ValueType result = init_result;

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for(int i=0; i<n ; i++) {
      ValueType tmp = init_result;
      op(i,tmp);
      join(result,tmp);
    }
    init_result = result;
  }


  /** \brief  Intra-thread vector parallel exclusive prefix sum. Executes op(iType i, ValueType & val, bool final)
   *          for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan operation is performed.
   * Depending on the target execution space the operator might be called twice: once with final=false
   * and once with final=true. When final==true val contains the prefix sum value. The contribution of this
   * "i" needs to be added to val no matter whether final==true or not. In a serial execution
   * (i.e. team_size==1) the operator is only called once with final==true. Scan_val will be set
   * to the final sum value over all vector lanes.
   * This functionality requires C++11 support.*/
  template< typename iType, class Operation, typename ValueType >
  KOKKOS_INLINE_FUNCTION  void vector_par_scan(const iType n, const Operation & op, ValueType& scan_val) const {

    scan_val = ValueType();

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for(int i=0; i<n ; i++) {
      op(i,scan_val,true);
    }
  }
#endif
  //----------------------------------------
  // Execution space specific:

  SerialTeamVectorMember( int arg_league_rank
                  , int arg_league_size
                  , int arg_shared_size
                  )
  : m_space( ((char *) SerialImpl::Sentinel::singleton().m_scratch) + SerialImpl::Sentinel::singleton().m_reduce_end
           , arg_shared_size )
  , m_league_rank( arg_league_rank )
  , m_league_size( arg_league_size )
  {}
};
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

template < class WorkArgTag >
class TeamPolicy< Kokkos::Serial , WorkArgTag > {
private:

  const int m_league_size ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Serial             execution_space ; ///< Execution space

  inline int team_size() const { return 1 ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size_request , int /* team_size_request */ )
    : m_league_size( league_size_request )
    { }

  TeamPolicy( int league_size_request , int /* team_size_request */ )
    : m_league_size( league_size_request )
    { }

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & ) { return 1 ; }

  typedef Impl::SerialTeamMember  member_type ;
};

template < unsigned int VectorLength, class WorkArgTag >
class TeamVectorPolicy< VectorLength, Kokkos::Serial , WorkArgTag > {
private:

  const int m_league_size ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Serial             execution_space ; ///< Execution space

  inline int team_size() const { return 1 ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamVectorPolicy( execution_space & , int league_size_request , int /* team_size_request */ )
    : m_league_size( league_size_request )
    { }

  TeamVectorPolicy( int league_size_request , int /* team_size_request */ )
    : m_league_size( league_size_request )
    { }

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & ) { return 1 ; }

  typedef Impl::SerialTeamVectorMember<VectorLength>  member_type ;
};

} /* namespace Kokkos */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 >
                 , Kokkos::Serial
                 >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 > Policy ;

public:
  // work tag is void
  template< class PType >
  inline
  ParallelFor( typename Impl::enable_if<
                 ( Impl::is_same< PType , Policy >::value &&
                   Impl::is_same< typename PType::work_tag , void >::value
                 ), const FunctorType & >::type functor
             , const PType & policy )
    {
      const typename PType::member_type e = policy.end();
      for ( typename PType::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( i );
      }
    }

  // work tag is non-void
  template< class PType >
  inline
  ParallelFor( typename Impl::enable_if<
                 ( Impl::is_same< PType , Policy >::value &&
                   ! Impl::is_same< typename PType::work_tag , void >::value
                 ), const FunctorType & >::type functor
             , const PType & policy )
    {
      const typename PType::member_type e = policy.end();
      for ( typename PType::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( typename PType::work_tag() , i );
      }
    }
};

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 >
                    , Kokkos::Serial
                    >
{
public:
  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 > Policy ;

  typedef ReduceAdapter< FunctorType >  Reduce ;
  typedef typename Reduce::pointer_type pointer_type ;

  // Work tag is void
  template< class ViewType , class PType >
  ParallelReduce( typename Impl::enable_if<
                    ( Impl::is_view< ViewType >::value &&
                      Impl::is_same< typename ViewType::memory_space , HostSpace >::value &&
                      Impl::is_same< PType , Policy >::value &&
                      Impl::is_same< typename PType::work_tag , void >::value
                    ), const FunctorType & >::type functor
                , const PType     & policy
                , const ViewType  & result
                )
    {
      pointer_type result_ptr = result.ptr_on_device();

      if ( ! result_ptr ) {
        result_ptr = (pointer_type)
          Kokkos::Serial::scratch_memory_resize( Reduce::value_size( functor ) , 0 );
      }

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      const typename PType::member_type e = policy.end();
      for ( typename PType::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( i , update );
      }

      Reduce::final( functor , result_ptr );
    }

  // Work tag is non-void
  template< class ViewType , class PType >
  ParallelReduce( typename Impl::enable_if<
                    ( Impl::is_view< ViewType >::value &&
                      Impl::is_same< typename ViewType::memory_space , HostSpace >::value &&
                      Impl::is_same< PType , Policy >::value &&
                      ! Impl::is_same< typename PType::work_tag , void >::value
                    ), const FunctorType & >::type functor
                , const PType     & policy
                , const ViewType  & result
                )
    {
      pointer_type result_ptr = result.ptr_on_device();

      if ( ! result_ptr ) {
        result_ptr = (pointer_type)
          Kokkos::Serial::scratch_memory_resize( Reduce::value_size( functor ) , 0 );
      }

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      const typename PType::member_type e = policy.end();
      for ( typename PType::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( typename PType::work_tag() , i , update );
      }

      Reduce::final( functor , result_ptr );
    }
};

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 >
                  , Kokkos::Serial
                  >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 > Policy ;
  typedef ReduceAdapter< FunctorType >  Reduce ;

public:

  typedef typename Reduce::pointer_type pointer_type ;

  // work tag is void
  template< class PType >
  inline
  ParallelScan( typename Impl::enable_if<
                 ( Impl::is_same< PType , Policy >::value &&
                   Impl::is_same< typename PType::work_tag , void >::value
                 ), const FunctorType & >::type functor
             , const PType & policy )
    {
      pointer_type result_ptr = (pointer_type)
        Kokkos::Serial::scratch_memory_resize( Reduce::value_size( functor ) , 0 );

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      const typename PType::member_type e = policy.end();
      for ( typename PType::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( i , update , true );
      }

      Reduce::final( functor , result_ptr );
    }

  // work tag is non-void
  template< class PType >
  inline
  ParallelScan( typename Impl::enable_if<
                 ( Impl::is_same< PType , Policy >::value &&
                   ! Impl::is_same< typename PType::work_tag , void >::value
                 ), const FunctorType & >::type functor
             , const PType & policy )
    {
      pointer_type result_ptr = (pointer_type)
        Kokkos::Serial::scratch_memory_resize( Reduce::value_size( functor ) , 0 );

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      const typename PType::member_type e = policy.end();
      for ( typename PType::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( typename PType::work_tag() , i , update , true );
      }

      Reduce::final( functor , result_ptr );
    }
};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , Kokkos::TeamPolicy< Kokkos::Serial , void > , Kokkos::Serial >
{
private:
  typedef Kokkos::TeamPolicy< Kokkos::Serial , void > Policy ;
public:

  ParallelFor( const FunctorType & functor
             , const Policy      & policy )
    {
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );

      Kokkos::Serial::scratch_memory_resize( 0 , shared_size );

      for ( int ileague = 0 ; ileague < policy.league_size() ; ++ileague ) {
        functor( typename Policy::member_type(ileague,policy.league_size(),shared_size) );
      }
    }
};

template< unsigned int VectorLength, class FunctorType >
class ParallelFor< FunctorType , Kokkos::TeamVectorPolicy< VectorLength, Kokkos::Serial , void > , Kokkos::Serial >
{
private:
  typedef Kokkos::TeamVectorPolicy< VectorLength,Kokkos::Serial , void > Policy ;
public:

  ParallelFor( const FunctorType & functor
             , const Policy      & policy )
    {
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );

      Kokkos::Serial::scratch_memory_resize( 0 , shared_size );

      for ( int ileague = 0 ; ileague < policy.league_size() ; ++ileague ) {
        functor( typename Policy::member_type(ileague,policy.league_size(),shared_size) );
      }
    }
};

template< class FunctorType >
class ParallelReduce< FunctorType , Kokkos::TeamPolicy< Kokkos::Serial , void > , Kokkos::Serial > {
private:
  typedef Kokkos::TeamPolicy< Kokkos::Serial , void > Policy ;
  typedef ReduceAdapter< FunctorType >  Reduce ;
public:
  typedef typename Reduce::pointer_type pointer_type ;

  template< class ViewType >
  ParallelReduce( const FunctorType  & functor
                , const Policy       & policy
                , const ViewType     & result
                )
    {
      const int reduce_size = Reduce::value_size( functor );
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );
      void * const scratch_reduce = Kokkos::Serial::scratch_memory_resize( reduce_size , shared_size );

      const pointer_type result_ptr =
        result.ptr_on_device() ? result.ptr_on_device() 
                               : (pointer_type) scratch_reduce ;

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      for ( int ileague = 0 ; ileague < policy.league_size() ; ++ileague ) {
        functor( typename Policy::member_type(ileague,policy.league_size(),shared_size) , update );
      }

      Reduce::final( functor , result_ptr );
    }
};

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_SERIAL_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

