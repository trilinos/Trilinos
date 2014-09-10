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

#ifndef KOKKOS_EXECPOLICY_HPP
#define KOKKOS_EXECPOLICY_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Execution policy for work over a range of an integral type.
 *
 * Valid template argument options:
 *
 *  With a specificed execution space:
 *    < ExecSpace , WorkTag , { IntConst | IntType } >
 *    < ExecSpace , WorkTag , void >
 *    < ExecSpace , { IntConst | IntType } , void >
 *    < ExecSpace , void , void >
 *
 *  With the default execution space:
 *    < WorkTag , { IntConst | IntType } , void >
 *    < WorkTag , void , void >
 *    < { IntConst | IntType } , void , void >
 *    < void , void , void >
 *
 *  IntType  is a fundamental integral type
 *  IntConst is an Impl::integral_constant< IntType , Blocking >
 *
 *  Blocking is the granularity of partitioning the range among threads.
 */
template< class Arg0 = void , class Arg1 = void , class Arg2 = void >
class RangePolicy {
private:

  // Default integral type and blocking factor:
  typedef int DefaultIntType ;
  enum { DefaultIntValue = 8 };

  enum { Arg0_Void = Impl::is_same< Arg0 , void >::value };
  enum { Arg1_Void = Impl::is_same< Arg1 , void >::value };
  enum { Arg2_Void = Impl::is_same< Arg2 , void >::value };

  enum { Arg0_ExecSpace = Impl::is_execution_space< Arg0 >::value };

  enum { Arg0_IntConst = Impl::is_integral_constant< Arg0 >::value };
  enum { Arg1_IntConst = Impl::is_integral_constant< Arg1 >::value };
  enum { Arg2_IntConst = Impl::is_integral_constant< Arg2 >::value };

  enum { Arg0_IntType = Impl::is_integral< Arg0 >::value };
  enum { Arg1_IntType = Impl::is_integral< Arg1 >::value };
  enum { Arg2_IntType = Impl::is_integral< Arg2 >::value };

  enum { Arg0_WorkTag = ! Arg0_ExecSpace && ! Arg0_IntConst && ! Arg0_IntType && ! Arg0_Void };
  enum { Arg1_WorkTag =   Arg0_ExecSpace && ! Arg1_IntConst && ! Arg1_IntType && ! Arg1_Void };


  enum { ArgOption_OK = Impl::StaticAssert< (
    ( Arg0_ExecSpace && Arg1_WorkTag && ( Arg2_IntConst || Arg2_IntType ) ) ||
    ( Arg0_ExecSpace && Arg1_WorkTag && Arg2_Void ) ||
    ( Arg0_ExecSpace && ( Arg1_IntConst || Arg2_IntType ) && Arg2_Void ) ||
    ( Arg0_ExecSpace && Arg1_Void && Arg2_Void ) ||
    ( Arg0_WorkTag && ( Arg1_IntConst || Arg2_IntType ) && Arg2_Void ) ||
    ( Arg0_WorkTag && Arg1_Void && Arg2_Void ) ||
    ( ( Arg0_IntConst || Arg0_IntType ) && Arg1_Void && Arg2_Void ) ||
    ( Arg0_Void && Arg1_Void && Arg2_Void )
    ) >::value };

  // The first argument is the execution space, otherwise the default execution space
  typedef typename Impl::if_c< Arg0_ExecSpace , Arg0 ,
          Kokkos::DefaultExecutionSpace
          >::type
    ExecSpace ;

  // The work argument tag is the first or second argument
  typedef typename Impl::if_c< Arg0_WorkTag , Arg0 ,
          typename Impl::if_c< Arg1_WorkTag , Arg1 , void
          >::type >::type
    WorkTag ;

  enum { Granularity = Arg0_IntConst ? unsigned(Impl::is_integral_constant<Arg0>::integral_value) : (
                       Arg1_IntConst ? unsigned(Impl::is_integral_constant<Arg1>::integral_value) : (
                       Arg2_IntConst ? unsigned(Impl::is_integral_constant<Arg2>::integral_value) : (
                                       unsigned(DefaultIntValue) ))) };

  // Only accept the integral type if the blocking is a power of two
  typedef typename Impl::enable_if< Impl::is_power_of_two< Granularity >::value ,
            typename Impl::if_c< Arg0_IntType , Arg0 ,
            typename Impl::if_c< Arg1_IntType , Arg1 ,
            typename Impl::if_c< Arg2_IntType , Arg2 ,
            typename Impl::if_c< Arg0_IntConst , typename Impl::is_integral_constant<Arg0>::integral_type ,
            typename Impl::if_c< Arg1_IntConst , typename Impl::is_integral_constant<Arg1>::integral_type ,
            typename Impl::if_c< Arg2_IntConst , typename Impl::is_integral_constant<Arg2>::integral_type ,
                                                 DefaultIntType
            >::type >::type >::type
            >::type >::type >::type
          >::type
    IntType ;

  enum { GranularityMask = IntType(Granularity) - 1 };

  IntType m_begin ;
  IntType m_end ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef ExecSpace                  execution_space ;
  typedef RangePolicy                execution_policy ;
  typedef IntType                    member_type ;
  typedef WorkTag                    work_tag ;

  KOKKOS_INLINE_FUNCTION member_type begin() const { return m_begin ; }
  KOKKOS_INLINE_FUNCTION member_type end()   const { return m_end ; }

  KOKKOS_INLINE_FUNCTION RangePolicy() : m_begin(0), m_end(0) {}

  /** \brief  Total range */
  KOKKOS_INLINE_FUNCTION
  RangePolicy( const member_type work_begin
             , const member_type work_end
             )
    : m_begin( work_begin < work_end ? work_begin : 0 )
    , m_end(   work_begin < work_end ? work_end : 0 )
    {}

  /** \brief  Subrange for a partition's rank and size.
   *
   *  Typically used to partition a range over a group of threads.
   */
  KOKKOS_INLINE_FUNCTION
  RangePolicy( const RangePolicy & range
             , const int part_rank
             , const int part_size
             )
    : m_begin(0), m_end(0)
    {
      if ( part_size ) {

        // Split evenly among partitions, then round up to the granularity.
        const member_type work_part =
          ( ( ( ( range.m_end - range.m_begin ) + ( part_size - 1 ) ) / part_size ) + GranularityMask ) & ~member_type(GranularityMask);

        m_begin = range.m_begin + work_part * part_rank ;
        m_end   = m_begin       + work_part ;

        if ( range.m_end < m_begin ) m_begin = range.m_end ;
        if ( range.m_end < m_end )   m_end   = range.m_end ;
      }
    }
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Execution policy for parallel work over a league of teams of threads.
 *
 *  The work functor is called for each thread of each team such that
 *  the team's member threads are guaranteed to be concurrent.
 *
 *  The team's threads have access to team shared scratch memory and
 *  team collective operations.
 *
 *  If the WorkTag is non-void then the first calling argument of the
 *  work functor's parentheses operator is 'const WorkTag &'.
 *  This allows a functor to have multiple work member functions.
 *
 *  template argument option with specified execution space:
 *    < ExecSpace , WorkTag >
 *    < ExecSpace , void >
 *
 *  template argument option with default execution space:
 *    < WorkTag , void >
 *    < void , void >
 */
template< class Arg0 = void , class Arg1 = void >
class TeamPolicy {
private:

  enum { Arg0_ExecSpace = Impl::is_execution_space< Arg0 >::value };
  enum { Arg0_Void      = Impl::is_same< Arg0 , void >::value };
  enum { Arg1_Void      = Impl::is_same< Arg1 , void >::value };

  enum { ArgOption_OK = Impl::StaticAssert< ( Arg0_ExecSpace || Arg1_Void ) >::value };

  typedef typename Impl::if_c< Arg0_ExecSpace , Arg0 , Kokkos::DefaultExecutionSpace >::type
    ExecSpace ;

  typedef typename Impl::if_c< Arg0_ExecSpace , Arg1 , Arg0 >::type
    WorkTag ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef ExecSpace                  execution_space ; ///< Execution space
  typedef WorkTag                    work_tag ;

  /** \brief  Query maximum team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   */
  template< class FunctorType >
  static int team_size_max( const FunctorType & );

  /** \brief  Construct policy with the given instance of the execution space */
  TeamPolicy( execution_space & , int league_size_request , int team_size_request );

  /** \brief  Construct policy with the default instance of the execution space */
  TeamPolicy( int league_size_request , int team_size_request );

  /** \brief  The actual league size (number of teams) of the policy.
   *
   *  This may be smaller than the requested league size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int league_size() const ;

  /** \brief  The actual team size (number of threads per team) of the policy.
   *
   *  This may be smaller than the requested team size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int team_size() const ;

  /** \brief  Parallel execution of a functor calls the functor once with
   *          each member of the execution policy.
   */
  struct member_type {

    /** \brief  Handle to the currently executing team shared scratch memory */
    KOKKOS_INLINE_FUNCTION
    typename execution_space::scratch_memory_space team_shmem() const ;

    /** \brief  Rank of this team within the league of teams */
    KOKKOS_INLINE_FUNCTION int league_rank() const ;

    /** \brief  Number of teams in the league */
    KOKKOS_INLINE_FUNCTION int league_size() const ;

    /** \brief  Rank of this thread within this team */
    KOKKOS_INLINE_FUNCTION int team_rank() const ;

    /** \brief  Number of threads in this team */
    KOKKOS_INLINE_FUNCTION int team_size() const ;

    /** \brief  Barrier among the threads of this team */
    KOKKOS_INLINE_FUNCTION void team_barrier() const ;

    /** \brief  Intra-team reduction. Returns join of all values of the team members. */
    template< class JoinOp >
    KOKKOS_INLINE_FUNCTION
    typename JoinOp::value_type team_reduce( const typename JoinOp::value_type
                                           , const JoinOp & ) const ;

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
     *
     *  The highest rank thread can compute the reduction total as
     *    reduction_total = dev.team_scan( value ) + value ;
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const ;

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
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const ;
  };
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Execution policy for parallel work over a league of teams of threads
 *          with explicit vector support.
 *
 *  The work functor is called for each thread of each team such that
 *  the team's member threads are guaranteed to be concurrent.
 *
 *  The team's threads have access to team shared scratch memory and
 *  team collective operations.
 *
 *  If the WorkArgTag is non-void then the first calling argument of the
 *  work functor's parentheses operator is 'const WorkArgTag &'.
 *  This allows a functor to have multiple work member functions.
 */
template< unsigned VectorLength
        , class ExecSpace  = DefaultExecutionSpace
        , class WorkArgTag = void >
class TeamVectorPolicy {
public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef ExecSpace                  execution_space ; ///< Execution space

  /** \brief  Query maximum team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   */
  template< class FunctorType >
  static int team_size_max( const FunctorType & );

  /** \brief  Construct policy with the given instance of the execution space */
  TeamVectorPolicy( execution_space & , int league_size_request , int team_size_request );

  /** \brief  Construct policy with the default instance of the execution space */
  TeamVectorPolicy( int league_size_request , int team_size_request );

  /** \brief  The actual league size (number of teams) of the policy.
   *
   *  This may be smaller than the requested league size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int league_size() const ;

  /** \brief  The actual team size (number of threads per team) of the policy.
   *
   *  This may be smaller than the requested team size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int team_size() const ;

  /** \brief  Parallel execution of a functor calls the functor once with
   *          each member of the execution policy.
   */
  struct member_type {

    /** \brief  Handle to the currently executing team shared scratch memory */
    KOKKOS_INLINE_FUNCTION
    typename execution_space::scratch_memory_space team_shmem() const ;

    /** \brief  Rank of this team within the league of teams */
    KOKKOS_INLINE_FUNCTION int league_rank() const ;

    /** \brief  Number of teams in the league */
    KOKKOS_INLINE_FUNCTION int league_size() const ;

    /** \brief  Rank of this thread within this team */
    KOKKOS_INLINE_FUNCTION int team_rank() const ;

    /** \brief  Number of threads in this team */
    KOKKOS_INLINE_FUNCTION int team_size() const ;

    /** \brief  Barrier among the threads of this team */
    KOKKOS_INLINE_FUNCTION void team_barrier() const ;

    /** \brief  Intra-team reduction. Returns join of all values of the team members. */
    template< class JoinOp >
    KOKKOS_INLINE_FUNCTION
    typename JoinOp::value_type team_reduce( const typename JoinOp::value_type
                                           , const JoinOp & ) const ;

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
     *
     *  The highest rank thread can compute the reduction total as
     *    reduction_total = dev.team_scan( value ) + value ;
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const ;

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
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const ;

#ifdef KOKKOS_HAVE_CXX11
    /** \brief  Intra-team parallel for. Executes op(const iType i) for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all threads in the team.
     * This functionality requires C++11 support.*/
    //template< typename iType, class Operation, typename ValueType >
    //KOKKOS_INLINE_FUNCTION void team_par_for(const iType n, const Operation & op) const ;

    /** \brief  Intra-team parallel sum reduction. Executes op(const iType i, ValueType & val) for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all threads in the team and an ordered summation over val of all threads
     * takes place. This functionality requires C++11 support.*/
    //template< typename iType, class Operation, typename ValueType >
    //KOKKOS_INLINE_FUNCTION void team_par_reduce(const iType n, const Operation & op, ValueType& result) const ;

    /** \brief  Intra-team parallel reduction. Executes op(const iType i, ValueType & val) for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all threads in the team and an ordered reduction over val of all threads
     * takes place where JoinType(ValueType& val, const ValueType& update). This functionality requires C++11 support.*/
    //template< typename iType, class Operation, typename ValueType, class JoinType >
    //KOKKOS_INLINE_FUNCTION void team_par_reduce(const iType n, const Operation & op, ValueType& result, const JoinType & join) const ;

    /** \brief  Intra-team parallel exclusive prefix sum. Executes op(const iType i, ValueType & val, bool final)
     *          for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all threads in the team and a scan operation is performed.
     * Depending on the target execution space the operator might be called twice: once with final=false
     * and once with final=true. When final==true val contains the prefix sum value. The contribution of this
     * "i" needs to be added to val no matter whether final==true or not. In a serial execution
     * (i.e. team_size==1) the operator is only called once with final==true. Scan_val will be set
     * to the final sum value over all threads.
     * This functionality requires C++11 support.*/
    //template< typename iType, class Operation, typename ValueType >
    //KOKKOS_INLINE_FUNCTION  void team_par_scan(const iType n, const Operation & op, ValueType& scan_val) const ;

    /** \brief  Guarantees execution of op() with only a single vector lane of this thread. */
    template< class Operation >
    KOKKOS_INLINE_FUNCTION void vector_single(const Operation & op) const ;

    /** \brief  Intra-thread vector parallel for. Executes op(const iType i) for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
     * This functionality requires C++11 support.*/
    template< typename iType, class Operation, typename ValueType >
    KOKKOS_INLINE_FUNCTION void vector_par_for(const iType n, const Operation & op) const ;

    /** \brief  Intra-thread vector parallel reduce. Executes op(const iType i, ValueType & val) for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
     * val is performed and put into result. This functionality requires C++11 support.*/
    template< typename iType, class Operation, typename ValueType >
    KOKKOS_INLINE_FUNCTION void vector_par_reduce(const iType n, const Operation & op, ValueType& result) const ;

    /** \brief  Intra-thread vector parallel reduce. Executes op(const iType i, ValueType & val) for each i=0..N-1.
     *
     * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
     * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
     * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
     * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
     * '1 for *'). This functionality requires C++11 support.*/
    template< typename iType, class Operation, typename ValueType, class JoinType >
    KOKKOS_INLINE_FUNCTION void vector_par_reduce(const iType n, const Operation & op, ValueType& init_result, const JoinType & join) const ;

    /** \brief  Intra-thread vector parallel exclusive prefix sum. Executes op(const iType i, ValueType & val, bool final)
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
    KOKKOS_INLINE_FUNCTION  void vector_par_scan(const iType n, const Operation & op, ValueType& scan_val) const ;
#endif
  };
};

} // namespace Kokkos

#endif /* #define KOKKOS_EXECPOLICY_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

