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
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

template< class ExecSpace , typename IntType , unsigned GranularityPowerOfTwo = 3 >
class ExecPolicyRange {
private:

  enum { Granularity     = IntType(1) << GranularityPowerOfTwo };
  enum { GranularityMask = IntType(Granularity) - 1 };

public:

  typedef IntType index_type ;

  index_type begin ;
  index_type end ;

  KOKKOS_INLINE_FUNCTION
  ExecPolicyRange() : begin(0), end(0) {}

  KOKKOS_INLINE_FUNCTION
  ExecPolicyRange( const int part_rank
                 , const int part_size
                 , const index_type work_begin
                 , const index_type work_end
                 )
    : begin(0), end(0)
    {
      if ( part_size && work_begin < work_end ) {

        // Split evenly among partitions, then round up to the granularity.
        index_type work_part = ( ( work_end - work_begin ) + ( part_size - 1 ) ) / part_size ;

        if ( GranularityMask ) { work_part = ( work_part + GranularityMask ) & ~GranularityMask ; }

        begin = work_begin + work_part * part_rank ;
        end   = begin      + work_part ;

        if ( work_end < begin ) begin = work_end ;
        if ( work_end < end )   end   = work_end ;
      }
    }
};

template< class ExecSpace >
class ExecPolicyTeam {
public:
  /**\brief Tag identifying this class' concept */
  typedef Impl::ExecutionPolicyTag   kokkos_tag ;
  typedef ExecSpace                  execution_space ;

  struct index_type {

    KOKKOS_INLINE_FUNCTION int league_rank() const ;
    KOKKOS_INLINE_FUNCTION int league_size() const ;
    KOKKOS_INLINE_FUNCTION int team_rank() const ;
    KOKKOS_INLINE_FUNCTION int team_size() const ;

    KOKKOS_INLINE_FUNCTION void barrier();

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
     *
     *  The highest rank thread can compute the reduction total as
     *    reduction_total = dev.team_scan( value ) + value ;
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type scan( const Type & value );

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
    KOKKOS_INLINE_FUNCTION Type scan( const Type & value , Type * const global_accum );

    KOKKOS_INLINE_FUNCTION void * get_shmem( const int size );
  };
};

} // namespace Kokkos

#endif /* #define KOKKOS_EXECPOLICY_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

