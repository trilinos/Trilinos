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

#ifndef KOKKOSARRAY_HOST_THREADDATA_HPP
#define KOKKOSARRAY_HOST_THREADDATA_HPP

#include <KokkosArray_HostSpace.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>

namespace KokkosArray {
namespace Impl {

class HostThreadSentinel ;

//----------------------------------------------------------------------------
/** \brief  A thread within the pool. */

class HostThread {
public:

  typedef HostSpace::size_type size_type ;

  inline size_type rank() const { return m_thread_rank ; }

  inline size_type gang_rank() const { return m_gang_rank ; }

  inline size_type worker_rank() const { return m_worker_rank ; }

  //----------------------------------------------------------------------
  /** \brief  Compute a range of work for this thread's rank */

  std::pair< size_type , size_type >
  inline
  work_range( const size_type work_count ) const
    {
      enum { work_align = HostSpace::WORK_ALIGNMENT };
      enum { work_shift = power_of_two< work_align >::value };
      enum { work_mask  = work_align - 1 };

      // unit_of_work_count = ( work_count + work_mask ) >> work_shift
      // unit_of_work_per_thread = ( unit_of_work_count + thread_count - 1 ) / thread_count
      // work_per_thread = unit_of_work_per_thread * work_align

      const size_type work_per_thread =
        ( ( ( ( work_count + work_mask ) >> work_shift ) + m_thread_count - 1 ) / m_thread_count ) << work_shift ;

      const size_type work_begin = std::min( m_thread_rank * work_per_thread , work_count );
      const size_type work_end   = std::min( work_begin + work_per_thread , work_count );

      return std::pair<size_type,size_type>( work_begin , work_end );
    }

  //----------------------------------------------------------------------

  inline
  void * reduce_data() const
    {
#if defined( __INTEL_COMPILER )
__assume_aligned(m_reduce,HostSpace::MEMORY_ALIGNMENT);
#endif
      return m_reduce ;
    }

  //----------------------------------------------------------------------

  inline
  size_type fan_count() { return m_fan_count ; }

  inline
  HostThread & fan( unsigned i ) const { return *m_fan[i]; }

  //----------------------------------------------------------------------

  inline static
  HostThread * get_thread( const unsigned global_rank )
    { return m_thread[ global_rank ]; }

  static
  void set_thread( const unsigned global_rank , HostThread * );

  static
  void clear_thread( const unsigned global_rank );

  ~HostThread();
  HostThread();

  //----------------------------------------------------------------------

  static const unsigned max_fan_count = 16 ;
  static const unsigned max_thread_count = 1 << max_fan_count ;

private:

  HostThread( const HostThread & );
  HostThread & operator = ( const HostThread & );

  static HostThread * m_thread[ max_thread_count ];

  HostThread  * m_fan[ max_fan_count ] ;
  size_type     m_fan_count ;
  size_type     m_thread_rank ;
  size_type     m_thread_count ;
  size_type     m_gang_rank ;
  size_type     m_gang_count ;
  size_type     m_worker_rank ;
  size_type     m_worker_count ;
  void        * m_reduce ;    ///< Reduction memory

public:

  int  volatile m_state ;     ///< Thread control flag

  /** \brief States of a worker thread */
  enum State { ThreadTerminating ///<  Exists, termination in progress
             , ThreadInactive    ///<  Exists, waiting for work
             , ThreadActive      ///<  Exists, performing work
             , ThreadRendezvous  ///<  Exists, waiting in a barrier or reduce
             };

  friend class HostThreadSentinel ;
  friend class HostInternal ;
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HOST_THREADDATA_HPP */

