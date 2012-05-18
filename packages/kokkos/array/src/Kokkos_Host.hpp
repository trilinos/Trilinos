/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_HOST_HPP
#define KOKKOS_HOST_HPP

#include <cstddef>
#include <impl/Kokkos_IndexMap.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/*--------------------------------------------------------------------------*/
/** \brief  Kokkos device for multicores in the host memory space
 *          with the prefered multidimensional array mapping.
 */
class Host {
public:
  /*--------------------------------*/
  /* required type declarations and functions for a device */

  typedef Host    memory_space ;
  typedef size_t  size_type ;

  template< unsigned Rank = 0 ,
            unsigned N1 = 0 , unsigned N2 = 0 , unsigned N3 = 0 ,
            unsigned N4 = 0 , unsigned N5 = 0 , unsigned N6 = 0 ,
            unsigned N7 = 0 >
  struct IndexMap {
    typedef Impl::IndexMapRight<memory_space,Rank,N1,N2,N3,N4,N5,N6,N7> type ;
  };

  /*--------------------------------*/
  /** \brief  Set the device in a 'sleep' state.
   *
   *  This function sets the device in a "sleep" state
   *  in which it is not ready for work and should
   *  consume less resources.
   *
   *  Return 'true' if the device is in the 'sleep' state.
   *  Return 'false' if the device is actively working and
   *  could not enter the 'sleep' state.
   */
  static bool sleep();

  /** \brief  Wake the device from the 'sleep' state so
   *          it is ready for work.
   *
   *  Return 'true' if the device is in the 'ready' state.
   *  Return 'false' if the device is actively working.
   */
  static bool wake();

  /** \brief  The parallel_for or parallel_reduce dispatch of a
   *          functor may return before the functor completes.
   *          Wait until all dispatched functors complete.
   */
  static void fence() {}

  /*--------------------------------*/
  /** \brief  Finalize: terminate spawned threads */
  static void finalize();

  /*--------------------------------*/
  /* -- Device specific functions --*/

  /** \brief  Initialize the device in the 'ready to work' state.
   *
   *  The device is initialized in a "ready to work" state
   *  which reduces latency / improves performance when
   *  dispatching work; however, resources are consumed
   *  even when no work is being done.
   *
   *  Initialize with number of NUMA nodes
   *  and number of threads on each NUMA node.
   */
  static void initialize( const size_type use_node_count ,
                          const size_type use_node_thread_count );

  /** \brief  Detect number of admissible NUMA nodes */
  static size_type detect_node_count();

  /** \brief  Detect number of cores per NUMA node */
  static size_type detect_node_core_count();

  /** \brief  An alignment size for large arrays */
  static size_type detect_memory_page_size();
};

/*--------------------------------------------------------------------------*/
/** \brief  Host memory space with another device's multi-index mapping. */
template< class Device >
struct HostMapped {
public:
  typedef HostMapped< Device >        type ;
  typedef typename Device::size_type  size_type ;
  typedef Host::memory_space          memory_space ;

  template< unsigned Rank = 0,
            unsigned N1 = 0, unsigned N2 = 0, unsigned N3 = 0,
            unsigned N4 = 0, unsigned N5 = 0, unsigned N6 = 0,
            unsigned N7 = 0 >
  struct IndexMap {
    typedef typename
      Device::template IndexMap<Rank,N1,N2,N3,N4,N5,N6,N7>::type type ;
  };
};

/** \brief  The host mapped onto the host is the host */
template<> struct HostMapped<Host> { typedef Host type ; };

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#include <Host/Kokkos_Host_MemoryManager.hpp>
#include <Host/Kokkos_Host_Parallel.hpp>
#include <Host/Kokkos_Host_ParallelFor.hpp>
#include <Host/Kokkos_Host_ParallelReduce.hpp>

#endif /* #define KOKKOS_HOST_HPP */

//----------------------------------------------------------------------------
/* Partial specializations for optional data structures */

#if   defined( KOKKOS_VALUE_HPP ) && \
    ! defined( KOKKOS_HOST_VALUE_HPP )
#include <Host/Kokkos_Host_Value.hpp>
#endif

#if   defined( KOKKOS_MULTIVECTOR_HPP ) && \
    ! defined( KOKKOS_HOST_MULTIVECTOR_HPP )
#include <Host/Kokkos_Host_MultiVector.hpp>
#endif

#if   defined( KOKKOS_ARRAY_HPP ) && \
    ! defined( KOKKOS_HOST_ARRAY_HPP )
#include <Host/Kokkos_Host_Array.hpp>
#endif

#if   defined( KOKKOS_MDARRAY_HPP ) && \
    ! defined( KOKKOS_HOST_MDARRAY_HPP )
#include <Host/Kokkos_Host_MDArray.hpp>
#endif

#if   defined( KOKKOS_PREFIXSUM_HPP ) && \
    ! defined( KOKKOS_HOST_PREFIXSUM_HPP )
#include <Host/Kokkos_Host_PrefixSum.hpp>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

