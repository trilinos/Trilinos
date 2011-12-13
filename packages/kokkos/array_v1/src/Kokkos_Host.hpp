/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
#include <impl/Kokkos_MDArrayIndexMap.hpp>

#define KOKKOS_HOST  Kokkos::Host

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

  typedef size_t                                     size_type ;
  typedef Host                                       memory_space ;
  typedef Impl::MDArrayIndexMapRight< memory_space > mdarray_map ;

  /*--------------------------------*/

  struct DetectAndUseAllCores {};

  struct SetThreadCount {
    const size_t thread_count ;
    SetThreadCount( size_t n ) : thread_count( n ) {}
  };

  /** \brief  Initialize the device in the 'ready to work' state.
   *
   *  The device is initialized in a "ready to work" state
   *  which reduces latency / improves performance when
   *  dispatching work; however, resources are consumed
   *  even when no work is being done.
   */
  static void initialize( const DetectAndUseAllCores );

  static void initialize( const SetThreadCount );

  static void finalize();

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
  /* -- Device specific functions --*/

  /** \brief  If core count detection functionality is available
   *          return the core count, otherwise return zero.
   */
  static size_type detect_core_count();
};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#include <Host/Kokkos_Host_MemoryManager.hpp>
#include <Host/Kokkos_Host_Parallel.hpp>
#include <Host/Kokkos_Host_ParallelFor.hpp>
#include <Host/Kokkos_Host_ParallelReduce.hpp>

#endif /* #define KOKKOS_HOST_HPP */

/* Partial specializations for optional data structures */
#include <Host/Kokkos_Host_Specialize.hpp>

