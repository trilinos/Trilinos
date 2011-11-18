/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_DEVICENUMA_HPP
#define KOKKOS_DEVICENUMA_HPP

#include <Kokkos_DeviceHost.hpp>

#define KOKKOS_DEVICE_NUMA  Kokkos::DeviceNUMA

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class DeviceNUMA {
public:

  /** \brief  On the NUMA device use size_t for indexing */
  typedef size_t                size_type ;

  /** \brief  The NUMA device uses the Host memory space */
  typedef HostMemory            memory_space ;

  /** \brief  Default mdarray map is index from right */
  typedef Impl::MDArrayIndexMapRight  mdarray_map ;

  /*--------------------------------*/

  /** \brief  The parallel_for or parallel_reduce dispactch of a
   *          functor may return before the functor completes.
   *          Wait until all dispatched functors complete.
   */
  static void wait_functor_completion() {}

  /*--------------------------------*/
  /** \brief  If core count detection functionality is available
   *          return the core count, otherwise return zero.
   */
  static size_type detect_core_count(); 

  enum UtilizationStrategy
    { DETECT_AND_USE_ALL_CORES   ///< Use every core of every node
    , DETECT_AND_USE_MOST_CORES  ///< Use all but one core of every node
    , MANUALLY_SET_THREAD_COUNT  ///< Manually specify number of threads
    };

  /** \brief  Create threads on the numa-host manycore device
   *          with the given strategy.
   */
  static void initialize( UtilizationStrategy ,
                          size_type manually_set_thread_count = 0 );

  /** \brief  Destroy created threads */
  static void finalize();

  /** \brief  Threads are "spinning" for
   *          optimal work-spawning performance.
   *          Block the threads so they do not consume
   *          resources on the cores.
   *
   *  Throw an exception if the device is active or already blocked.
   */
  static void block();

  /** \brief  Unblock threads and return them to the
   *          "spinning" and ready-to-work state.
   *
   *  Throw an exception if the device is not blocked.
   */
  static void unblock();
};

} // namespace Kokkos

#include <DeviceNUMA/Kokkos_DeviceNUMA_Impl.hpp>

#endif /* #define KOKKOS_DEVICENUMA_HPP */

