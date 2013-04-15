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

#ifndef KOKKOSARRAY_HWLOC_HPP
#define KOKKOSARRAY_HWLOC_HPP

#include <utility>

namespace KokkosArray {

/** \brief  Minimal subset of logical 'hwloc' functionality available
 *          from http://www.open-mpi.org/projects/hwloc/.
 *
 *  The calls are NOT thread safe in order to avoid mutexes,
 *  memory allocations, or other actions which could give the
 *  runtime system an opportunity to migrate the threads or
 *  touch allocated memory during the function calls.
 *
 *  All calls to these functions should be perform by a thread
 *  when it has guaranteed exclusive access; e.g., for OpenMP
 *  within a 'critical' region.
 */
struct hwloc {

  /** \brief  Query the core topology of ( NUMA x Core/NUMA ).
   *
   *  The topology is limited by the process binding,
   *  which may have been set by MPI.  NUMA rank #0
   *  contains the core on which the process / master thread
   *  is running.  The master thread should only be bound
   *  to its original NUMA rank - because moving it to
   *  a different NUMA rank will displace it from all of
   *  the memory which it has already touched.
   */
  static std::pair<unsigned,unsigned> get_core_topology();

  /** \brief  Number of concurrent threads per core.
   *
   *  This typically reflects the number of hyperthreads
   *  the core can support.
   */
  static unsigned get_core_capacity();

  /** \brief  Query core-coordinate of the current thread
   *          with respect to the core_topology.
   *
   *  As long as the thread is running within the 
   *  process binding the following condition holds.
   *
   *  core_coordinate.first  < core_topology.first
   *  core_coordinate.second < core_topology.second
   */
  static std::pair<unsigned,unsigned> get_this_thread_coordinate();

  /** \brief  Bind the current thread to a core. */
  static bool bind_this_thread( const std::pair<unsigned,unsigned> );

  /** \brief  Unbind the current thread back to the original process binding */
  static bool unbind_this_thread();

private:

  hwloc();
  ~hwloc();

  static void sentinel();
};

} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HWLOC_HPP */

