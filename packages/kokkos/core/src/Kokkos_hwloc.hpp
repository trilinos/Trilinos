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

#ifndef KOKKOS_HWLOC_HPP
#define KOKKOS_HWLOC_HPP

#include <utility>

namespace Kokkos {

/** \brief  Minimal subset of logical 'hwloc' functionality available
 *          from http://www.open-mpi.org/projects/hwloc/.
 *
 *  The calls are NOT thread safe in order to avoid mutexes,
 *  memory allocations, or other actions which could give the
 *  runtime system an opportunity to migrate the threads or
 *  touch allocated memory during the function calls.
 *
 *  All calls to these functions should be performed by a thread
 *  when it has guaranteed exclusive access; e.g., for OpenMP
 *  within a 'critical' region.
 */
namespace hwloc {

/** \brief  Query if hwloc is available */
bool available();

/** \brief  Query number of available NUMA regions.
 *          This will be less than the hardware capacity
 *          if the MPI process is pinned to a NUMA region.
 */
unsigned get_available_numa_count();

/** \brief  Query number of available cores per NUMA regions.
 *          This will be less than the hardware capacity
 *          if the MPI process is pinned to a set of cores.
 */
unsigned get_available_cores_per_numa();

/** \brief  Query number of available "hard" threads per core; i.e., hyperthreads */
unsigned get_available_threads_per_core();


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
std::pair<unsigned,unsigned> get_core_topology();

/** \brief  Number of concurrent threads per core.
 *
 *  This typically reflects the number of hyperthreads
 *  the core can support.
 */
unsigned get_core_capacity();

} /* namespace hwloc */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Internal functions for binding persistent spawned threads.

namespace Kokkos {
namespace hwloc {

/** \brief  Determine best use of cores for a given thread count */
std::pair<unsigned,unsigned> use_core_topology( const unsigned thread_count );

/** \brief  Query core-coordinate of the current thread
 *          with respect to the core_topology.
 *
 *  As long as the thread is running within the 
 *  process binding the following condition holds.
 *
 *  core_coordinate.first  < core_topology.first
 *  core_coordinate.second < core_topology.second
 */
std::pair<unsigned,unsigned> get_this_thread_coordinate();

/** \brief  Bind the current thread to a core. */
bool bind_this_thread( const std::pair<unsigned,unsigned> );

/** \brief  Bind the current thread to one of the cores in the list.
 *          Set that entry to (~0,~0) and return the index.
 *          If binding fails return ~0.
 */
unsigned bind_this_thread( const unsigned               coordinate_count ,
                           std::pair<unsigned,unsigned> coordinate[] );

/** \brief  Unbind the current thread back to the original process binding */
bool unbind_this_thread();

void thread_mapping( const std::pair<unsigned,unsigned> team_topo ,
                     const std::pair<unsigned,unsigned> core_use ,
                     const std::pair<unsigned,unsigned> core_topo ,
                           std::pair<unsigned,unsigned> thread_coord[] );

void thread_mapping( const std::pair<unsigned,unsigned> team_topo ,
                     const std::pair<unsigned,unsigned> core_use ,
                     const std::pair<unsigned,unsigned> core_topo ,
                     const std::pair<unsigned,unsigned> master_coord ,
                           std::pair<unsigned,unsigned> thread_coord[] );

} /* namespace hwloc */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void host_thread_mapping( const std::pair<unsigned,unsigned> team_topo ,
                          const std::pair<unsigned,unsigned> core_use ,
                          const std::pair<unsigned,unsigned> core_topo ,
                                std::pair<unsigned,unsigned> thread_coord[] );

void host_thread_mapping( const std::pair<unsigned,unsigned> team_topo ,
                          const std::pair<unsigned,unsigned> core_use ,
                          const std::pair<unsigned,unsigned> core_topo ,
                          const std::pair<unsigned,unsigned> master_coord ,
                                std::pair<unsigned,unsigned> thread_coord[] );

}
}

#endif /* #define KOKKOS_HWLOC_HPP */

