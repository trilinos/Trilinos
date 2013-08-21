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

#ifndef KOKKOS_THREADS_HPP
#define KOKKOS_THREADS_HPP

#include <cstddef>
#include <iosfwd>
#include <Kokkos_Layout.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Host.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class ThreadsExec ;
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Device for a pool of Pthreads or C11 threads on a CPU. */
class Threads {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  typedef Threads                  device_type ;
  typedef Kokkos::HostSpace        memory_space ;
  typedef memory_space::size_type  size_type ;
  typedef Kokkos::LayoutRight      array_layout ;
  typedef Kokkos::Threads          host_mirror_device_type ;

  //@}
  /*------------------------------------------------------------------------*/
  //! \name Static functions that all Kokkos devices must implement.
  //@{

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
  static void fence();

  /// \brief Free any resources being consumed by the device.
  ///
  /// For the Host device, this terminates spawned worker threads.
  static void finalize();

  /** \brief  Print configuration information */
  static void print_configuration( std::ostream & , bool detail = false );

  //@}
  /*------------------------------------------------------------------------*/
  /** \name Function for the functor device interface */
  /**@{ */

  inline int league_rank() const ;
  inline int league_size() const ;
  inline int team_rank() const ;
  inline int team_size() const ;

  inline void team_barrier();

  template< typename T >
  inline T * get_shmem( const int count );

  inline std::pair<size_t,size_t> work_range( size_t ) const ;

  explicit inline Threads( Impl::ThreadsExec & );

  /**@} */
  /*------------------------------------------------------------------------*/
  //! \name Device-specific functions
  //@{

  /** \brief Initialize the device in the "ready to work" state.
   *
   *  The device is initialized in a "ready to work" or "awake" state.
   *  This state reduces latency and thus improves performance when
   *  dispatching work.  However, the "awake" state consumes resources
   *  even when no work is being done.  You may call sleep() to put
   *  the device in a "sleeping" state that does not consume as many
   *  resources, but it will take time (latency) to awaken the device
   *  again (via the wake()) method so that it is ready for work.
   *
   *  Initialize with ( league_size , team_size ).
   *  All worker threads of a team will occupy the same NUMA node.
   *
   *  The core topology must be less than or equal to hwloc::get_core_topology().
   *  If core topology is not input then the full hwloc::get_core_topology() is used.
   */
  static void initialize( const std::pair<unsigned,unsigned> league_team ,
                          const std::pair<unsigned,unsigned> hardware_topology =
                                std::pair<unsigned,unsigned>(0u,0u) );

  //@}
  /*------------------------------------------------------------------------*/

private:

  friend class Impl::ThreadsExec ;

  Impl::ThreadsExec & m_exec ;
};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#include <Kokkos_Parallel.hpp>
#include <Threads/Kokkos_ThreadsExec.hpp>
#include <Threads/Kokkos_Threads_Parallel.hpp>

#endif /* #define KOKKOS_THREADS_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

