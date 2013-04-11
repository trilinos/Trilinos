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

#ifndef KOKKOSARRAY_HOST_INTERNAL_HPP
#define KOKKOSARRAY_HOST_INTERNAL_HPP

#include <KokkosArray_Host.hpp>
#include <KokkosArray_HostSpace.hpp>
#include <Host/KokkosArray_Host_Parallel.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

/**
 * \class HostInternal
 * \brief Internal implementation of intraprocess parallelism on the host.
 *
 * Hardware model
 * ==============
 *
 * - The Host process is running within a NUMA multiprocessor environment.
 * - The hardware locality (hwloc) library defines a 'node' as a collection
 *   of processing units associated with a NUMA region.
 * - If the Host process is pinned to a particular NUMA node we assume
 *   that the threads of the Host process are also restricted to that node.
 *
 * Be aware that "node" here means "CPU processing units associated
 * with a NUMA affinity region," and does not have its traditional
 * high-performance computing hardware meaning.
 */
class HostInternal {
protected:

  unsigned         m_thread_count ;  // Number of threads
  HostThread       m_master_thread ;

  const HostThreadWorker * volatile m_worker ;

  ~HostInternal();

  HostInternal();

private:
  /// \brief Spawn the worker threads, and set up interthread communication.
  ///
  /// \param use_node_count [in] The number of NUMA regions ("nodes")
  ///   to use.
  /// \param use_node_thread_count [in] The number of worker threads
  ///   to use per NUMA region.
  ///
  /// \return Whether the threads successfully spawned.
  ///
  /// The calling thread also gets bound as a worker thread.  This has
  /// implications for parallel kernels: in particular, they are not
  /// asynchronous.
  bool spawn_threads( const unsigned gang_count ,
                      const unsigned worker_count);

  bool spawn();

  unsigned bind_host_thread();

  void activate_threads();

public:
  /// \brief Assert at run time that the calling worker thread is inactive.
  ///
  /// \param method [in] Name of the method invoking the assertion.
  ///   Used only for constructing the exception message in case the
  ///   assertion fails.
  void verify_inactive( const char * const method ) const ;

  /// \brief Initialize the worker threads for parallel kernels.
  ///
  /// \param use_node_count [in] The number of NUMA regions ("nodes")
  ///   to use.
  /// \param use_node_thread_count [in] The number of worker threads
  ///   to use per NUMA region.
  ///
  /// \exception std::runtime_error if initialization failed.
  ///
  /// The calling thread, which is the master thread, also gets bound
  /// as a worker thread.  This has implications for parallel kernels:
  /// in particular, they are not asynchronous.  Tasks get assigned to
  /// the master thread as well as to the other worker threads.
  void initialize( const unsigned use_node_count ,
                   const unsigned use_node_thread_count );

  void finalize();

  void print_configuration( std::ostream & ) const ;

  void execute_serial( const HostThreadWorker & worker );

  inline void execute( const HostThreadWorker & worker );

  static void driver();

  bool is_master_thread() const ;

  void * reduce_scratch() const ;

  //! Access the one HostInternal instance.
  static HostInternal & singleton();

  friend class KokkosArray::Host ;
  friend class KokkosArray::HostSpace ;
};

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_HOST_INTERNAL_HPP */

