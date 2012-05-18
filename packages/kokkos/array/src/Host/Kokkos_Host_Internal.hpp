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

#ifndef KOKKOS_HOST_INTERNAL_HPP
#define KOKKOS_HOST_INTERNAL_HPP

#include <Kokkos_Host.hpp>
#include <Host/Kokkos_Host_Parallel.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class HostWorkerBlock : public HostThreadWorker<void> {
public:
  void execute_on_thread( HostThread & ) const ;

  HostWorkerBlock()  {}
  ~HostWorkerBlock() {}
};

//----------------------------------------------------------------------------

/**
 *  Model:
 *  The Host process is running within a NUMA multiprocessor environment.
 *  The hardware locality (hwloc) library defines a 'node' as a collection
 *  of processing units associated with a NUMA region.
 *  If the Host process is pinned to a particular NUMA node we assume
 *  that the threads of the Host process are also restricted to that node.
 */
class HostInternal {
protected:

  enum { THREAD_COUNT_MAX = 1023 };

  typedef Host::size_type size_type ;

  HostWorkerBlock  m_worker_block ;
  size_type        m_node_rank ;       //
  size_type        m_node_count ;      //
  size_type        m_node_pu_count ;   // Assuming all nodes are equivalent
  size_type        m_page_size ;       // 
  size_type        m_thread_count ;    // 
  size_type        m_node_thread_count ; // 
  HostThread       m_thread[ THREAD_COUNT_MAX + 1 ];

  const HostThreadWorker<void> * volatile m_worker ;

  virtual ~HostInternal();
  HostInternal();

  bool spawn_threads( const size_type use_node_count ,
                      const size_type use_node_thread_count );

private:

  bool spawn( HostThread * thread );

public:

  void verify_inactive( const char * const method ) const ;

  void initialize( const size_type use_node_count ,
                   const size_type use_node_thread_count );

  void finalize();

  virtual bool bind_to_node( const HostThread & ) const ;

  inline void execute( const HostThreadWorker<void> & worker );

  inline void execute( HostThread & thread ) const
    { m_worker->execute_on_thread( thread ); }

  static HostInternal & singleton();

  friend class Kokkos::Host ;
};

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_HOST_INTERNAL_HPP */

