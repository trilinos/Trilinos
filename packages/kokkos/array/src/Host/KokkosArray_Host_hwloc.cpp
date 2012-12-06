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

/*--------------------------------------------------------------------------*/

#define DEBUG_PRINT 0

#include <iostream>
#include <limits>

/* KokkosArray interfaces */

#include <KokkosArray_Host.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Third Party Libraries */

/* Hardware locality library: http://www.open-mpi.org/projects/hwloc/ */
#include <hwloc.h>

#define  REQUIRED_HWLOC_API_VERSION  0x000010300

#if HWLOC_API_VERSION < REQUIRED_HWLOC_API_VERSION
#error "Requires  http://www.open-mpi.org/projects/hwloc/  Version 1.3 or greater"
#endif

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

class HostInternalHWLOC : public HostInternal {
private:

  hwloc_topology_t m_host_topology ;
  hwloc_obj_type_t m_node_type ;          // hwloc type for a "node"
  unsigned         m_node_core_count ;    // Cores per node
  unsigned         m_node_core_pu_count ; // Processing units per core per node

public:

  ~HostInternalHWLOC();
  HostInternalHWLOC();

  bool bind_thread( const unsigned thread_rank ) const ;
};

//----------------------------------------------------------------------------

HostInternal & HostInternal::singleton()
{
  static HostInternalHWLOC self ; return self ;
}

//----------------------------------------------------------------------------

bool HostInternalHWLOC::bind_thread( const unsigned thread_rank ) const
{
  bool result = true ;

  // Can only safely bind threads if
  // (1) the number of cores is known and
  // (2) the process is bound to a NUMA node.

  if ( 0 < m_node_core_count && 0 <= m_node_rank) {

    // How many cores will be used:
    const unsigned max_worker_per_core =
      ( HostInternal::m_worker_count + m_node_core_count - 1 ) / m_node_core_count ;

    const unsigned min_worker_per_core =
      1 == max_worker_per_core ? 1 : max_worker_per_core - 1 ;

    const unsigned core_base =
      m_node_core_count * max_worker_per_core - HostInternal::m_worker_count ;

    const unsigned core_base_worker_count = core_base * min_worker_per_core ;

    // Which node -> core -> processing unit

    const unsigned gang_rank   = thread_rank / HostInternal::m_worker_count ;
    const unsigned worker_rank = thread_rank % HostInternal::m_worker_count ;

    const unsigned node_rank =
      ( gang_rank + HostInternal::m_node_rank ) % HostInternal::m_node_count ;

    const unsigned core_rank = 
      worker_rank < core_base_worker_count ?
      worker_rank / min_worker_per_core :
      core_base + ( worker_rank - core_base_worker_count ) / max_worker_per_core ;

    const unsigned pu_rank = worker_rank % m_node_core_pu_count ;

    const hwloc_obj_t node =
      hwloc_get_obj_by_type( m_host_topology, m_node_type, node_rank );

    const hwloc_obj_t core =
      hwloc_get_obj_inside_cpuset_by_type( m_host_topology ,
                                           node->allowed_cpuset ,
                                           HWLOC_OBJ_CORE ,
                                           core_rank );

    const hwloc_obj_t pu =
      hwloc_get_obj_inside_cpuset_by_type( m_host_topology ,
                                           core->allowed_cpuset ,
                                           HWLOC_OBJ_PU ,
                                           pu_rank );

    result = 0 == hwloc_set_cpubind( m_host_topology ,
                                     pu->allowed_cpuset ,
                                     HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );

    if ( result ) {

      hwloc_cpuset_t thread_cpuset = hwloc_bitmap_alloc();

      hwloc_get_cpubind( m_host_topology, thread_cpuset, HWLOC_CPUBIND_THREAD );

      result = hwloc_bitmap_isequal( thread_cpuset , pu->allowed_cpuset );

      hwloc_bitmap_free( thread_cpuset );
    }

#if DEBUG_PRINT
    std::cout << ( result ? "SUCCESS " : "FAILED " )
              << "HWLOC::bind_thread thread[ "
              << thread_rank
              << " @ " << gang_rank
              << "." << worker_rank
              << " ] to node[" << node_rank
              << "].core[" << core_rank
              << "].pu[" << pu_rank
              << "]"
              << std::endl ;
#endif

  }

  return result ;
}

//----------------------------------------------------------------------------

HostInternalHWLOC::HostInternalHWLOC()
  : HostInternal()
  , m_host_topology()
  , m_node_type( HWLOC_OBJ_TYPE_MAX )
  , m_node_core_count( 0 )
  , m_node_core_pu_count( 0 )
{
  hwloc_topology_init( & m_host_topology );
  hwloc_topology_load( m_host_topology );

  size_t node_count = 0 ;

  //------------------------------------
  // Try for NUMA node knowledge
  m_node_type = HWLOC_OBJ_NODE ;
  node_count  = hwloc_get_nbobjs_by_type( m_host_topology , m_node_type );

  if ( 0 == node_count ) {
    // No knowledge of NUMA, try for SOCKET knowledge
    m_node_type = HWLOC_OBJ_SOCKET ;
    node_count  = hwloc_get_nbobjs_by_type( m_host_topology , m_node_type );
  }

  if ( 0 == node_count ) {
    // No knowledge of NUMA or SOCKET, try for MACHINE
    m_node_type = HWLOC_OBJ_MACHINE ;
    node_count  = hwloc_get_nbobjs_by_type( m_host_topology , m_node_type );
  }
  //------------------------------------

#if DEBUG_PRINT
  std::cout << "HWLOC node_count = " << node_count << std::endl ;
#endif

  if ( node_count ) {
    // Get cpuset binding of this process.
    // This may have been bound by 'mpirun'.

    bool node_symmetry = true ;
    bool page_symmetry = true ;

    size_t cache_line_size = 0 ;
    size_t page_size = 0 ;
    size_t node_core_count = 0 ;
    size_t node_core_pu_count = 0 ;
    size_t node_rank = std::numeric_limits<size_t>::max();

    hwloc_cpuset_t proc_cpuset = hwloc_bitmap_alloc();

    hwloc_get_cpubind( m_host_topology , proc_cpuset , 0 );

    // Is this process bound to a particular NUMA node?
    // This may have been done by 'mpirun'.
    // If so then restrict this device to that NUMA node.

    for ( size_t i = 0 ; i < node_count && node_count < node_rank ; ++i ) {

      const hwloc_obj_t node =
        hwloc_get_obj_by_type( m_host_topology , m_node_type , i );

      // If the process' cpu set is included in the node cpuset
      // then assumed pinned to that node.
      if ( hwloc_bitmap_isincluded( proc_cpuset , node->allowed_cpuset ) ) {
        node_rank = i ;

#if DEBUG_PRINT
        std::cout << "HWLOC: process is bound to node[" << i << "]"
                  << std::endl ;
#endif

      }
    }

    const bool bound_proc = node_rank < node_count ;

    // If master process is not bound, choose NUMA region #0 to bind.
    // This region will always exist.

    HostInternal::m_node_rank  = bound_proc ? node_rank : 0 ;
    HostInternal::m_node_count = node_count ;

    for ( unsigned i = 0 ; i < node_count ; ++i ) {

      const hwloc_obj_t node =
        hwloc_get_obj_by_type( m_host_topology , m_node_type , i );

      const size_t core_count =
        hwloc_get_nbobjs_inside_cpuset_by_type( m_host_topology ,
                                                node->allowed_cpuset ,
                                                HWLOC_OBJ_CORE );

      if ( 0 == node_core_count ) { node_core_count = core_count ; }

      if ( core_count != node_core_count ) { node_symmetry = false ; }

      for ( size_t j = 0 ; j < core_count ; ++j ) {

        const hwloc_obj_t core =
          hwloc_get_obj_inside_cpuset_by_type( m_host_topology ,
                                               node->allowed_cpuset ,
                                               HWLOC_OBJ_CORE , j );

#if DEBUG_PRINT
        if ( hwloc_bitmap_isincluded( proc_cpuset , core->allowed_cpuset ) ) {
          std::cout << "HWLOC: process is bound to node[" << i << "]"
                    << ".core[" << j << "]"
                    << std::endl ;
        }
#endif

        const size_t pu_count =
          hwloc_get_nbobjs_inside_cpuset_by_type( m_host_topology ,
                                                  core->allowed_cpuset ,
                                                  HWLOC_OBJ_PU );

        if ( 0 == node_core_pu_count ) { node_core_pu_count = pu_count ; }

        if ( pu_count != node_core_pu_count ) { node_symmetry = false ; }

        // Use the largest cache line size
        // assuming the largest will be a multiple of the smallest...

        const hwloc_obj_t core_cache_info =
          hwloc_get_shared_cache_covering_obj( m_host_topology , core );

        if ( core_cache_info && core_cache_info->attr ) {
          if ( cache_line_size < core_cache_info->attr->cache.linesize ) {
            cache_line_size = core_cache_info->attr->cache.linesize ;
          }
        }
      }

      for ( unsigned j = 0 ; j < node->memory.page_types_len ; ++j ) {
        if ( node->memory.page_types[j].count ) {
          if ( 0 == page_size ) {
            page_size = node->memory.page_types[j].size ;
          }
          page_symmetry = node->memory.page_types[j].size == page_size ;
        }
      }
    }

    hwloc_bitmap_free( proc_cpuset );

    if ( node_symmetry && node_core_count ) {

      m_node_core_count    = node_core_count ;
      m_node_core_pu_count = node_core_pu_count ;

      HostInternal::m_node_pu_count = node_core_count * node_core_pu_count ;
    }
    else {
      std::cerr << "KokkosArray::Host WARNING: multicore CPUs are not symmetric"
                << std::endl ;
    }

    if ( page_symmetry && page_size ) {
      HostInternal::m_page_size = page_size ;
    }

    if ( cache_line_size ) {
      HostInternal::m_cache_line_size = cache_line_size ;
      HostInternal::m_work_chunk      = cache_line_size / sizeof(void*);
    }
  }
}

HostInternalHWLOC::~HostInternalHWLOC()
{
  hwloc_topology_destroy( m_host_topology );
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} /* namespace Impl */
} /* namespace KokkosArray */

