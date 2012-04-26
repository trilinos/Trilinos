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

/*--------------------------------------------------------------------------*/

#include <iostream>

/* Kokkos interfaces */

#include <Kokkos_Host.hpp>
#include <Host/Kokkos_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Third Party Libraries */

/* Hardware locality library: http://www.open-mpi.org/projects/hwloc/ */
#include <hwloc.h>

#define  REQUIRED_HWLOC_API_VERSION  0x000010300

#if HWLOC_API_VERSION < REQUIRED_HWLOC_API_VERSION
#error "Requires  http://www.open-mpi.org/projects/hwloc/  Version 1.3 or greater"
#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

class HostInternalHWLOC : public HostInternal {
private:

  hwloc_topology_t  m_host_topology ;

public:

  ~HostInternalHWLOC();
  HostInternalHWLOC();

  bool bind_to_node( const HostThread & thread ) const ;
};

//----------------------------------------------------------------------------

HostInternal & HostInternal::singleton()
{
  static HostInternalHWLOC self ; return self ;
}

//----------------------------------------------------------------------------

bool HostInternalHWLOC::bind_to_node( const HostThread & thread ) const
{
  bool result = true ;

  if ( HostInternal::m_node_count ) {

    const hwloc_obj_t node =
      hwloc_get_obj_by_type( m_host_topology, HWLOC_OBJ_NODE, thread.node() );

    result = 0 == hwloc_set_cpubind( m_host_topology ,
                                     node->allowed_cpuset ,
                                     HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );

#if 0
    std::cout << ( result ? "SUCCESS " : "FAILED " )
              << "HWLOC::bind_to_node thread["
              << thread.rank()
              << "] to node["
              << thread.node()
              << "] with PU_count["
              << hwloc_bitmap_weight( node->allowed_cpuset )
              << "]"
              << std::endl ;
#endif

  }
  return result ;
};

//----------------------------------------------------------------------------

HostInternalHWLOC::HostInternalHWLOC()
  : HostInternal()
{
  hwloc_topology_init( & m_host_topology );
  hwloc_topology_load( m_host_topology );

  const size_t node_count =
    hwloc_get_nbobjs_by_type( m_host_topology , HWLOC_OBJ_NODE );

  if ( node_count ) {
    // Get cpuset binding of this process.
    // This may have been bound by 'mpirun'.

    bool node_symmetry = true ;
    bool page_symmetry = true ;
    size_t page_size = 0 ;
    size_t node_core_count = 0 ;

    hwloc_cpuset_t proc_cpuset = hwloc_bitmap_alloc();

    hwloc_get_cpubind( m_host_topology , proc_cpuset , 0 );

    for ( size_t i = 0 ; i < node_count ; ++i ) {

      const hwloc_obj_t node =
        hwloc_get_obj_by_type( m_host_topology , HWLOC_OBJ_NODE , i );

      // If the process' cpu set is included in the node cpuset
      // then assumed pinned to that node.
      if ( hwloc_bitmap_isincluded( proc_cpuset , node->allowed_cpuset ) ) {
        m_node_rank = i ;
      }

      const unsigned count =
        hwloc_get_nbobjs_inside_cpuset_by_type( m_host_topology ,
                                                node->allowed_cpuset ,
                                                HWLOC_OBJ_PU );
                                                // HWLOC_OBJ_CORE );

      if ( 0 == node_core_count ) {
        node_core_count = count ;
      }

      node_symmetry = count == node_core_count ;

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
      HostInternal::m_node_core_count = node_core_count ;
    }

    if ( page_symmetry && page_size ) {
      HostInternal::m_page_size = page_size ;
    }

    HostInternal::m_node_count = node_count ;
  }
}

HostInternalHWLOC::~HostInternalHWLOC()
{
  hwloc_topology_destroy( m_host_topology );
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} /* namespace Impl */
} /* namespace Kokkos */

