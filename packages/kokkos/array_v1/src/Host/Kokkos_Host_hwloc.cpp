/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2011 Sandia Corporation
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

/*--------------------------------------------------------------------------*/
/* Kokkos interfaces */

#include <Kokkos_Host.hpp>
#include <Host/Kokkos_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Third Party Libraries */

/* Hardware locality library: http://www.open-mpi.org/projects/hwloc/ */
#include <hwloc.h>

#define  REQURED_HWLOC_API_VERSION  0x000010300

#if HWLOC_API_VERSION < REQUIRED_HWLOC_API_VERSION
#error "Requires  http://www.open-mpi.org/projects/hwloc/  Version 1.3 or greater"
#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

namespace {

class Host_hwloc {
private:

  hwloc_topology_t  m_host_topology ;
  unsigned          m_node_count ;
  unsigned          m_core_per_node ;
  unsigned          m_page_size ;

  ~Host_hwloc();
  Host_hwloc();
  Host_hwloc( const Host_hwloc & );
  Host_hwloc & operator = ( const Host_hwloc & );

public:

  static const Host_hwloc & singleton();

  inline
  unsigned page_size() const
  { return m_page_size ; }

  inline
  unsigned node_count() const
  { return m_node_count ; }

  inline
  unsigned core_per_node() const
  { return m_core_per_node ; }

  inline
  bool cpubind( unsigned node_rank ) const
  {
    const hwloc_obj_t node =
      hwloc_get_obj_by_type( m_host_topology, HWLOC_OBJ_NODE, node_rank );

    return -1 != hwloc_set_cpubind( m_host_topology , node->allowed_cpuset ,
                                    HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );
  }
};

Host_hwloc::Host_hwloc()
{
  m_node_count    = 0 ;
  m_core_per_node = 0 ;
  m_page_size     = 0 ;

  hwloc_topology_init( & m_host_topology );
  hwloc_topology_load( m_host_topology );

  const int hwloc_depth =
    hwloc_get_type_depth( m_host_topology , HWLOC_OBJ_NODE );

  if ( HWLOC_TYPE_DEPTH_UNKNOWN != hwloc_depth ) {

    m_node_count =
      hwloc_get_nbobjs_by_depth( m_host_topology , hwloc_depth );

    bool ok = 0 < m_node_count ;

    for ( unsigned i = 0 ; i < m_node_count && ok ; ++i ) {

      const hwloc_obj_t node =
        hwloc_get_obj_by_type( m_host_topology , HWLOC_OBJ_NODE , i );

      const unsigned count = hwloc_bitmap_weight( node->allowed_cpuset );

      if ( 0 == m_core_per_node ) { m_core_per_node = count ; }

      ok = count == m_core_per_node ;

      bool homogeneous_page_size = true ;

      unsigned page_size = 0 ;

      for ( unsigned j = 0 ; homogeneous_page_size &&
                             j < node->memory.page_types_len ; ++j ) {
        if ( node->memory.page_types[j].count ) {
          if ( 0 == page_size ) {
            page_size = node->memory.page_types[j].size ;
          }
          homogeneous_page_size = node->memory.page_types[j].size == page_size ;
        }
      }

      if ( homogeneous_page_size ) { m_page_size = page_size ; }
    }

    if ( ! ok ) {
      m_core_per_node = 0 ;
      m_node_count    = 0 ;
      m_page_size     = 0 ;
    }
  }
}

Host_hwloc::~Host_hwloc()
{
  hwloc_topology_destroy( m_host_topology );
}

const Host_hwloc & Host_hwloc::singleton()
{
  static Host_hwloc self ;
  return self ;
}

} // namespace <>
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

unsigned host_internal_page_size()
{ return Host_hwloc::singleton().page_size(); }

unsigned host_internal_node_count()
{ return Host_hwloc::singleton().node_count(); }

unsigned host_internal_core_per_node()
{ return Host_hwloc::singleton().core_per_node(); }

bool host_internal_bind_this_thread_to_node( unsigned node_rank )
{ return Host_hwloc::singleton().cpubind( node_rank ); }

} // namespace Impl
} // namespace Kokkos

