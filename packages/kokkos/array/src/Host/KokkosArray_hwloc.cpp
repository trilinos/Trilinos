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

#include <sstream>
#include <iostream>
#include <limits>
#include <utility>

/* KokkosArray interfaces */

#include <Host/KokkosArray_hwloc.hpp>

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
namespace {

void print_bitmap( std::ostream & s , const hwloc_bitmap_t bitmap )
{
  s << "{" ;
  for ( int i = hwloc_bitmap_first( bitmap ) ;
        -1 != i ; i = hwloc_bitmap_next( bitmap , i ) ) {
    s << " " << i ;
  }
  s << " }" ;
}

struct HWLOC_Singleton {
  enum { MAX_ROOT_NODE = 1024 };

  hwloc_topology_t m_topology ;
  hwloc_obj_type_t m_root_type ;

  unsigned m_root_rank[ MAX_ROOT_NODE ];
  unsigned m_capacity[ hwloc::max_depth ];
  unsigned m_capacity_depth ;

  static HWLOC_Singleton & singleton();
  HWLOC_Singleton();
  ~HWLOC_Singleton();
};

HWLOC_Singleton &
HWLOC_Singleton::singleton()
{
  static HWLOC_Singleton self ;
  return self ;
}

HWLOC_Singleton::HWLOC_Singleton()
{
  m_root_type = HWLOC_OBJ_TYPE_MAX ;
  m_capacity_depth = 0 ;

  for ( unsigned i = 0 ; i < hwloc::max_depth ; ++i ) {
    m_capacity[i] = 0 ;
  }

  hwloc_topology_init( & m_topology );
  hwloc_topology_load( m_topology );

  { // Choose a hwloc type for a 'node' from, in search order, the following:
    static const hwloc_obj_type_t candidate_node_type[] =
      { HWLOC_OBJ_NODE   /* NUMA region     */
      , HWLOC_OBJ_SOCKET /* hardware socket */
      , HWLOC_OBJ_CORE   /* hardware core   */
      };

    enum { CANDIDATE_NODE_TYPE_COUNT =
             sizeof(candidate_node_type) / sizeof(hwloc_obj_type_t) };

    for ( int k = 0 ; k < CANDIDATE_NODE_TYPE_COUNT && HWLOC_OBJ_TYPE_MAX == m_root_type ; ++k ) {
      if ( 1 < hwloc_get_nbobjs_by_type( m_topology , candidate_node_type[k] ) ) {
        m_root_type = candidate_node_type[k] ;
      }
    }
  }

  {
    // Determine which of these 'node' types are available to this process.
    // The process may have been bound (e.g., by MPI) to a subset of these node types.

    hwloc_bitmap_t proc_cpuset = hwloc_bitmap_alloc();

    hwloc_get_cpubind( m_topology , proc_cpuset , HWLOC_CPUBIND_PROCESS );

    int node_count    = 0 ;
    int core_per_node = 0 ;
    int pu_per_core   = 0 ;
    bool symmetry     = true ;

    const int max_count = hwloc_get_nbobjs_by_type( m_topology , m_root_type );

    for ( int i = 0 ; i < max_count ; ++i ) {

      // Candidate node:
      const hwloc_obj_t node = hwloc_get_obj_by_type( m_topology , m_root_type , i );

      // Cores within this node:
      const int core_count =
        hwloc_get_nbobjs_inside_cpuset_by_type( m_topology ,
                                                node->allowed_cpuset ,
                                                HWLOC_OBJ_CORE );

      int core_available = 0 ;
      int pu_available   = 0 ;

      // All cores within this 'node' must be available
      // to add this node to the list.

      for ( int j = 0 ; j < core_count ; ++j ) {

        const hwloc_obj_t core =
          hwloc_get_obj_inside_cpuset_by_type( m_topology ,
                                               node->allowed_cpuset ,
                                               HWLOC_OBJ_CORE , j );

        // If process' cpuset intersects core's cpuset
        // then process can access this core.
        // Must use intersection instead of inclusion because the Intel-Phi
        // MPI may bind the process to only one of the core's hyperthreads.
        //
        // Assumption: if the process can access any hyperthread of the core
        // then it has ownership of the entire core.
        // This assumes that it would be performance-detrimental
        // to spawn more than one MPI process per core and use nested threading.

        if ( hwloc_bitmap_intersects( proc_cpuset , core->allowed_cpuset ) ) {

          ++core_available ;

          const int pu_count =
            hwloc_get_nbobjs_inside_cpuset_by_type( m_topology ,
                                                    core->allowed_cpuset ,
                                                    HWLOC_OBJ_PU );

          if ( 0 == pu_available ) pu_available = pu_count ;

          if ( pu_count != pu_available ) { symmetry = false ; }
        }
      }

      if ( core_count && core_count == core_available ) {

        if ( 0 == core_per_node ) core_per_node = core_count ;
        if ( 0 == pu_per_core )   pu_per_core   = pu_available ;

        if ( core_count != core_per_node ) {
          symmetry      = false ;
          core_per_node = std::min( core_count , core_per_node );
        }
        if ( pu_available != pu_per_core ) {

          symmetry    = false ;
          pu_per_core = std::min( pu_available , pu_per_core );
        }

        m_root_rank[ node_count++ ] = i ;
      }
    }

    if ( node_count ) {
      m_capacity[0] = node_count ;
      m_capacity_depth = 1 ;

      if ( 1 < core_per_node ) {
        m_capacity[1] = core_per_node ;
        m_capacity_depth = 2 ;
      }

      if ( 1 < pu_per_core ) {
        m_capacity[ m_capacity_depth ] = pu_per_core ;
        ++m_capacity_depth ;
      }
    }

    hwloc_bitmap_free( proc_cpuset );

    if ( ! symmetry ) {
      std::cerr << "KokkosArray::Impl::hwloc WARNING Asymmetric hardware core hierarchy" << std::endl ;
    }
  }
}

HWLOC_Singleton::~HWLOC_Singleton()
{
  hwloc_topology_destroy( m_topology );
}

} // namespace
} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

void hwloc::print_thread_capacity( std::ostream & s )
{
  const HWLOC_Singleton & h = HWLOC_Singleton::singleton();

  static const char name_numa[] = "NUMA" ;
  static const char name_socket[] = "SOCKET" ;
  static const char name_core[] = "CORE" ;
  static const char name_pu[] = "PU" ;

  const int max_count = hwloc_get_nbobjs_by_type( h.m_topology , h.m_root_type );

  const char * name_level_0 = 0 ;
  const char * name_level_1 = 0 ;
  const char * name_level_2 = 0 ;

  switch( h.m_root_type ) {
  case HWLOC_OBJ_NODE :
    name_level_0 = name_numa ;
    if ( 1 < h.m_capacity_depth ) name_level_1 = name_core ;
    if ( 2 < h.m_capacity_depth ) name_level_2 = name_pu ;
    break ;

  case HWLOC_OBJ_SOCKET :
    name_level_0 = name_socket ;
    if ( 1 < h.m_capacity_depth ) name_level_1 = name_core ;
    if ( 2 < h.m_capacity_depth ) name_level_2 = name_pu ;
    break ;

  case HWLOC_OBJ_CORE :
    name_level_0 = name_core ;
    if ( 1 < h.m_capacity_depth ) name_level_1 = name_pu ;
    break ;

  case HWLOC_OBJ_PU :
    name_level_0 = name_pu ;
    break ;

  default:
    break ;
  }

  s << "KokkosArray::Impl::hwloc{" ;
  if ( name_level_0 ) {

    s << " " << name_level_0 << "[" ;

    if ( h.m_capacity[0] < (unsigned) max_count ) {
      s << "(" ;
      for ( unsigned i = 0 ; i < h.m_capacity[0] ; ++i ) {
        s << " " << h.m_root_rank[i] ;
      }
      s << " ) / " ;
    }
    s << max_count << "]" ;
  }
  if ( name_level_1 ) s << " x " << name_level_1 << "[" << h.m_capacity[1] << "]" ;
  if ( name_level_2 ) s << " x " << name_level_2 << "[" << h.m_capacity[2] << "]" ;
  s << " }" ;
}

unsigned hwloc::get_thread_capacity_depth()
{
  HWLOC_Singleton & h = HWLOC_Singleton::singleton();

  return h.m_capacity_depth ;
}

void hwloc::get_thread_capacity( unsigned capacity[] )
{
  HWLOC_Singleton & h = HWLOC_Singleton::singleton();

  for ( unsigned i = 0 ; i < hwloc::max_depth ; ++i ) {
    capacity[i] = h.m_capacity[i] ;
  }
}

//----------------------------------------------------------------------------

void hwloc::map_thread( const hwloc::BindingPolicy policy ,
                        const unsigned rank ,
                        const unsigned count ,
                        unsigned coordinate[] )
{
  HWLOC_Singleton & h = HWLOC_Singleton::singleton();

  hwloc::map_thread( policy , h.m_capacity_depth , h.m_capacity ,
                     rank , count , coordinate );
}

//----------------------------------------------------------------------------

void hwloc::map_thread( const hwloc::BindingPolicy policy ,
                        const unsigned capacity_depth ,
                        const unsigned capacity[] ,
                        const unsigned rank ,
                        const unsigned count ,
                        unsigned coordinate[] )
{
  unsigned capacity_count = 1 ;

  for ( unsigned i = 0 ; i < hwloc::max_depth ; ++i ) {
    capacity_count *= capacity[i] ;
    coordinate[i] = 0 ;
  }

  if ( SPREAD == policy || capacity_count <= count ) {
    // Spread threads across higher levels of the topology
    // and use higher ranking hardware resource.

    unsigned n = count ;
    unsigned r = rank ;
    for ( unsigned i = 0 ; i < capacity_depth ; ++i ) {
      // n = k * bin + ( capacity[i] - k ) * ( bin + 1 )
      const unsigned bin  = n / capacity[i] ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = capacity[i] * bin1 - n ;
      const unsigned part = k * bin ;

      if ( r < part ) {
        coordinate[i]  = r / bin ;
        r = r - bin * coordinate[i] ;
        n = bin ;
      }
      else {
        const unsigned r1 = r - part ;
        const unsigned c1 = r1 / bin1 ;

        coordinate[i]  = c1 + k ;
        r = r1 - c1 * bin1 ;
        n = bin1 ;
      }
    }
  }
  else if ( PACK == policy ) {
    unsigned r0 = capacity_count - count ;
    unsigned r = rank + r0 ;

    for ( unsigned i = 0 ; i < capacity_depth ; ++i ) {
      capacity_count /= capacity[i] ;

      const unsigned c0 = r0 / capacity_count ;
      coordinate[i]  = r / capacity_count ;

      r  = r  % capacity_count ;
      r0 = ( c0 < coordinate[i] ) ? 0 : r0 % capacity_count ;
    }
  }
}

//----------------------------------------------------------------------------

bool hwloc::bind_this_thread( const unsigned coordinate[] )
{
  HWLOC_Singleton & h = HWLOC_Singleton::singleton();

  bool result = true ;

  for ( unsigned i = 0 ; i < h.m_capacity_depth ; ++i ) {
    if ( h.m_capacity[i] <= coordinate[i] ) {
      result = false ;
    }
  }

  if ( result ) {

    const int node_rank = h.m_root_rank[ coordinate[0] ] ;
    int core_rank = 0 ;
    int pu_rank = 0 ;

    switch( h.m_root_type ) {
    case HWLOC_OBJ_NODE :
    case HWLOC_OBJ_SOCKET :
      if ( 1 < h.m_capacity_depth ) core_rank = coordinate[1] ;
      if ( 2 < h.m_capacity_depth ) pu_rank   = coordinate[2] ;
      break ;

    case HWLOC_OBJ_CORE :
      if ( 1 < h.m_capacity_depth ) pu_rank   = coordinate[1] ;
      break ;

    default:
      break ;
    }

    const hwloc_obj_t node =
      hwloc_get_obj_by_type( h.m_topology ,
                             h.m_root_type ,
                             node_rank );

    const hwloc_obj_t core =
      hwloc_get_obj_inside_cpuset_by_type( h.m_topology ,
                                           node->allowed_cpuset ,
                                           HWLOC_OBJ_CORE ,
                                           core_rank );

    const hwloc_obj_t pu =
      hwloc_get_obj_inside_cpuset_by_type( h.m_topology ,
                                           core->allowed_cpuset ,
                                           HWLOC_OBJ_PU ,
                                           pu_rank );

    result = 0 == hwloc_set_cpubind( h.m_topology ,
                                     pu->allowed_cpuset ,
                                     HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );

    if ( result ) {

      hwloc_cpuset_t thread_cpuset = hwloc_bitmap_alloc();

      hwloc_get_cpubind( h.m_topology, thread_cpuset, HWLOC_CPUBIND_THREAD );

      result = hwloc_bitmap_isequal( thread_cpuset , pu->allowed_cpuset );

      hwloc_bitmap_free( thread_cpuset );
    }

#if 0
  std::ostringstream msg ;
  msg << "KokkosArray::Impl::hwloc::bind_this_thread("
      << coordinate[0] << ","
      << coordinate[1] << ","
      << coordinate[2] << ","
      << coordinate[3] << ") " ;
  msg << "node" ;  print_bitmap(msg,node->allowed_cpuset);
  msg << " core" ; print_bitmap(msg,core->allowed_cpuset);
  msg << " pu" ;   print_bitmap(msg,pu->allowed_cpuset);
  msg << std::endl ;
  std::cout << msg.str();
#endif

  }

  return result ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

