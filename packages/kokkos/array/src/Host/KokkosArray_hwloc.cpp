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

#include <iostream>
#include <sstream>
#include <stdexcept>

/* KokkosArray interfaces */

#include <KokkosArray_hwloc.hpp>

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
namespace {

enum { MAX_CORE = 1024 };

std::pair<unsigned,unsigned> s_core_topology(0,0);
unsigned                     s_core_capacity(0);
hwloc_topology_t             s_hwloc_topology(0);
hwloc_bitmap_t               s_hwloc_location(0);
hwloc_bitmap_t               s_process_binding(0);
hwloc_bitmap_t               s_core[ MAX_CORE ];

}

inline
void print_bitmap( std::ostream & s , const hwloc_const_bitmap_t bitmap )
{
  s << "{" ;
  for ( int i = hwloc_bitmap_first( bitmap ) ;
        -1 != i ; i = hwloc_bitmap_next( bitmap , i ) ) {
    s << " " << i ;
  }
  s << " }" ;
}

//----------------------------------------------------------------------------

bool hwloc::bind_this_thread( const std::pair<unsigned,unsigned> coord )
{
  // As safe and fast as possible.
  // Fast-lookup by caching the coordinate -> hwloc cpuset mapping in 's_core'.
  return coord.first  < s_core_topology.first &&
         coord.second < s_core_topology.second &&
         0 == hwloc_set_cpubind( s_hwloc_topology ,
                                 s_core[ coord.second + coord.first * s_core_topology.second ] ,
                                 HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );
}

bool hwloc::unbind_this_thread()
{
  return s_hwloc_topology &&
         0 == hwloc_set_cpubind( s_hwloc_topology ,
                                 s_process_binding ,
                                 HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );
}

//----------------------------------------------------------------------------

std::pair<unsigned,unsigned> hwloc::get_this_thread_coordinate()
{
  const unsigned n = s_core_topology.first * s_core_topology.second ;

  std::pair<unsigned,unsigned> coord(0,0);

  // Using the pre-allocated 's_hwloc_location' to avoid memory
  // allocation by this thread.  This call is NOT thread-safe.
  hwloc_get_last_cpu_location( s_hwloc_topology ,
                               s_hwloc_location , HWLOC_CPUBIND_THREAD );

  unsigned i = 0 ;

  while ( i < n && ! hwloc_bitmap_intersects( s_hwloc_location , s_core[ i ] ) ) ++i ;

  if ( i < n ) {
    coord.first  = i / s_core_topology.second ;
    coord.second = i % s_core_topology.second ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::hwloc::get_this_thread_coordinate() FAILED :" ;

    if ( 0 != s_process_binding && 0 != s_hwloc_location ) {
      msg << " cpu_location" ;
      print_bitmap( msg , s_hwloc_location );
      msg << " is not a member of the process_cpu_set" ;
      print_bitmap( msg , s_process_binding );
    }
    else {
      msg << " not initialized" ;
    }
    throw std::runtime_error( msg.str() );
  }
  return coord ;
}

//----------------------------------------------------------------------------

std::pair<unsigned,unsigned>
hwloc::get_core_topology()
{ sentinel(); return s_core_topology ; }

unsigned
hwloc::get_core_capacity()
{ sentinel(); return s_core_capacity ; }

void hwloc::sentinel()
{ static hwloc self ; }

//----------------------------------------------------------------------------

hwloc::~hwloc()
{
  hwloc_topology_destroy( s_hwloc_topology );
  hwloc_bitmap_free( s_process_binding );
  hwloc_bitmap_free( s_hwloc_location );
}

hwloc::hwloc()
{
  s_core_topology   = std::pair<unsigned,unsigned>(0,0);
  s_core_capacity   = 0 ;
  s_hwloc_topology  = 0 ;
  s_hwloc_location  = 0 ;
  s_process_binding = 0 ;

  for ( unsigned i = 0 ; i < MAX_CORE ; ++i ) s_core[i] = 0 ;

  hwloc_topology_init( & s_hwloc_topology );
  hwloc_topology_load( s_hwloc_topology );

  s_hwloc_location  = hwloc_bitmap_alloc();
  s_process_binding = hwloc_bitmap_alloc();

  hwloc_get_cpubind( s_hwloc_topology , s_process_binding ,  HWLOC_CPUBIND_PROCESS );

  // Choose a hwloc object type for the NUMA level, which may not exist.

  hwloc_obj_type_t root_type = HWLOC_OBJ_TYPE_MAX ;

  {
    // Object types to search, in order.
    static const hwloc_obj_type_t candidate_root_type[] =
      { HWLOC_OBJ_NODE     /* NUMA region     */
      , HWLOC_OBJ_SOCKET   /* hardware socket */
      , HWLOC_OBJ_MACHINE  /* local machine   */
      };

    enum { CANDIDATE_ROOT_TYPE_COUNT =
             sizeof(candidate_root_type) / sizeof(hwloc_obj_type_t) };

    for ( int k = 0 ; k < CANDIDATE_ROOT_TYPE_COUNT && HWLOC_OBJ_TYPE_MAX == root_type ; ++k ) {
      if ( 0 < hwloc_get_nbobjs_by_type( s_hwloc_topology , candidate_root_type[k] ) ) {
        root_type = candidate_root_type[k] ;
      }
    }
  }

  // Determine which of these 'root' types are available to this process.
  // The process may have been bound (e.g., by MPI) to a subset of these root types.
  // Determine current location of the master (calling) process>

  hwloc_bitmap_t proc_cpuset_location = hwloc_bitmap_alloc();

  hwloc_get_last_cpu_location( s_hwloc_topology , proc_cpuset_location , HWLOC_CPUBIND_THREAD );

  const unsigned max_root = hwloc_get_nbobjs_by_type( s_hwloc_topology , root_type );

  unsigned root_base     = max_root ;
  unsigned root_count    = 0 ;
  unsigned core_per_root = 0 ;
  unsigned pu_per_core   = 0 ;
  bool     symmetric     = true ;

  for ( unsigned i = 0 ; i < max_root ; ++i ) {

    const hwloc_obj_t root = hwloc_get_obj_by_type( s_hwloc_topology , root_type , i );

    if ( hwloc_bitmap_intersects( s_process_binding , root->allowed_cpuset ) ) {

      ++root_count ;

      // Remember which root (NUMA) object the master thread is running on.
      // This will be logical NUMA rank #0 for this process.

      if ( hwloc_bitmap_intersects( proc_cpuset_location, root->allowed_cpuset ) ) {
        root_base = i ;
      }

      // Count available cores:

      const unsigned max_core =
        hwloc_get_nbobjs_inside_cpuset_by_type( s_hwloc_topology ,
                                                root->allowed_cpuset ,
                                                HWLOC_OBJ_CORE );

      unsigned core_count = 0 ;

      for ( unsigned j = 0 ; j < max_core ; ++j ) {

        const hwloc_obj_t core =
          hwloc_get_obj_inside_cpuset_by_type( s_hwloc_topology ,
                                               root->allowed_cpuset ,
                                               HWLOC_OBJ_CORE , j );

        // If process' cpuset intersects core's cpuset then process can access this core.
        // Must use intersection instead of inclusion because the Intel-Phi
        // MPI may bind the process to only one of the core's hyperthreads.
        //
        // Assumption: if the process can access any hyperthread of the core
        // then it has ownership of the entire core.
        // This assumes that it would be performance-detrimental
        // to spawn more than one MPI process per core and use nested threading.

        if ( hwloc_bitmap_intersects( s_process_binding , core->allowed_cpuset ) ) {

          ++core_count ;

          const unsigned pu_count =
            hwloc_get_nbobjs_inside_cpuset_by_type( s_hwloc_topology ,
                                                    core->allowed_cpuset ,
                                                    HWLOC_OBJ_PU );

          if ( pu_per_core == 0 ) pu_per_core = pu_count ;

          // Enforce symmetry by taking the minimum:

          pu_per_core = std::min( pu_per_core , pu_count );

          if ( pu_count != pu_per_core ) symmetric = false ;
        }
      }

      if ( 0 == core_per_root ) core_per_root = core_count ;

      // Enforce symmetry by taking the minimum:

      core_per_root = std::min( core_per_root , core_count );

      if ( core_count != core_per_root ) symmetric = false ;
    }
  }

  s_core_topology.first  = root_count ;
  s_core_topology.second = core_per_root ;
  s_core_capacity        = pu_per_core ;

  // Fill the 's_core' array for fast mapping from a core coordinate to the
  // hwloc cpuset object required for thread location querying and binding.

  for ( unsigned i = 0 ; i < max_root ; ++i ) {

    const unsigned root_rank = ( i + root_base ) % max_root ;

    const hwloc_obj_t root = hwloc_get_obj_by_type( s_hwloc_topology , root_type , root_rank );

    if ( hwloc_bitmap_intersects( s_process_binding , root->allowed_cpuset ) ) {

      const unsigned max_core =
        hwloc_get_nbobjs_inside_cpuset_by_type( s_hwloc_topology ,
                                                root->allowed_cpuset ,
                                                HWLOC_OBJ_CORE );

      unsigned core_count = 0 ;

      for ( unsigned j = 0 ; j < max_core && core_count < core_per_root ; ++j ) {

        const hwloc_obj_t core =
          hwloc_get_obj_inside_cpuset_by_type( s_hwloc_topology ,
                                               root->allowed_cpuset ,
                                               HWLOC_OBJ_CORE , j );

        if ( hwloc_bitmap_intersects( s_process_binding , core->allowed_cpuset ) ) {

          s_core[ core_count + core_per_root * i ] = core->allowed_cpuset ;

          ++core_count ;
        }
      }
    }
  }

  hwloc_bitmap_free( proc_cpuset_location );

  if ( ! symmetric ) {
    std::cout << "KokkosArray::hwloc WARNING: Using a symmetric subset of a non-symmetric core topology."
              << std::endl ;
  }
}

//----------------------------------------------------------------------------

} /* namespace KokkosArray */


