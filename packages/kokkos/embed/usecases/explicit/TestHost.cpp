
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <Kokkos_Host.hpp>
#include <Kokkos_Array.hpp>
#include <impl/Kokkos_ArrayAnalyzeShape.hpp>
#include <impl/Kokkos_ArrayViewDefault.hpp>
#include <impl/Kokkos_Timer.hpp>

#include <TestHexGrad.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_explicit( size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t elem_count ,
                         size_t iter_count )
{
  Kokkos::Host::initialize( numa_node_count , numa_node_thread_count );

  Explicit::test< Kokkos::Host >( "Host" , elem_count , iter_count );

  Kokkos::Host::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

