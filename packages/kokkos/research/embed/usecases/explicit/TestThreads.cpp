
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <Kokkos_Core.hpp>
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
  Kokkos::Threads::initialize( numa_node_count * numa_node_thread_count , numa_node_count );

  Explicit::test< Kokkos::Threads >( "Threads" , elem_count , iter_count );

  Kokkos::Threads::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

