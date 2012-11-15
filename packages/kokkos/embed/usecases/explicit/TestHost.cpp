
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <KokkosArray_Host.hpp>
#include <KokkosArray_Array.hpp>
#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewOperLeft.hpp>
#include <impl/KokkosArray_ArrayViewOperRight.hpp>
#include <impl/KokkosArray_Timer.hpp>

#include <TestHexGrad.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_explicit( size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t elem_count ,
                         size_t iter_count )
{
  KokkosArray::Host::initialize( numa_node_count , numa_node_thread_count );

  Explicit::test< KokkosArray::Host >( "Host" , elem_count , iter_count );

  KokkosArray::Host::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

