//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

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

void test_cuda_explicit( size_t elem_count ,
                         size_t iter_count )
{
#if __CUDACC__

  Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(1) );

  Explicit::test< Kokkos::Cuda >( "Cuda" , elem_count , iter_count );

  Kokkos::Cuda::finalize();
#endif
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

