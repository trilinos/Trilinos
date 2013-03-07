//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Array.hpp>
#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewLeft.hpp>
#include <impl/KokkosArray_ArrayViewRight.hpp>
#include <impl/KokkosArray_Timer.hpp>

#include <TestHexGrad.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_cuda_explicit( size_t elem_count ,
                         size_t iter_count )
{
#if __CUDACC__

  KokkosArray::Cuda::initialize( KokkosArray::Cuda::SelectDevice(1) );

  Explicit::test< KokkosArray::Cuda >( "Cuda" , elem_count , iter_count );

  KokkosArray::Cuda::finalize();
#endif
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

