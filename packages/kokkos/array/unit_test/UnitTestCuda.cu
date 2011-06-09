
#include <iostream>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>

#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_DeviceCuda_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestMDArrayIndexMap.hpp>
#include <UnitTestReduce.hpp>

namespace Test {

void test_device_cuda()
{
#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )
  try {
    Kokkos::DeviceCuda::initialize();

    UnitTestDeviceMemoryManagement< Kokkos::DeviceCuda >();
    UnitTestValueView<       Kokkos::DeviceCuda >();
    UnitTestMultiVectorView< Kokkos::DeviceCuda >();
    UnitTestMDArrayView<     Kokkos::DeviceCuda >();
    UnitTestMDArrayDeepCopy< Kokkos::DeviceCuda >();

    Test::UnitTestMDArrayIndexMap< Kokkos::DeviceCuda >();

    UnitTestReduce< long ,   Kokkos::DeviceCuda >( 1000000 );
    UnitTestReduce< double , Kokkos::DeviceCuda >( 1000000 );

    std::cout << "PASSED : UnitTestCuda" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : UnitTestCuda : " << x.what() << std::endl ;
  }
#else
  std::cout << "PASSED : SKIPPED UnitTestCuda - NO DEVICE CUDA" << std::endl ;
#endif
}

}

#include <Kokkos_DeviceClear_macros.hpp>

//----------------------------------------------------------------------------

