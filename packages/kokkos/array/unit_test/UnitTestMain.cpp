
#include <Kokkos_DeviceHost.hpp>

#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>

#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_DeviceHost_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestMDArrayIndexMap.hpp>
#include <UnitTestReduce.hpp>

#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

void test_device_host()
{
  try {
    UnitTestDeviceMemoryManagement< Kokkos::DeviceHost >();
    UnitTestValueView<       Kokkos::DeviceHost >();
    UnitTestMultiVectorView< Kokkos::DeviceHost >();
    UnitTestMDArrayView<     Kokkos::DeviceHost >();
    UnitTestMDArrayDeepCopy< Kokkos::DeviceHost >();

    Test::UnitTestMDArrayIndexMap< Kokkos::DeviceHost >();

    UnitTestReduce< long ,   Kokkos::DeviceHost >( 1000000 );
    UnitTestReduce< double , Kokkos::DeviceHost >( 1000000 );

    std::cout << "PASSED : UnitTestHost" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : UnitTestHost : " << x.what() << std::endl ;
  }
}

}

//----------------------------------------------------------------------------

namespace Test {
void test_device_tpi();
void test_device_cuda();
}

//----------------------------------------------------------------------------

int main()
{
  Test::test_device_host();
  Test::test_device_tpi();
  Test::test_device_cuda();

  return 0 ;
}

