
#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>

#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <impl/Kokkos_DeviceHost_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestReduce.hpp>

#include <impl/Kokkos_DeviceClear_macros.hpp>

namespace {

void test_device_host()
{
  std::cout << "test_device_host()" << std::endl ;

  UnitTestDeviceMemoryManagement< Kokkos::DeviceHost >();
  UnitTestValueView<       Kokkos::DeviceHost >();
  UnitTestMultiVectorView< Kokkos::DeviceHost >();
  UnitTestMDArrayView<     Kokkos::DeviceHost >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceHost >();

  UnitTestReduce< long ,   Kokkos::DeviceHost >( 1000000 );
  UnitTestReduce< double , Kokkos::DeviceHost >( 1000000 );
}

}

//----------------------------------------------------------------------------

#include <impl/Kokkos_DeviceCuda_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
// #include <UnitTestReduce.hpp>

namespace {

void test_device_cuda()
{
  std::cout << "test_device_cuda()" << std::endl ;

  Kokkos::DeviceCuda::initialize();

  UnitTestDeviceMemoryManagement< Kokkos::DeviceCuda >();
  UnitTestValueView<       Kokkos::DeviceCuda >();
  UnitTestMultiVectorView< Kokkos::DeviceCuda >();
  UnitTestMDArrayView<     Kokkos::DeviceCuda >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceCuda >();

  // UnitTestReduce< long ,   Kokkos::DeviceCuda >( 1000000 );
  // UnitTestReduce< double , Kokkos::DeviceCuda >( 1000000 );
}

}

#include <impl/Kokkos_DeviceClear_macros.hpp>

//----------------------------------------------------------------------------

int main()
{
  test_device_host();
  test_device_cuda();

  return 0 ;
}

