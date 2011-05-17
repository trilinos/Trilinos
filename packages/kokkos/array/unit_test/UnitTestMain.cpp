
#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_MDArrayDeepCopy.hpp>
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

//----------------------------------------------------------------------------

#include <impl/Kokkos_DeviceTPI_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestReduce.hpp>

#include <impl/Kokkos_DeviceClear_macros.hpp>

//----------------------------------------------------------------------------

namespace {

void test_device_host()
{
  UnitTestDeviceMemoryManagement< Kokkos::DeviceHost >();
  UnitTestValueView<       Kokkos::DeviceHost >();
  UnitTestMultiVectorView< Kokkos::DeviceHost >();
  UnitTestMDArrayView<     Kokkos::DeviceHost >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceHost >();

  UnitTestReduce< long ,   Kokkos::DeviceHost >( 1000000 );
  UnitTestReduce< double , Kokkos::DeviceHost >( 1000000 );
}

void test_device_tpi()
{
  Kokkos::DeviceTPI::initialize( 4 );

  UnitTestDeviceMemoryManagement< Kokkos::DeviceTPI >();
  UnitTestValueView<       Kokkos::DeviceTPI >();
  UnitTestMultiVectorView< Kokkos::DeviceTPI >();
  UnitTestMDArrayView<     Kokkos::DeviceTPI >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceTPI >();

  UnitTestReduce< long ,   Kokkos::DeviceTPI >( 1000000 );
  UnitTestReduce< double , Kokkos::DeviceTPI >( 1000000 );
}

}

int main()
{
  test_device_host();
  test_device_tpi();

  return 0 ;
}

