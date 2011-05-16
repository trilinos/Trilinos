
#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_MDArrayDeepCopy.hpp>

#define KOKKOS_MACRO_DEVICE DeviceHost
#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#undef KOKKOS_MACRO_DEVICE

#define KOKKOS_MACRO_DEVICE DeviceTPI
#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
// #include <UnitTestMDArrayDeepCopy.hpp>
#undef KOKKOS_MACRO_DEVICE

int main()
{
  UnitTestDeviceMemoryManagement< Kokkos::DeviceHost >();
  UnitTestDeviceMemoryManagement< Kokkos::DeviceTPI >();
  UnitTestValueView< Kokkos::DeviceHost >();
  UnitTestValueView< Kokkos::DeviceTPI >();
  UnitTestMultiVectorView< Kokkos::DeviceHost >();
  UnitTestMultiVectorView< Kokkos::DeviceTPI >();
  UnitTestMDArrayView< Kokkos::DeviceHost >();
  UnitTestMDArrayView< Kokkos::DeviceTPI >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceHost >();
  // UnitTestMDArrayDeepCopy< Kokkos::DeviceTPI >();

  return 0 ;
}

