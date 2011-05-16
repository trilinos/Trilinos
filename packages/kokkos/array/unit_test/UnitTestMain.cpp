
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

int main()
{
  Kokkos::DeviceTPI::initialize( 4 );

  UnitTestDeviceMemoryManagement< Kokkos::DeviceHost >();
  UnitTestDeviceMemoryManagement< Kokkos::DeviceTPI >();
  UnitTestValueView< Kokkos::DeviceHost >();
  UnitTestValueView< Kokkos::DeviceTPI >();
  UnitTestMultiVectorView< Kokkos::DeviceHost >();
  UnitTestMultiVectorView< Kokkos::DeviceTPI >();
  UnitTestMDArrayView< Kokkos::DeviceHost >();
  UnitTestMDArrayView< Kokkos::DeviceTPI >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceHost >();
  UnitTestMDArrayDeepCopy< Kokkos::DeviceTPI >();

  UnitTestReduce< long , Kokkos::DeviceHost >( 1000000 );
  UnitTestReduce< double , Kokkos::DeviceHost >( 1000000 );
  UnitTestReduce< long , Kokkos::DeviceTPI >( 1000000 );
  UnitTestReduce< double , Kokkos::DeviceTPI >( 1000000 );

  return 0 ;
}

