
#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceTPI.hpp>

#include <Kokkos_ValueView.hpp>
#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_MDArrayView.hpp>

#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_DeviceTPI_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestReduce.hpp>

#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

void test_device_tpi()
{
  try {
    Kokkos::DeviceTPI::initialize( 4 );

    UnitTestDeviceMemoryManagement< Kokkos::DeviceTPI >();
    UnitTestValueView<       Kokkos::DeviceTPI >();
    UnitTestMultiVectorView< Kokkos::DeviceTPI >();
    UnitTestMDArrayView<     Kokkos::DeviceTPI >();
    UnitTestMDArrayDeepCopy< Kokkos::DeviceTPI >();

    UnitTestReduce< long ,   Kokkos::DeviceTPI >( 1000000 );
    UnitTestReduce< double , Kokkos::DeviceTPI >( 1000000 );

    std::cout << "PASSED : UnitTestTPI" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : UnitTestTPI : " << x.what() << std::endl ;
  }

  Kokkos::DeviceTPI::finalize();
}

}

