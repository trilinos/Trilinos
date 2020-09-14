#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

static const std::string MM_TEST_FILE="test_dcomplex";

#define TEST_BEGIN
#define TEST_END  
//#define TEST_BEGIN Kokkos::initialize()
//#define TEST_END   Kokkos::finalize()

#define __TACHO_TEST_CUDA__
#include "Tacho_config.h"
#include "Tacho_Util.hpp"

typedef typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type HostDeviceType;
typedef typename Tacho::UseThisDevice<Kokkos::Cuda>::device_type DeviceType;

typedef Kokkos::complex<double> ValueType;
typedef double MagnitudeType;

#include "Tacho_Test.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  TEST_BEGIN;

  const bool detail = false;
  printExecSpaceConfiguration<typename DeviceType::execution_space>("DeviceSpace", detail);
  printExecSpaceConfiguration<typename HostDeviceType::execution_space>("HostSpace",   detail);

  TEST_END;

  ::testing::InitGoogleTest(&argc, argv);
  const int r_val = RUN_ALL_TESTS();

  Kokkos::finalize();

  return r_val;
}
