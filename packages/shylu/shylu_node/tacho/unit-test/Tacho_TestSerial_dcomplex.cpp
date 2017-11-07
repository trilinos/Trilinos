#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
typedef Kokkos::Serial DeviceSpaceType;
typedef Kokkos::complex<double> ValueType;

static const std::string MM_TEST_FILE="test_dcomplex";

#define __TACHO_TEST_SERIAL__
#include "ShyLU_NodeTacho_config.h"
#include "Tacho_Test.hpp"

using namespace Tacho::Experimental;

int main (int argc, char *argv[]) {

  const bool detail = false;
  printExecSpaceConfiguration<DeviceSpaceType>("DeviceSpace", detail);
  printExecSpaceConfiguration<HostSpaceType>  ("HostSpace",   detail);
  
  Kokkos::initialize(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  const int r_val = RUN_ALL_TESTS();

  Kokkos::finalize();

  return r_val;
}
