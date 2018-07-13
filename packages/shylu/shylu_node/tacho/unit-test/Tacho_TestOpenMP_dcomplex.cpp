#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
typedef Kokkos::OpenMP                    DeviceSpaceType;
typedef Kokkos::complex<double> ValueType;
typedef double MagnitudeType;

static const std::string MM_TEST_FILE="test_dcomplex";

#define TEST_BEGIN
#define TEST_END

#define __TACHO_TEST_OPENMP__
#include "ShyLU_NodeTacho_config.h"
#include "Tacho_Test.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  const bool detail = false;
  printExecSpaceConfiguration<DeviceSpaceType>("DeviceSpace", detail);
  printExecSpaceConfiguration<HostSpaceType>  ("HostSpace",   detail);

  ::testing::InitGoogleTest(&argc, argv);
  const int r_val = RUN_ALL_TESTS();

  Kokkos::finalize();

  return r_val;
}
