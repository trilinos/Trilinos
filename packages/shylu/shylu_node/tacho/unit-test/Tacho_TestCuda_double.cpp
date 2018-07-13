#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
typedef Kokkos::Cuda DeviceSpaceType;
typedef double ValueType;
typedef double MagnitudeType;

static const std::string MM_TEST_FILE="test_double";

#define TEST_BEGIN 
#define TEST_END   
//#define TEST_BEGIN Kokkos::initialize()
//#define TEST_END   Kokkos::finalize()

#define __TACHO_TEST_CUDA__
#include "ShyLU_NodeTacho_config.h"
#include "Tacho_Test.hpp"

using namespace Tacho;

int main (int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  TEST_BEGIN;

  const bool detail = false;
  printExecSpaceConfiguration<DeviceSpaceType>("DeviceSpace", detail);
  printExecSpaceConfiguration<HostSpaceType>  ("HostSpace",   detail);

  TEST_END;

  ::testing::InitGoogleTest(&argc, argv);
  const int r_val = RUN_ALL_TESTS();

  Kokkos::finalize();

  return r_val;
}
