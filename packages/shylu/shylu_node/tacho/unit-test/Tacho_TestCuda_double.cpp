#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
typedef Kokkos::Cuda DeviceSpaceType;
typedef double ValueType;
typedef double MagnitudeType;

static const std::string MM_TEST_FILE="test_double";

#define __TACHO_TEST_CUDA__
#include "ShyLU_NodeTacho_config.h"
//#include "Tacho_Test.hpp"

//#include "Tacho_TestCrsMatrixBase.hpp"
//#include "Tacho_TestGraph.hpp"
//#include "Tacho_TestSymbolic.hpp"
//#include "Tacho_TestNumeric.hpp"
//#include "Tacho_TestTaskFunctor.hpp"

//#include "Tacho_TestDenseMatrixView.hpp"
//#include "Tacho_TestDenseByBlocks.hpp"

#include "Tacho_TestDenseLinearAlgebra.hpp"

using namespace Tacho;

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
