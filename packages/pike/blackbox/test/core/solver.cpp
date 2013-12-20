#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include <iostream>

namespace pike {

  TEUCHOS_UNIT_TEST(solver, ostream_overload)
  {
    pike::BlockGaussSeidel solver;

    std::cout << solver << std::endl;
  }

}
