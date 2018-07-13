// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include "Tempus_Version.hpp"

namespace Tempus_Test {

  TEUCHOS_UNIT_TEST(version, default)
  {
    Tempus::version();
  }

}
