// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_BackwardEuler_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Combined_FSA)
{
  test_sincos_fsa(true, false, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Combined_FSA_Tangent)
{
  test_sincos_fsa(true, true, out, success);
}

} // namespace Tempus_Test
