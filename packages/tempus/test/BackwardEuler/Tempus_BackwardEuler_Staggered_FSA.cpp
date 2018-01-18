// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_BackwardEuler_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Staggered_FSA)
{
  test_sincos_fsa(false, false, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Staggered_FSA_Tangent)
{
  test_sincos_fsa(false, true, out, success);
}

} // namespace Tempus_Test
