// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_BackwardEuler_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(BackwardEuler, SteadyQuadratic_PseudoTransient_FSA)
{
  test_pseudotransient_fsa(false, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SteadyQuadratic_PseudoTransient_FSA_Tangent)
{
  test_pseudotransient_fsa(true, out, success);
}

} // namespace Tempus_Test
