//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

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

// NOX row sum scaling: the state and sensitivity steppers are created from the
// same parameters, in which the NOX solver of the state stepper stores its row
// sum scaling objects (in the state residual space).  The sensitivity stepper
// must not use them.
TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Staggered_FSA_RowSumScaling)
{
  test_sincos_fsa(false, false, out, success, true);
}

}  // namespace Tempus_Test
