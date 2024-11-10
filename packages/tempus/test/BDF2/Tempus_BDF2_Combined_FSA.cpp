//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_BDF2_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(BDF2, SinCos_Combined_FSA)
{
  test_sincos_fsa(true, false, out, success);
}

TEUCHOS_UNIT_TEST(BDF2, SinCos_Combined_FSA_Tangent)
{
  test_sincos_fsa(true, true, out, success);
}

}  // namespace Tempus_Test
