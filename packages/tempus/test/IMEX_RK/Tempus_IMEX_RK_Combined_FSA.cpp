//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_IMEX_RK_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(IMEX_RK, VanDerPol_Combined_FSA)
{
  test_vdp_fsa(true, false, out, success);
}

}  // namespace Tempus_Test
