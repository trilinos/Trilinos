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

TEUCHOS_UNIT_TEST(IMEX_RK, VanDerPol_Staggered_FSA_Tangent)
{
  test_vdp_fsa(false, true, out, success);
}

}  // namespace Tempus_Test
