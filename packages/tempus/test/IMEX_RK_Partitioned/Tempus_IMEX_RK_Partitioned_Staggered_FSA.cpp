// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_IMEX_RK_Partitioned_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(IMEX_RK_Partitioned, VanDerPol_Staggered_FSA)
{
  test_vdp_fsa(false, false, out, success);
}

TEUCHOS_UNIT_TEST(IMEX_RK_Partitioned, VanDerPol_Staggered_FSA_Tangent)
{
  test_vdp_fsa(false, true, out, success);
}

} // namespace Tempus_Test
