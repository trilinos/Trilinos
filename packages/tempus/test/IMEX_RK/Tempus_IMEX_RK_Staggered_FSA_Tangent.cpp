// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_IMEX_RK_FSA.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(IMEX_RK, VanDerPol_Staggered_FSA_Tangent)
{
  test_vdp_fsa(false, true, out, success);
}

} // namespace Tempus_Test
