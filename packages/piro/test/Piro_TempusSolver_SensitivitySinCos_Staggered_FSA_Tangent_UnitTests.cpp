// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_TempusSolver_SensitivitySinCos_FSA.hpp"

#ifdef HAVE_PIRO_TEMPUS

TEUCHOS_UNIT_TEST(Piro_TempusSolver, SinCos_Staggered_FSA_Tangent)
{
  test_sincos_fsa(false, true, out, success);
}

#endif /* HAVE_PIRO_TEMPUS */
