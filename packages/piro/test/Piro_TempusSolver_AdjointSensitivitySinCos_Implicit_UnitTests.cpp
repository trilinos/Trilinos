// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_TempusSolver_AdjointSensitivitySinCosUnitTests.hpp"

#ifdef HAVE_PIRO_TEMPUS

TEUCHOS_UNIT_TEST(Piro_TempusSolver, SinCos_AdjointSensitivities_ImplicitAdjointModel)
{
  test_sincos_asa(out, success, false);
}

#endif /* HAVE_PIRO_TEMPUS */
