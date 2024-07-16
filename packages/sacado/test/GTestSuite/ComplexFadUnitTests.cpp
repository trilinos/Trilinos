// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_SACADO_COMPLEX

#include "FadUnitTests2.hpp"

#include "Sacado_Fad_SimpleFad.hpp"

typedef ::testing::Types<
  Sacado::Fad::DFad<std::complex<double>>,
  Sacado::Fad::SFad<std::complex<double>,5>,
  Sacado::Fad::SLFad<std::complex<double>,5>,
  Sacado::Fad::SimpleFad<std::complex<double>>,
  Sacado::Fad::DVFad<std::complex<double>>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(ComplexFad, FadOpsUnitTest2, FadTypes);

#endif
