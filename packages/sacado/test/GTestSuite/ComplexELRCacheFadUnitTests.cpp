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

typedef ::testing::Types<
  Sacado::ELRCacheFad::DFad<std::complex<double>>,
  Sacado::ELRCacheFad::SFad<std::complex<double>,5>,
  Sacado::ELRCacheFad::SLFad<std::complex<double>,5>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(ComplexELRCacheFad, FadOpsUnitTest2, FadTypes);

#endif
