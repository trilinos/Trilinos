// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "FadUnitTests2.hpp"

typedef ::testing::Types<
  Sacado::ELRCacheFad::DFad<double>,
  Sacado::ELRCacheFad::SFad<double,5>,
  Sacado::ELRCacheFad::SLFad<double,5>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(RealELRCacheFad, FadOpsUnitTest2, FadTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(RealELRCacheFad, RealFadOpsUnitTest2, FadTypes);
