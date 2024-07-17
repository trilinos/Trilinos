// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "FadUnitTests.hpp"

typedef ::testing::Types<
  Sacado::CacheFad::DFad<double>,
  Sacado::CacheFad::SFad<double,5>,
  Sacado::CacheFad::SLFad<double,5>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(CacheFad, FadOpsUnitTest, FadTypes);
