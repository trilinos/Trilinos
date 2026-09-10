// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "NestedFadUnitTests.hpp"

typedef ::testing::Types<

  Sacado::Fad::DFad<Sacado::Fad::DFad<double> >,
  Sacado::Fad::SFad<Sacado::Fad::SFad<double,3>,5>,
  Sacado::Fad::SLFad<Sacado::Fad::SLFad<double,3>,5>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(FadFad, FadFadOpsUnitTest, FadTypes);
