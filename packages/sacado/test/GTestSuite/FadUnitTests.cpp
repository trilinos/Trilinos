// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "FadUnitTests.hpp"

#include "Sacado_Fad_SimpleFad.hpp"

typedef ::testing::Types<
  Sacado::Fad::DFad<double>,
  Sacado::Fad::SFad<double,5>,
  Sacado::Fad::SLFad<double,5>,

  Sacado::Fad::SimpleFad<double>,
  Sacado::Fad::DVFad<double>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(Fad, FadOpsUnitTest, FadTypes);
