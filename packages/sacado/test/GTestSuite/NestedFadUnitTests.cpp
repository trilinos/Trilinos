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
#if defined(SACADO_ENABLE_NEW_DESIGN) && !defined(SACADO_NEW_FAD_DESIGN_IS_DEFAULT)
  Sacado::Fad::Exp::DFad<Sacado::Fad::Exp::DFad<double> >,
  Sacado::Fad::Exp::SFad<Sacado::Fad::Exp::SFad<double,3>,5>,
  Sacado::Fad::Exp::SLFad<Sacado::Fad::Exp::SLFad<double,3>,5>,
#endif
  Sacado::Fad::DFad<Sacado::Fad::DFad<double> >,
  Sacado::Fad::SFad<Sacado::Fad::SFad<double,3>,5>,
  Sacado::Fad::SLFad<Sacado::Fad::SLFad<double,3>,5>
  > FadTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(FadFad, FadFadOpsUnitTest, FadTypes);
