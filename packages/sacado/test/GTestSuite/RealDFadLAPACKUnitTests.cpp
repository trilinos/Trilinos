// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "FadLAPACKUnitTests.hpp"

typedef ::testing::Types< Sacado::Fad::DFad<double> > FadTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(DFad, FadLAPACKUnitTests, FadTypes);
