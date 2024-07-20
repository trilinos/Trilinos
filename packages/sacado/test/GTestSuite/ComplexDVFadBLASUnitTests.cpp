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

#include "FadBLASUnitTests.hpp"

typedef ::testing::Types< Sacado::Fad::DVFad< std::complex<double> > > FadTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(ComplexDVFad, FadBLASUnitTests, FadTypes);

#endif
