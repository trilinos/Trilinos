// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_ADOLC

#include "TayUnitTests.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"

typedef ::testing::Types<
  Sacado::Tay::Taylor<double>,
  Sacado::Tay::CacheTaylor<double>
  > TayTypes;

typedef ::testing::Types<
  Sacado::Tay::Taylor<double>
  > MaxMinTayTypes;

INSTANTIATE_TYPED_TEST_SUITE_P(Taylor, TaylorOpsUnitTest, TayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MaxMinTaylor, TaylorMaxMinUnitTest, MaxMinTayTypes);

#endif
