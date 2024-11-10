// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Array_UnitTest_helpers.hpp"


int ArrayUnitTestHelpers::n = 4;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "n", &ArrayUnitTestHelpers::n, "Number of elements in the array" );
}
