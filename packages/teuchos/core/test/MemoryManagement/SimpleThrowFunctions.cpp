// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "SimpleThrowFunctions.hpp"
#include "Teuchos_Assert.hpp"


void func_a()
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "This is an exception I throw!");
}


void func_b()
{
  func_a();
}
