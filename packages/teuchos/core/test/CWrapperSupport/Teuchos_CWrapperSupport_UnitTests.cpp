// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CWrapperSupport_Cpp.hpp"
#include "someCFunc.h"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


void unitTestSetup(Teuchos::FancyOStream &out)
{
  Teuchos::CWrapperErrorHandling::setPrintErrorOStream(
    Teuchos::rcpFromRef(out));
}


TEUCHOS_UNIT_TEST( CWrapperErrorHandling, normalReturn )
{
  unitTestSetup(out);
  int ierr = 0;
  ECHO(const int output = someCFunc(1, &ierr));
  TEST_EQUALITY_CONST(output, 1);
  TEST_EQUALITY_CONST(ierr, 0);
}


TEUCHOS_UNIT_TEST( CWrapperErrorHandling, throwReturnErrorCode )
{
  unitTestSetup(out);
  int ierr = 0;
  ECHO(const int output = someCFunc(-3, &ierr));
  TEST_EQUALITY_CONST(output, -1);
  TEST_EQUALITY_CONST(ierr, -1);
}


TEUCHOS_UNIT_TEST( CWrapperErrorHandling, nonthrowErrorReturn )
{
  unitTestSetup(out);
  int ierr = 0;
  ECHO(const int output = someCFunc(11, &ierr));
  TEST_EQUALITY_CONST(output, -1);
  TEST_EQUALITY_CONST(ierr, -2);
}


} // namespace
