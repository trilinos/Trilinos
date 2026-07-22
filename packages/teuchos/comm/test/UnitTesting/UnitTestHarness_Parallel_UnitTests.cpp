// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( UnitTestHarness, nonRootFails ) {
  out << "Pass on even procs but fail on other procs!\n";
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  const int procRank = comm->getRank();
  TEST_EQUALITY_CONST(procRank%2, 0);
}


TEUCHOS_UNIT_TEST( UnitTestHarness, nonRootThrowsTeuchosExcept ) {
  out << "Pass on even procs but throws Teuchos exception on other processes!\n";
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  const int procRank = comm->getRank();
  TEUCHOS_ASSERT_EQUALITY(procRank%2, 0); // Throws on non-root processes
  int myval = 1;
  TEST_EQUALITY_CONST(myval, 1);
}


TEUCHOS_UNIT_TEST( UnitTestHarness, nonRootThrowsIntExcept ) {
  out << "Pass on even procs but throws int exception on other processes!\n";
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  const int procRank = comm->getRank();
  if (procRank%2 != 0) {
    throw procRank;
  }
  int myval = 1;
  TEST_EQUALITY_CONST(myval, 1);
}


} // namespace Teuchos



