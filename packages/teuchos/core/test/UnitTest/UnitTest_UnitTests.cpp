// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_dyn_cast.hpp"


namespace {


TEUCHOS_UNIT_TEST( UnitTest, verbose )
{
  // This test checks to see that the 'verbose' bool is returned correctly.
  // This test uses knowlege of the internals of
  // Teuchos::UnitTestRepository::runUnitTests(...) to determine if this is
  // set correctly according to the --details input option.  This test is
  // *very* closely tied to the interneral implemetation of
  // Teuchos::UnitTestRepository.
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::FancyOStream;
  bool verbose = Teuchos::UnitTestRepository::verboseUnitTests();
  const Teuchos::RCP<std::ostream> wrappedOut = out.getOStream();
  if (verbose) {
    TEST_THROW(dyn_cast<std::ostringstream>(*wrappedOut), Teuchos::m_bad_cast);
  }
  else {
    TEST_NOTHROW(dyn_cast<std::ostringstream>(*wrappedOut));
  }
}


} // namespace
