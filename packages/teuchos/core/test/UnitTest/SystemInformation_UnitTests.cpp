// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_SystemInformation.hpp"
#include <stdlib.h>

namespace {


TEUCHOS_UNIT_TEST( SystemInformation, Commands )
{
  if (Teuchos::SystemInformation::commandIsAvailable("ls")) {
    auto output = Teuchos::SystemInformation::runCommandAndCaptureOutput("ls");
    TEST_ASSERT(output.find("TeuchosCore_SystemInformation_UnitTests.exe") != std::string::npos);
  }
}


TEUCHOS_UNIT_TEST( SystemInformation, EnvVariables ) {

  Teuchos::SystemInformation::registerEnvironmentVariable("BLAH_BLAH");
  {
    auto values = Teuchos::SystemInformation::collectSystemInformation();
    TEST_ASSERT(values.find("BLAH_BLAH") != values.end());
    TEST_EQUALITY_CONST(values["BLAH_BLAH"], "NOT SET");
  }

  setenv("BLAH_BLAH", "test", 1);
  {
    auto values = Teuchos::SystemInformation::collectSystemInformation();
    TEST_ASSERT(values.find("BLAH_BLAH") != values.end());
    TEST_EQUALITY_CONST(values["BLAH_BLAH"], "test");
  }
  unsetenv("BLAH_BLAH");
}


TEUCHOS_UNIT_TEST( SystemInformation, CommonlyUsed ) {
  Teuchos::SystemInformation::initializeCollection();
  auto values = Teuchos::SystemInformation::collectSystemInformation();
  TEST_ASSERT(values.find("lscpu") != values.end());
  TEST_ASSERT(values.find("sensors") != values.end());
}

} // namespace
