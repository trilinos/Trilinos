// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cstdlib>
#include <cstring>
#include <utility>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_SystemInformation.hpp"

namespace {

// setenv/unsetenv are POSIX; MSVC CRT provides _putenv_s instead.
// Other Windows toolchains (e.g. MinGW) may still expose setenv/unsetenv.
#if defined(_MSC_VER)
void teuchosTestSetEnv(const char* name, const char* value, int overwrite)
{
  if(name == nullptr || value == nullptr || std::strchr(name, '=') != nullptr) {
    return;
  }
  if(overwrite == 0) {
    char* buf{};
    size_t bufSize{};
    if(_dupenv_s(std::addressof(buf), std::addressof(bufSize), name) == 0 && buf != nullptr) {
      std::free(buf);
      return;
    }
  }
  (void)_putenv_s(name, value);
}

void teuchosTestUnsetEnv(const char* name)
{
  if(name == nullptr || std::strchr(name, '=') != nullptr) {
    return;
  }
  (void)_putenv_s(name, "");
}
#else  // !defined(_MSC_VER)
void teuchosTestSetEnv(const char* name, const char* value, int overwrite)
{
  if(name == nullptr || value == nullptr || std::strchr(name, '=') != nullptr) {
    return;
  }
  (void)setenv(name, value, overwrite);
}

void teuchosTestUnsetEnv(const char* name)
{
  if(name == nullptr || std::strchr(name, '=') != nullptr) {
    return;
  }
  (void)unsetenv(name);
}
#endif


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

  teuchosTestSetEnv("BLAH_BLAH", "test", 1);
  {
    auto values = Teuchos::SystemInformation::collectSystemInformation();
    TEST_ASSERT(values.find("BLAH_BLAH") != values.end());
    TEST_EQUALITY_CONST(values["BLAH_BLAH"], "test");
  }
  teuchosTestUnsetEnv("BLAH_BLAH");
}


TEUCHOS_UNIT_TEST( SystemInformation, CommonlyUsed ) {
  Teuchos::SystemInformation::initializeCollection();
  auto values = Teuchos::SystemInformation::collectSystemInformation();
  TEST_ASSERT(values.find("lscpu") != values.end());
  TEST_ASSERT(values.find("sensors") != values.end());
}

} // namespace
