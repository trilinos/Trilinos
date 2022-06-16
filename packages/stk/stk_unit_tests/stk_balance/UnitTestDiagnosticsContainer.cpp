// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "gtest/gtest.h"
#include "stk_balance/internal/Diagnostics.hpp"
#include "stk_balance/internal/DiagnosticsContainer.hpp"

namespace {

class UnitTestDiagnosticsContainer : public ::testing::Test {};

class UnsignedDiagnosticTester : public stk::balance::UnsignedDiagnostic
{
public:
  virtual std::string print_header1(unsigned ) override { return "Test"; }
  virtual std::string print_header2(unsigned ) override { return "Diagnostic"; }
};

class UnsignedDiagnosticTester2 : public stk::balance::UnsignedDiagnostic
{
public:
  virtual std::string print_header1(unsigned ) override { return "Test"; }
  virtual std::string print_header2(unsigned ) override { return "Diag2"; }
};


TEST_F(UnitTestDiagnosticsContainer, queryDiagnostic_notFound)
{
  stk::balance::DiagnosticsContainer diagContainer;
  EXPECT_EQ(diagContainer.get<UnsignedDiagnosticTester>(), nullptr);
}

TEST_F(UnitTestDiagnosticsContainer, registration)
{
  stk::balance::DiagnosticsContainer diagContainer;
  diagContainer.register_diagnostic<UnsignedDiagnosticTester>();
  diagContainer.register_diagnostic<UnsignedDiagnosticTester2>();

  auto retrievedUnsignedDiagnosticTester = diagContainer.get<UnsignedDiagnosticTester>();
  auto retrievedUnsignedDiagnosticTester2 = diagContainer.get<UnsignedDiagnosticTester2>();

  ASSERT_NE(retrievedUnsignedDiagnosticTester, nullptr);
  ASSERT_NE(retrievedUnsignedDiagnosticTester2, nullptr);
  EXPECT_EQ(typeid(*retrievedUnsignedDiagnosticTester), typeid(UnsignedDiagnosticTester));
  EXPECT_EQ(typeid(*retrievedUnsignedDiagnosticTester2), typeid(UnsignedDiagnosticTester2));
}

TEST_F(UnitTestDiagnosticsContainer, duplicateRegistration)
{
  stk::balance::DiagnosticsContainer diagContainer;
  diagContainer.register_diagnostic<UnsignedDiagnosticTester>();
  EXPECT_ANY_THROW(diagContainer.register_diagnostic<UnsignedDiagnosticTester>());
}

TEST_F(UnitTestDiagnosticsContainer, iterateDiagnostics)
{
  stk::balance::DiagnosticsContainer diagContainer;
  diagContainer.register_diagnostic<UnsignedDiagnosticTester>();
  diagContainer.register_diagnostic<UnsignedDiagnosticTester2>();

  std::vector<const std::type_info *> expectedTypes { &typeid(UnsignedDiagnosticTester),
                                                      &typeid(UnsignedDiagnosticTester2) };

  unsigned index = 0;
  for (auto it = diagContainer.begin(); it != diagContainer.end(); ++it, ++index) {
    auto diagContainerEntry = *it;
    EXPECT_EQ(typeid(*diagContainerEntry), *expectedTypes[index]);
  }
}

}
