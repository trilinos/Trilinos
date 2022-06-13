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
#include "stk_util/parallel/Parallel.hpp"
#include "stk_balance/internal/Diagnostics.hpp"
#include "stk_balance/internal/DiagnosticsPrinter.hpp"
#include "stk_balance/internal/DiagnosticsContainer.hpp"

namespace {

class UnitTestDiagnosticsPrinter : public ::testing::Test
{
public:
  UnitTestDiagnosticsPrinter() = default;
  ~UnitTestDiagnosticsPrinter() override = default;

  virtual void SetUp() override {
    stk::balance::impl::g_diagnosticsContainer.clear();
  }

  virtual void TearDown() override {
    stk::balance::impl::g_diagnosticsContainer.clear();
  }
};

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
  virtual std::string print_header2(unsigned ) override { return "Diag"; }
};

class UnsignedWithPercentDiagnosticTester : public stk::balance::UnsignedWithPercentDiagnostic
{
public:
  virtual std::string print_header1(unsigned ) override { return "Test"; }
  virtual std::string print_header2(unsigned ) override { return "Diagnostic"; }
};

class DoubleWithPercentDiagnosticTester : public stk::balance::DoubleWithPercentDiagnostic
{
public:
  virtual std::string print_header1(unsigned ) override { return "Test"; }
  virtual std::string print_header2(unsigned ) override { return "Diagnostic"; }
};

class MultiUnsignedDiagnosticTester : public stk::balance::MultiUnsignedDiagnostic
{
public:
  MultiUnsignedDiagnosticTester(unsigned numColumns)
    : MultiUnsignedDiagnostic(numColumns)
  {}
  virtual std::string print_header1(unsigned ) override { return "Test"; }
  virtual std::string print_header2(unsigned column) override { return "Diag " + std::to_string(column + 1); }
};

class MultiUnsignedWithPercentDiagnosticTester : public stk::balance::MultiUnsignedWithPercentDiagnostic
{
public:
  MultiUnsignedWithPercentDiagnosticTester(unsigned numColumns)
    : MultiUnsignedWithPercentDiagnostic(numColumns)
  {}
  virtual std::string print_header1(unsigned ) override { return "Test"; }
  virtual std::string print_header2(unsigned column) override { return "Diag " + std::to_string(column + 1); }
};

TEST_F(UnitTestDiagnosticsPrinter, empty)
{
  const int numLogicalRanks = 2;
  std::ostringstream os;

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  EXPECT_EQ(os.str(), "");
}

TEST_F(UnitTestDiagnosticsPrinter, unsignedDiagnostic)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 2;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedDiagnosticTester>();

  auto * diag = stk::balance::get_diagnostic<UnsignedDiagnosticTester>();
  if (parallelRank == 0) {
    diag->store_value(parallelRank, 10);
  }
  if (parallelRank == 1) {
    diag->store_value(parallelRank, 100);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "===================\n"
                                     "  MPI |    Test    \n"
                                     " Rank | Diagnostic \n"
                                     "===================\n"
                                     "   0  |      10    \n"
                                     "   1  |     100    \n"
                                     "===================\n"
                                     "\n"
                                     "===================\n"
                                     "      |    Test    \n"
                                     "  Qty | Diagnostic \n"
                                     "===================\n"
                                     "  Min |      10    \n"
                                     "  Max |     100    \n"
                                     "  Avg |      55    \n"
                                     "===================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, twoUnsignedDiagnostics)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 2;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedDiagnosticTester>();
  stk::balance::register_diagnostic<UnsignedDiagnosticTester2>();

  auto * diag1 = stk::balance::get_diagnostic<UnsignedDiagnosticTester>();
  auto * diag2 = stk::balance::get_diagnostic<UnsignedDiagnosticTester2>();
  if (parallelRank == 0) {
    diag1->store_value(parallelRank, 10);
    diag2->store_value(parallelRank, 1000000000);
  }
  if (parallelRank == 1) {
    diag1->store_value(parallelRank, 100);
    diag2->store_value(parallelRank, 10);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "================================\n"
                                     "  MPI |    Test    |    Test    \n"
                                     " Rank | Diagnostic |    Diag    \n"
                                     "================================\n"
                                     "   0  |      10    | 1000000000 \n"
                                     "   1  |     100    |         10 \n"
                                     "================================\n"
                                     "\n"
                                     "================================\n"
                                     "      |    Test    |    Test    \n"
                                     "  Qty | Diagnostic |    Diag    \n"
                                     "================================\n"
                                     "  Min |      10    |         10 \n"
                                     "  Max |     100    | 1000000000 \n"
                                     "  Avg |      55    |  500000005 \n"
                                     "================================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, twoUnsignedDiagnostics_oneLogicalRank)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 1;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedDiagnosticTester>();
  stk::balance::register_diagnostic<UnsignedDiagnosticTester2>();

  auto * diag1 = stk::balance::get_diagnostic<UnsignedDiagnosticTester>();
  auto * diag2 = stk::balance::get_diagnostic<UnsignedDiagnosticTester2>();
  if (parallelRank == 1) {
    diag1->store_value(0, 10);
    diag2->store_value(0, 1000000000);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "================================\n"
                                     "  MPI |    Test    |    Test    \n"
                                     " Rank | Diagnostic |    Diag    \n"
                                     "================================\n"
                                     "   0  |     10     | 1000000000 \n"
                                     "================================\n"
                                     "\n"
                                     "================================\n"
                                     "      |    Test    |    Test    \n"
                                     "  Qty | Diagnostic |    Diag    \n"
                                     "================================\n"
                                     "  Min |     10     | 1000000000 \n"
                                     "  Max |     10     | 1000000000 \n"
                                     "  Avg |     10     | 1000000000 \n"
                                     "================================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, twoUnsignedDiagnostics_threeLogicalRanks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 3;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedDiagnosticTester>();
  stk::balance::register_diagnostic<UnsignedDiagnosticTester2>();

  auto * diag1 = stk::balance::get_diagnostic<UnsignedDiagnosticTester>();
  auto * diag2 = stk::balance::get_diagnostic<UnsignedDiagnosticTester2>();
  if (parallelRank == 0) {
    diag1->store_value(0, 10);
    diag1->store_value(2, 55);
    diag2->store_value(0, 1000000000);
  }
  if (parallelRank == 1) {
    diag1->store_value(1, 100);
    diag2->store_value(1, 10);
    diag2->store_value(2, 500000005);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "================================\n"
                                     "  MPI |    Test    |    Test    \n"
                                     " Rank | Diagnostic |    Diag    \n"
                                     "================================\n"
                                     "   0  |      10    | 1000000000 \n"
                                     "   1  |     100    |         10 \n"
                                     "   2  |      55    |  500000005 \n"
                                     "================================\n"
                                     "\n"
                                     "================================\n"
                                     "      |    Test    |    Test    \n"
                                     "  Qty | Diagnostic |    Diag    \n"
                                     "================================\n"
                                     "  Min |      10    |         10 \n"
                                     "  Max |     100    | 1000000000 \n"
                                     "  Avg |      55    |  500000005 \n"
                                     "================================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, unsignedWithPercentDiagnostics)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 4;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedWithPercentDiagnosticTester>();

  auto * diag = stk::balance::get_diagnostic<UnsignedWithPercentDiagnosticTester>();
  if (parallelRank == 0) {
    diag->store_value(0, 10);
    diag->store_value(3, 170);
  }
  if (parallelRank == 1) {
    diag->store_value(1, 400);
    diag->store_value(2, 100);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "======================\n"
                                     "  MPI |      Test     \n"
                                     " Rank |   Diagnostic  \n"
                                     "======================\n"
                                     "   0  |  10 ( -94.1%) \n"
                                     "   1  | 400 (+135.3%) \n"
                                     "   2  | 100 ( -41.2%) \n"
                                     "   3  | 170 (  +0.0%) \n"
                                     "======================\n"
                                     "\n"
                                     "======================\n"
                                     "      |      Test     \n"
                                     "  Qty |   Diagnostic  \n"
                                     "======================\n"
                                     "  Min |  10 ( -94.1%) \n"
                                     "  Max | 400 (+135.3%) \n"
                                     "  Avg | 170           \n"
                                     "======================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, doubleWithPercentDiagnostics)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 3;
  std::ostringstream os;

  stk::balance::register_diagnostic<DoubleWithPercentDiagnosticTester>();

  auto * diag = stk::balance::get_diagnostic<DoubleWithPercentDiagnosticTester>();
  if (parallelRank == 0) {
    diag->store_value(0, -0.2);
    diag->store_value(2, 0.2);
  }
  if (parallelRank == 1) {
    diag->store_value(1, 0.6);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "=========================\n"
                                     "  MPI |       Test       \n"
                                     " Rank |    Diagnostic    \n"
                                     "=========================\n"
                                     "   0  | -0.200 (-200.0%) \n"
                                     "   1  |  0.600 (+200.0%) \n"
                                     "   2  |  0.200 (  +0.0%) \n"
                                     "=========================\n"
                                     "\n"
                                     "=========================\n"
                                     "      |       Test       \n"
                                     "  Qty |    Diagnostic    \n"
                                     "=========================\n"
                                     "  Min | -0.200 (-200.0%) \n"
                                     "  Max |  0.600 (+200.0%) \n"
                                     "  Avg |  0.200           \n"
                                     "=========================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, oneUnsignedDiagnostic_oneMultiUnsignedDiagnostic_threeLogicalRanks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 3;
  const unsigned numColumns = 2;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedDiagnosticTester>();
  stk::balance::register_diagnostic<MultiUnsignedDiagnosticTester>(numColumns);

  auto * diag1 = stk::balance::get_diagnostic<UnsignedDiagnosticTester>();
  auto * diag2 = stk::balance::get_diagnostic<MultiUnsignedDiagnosticTester>();
  if (parallelRank == 0) {
    diag1->store_value(0, 10);
    diag1->store_value(2, 55);

    diag2->store_value(0, 0, 10000000);
    diag2->store_value(1, 0, 1000000000);
  }
  if (parallelRank == 1) {
    diag1->store_value(1, 100);

    diag2->store_value(0, 1, 10);
    diag2->store_value(0, 2, 5000005);
    diag2->store_value(1, 1, 10);
    diag2->store_value(1, 2, 500000005);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "===========================================\n"
                                     "  MPI |    Test    |   Test   |    Test    \n"
                                     " Rank | Diagnostic |  Diag 1  |   Diag 2   \n"
                                     "===========================================\n"
                                     "   0  |      10    | 10000000 | 1000000000 \n"
                                     "   1  |     100    |       10 |         10 \n"
                                     "   2  |      55    |  5000005 |  500000005 \n"
                                     "===========================================\n"
                                     "\n"
                                     "===========================================\n"
                                     "      |    Test    |   Test   |    Test    \n"
                                     "  Qty | Diagnostic |  Diag 1  |   Diag 2   \n"
                                     "===========================================\n"
                                     "  Min |      10    |       10 |         10 \n"
                                     "  Max |     100    | 10000000 | 1000000000 \n"
                                     "  Avg |      55    |  5000005 |  500000005 \n"
                                     "===========================================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

TEST_F(UnitTestDiagnosticsPrinter, oneUnsignedDiagnostic_oneMultiUnsignedWithPercentDiagnostic_fourLogicalRanks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 4;
  const unsigned numColumns = 2;
  std::ostringstream os;

  stk::balance::register_diagnostic<UnsignedDiagnosticTester>();
  stk::balance::register_diagnostic<MultiUnsignedWithPercentDiagnosticTester>(numColumns);

  auto * diag1 = stk::balance::get_diagnostic<UnsignedDiagnosticTester>();
  auto * diag2 = stk::balance::get_diagnostic<MultiUnsignedWithPercentDiagnosticTester>();
  if (parallelRank == 0) {
    diag1->store_value(0, 10);
    diag1->store_value(3, 50);

    diag2->store_value(0, 0, 10);
    diag2->store_value(0, 3, 170);
    diag2->store_value(1, 0, 100);
    diag2->store_value(1, 3, 1700);
  }
  if (parallelRank == 1) {
    diag1->store_value(1, 100);
    diag1->store_value(2, 40);

    diag2->store_value(0, 1, 400);
    diag2->store_value(0, 2, 100);
    diag2->store_value(1, 1, 4000);
    diag2->store_value(1, 2, 1000);
  }

  stk::balance::DiagnosticsPrinter diagPrinter(MPI_COMM_WORLD, numLogicalRanks);
  diagPrinter.print(os);

  const std::string expectedResult = "\n"
                                     "====================================================\n"
                                     "  MPI |    Test    |      Test     |      Test      \n"
                                     " Rank | Diagnostic |     Diag 1    |     Diag 2     \n"
                                     "====================================================\n"
                                     "   0  |      10    |  10 ( -94.1%) |  100 ( -94.1%) \n"
                                     "   1  |     100    | 400 (+135.3%) | 4000 (+135.3%) \n"
                                     "   2  |      40    | 100 ( -41.2%) | 1000 ( -41.2%) \n"
                                     "   3  |      50    | 170 (  +0.0%) | 1700 (  +0.0%) \n"
                                     "====================================================\n"
                                     "\n"
                                     "====================================================\n"
                                     "      |    Test    |      Test     |      Test      \n"
                                     "  Qty | Diagnostic |     Diag 1    |     Diag 2     \n"
                                     "====================================================\n"
                                     "  Min |      10    |  10 ( -94.1%) |  100 ( -94.1%) \n"
                                     "  Max |     100    | 400 (+135.3%) | 4000 (+135.3%) \n"
                                     "  Avg |      50    | 170           | 1700           \n"
                                     "====================================================\n";
  if (parallelRank == 0) {
    EXPECT_EQ(os.str(), expectedResult);
  }
}

}
