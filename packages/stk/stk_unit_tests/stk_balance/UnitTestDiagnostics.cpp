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

namespace {

class UnitTestDiagnostics : public ::testing::Test {};

class UnsignedDiagnosticTester : public stk::balance::UnsignedDiagnostic
{
public:
  virtual std::string print_header1() override { return "Test"; }
  virtual std::string print_header2() override { return "Diagnostic"; }
};

class UnsignedDiagnosticTester2 : public stk::balance::UnsignedDiagnostic
{
public:
  virtual std::string print_header1() override { return "Test"; }
  virtual std::string print_header2() override { return "Diag"; }
};

class UnsignedWithPercentDiagnosticTester : public stk::balance::UnsignedWithPercentDiagnostic
{
public:
  virtual std::string print_header1() override { return "Test"; }
  virtual std::string print_header2() override { return "Diagnostic"; }
};

class DoubleWithPercentDiagnosticTester : public stk::balance::DoubleWithPercentDiagnostic
{
public:
  virtual std::string print_header1() override { return "Test"; }
  virtual std::string print_header2() override { return "Diagnostic"; }
};


TEST_F(UnitTestDiagnostics, oneRealRank_oneLogicalRank_tooManyValues)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const int numLogicalRanks = 1;
  UnsignedDiagnosticTester diagTester;

  diagTester.store_value(0, 123);
  diagTester.store_value(1, 234);

  EXPECT_ANY_THROW(diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks));
}

TEST_F(UnitTestDiagnostics, oneRealRank_oneLogicalRank)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const int numLogicalRanks = 1;
  UnsignedDiagnosticTester diagTester;

  diagTester.store_value(0, 123);

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  EXPECT_EQ(diagTester.print_header1(), "Test");
  EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

  EXPECT_EQ(diagTester.print_rank_value(0), "123");

  EXPECT_EQ(diagTester.print_min(), "123");
  EXPECT_EQ(diagTester.print_max(), "123");
  EXPECT_EQ(diagTester.print_avg(), "123");
}

TEST_F(UnitTestDiagnostics, oneRealRank_twoLogicalRanks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const int numLogicalRanks = 2;
  UnsignedDiagnosticTester diagTester;

  diagTester.store_value(0, 10);
  diagTester.store_value(1, 20);

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  EXPECT_EQ(diagTester.print_header1(), "Test");
  EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

  EXPECT_EQ(diagTester.print_rank_value(0), "10");
  EXPECT_EQ(diagTester.print_rank_value(1), "20");

  EXPECT_EQ(diagTester.print_min(), "10");
  EXPECT_EQ(diagTester.print_max(), "20");
  EXPECT_EQ(diagTester.print_avg(), "15");
}

TEST_F(UnitTestDiagnostics, twoRealRanks_oneLogicalRank)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 1;

  UnsignedDiagnosticTester diagTester;

  if (parallelRank == 1) {
    diagTester.store_value(0, 10);
  }

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  if (parallelRank == 0) {
    EXPECT_EQ(diagTester.print_header1(), "Test");
    EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

    EXPECT_EQ(diagTester.print_rank_value(0), "10");

    EXPECT_EQ(diagTester.print_min(), "10");
    EXPECT_EQ(diagTester.print_max(), "10");
    EXPECT_EQ(diagTester.print_avg(), "10");
  }
}

TEST_F(UnitTestDiagnostics, twoRealRanks_twoLogicalRanks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 2;

  UnsignedDiagnosticTester diagTester;

  if (parallelRank == 0) {
    diagTester.store_value(0, 10);
  }
  if (parallelRank == 1) {
    diagTester.store_value(1, 20);
  }

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  if (parallelRank == 0) {
    EXPECT_EQ(diagTester.print_header1(), "Test");
    EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

    EXPECT_EQ(diagTester.print_rank_value(0), "10");
    EXPECT_EQ(diagTester.print_rank_value(1), "20");

    EXPECT_EQ(diagTester.print_min(), "10");
    EXPECT_EQ(diagTester.print_max(), "20");
    EXPECT_EQ(diagTester.print_avg(), "15");
  }
}

TEST_F(UnitTestDiagnostics, twoRealRanks_threeLogicalRanks)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 3;

  UnsignedDiagnosticTester diagTester;

  if (parallelRank == 0) {
    diagTester.store_value(1, 15);
  }
  if (parallelRank == 1) {
    diagTester.store_value(0, 10);
    diagTester.store_value(2, 20);
  }

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  if (parallelRank == 0) {
    EXPECT_EQ(diagTester.print_header1(), "Test");
    EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

    EXPECT_EQ(diagTester.print_rank_value(0), "10");
    EXPECT_EQ(diagTester.print_rank_value(1), "15");
    EXPECT_EQ(diagTester.print_rank_value(2), "20");

    EXPECT_EQ(diagTester.print_min(), "10");
    EXPECT_EQ(diagTester.print_max(), "20");
    EXPECT_EQ(diagTester.print_avg(), "15");
  }
}

TEST_F(UnitTestDiagnostics, twoRealRanks_threeLogicalRanks_unsignedWithPercent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 3;

  UnsignedWithPercentDiagnosticTester diagTester;

  if (parallelRank == 0) {
    diagTester.store_value(0, 10);
    diagTester.store_value(2, 400);
  }
  if (parallelRank == 1) {
    diagTester.store_value(1, 100);
  }

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  if (parallelRank == 0) {
    EXPECT_EQ(diagTester.print_header1(), "Test");
    EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

    EXPECT_EQ(diagTester.print_rank_value(0),  "10 ( -94.1%)");
    EXPECT_EQ(diagTester.print_rank_value(1), "100 ( -41.2%)");
    EXPECT_EQ(diagTester.print_rank_value(2), "400 (+135.3%)");

    EXPECT_EQ(diagTester.print_min(),  "10 ( -94.1%)");
    EXPECT_EQ(diagTester.print_max(), "400 (+135.3%)");
    EXPECT_EQ(diagTester.print_avg(), "170          ");
  }
}

TEST_F(UnitTestDiagnostics, twoRealRanks_threeLogicalRanks_doubleWithPercent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;
  const int parallelRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const int numLogicalRanks = 3;

  DoubleWithPercentDiagnosticTester diagTester;

  if (parallelRank == 0) {
    diagTester.store_value(0, -0.2);
    diagTester.store_value(2, 0.2);
  }
  if (parallelRank == 1) {
    diagTester.store_value(1, 0.6);
  }

  diagTester.collect_data(MPI_COMM_WORLD, numLogicalRanks);
  diagTester.process_data(MPI_COMM_WORLD);

  if (parallelRank == 0) {
    EXPECT_EQ(diagTester.print_header1(), "Test");
    EXPECT_EQ(diagTester.print_header2(), "Diagnostic");

    EXPECT_EQ(diagTester.print_rank_value(0), "-0.200 (-200.0%)");
    EXPECT_EQ(diagTester.print_rank_value(1),  "0.600 (+200.0%)");
    EXPECT_EQ(diagTester.print_rank_value(2),  "0.200 (  +0.0%)");

    EXPECT_EQ(diagTester.print_min(), "-0.200 (-200.0%)");
    EXPECT_EQ(diagTester.print_max(),  "0.600 (+200.0%)");
    EXPECT_EQ(diagTester.print_avg(),  "0.200          ");
  }
}

}
