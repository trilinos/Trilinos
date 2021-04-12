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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Version.hpp>
#include <stk_coupling/impl_VersionUtils.hpp>
#include "TestCompatibilityMode.hpp"

namespace {

class TestCouplingVersion : public stk::coupling::impl::OfficialCouplingVersion
{
public:
  TestCouplingVersion(int inputVersion)
  {
    m_couplingVersion = inputVersion;
  }

  ~TestCouplingVersion()
  {
    stk::coupling::impl::OfficialCouplingVersion().set_version();
  }
};

TEST(UnitTestVersionCompatibility, get_coupling_compatibility_mode_official_version)
{
  unset_compatibility_mode();
  EXPECT_EQ(stk::coupling::impl::Current, stk::coupling::impl::get_coupling_compatibility_mode(MPI_COMM_WORLD));
}

TEST(UnitTestVersionCompatibility, get_coupling_compatibility_mode_backwards_compatible)
{
  TestCompatibilityMode testMode(stk::coupling::impl::BackwardsCompatible);
  EXPECT_EQ(stk::coupling::impl::BackwardsCompatible, stk::coupling::impl::get_coupling_compatibility_mode(MPI_COMM_WORLD));
}

TEST(UnitTestVersionCompatibility, get_coupling_compatibility_mode_deprecated)
{
  TestCompatibilityMode testMode(stk::coupling::impl::Deprecated);
  EXPECT_EQ(stk::coupling::impl::Deprecated, stk::coupling::impl::get_coupling_compatibility_mode(MPI_COMM_WORLD));
}

TEST(UnitTestVersionCompatibility, coupling_compatibility_mode_incompatible)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  int myVersion = (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) ? 1 : 3;
  TestCouplingVersion testCouplingVersion(myVersion);
  testCouplingVersion.set_version();
  EXPECT_EQ(myVersion, stk::coupling::version());

  unset_compatibility_mode();
  EXPECT_EQ(stk::coupling::impl::Incompatible, stk::coupling::impl::get_coupling_compatibility_mode(MPI_COMM_WORLD));
}

TEST(UnitTestVersionCompatibility, coupling_compatibility_mode_version_off_by_one)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  int myVersion = (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) ? 2 : 3;
  TestCouplingVersion testCouplingVersion(myVersion);
  testCouplingVersion.set_version();
  EXPECT_EQ(myVersion, stk::coupling::version());

  stk::coupling::impl::CouplingCompatibilityMode myMode = (myVersion == 2) ? stk::coupling::impl::Deprecated : stk::coupling::impl::BackwardsCompatible;
  unset_compatibility_mode();
  EXPECT_EQ(myMode, stk::coupling::impl::get_coupling_compatibility_mode(MPI_COMM_WORLD));
}

TEST(UnitTestVersionCompatibility, coupling_compatibility_mode_current)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  int myVersion = 52;
  TestCouplingVersion testCouplingVersion(myVersion);
  testCouplingVersion.set_version();
  EXPECT_EQ(myVersion, stk::coupling::version());

  unset_compatibility_mode();
  EXPECT_EQ(stk::coupling::impl::Current, stk::coupling::impl::get_coupling_compatibility_mode(MPI_COMM_WORLD));
}

}
