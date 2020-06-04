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

#include "gtest/gtest.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_balance/balanceUtils.hpp"

class BalanceSettingsTester : public stk::unit_test_util::MeshFixture
{
public:
  stk::balance::BalanceSettings& get_balance_settings()
  {
    return m_settings;
  }

private:
  stk::balance::StkBalanceSettings m_settings;
};

TEST_F(BalanceSettingsTester, defaultInitialDecompMethod)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  EXPECT_EQ(settings.getInitialDecompMethod(), "RIB");
}

TEST_F(BalanceSettingsTester, inputFilename)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  settings.set_input_filename("input.g");
  EXPECT_EQ(settings.get_input_filename(), "input.g");
  settings.set_input_filename("cube.g");
  EXPECT_EQ(settings.get_input_filename(), "cube.g");
}

TEST_F(BalanceSettingsTester, outputFilename)
{
  stk::balance::BalanceSettings& settings = get_balance_settings();
  settings.set_output_filename("output.g");
  EXPECT_EQ(settings.get_output_filename(), "output.g");
  settings.set_output_filename("output/cube.g");
  EXPECT_EQ(settings.get_output_filename(), "output/cube.g");
}

