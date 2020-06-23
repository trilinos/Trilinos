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
#include "StkBalanceRunner.hpp"
#include "stk_balance/setup/LifeCycle.hpp"

namespace stk {
namespace integration_test_utils {

StkBalanceRunner::StkBalanceRunner(MPI_Comm c)
  : m_comm(c),
    m_execName("stk_balance"),
    m_inFile(""),
    m_outputDirectory(""),
    m_appTypeDefaultArg(""),
    m_decompMethodArg("")
{ }

void StkBalanceRunner::set_filename(const std::string& name)
{
  m_inFile = name;
}

void StkBalanceRunner::set_output_dir(const std::string& name)
{
  m_outputDirectory = name;
}

void StkBalanceRunner::set_app_type_defaults(const std::string& defaults)
{
  ThrowRequireMsg(defaults == "sm" || defaults == "sd",
      "Passed the invalid app type default specification " << defaults
      << " to StkBalanceRunner");

  m_appTypeDefaultArg = "--" + defaults;
}

void StkBalanceRunner::set_decomp_method(const std::string& method)
{
  m_decompMethodArg = "--decomp-method=" + method;
}

void StkBalanceRunner::run_end_to_end() const
{
  std::vector<const char*> args = assemble_args();
  stk::balance::LifeCycle balance(m_comm, args.size(), args.data());
  balance.run();
  EXPECT_EQ(balance.exit_code(), 0);
}

std::vector<const char*> StkBalanceRunner::assemble_args() const
{
  std::vector<const char*> args = {m_execName.c_str(), m_inFile.c_str(), m_outputDirectory.c_str()};

  if (m_appTypeDefaultArg != "") args.push_back(m_appTypeDefaultArg.c_str());
  if (m_decompMethodArg != "") args.push_back(m_decompMethodArg.c_str());

  return args;
}

} }
