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

#ifndef STK_BALANCE_LIFE_CYCLE_HPP
#define STK_BALANCE_LIFE_CYCLE_HPP

#include "stk_balance/balanceUtils.hpp"
#include "stk_balance/setup/Parser.hpp"
#include "stk_balance/setup/FileValidator.hpp"

namespace stk {
namespace balance {

class LifeCycle {
public:
  LifeCycle(MPI_Comm c, int argc, const char** argv);

  void run();
  int exit_code() const;

private:
  void parse();
  bool serial_no_op() const;
  void balance();

  void set_output_streams();
  void print_parse_error(const char* what) const;
  void print_balance_error(const char* what) const;
  void print_no_op_message() const;
  void print_running_message() const;

  MPI_Comm m_comm;
  const int m_argc;
  const char** m_argv;
  int m_exitCode;
  const bool m_isProc0;

  const FileValidator m_validator;
  StkBalanceSettings m_settings;
  Parser m_parser;
};

} }

#endif
