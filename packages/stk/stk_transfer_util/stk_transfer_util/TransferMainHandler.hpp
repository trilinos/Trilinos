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

#ifndef STK_TRANSFER_MAIN_HANDLER_HPP
#define STK_TRANSFER_MAIN_HANDLER_HPP

#include "stk_io/FileValidator.hpp"
#include "stk_transfer_util/TransferMainSettings.hpp"
#include "stk_transfer_util/TransferMainParser.hpp"
#include <stk_util/environment/LogWithTimeAndMemory.hpp>
#include <stk_util/Version.hpp>

namespace stk {
namespace transfer_util {

enum class TransferMainStatus {
  SUCCESS                    = 0,
  PARSE_ONLY                 = 1,
  PARSE_ERROR                = 2,
  EXECUTION_ERROR            = 3
};


class TransferMainHandler {
public:
  TransferMainHandler(MPI_Comm c, int argc, const char** argv);

  void run();
  TransferMainStatus exit_code() const;

  const TransferMainSettings & get_transfer_settings() const { return m_settings; }

private:
  void parse();
  void print_parse_error(const char* what) const;

  void print_running_message() const;

  bool is_no_op() const;
  void print_no_op_message() const;

  void transfer();
  void print_transfer_error(const char* what) const;

  MPI_Comm m_comm;
  int m_argc;
  const char** m_argv;
  TransferMainStatus m_exitCode;
  const bool m_isProc0;

  const stk::io::FileValidator m_validator;
  TransferMainSettings m_settings;
  TransferMainParser m_parser;
};

} }

#endif  //STK_TRANSFER_MAIN_HANDLER_HPP
