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

#ifndef M2NPARSER_HPP
#define M2NPARSER_HPP

#include "stk_util/command_line/CommandLineParserParallel.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace stk {
namespace balance {

class M2NBalanceSettings;

struct M2NOptionNames
{
  const std::string infile = "infile";
  const std::string logfile = "logfile";
  const std::string nprocs = "nprocs";
  const std::string useNestedDecomp = "use-nested-decomp";
};

class M2NExamples
{
public:
  M2NExamples() {}

  void set_exec_name(const std::string& name) { m_execName = name; }

  std::string get_quick_example();
  std::string get_long_examples();

private:
  std::string m_execName;
  const M2NOptionNames m_optionNames;
};

class M2NParser {

public:
  M2NParser(MPI_Comm comm);

  void parse_command_line_options(int argc, const char** argv, M2NBalanceSettings& settings);
  std::string get_quick_error() const;

private:
  void add_options_to_parser();
  void setup_messages(const char** argv);

  void set_filename(M2NBalanceSettings& settings) const;
  void set_logfile(M2NBalanceSettings& settings) const;
  void set_num_procs(M2NBalanceSettings& settings) const;
  void set_use_nested_decomp(M2NBalanceSettings& settings) const;

  const MPI_Comm m_comm;
  const M2NOptionNames m_optionNames;
  stk::CommandLineParserParallel m_commandLineParser;
  M2NExamples m_examples;

  std::string m_execName;
  std::string m_quickExample;
  std::string m_longExamples;
  std::string m_quickError;
};

} }

#endif // M2NPARSER_HPP
