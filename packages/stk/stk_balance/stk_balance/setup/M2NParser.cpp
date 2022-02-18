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

#include "M2NParser.hpp"
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_util/util/string_utils.hpp"
#include "stk_balance/balanceUtils.hpp"

namespace stk {
namespace balance {

std::string M2NExamples::get_quick_example()
{
  std::string mpiCmd = "  > mpirun -n 16 " + m_execName + " ";
  std::string usage = "Usage:\n" + mpiCmd + stk::angle_it(m_optionNames.infile) + " " + stk::angle_it(m_optionNames.nprocs)
      + "\n" + mpiCmd + stk::dash_it(m_optionNames.infile) + " " + stk::angle_it(m_optionNames.infile)
      + " " + stk::dash_it(m_optionNames.nprocs) + " " + stk::angle_it(m_optionNames.nprocs)
      + "\n\n";
  return usage;
}

std::string M2NExamples::get_long_examples()
{
  std::string examples = "Examples:\n\n";
  std::string tab = "  ";
  examples += tab + "To change from 16 processor decomposition to 128 processor decomposition:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " file.exo 128\n";
  examples += "\n";
  examples += tab + "To change from 16 processor decomposition to 8 processor decomposition:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " file.exo 8\n";
  return examples;
}

M2NParser::M2NParser(MPI_Comm comm)
  : m_comm(comm),
    m_commandLineParser(m_comm)
{
  add_options_to_parser();
}

void M2NParser::parse_command_line_options(int argc, const char** argv, M2NBalanceSettings& settings)
{
  setup_messages(argv);

  stk::parse_command_line(argc, argv, m_quickExample, m_longExamples, m_commandLineParser, m_comm);

  set_filename(settings);
  set_logfile(settings);
  set_num_procs(settings);
  set_use_nested_decomp(settings);
}

std::string M2NParser::get_quick_error() const
{
  return m_quickError;
}

void M2NParser::add_options_to_parser()
{
  stk::CommandLineOption infile{m_optionNames.infile, "i", "input file decomposed for old number of processors"};
  stk::CommandLineOption nprocs{m_optionNames.nprocs, "n", "new number of processors"};
  stk::CommandLineOption logfile{m_optionNames.logfile, "l",
                                 "Output log file path, one of: 'cout', 'cerr', or a file path."};
  stk::CommandLineOption useNested{m_optionNames.useNestedDecomp, "",
        "Nest the new decomposition completely within the boundaries "
        "of the input decomposition.  The new number of processors "
        "must be an integer multiple of the input processors."};

  m_commandLineParser.add_required_positional<std::string>(infile);
  m_commandLineParser.add_required_positional<int>(nprocs);
  m_commandLineParser.add_optional<std::string>(logfile, "stk_balance_m2n.log");
  m_commandLineParser.add_flag(useNested);
}

void M2NParser::setup_messages(const char** argv)
{
  m_execName = stk::tailname(argv[0]);
  m_examples.set_exec_name(m_execName);

  m_quickExample = m_examples.get_quick_example();
  m_longExamples = m_examples.get_long_examples();

  m_quickError = stk::get_quick_error(m_execName, m_quickExample);
}

void M2NParser::set_filename(M2NBalanceSettings& settings) const
{
  settings.set_input_filename(m_commandLineParser.get_option_value<std::string>(m_optionNames.infile));
}

void M2NParser::set_num_procs(M2NBalanceSettings& settings) const
{
  const unsigned numOutputProcs = m_commandLineParser.get_option_value<unsigned>(m_optionNames.nprocs);
  ThrowRequireMsg(numOutputProcs > 0, "Please specify a valid target processor count.");
  settings.set_num_output_processors(numOutputProcs);
}

void M2NParser::set_logfile(M2NBalanceSettings& settings) const
{
  ThrowRequire(m_commandLineParser.is_option_provided(m_optionNames.logfile));
  settings.set_log_filename(m_commandLineParser.get_option_value<std::string>(m_optionNames.logfile));
}

void M2NParser::set_use_nested_decomp(M2NBalanceSettings& settings) const
{
  const bool useNestedDecomp = m_commandLineParser.is_option_provided(m_optionNames.useNestedDecomp);
  settings.set_use_nested_decomp(useNestedDecomp);
  if (useNestedDecomp) {
    const int initialNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    const int finalNumProcs = settings.get_num_output_processors();
    const bool isValidProcCount = (finalNumProcs % initialNumProcs) == 0;
    ThrowRequireMsg(isValidProcCount, "Final number of processors (" << finalNumProcs << ") must be an integer "
                    << "multiple of initial processors (" << initialNumProcs << ") to use a nested decomposition.");
  }
}

} }
