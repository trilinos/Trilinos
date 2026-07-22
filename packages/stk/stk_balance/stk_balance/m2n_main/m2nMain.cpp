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

#include <stk_util/parallel/Parallel.hpp>
#include <stk_balance/balanceUtils.hpp>
#include "stk_balance/internal/LogUtils.hpp"

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <sstream>

static std::string deprecationMsg = "stk_balance_m2n is deprecated. Please use stk_balance instead.";
static std::string infileOpt = "infile";
static std::string logfileOpt = "logfile";
static std::string nprocsOpt = "nprocs";
static std::string useNestedDecompOpt = "use-nested-decomp";

void parse_args(int argc, const char**argv, stk::CommandLineParserParallel& commandLineParser)
{
  stk::CommandLineOption infile{infileOpt, "i", "input file decomposed for old number of processors"};
  stk::CommandLineOption nprocs{nprocsOpt, "n", "new number of processors"};
  stk::CommandLineOption logfile{logfileOpt, "l", "Output log file path, one of: 'cout', 'cerr', or a file path."};
  stk::CommandLineOption useNested{useNestedDecompOpt, "",
        "Nest the new decomposition completely within the boundaries "
        "of the input decomposition.  The new number of processors "
        "must be an integer multiple of the input processors."};

  commandLineParser.add_required_positional<std::string>(infile);
  commandLineParser.add_required_positional<int>(nprocs);
  commandLineParser.add_optional<std::string>(logfile, "");
  commandLineParser.add_flag(useNested);

  stk::parse_command_line(argc, argv, deprecationMsg, deprecationMsg, commandLineParser, stk::parallel_machine_world());
}

std::string get_deprecation_msg(int argc, const char**argv)
{
  int commSize;
  MPI_Comm_size(stk::parallel_machine_world(), &commSize);

  stk::CommandLineParserParallel commandLineParser(stk::parallel_machine_world());
  parse_args(argc, argv, commandLineParser);

  std::ostringstream mpirunPrefix;
  if (commSize > 1) {
    mpirunPrefix << "mpirun -n " << commSize << " ";
  }

  std::ostringstream originalCmd;
  originalCmd << "\t> " << mpirunPrefix.str();
  originalCmd << stk::tailname(argv[0]) << " ";
  for (int i = 1; i < argc; ++i) {
    originalCmd << argv[i] << " ";
  }

  std::ostringstream convertToString;
  convertToString << "\t> " << mpirunPrefix.str();
  convertToString << "stk_balance ";

  convertToString << "--rebalance-to " << commandLineParser.get_option_value<unsigned>(nprocsOpt) << " ";
  if (commandLineParser.is_option_parsed(logfileOpt)) {
    convertToString << "--logfile " << commandLineParser.get_option_value<std::string>(logfileOpt) << " ";
  }
  if (commandLineParser.is_option_parsed(useNestedDecompOpt)) {
    convertToString << "--use-nested-decomp" << " ";
  }
  convertToString << commandLineParser.get_option_value<std::string>(infileOpt) << " ";

  if (commSize == 1) {
    convertToString << "\n\t\tOR\n"; 
    convertToString << "\t> mpirun -n " << commandLineParser.get_option_value<unsigned>(nprocsOpt) << " ";
    convertToString << "stk_balance ";
    if (commandLineParser.is_option_parsed(logfileOpt)) {
      convertToString << "--logfile " << commandLineParser.get_option_value<std::string>(logfileOpt) << " ";
    }
    if (commandLineParser.is_option_parsed(useNestedDecompOpt)) {
      convertToString << "--use-nested-decomp" << " ";
    }
    convertToString << commandLineParser.get_option_value<std::string>(infileOpt) << " ";
  }

  std::ostringstream outputMsg;
  outputMsg << deprecationMsg << "\n\n\n";
  outputMsg << "  Input stk_balance_m2n command:\n";
  outputMsg << originalCmd.str() << "\n\n";
  outputMsg << "  can be converted to a stk_balance command:\n";
  outputMsg << convertToString.str() << "\n";

  return outputMsg.str();
}

int main(int argc, const char**argv)
{
  MPI_Init(&argc, const_cast<char***>(&argv));

  int rank;
  MPI_Comm_rank(stk::parallel_machine_world(), &rank);

  if (rank == 0) {
    std::cout << get_deprecation_msg(argc, argv) << std::endl;
  }

  MPI_Finalize();
  std::abort();
}
