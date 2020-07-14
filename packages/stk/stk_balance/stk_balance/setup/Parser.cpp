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

#include "Parser.hpp"
#include "stk_balance/search_tolerance_algs/SecondShortestEdgeFaceSearchTolerance.hpp"
#include "stk_util/environment/Env.hpp"  // for deprecatedOutputDirectory warnings

namespace stk {
namespace balance {

std::string construct_output_file_name(const std::string& outputDirectory, const std::string& inputFile)
{
  std::size_t found = inputFile.find_last_of("/");
  std::string filename = inputFile;
  if (found != std::string::npos) {
    filename = inputFile.substr(found + 1);
  }
  return outputDirectory + "/" + filename;
}

std::string Examples::get_quick_example()
{
  std::string mpiCmd = "  > mpirun -n <numProcsDecomp> " + m_execName + " ";
  std::string usage = "Usage:\n"
                    + mpiCmd + stk::angle_it(m_optionNames.infile) + " " + stk::bracket_it("output directory") + " "
                             + stk::bracket_it("args") + "\n"
                    + mpiCmd + stk::dash_it(m_optionNames.infile) + " " + stk::angle_it(m_optionNames.infile) + " "
                             + stk::dash_it(m_optionNames.outputDirectory) + " " + stk::bracket_it("output directory") + " "
                             + stk::bracket_it("args") + "\n"
                    + "\n";
  return usage;
}

std::string Examples::get_long_examples()
{
  std::string examples = "Examples:\n\n";
  std::string tab = "  ";
  examples += tab + "To decompose for 16 processors:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " file.exo\n";
  examples += "\n";
  examples += tab + "To decompose for 512 processors and put the decomposition into a directory named 'temp1':\n";
  examples += tab + tab + "> mpirun -n 512 " + m_execName + " file.exo temp1\n";
  examples += "\n";
  examples += tab + "To decompose for 64 processors and use settings suitable for solving Solid Mechanics problems:\n";
  examples += tab + tab + "> mpirun -n 64 " + m_execName + " file.exo " + stk::dash_it(m_optionNames.smDefaults) + "\n";
  examples += "\n";
  examples += tab + "To decompose for 16 processors and use the default relative contact search tolerance:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " file.exo " + stk::dash_it(m_optionNames.faceSearchRelTol) + "\n";
  examples += "\n";
  examples += tab + "To decompose for 16 processors and use a relative contact search tolerance of 0.05:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " file.exo " + stk::dash_it(m_optionNames.faceSearchRelTol) + "=0.05\n";
  examples += "\n";
  examples += tab + "To decompose for 16 processors with the RCB decomposition method:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " file.exo " + stk::dash_it(m_optionNames.decompMethod) + "=rcb\n";
  return examples;
}

Parser::Parser(MPI_Comm comm)
  : m_comm(comm),
    m_commandLineParser(m_comm)
{
  add_options_to_parser();
}

void Parser::parse_command_line_options(int argc, const char** argv, BalanceSettings& settings)
{
  setup_messages(argv);

  stk::parse_command_line(argc, argv, m_quickExample, m_longExamples, m_commandLineParser, m_comm);

  set_filenames(settings);
  set_app_type_defaults(settings);
  set_contact_search(settings);
  set_contact_search_tolerance(settings);
  set_decomp_method(settings);
}

std::string Parser::get_quick_error() const
{
  return m_quickError;
}

void Parser::add_options_to_parser()
{
  stk::CommandLineOption infile{m_optionNames.infile, "i",
                                "undecomposed serial input mesh file"};
  stk::CommandLineOption outputDirectory{m_optionNames.outputDirectory, "o",
                                         "output directory for decomposition"};
  stk::CommandLineOption deprecatedOutputDirectory{m_optionNames.deprecatedOutputDirectory, "",
                                         "DEPRECATED: output directory for decomposition"};

  std::ostringstream smStream;
  smStream << "Use settings suitable for solving Solid Mechanics problems. "
           << "This flag implies:" << std::endl
           << "    " << stk::dash_it(m_optionNames.faceSearchRelTol) << "=" << m_defaults.faceSearchRelTol << std::endl
           << "    Face search graph vertex weight multiplier = " << m_defaults.smFaceSearchVertexMultiplier << std::endl
           << "    Face search graph edge weight = " << m_defaults.smFaceSearchEdgeWeight;
  stk::CommandLineOption smDefaults{m_optionNames.smDefaults, "", smStream.str()};

  std::ostringstream sdStream;
  sdStream << "Use settings suitable for solving Structural Dynamics problems. "
           << "This flag implies:" << std::endl
           << "    " << stk::dash_it(m_optionNames.faceSearchAbsTol) << "=" << m_defaults.faceSearchAbsTol << std::endl
           << "    Face search graph vertex weight multiplier = " << m_defaults.faceSearchVertexMultiplier << std::endl
           << "    Face search graph edge weight = " << m_defaults.faceSearchEdgeWeight << std::endl
           << "    Handle spider elements";
  stk::CommandLineOption sdDefaults{m_optionNames.sdDefaults, "", sdStream.str()};

  stk::CommandLineOption faceSearchAbsTol{m_optionNames.faceSearchAbsTol, "",
                                          "Use an absolute tolerance for face contact search. "
                                          "Optionally provide a numeric tolerance value."};
  stk::CommandLineOption faceSearchRelTol{m_optionNames.faceSearchRelTol, "",
                                          "Use a tolerance relative to the face size for face contact search. "
                                          "Optionally provide a numeric tolerance value."};
  stk::CommandLineOption contactSearch{m_optionNames.contactSearch, "",
                                       "Use proximity search for contact [on|off]"};
  stk::CommandLineOption decompMethod{m_optionNames.decompMethod, "",
                                      "Use this geometric decomposition method [rcb|rib|multijagged] "
                                      "or graph-based decomposition method [parmetis].\n"
                                      "Note that geometric methods do not use contact search and "
                                      "ignore all search-related options, as well as ignoring spider elements."};

  m_commandLineParser.add_required_positional<std::string>(infile);
  m_commandLineParser.add_optional_positional<std::string>(outputDirectory, ".");
  m_commandLineParser.add_optional_positional<std::string>(deprecatedOutputDirectory, "");
  m_commandLineParser.add_flag(smDefaults);
  m_commandLineParser.add_flag(sdDefaults);
  m_commandLineParser.add_optional_implicit(faceSearchAbsTol, m_defaults.faceSearchAbsTol);
  m_commandLineParser.add_optional_implicit(faceSearchRelTol, m_defaults.faceSearchRelTol);
  m_commandLineParser.add_optional(contactSearch, m_defaults.contactSearch);
  m_commandLineParser.add_optional(decompMethod, m_defaults.decompMethod);
}

void Parser::setup_messages(const char** argv)
{
  m_execName = stk::tailname(argv[0]);
  m_examples.set_exec_name(m_execName);

  m_quickExample = m_examples.get_quick_example();
  m_longExamples = m_examples.get_long_examples();

  m_quickError = stk::get_quick_error(m_execName, m_quickExample);
}

void Parser::set_filenames(BalanceSettings& settings) const
{
  std::string outputDirectory = m_commandLineParser.get_option_value<std::string>(m_optionNames.outputDirectory);
  if (m_commandLineParser.is_option_parsed(m_optionNames.deprecatedOutputDirectory)) {
    sierra::Env::outputP0() << "Warning:  The option " << stk::dash_it(m_optionNames.deprecatedOutputDirectory)
                            << " has been deprecated in favor of " << stk::dash_it(m_optionNames.outputDirectory)
                            << " and will be removed in a future release." << std::endl;
    outputDirectory = m_commandLineParser.get_option_value<std::string>(m_optionNames.deprecatedOutputDirectory);
  }

  const std::string inputFilename = m_commandLineParser.get_option_value<std::string>(m_optionNames.infile);
  const std::string outputFilename = construct_output_file_name(outputDirectory, inputFilename);

  settings.set_input_filename(inputFilename);
  settings.set_output_filename(outputFilename);
}

void Parser::set_app_type_defaults(BalanceSettings& settings) const
{
  bool useSM = m_commandLineParser.is_option_provided(m_optionNames.smDefaults);
  bool useSD = m_commandLineParser.is_option_provided(m_optionNames.sdDefaults);

  ThrowRequireMsg( !(useSM && useSD), "Can't set default settings for multiple apps at the same time");

  if (useSM) {
    settings.setEdgeWeightForSearch(m_defaults.smFaceSearchEdgeWeight);
    settings.setVertexWeightMultiplierForVertexInSearch(m_defaults.smFaceSearchVertexMultiplier);
    settings.setToleranceFunctionForFaceSearch(
        std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>(
          m_defaults.faceSearchRelTol
          ));
  }

  if (useSD) {
    settings.setShouldFixSpiders(true);
  }
}

void Parser::set_contact_search(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.contactSearch)) {
    std::string searchOption = m_commandLineParser.get_option_value<std::string>(m_optionNames.contactSearch);
    std::transform(searchOption.begin(), searchOption.end(), searchOption.begin(), ::tolower);

    ThrowRequireMsg(searchOption == "on" || searchOption == "off",
        "Invalid contact search type (" + searchOption + ").  Must be one of: [on|off]");

    settings.setIncludeSearchResultsInGraph(searchOption == "on");
  }
}

void Parser::set_contact_search_tolerance(BalanceSettings& settings) const
{
  bool useAbsTol = m_commandLineParser.is_option_provided(m_optionNames.faceSearchAbsTol);
  bool useRelTol = m_commandLineParser.is_option_provided(m_optionNames.faceSearchRelTol);

  ThrowRequireMsg( !(useAbsTol && useRelTol), "Must not specify both an absolute and relative tolerance");

  if (useAbsTol) {
    settings.setToleranceForFaceSearch(m_commandLineParser.get_option_value<double>(m_optionNames.faceSearchAbsTol));
  }

  if (useRelTol) {
    settings.setToleranceFunctionForFaceSearch(
        std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>(
          m_commandLineParser.get_option_value<double>(m_optionNames.faceSearchRelTol)));
  }
}

void Parser::set_decomp_method(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.decompMethod)) {
    settings.setDecompMethod(m_commandLineParser.get_option_value<std::string>(m_optionNames.decompMethod));
  }
}

} }
