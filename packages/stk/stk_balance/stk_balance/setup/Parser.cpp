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

namespace stk {
namespace balance {

Parser::Parser(MPI_Comm comm)
  : m_comm(comm)
{ }

void Parser::parse_command_line_options(int argc, const char** argv, BalanceSettings& settings)
{
  m_execName = stk::tailname(argv[0]);
  stk::balance::CommandLineOptions cmdLineOptions;
  m_quickExample = stk::balance::get_quick_example(m_execName, cmdLineOptions.infile.name, cmdLineOptions.outputDirectory.name, m_comm);

  m_options = stk::balance::parse_balance_command_line(argc, argv, m_execName, m_comm);
  generate_settings(m_options, settings);

  stk::balance::print_running_msg(m_execName, m_options, m_comm);
}

std::string Parser::get_quick_error() const
{
  return stk::get_quick_error(m_execName, m_quickExample);
}

void Parser::generate_settings(const ParsedOptions& options, BalanceSettings& settings)
{
  const std::string& inputFilename = m_options.m_inFile;

  std::string outputFilename = construct_output_file_name(m_options.outputDirectory, inputFilename);
  settings.set_input_filename(inputFilename);
  settings.set_output_filename(outputFilename);

  SearchToleranceType searchToleranceType = ABSOLUTE;

  if (m_options.is_option_provided(stk::balance::ParsedOptions::APP_TYPE)) {
    if (m_options.appTypeDefaults == stk::balance::SD_DEFAULTS) {
      settings.setShouldFixSpiders(true);
    }
    else if (m_options.appTypeDefaults == stk::balance::SM_DEFAULTS) {
      settings.setEdgeWeightForSearch(3.0);
      settings.setVertexWeightMultiplierForVertexInSearch(10.0);
      settings.setToleranceFunctionForFaceSearch(
            std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>());
      searchToleranceType = RELATIVE;
    }
  }

  if (m_options.is_option_provided(stk::balance::ParsedOptions::CONTACT_SEARCH)) {
    settings.setIncludeSearchResultsInGraph(m_options.useContactSearch);
  }

  if (m_options.is_option_provided(stk::balance::ParsedOptions::FACE_SEARCH_ABS_TOL)) searchToleranceType = ABSOLUTE;
  if (m_options.is_option_provided(stk::balance::ParsedOptions::FACE_SEARCH_REL_TOL)) searchToleranceType = RELATIVE;

  if (searchToleranceType == ABSOLUTE) {
    const double tolerance = m_options.is_option_provided(stk::balance::ParsedOptions::FACE_SEARCH_ABS_TOL) ?
                             m_options.faceSearchAbsTol : stk::balance::defaultFaceSearchTolerance;
    settings.setToleranceForFaceSearch(tolerance);
  }
  else if (searchToleranceType == RELATIVE) {
    if (m_options.is_option_provided(stk::balance::ParsedOptions::FACE_SEARCH_REL_TOL)) {
      settings.setToleranceFunctionForFaceSearch(
            std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>(m_options.faceSearchRelTol));
    }
    else {
      settings.setToleranceFunctionForFaceSearch(
            std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>());
    }
  }

  if (m_options.is_option_provided(stk::balance::ParsedOptions::DECOMP_METHOD)) {
    settings.setDecompMethod(m_options.decompMethod);
  }
}

} }
