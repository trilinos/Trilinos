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

#ifndef STK_BALANCE_PARSER_HPP
#define STK_BALANCE_PARSER_HPP

#include "mpi.h"

#include "stk_util/command_line/CommandLineParserParallel.hpp"
#include "stk_balance/balanceUtils.hpp"
#include "stk_balance/setup/DefaultSettings.hpp"

namespace stk {
namespace balance {


struct OptionNames
{
  const std::string infile = "infile";
  const std::string outputDirectory = "output-directory";
  const std::string rebalanceTo = "rebalance-to";
  const std::string logfile = "logfile";
  const std::string printDiagnostics = "print-diagnostics";
  const std::string smDefaults = "sm";
  const std::string sdDefaults = "sd";
  const std::string faceSearchAbsTol = "face-search-abs-tol";
  const std::string faceSearchRelTol = "face-search-rel-tol";
  const std::string contactSearch = "contact-search";
  const std::string fixSpiders = "fix-spiders";
  const std::string fixMechanisms = "fix-mechanisms";
  const std::string decompMethod = "decomp-method";
  const std::string vertexWeightBlockMultiplier = "block-weights";
  const std::string cohesiveElements = "cohesive-elements";
  const std::string useNestedDecomp = "use-nested-decomp";

  const std::string vertexWeightMethod = "vertex-weight-method";
  const std::string contactSearchEdgeWeight = "contact-search-edge-weight";
  const std::string contactSearchVertexWeightMultiplier = "contact-search-vertex-weight-mult";
  const std::string edgeWeightMultiplier = "edge-weight-mult";
  const std::string vertexWeightFieldName = "vertex-weight-field-name";
};

class Examples
{
public:
  Examples() {}

  void set_exec_name(const std::string& name) { m_execName = name; }

  std::string get_quick_example();
  std::string get_long_examples();

private:
  std::string m_execName;
  const OptionNames m_optionNames;
};

class Parser {

public:
  Parser(MPI_Comm comm);

  void parse_command_line_options(int argc, const char** argv, BalanceSettings& settings);
  std::string get_quick_error() const;

private:
  void add_options_to_parser();
  void setup_messages(const char** argv);

  void set_filenames(BalanceSettings& settings) const;
  void set_logfile(BalanceSettings& settings) const;
  void set_processors(BalanceSettings& settings) const;
  void set_use_nested_decomp(BalanceSettings& settings) const;
  void set_app_type_defaults(BalanceSettings& settings) const;
  void set_contact_search(BalanceSettings& settings) const;
  void set_fix_spiders(BalanceSettings& settings) const;
  void set_fix_mechanisms(BalanceSettings& settings) const;
  void set_contact_search_tolerance(BalanceSettings& settings) const;
  void set_decomp_method(BalanceSettings& settings) const;
  void set_vertex_weight_block_multiplier(BalanceSettings& settings) const;
  void set_cohesive_elements(BalanceSettings& settings) const;
  void set_print_diagnostics(BalanceSettings& settings) const;

  void set_vertex_weight_method(BalanceSettings& settings) const;
  void set_vertex_weight_field_name(BalanceSettings& settings) const;
  void set_contact_search_edge_weight(BalanceSettings& settings) const;
  void set_contact_search_vertex_weight_multiplier(BalanceSettings& settings) const;
  void set_edge_weight_multiplier(BalanceSettings& settings) const;

  const MPI_Comm m_comm;
  const OptionNames m_optionNames;
  stk::CommandLineParserParallel m_commandLineParser;
  Examples m_examples;

  std::string m_execName;
  std::string m_quickExample;
  std::string m_longExamples;
  std::string m_quickError;
};

} }

#endif
