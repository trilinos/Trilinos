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
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_util/util/string_utils.hpp"
#include "stk_io/FileValidator.hpp"

namespace stk {
namespace balance {

std::string Examples::get_quick_example()
{
  std::string mpiCmd = "  > mpirun -n ";
  std::string usage = "Usage:\n"
                    + mpiCmd + stk::angle_it("outProcs") + " " + m_execName + " " + stk::angle_it(m_optionNames.infile) + " "
                             + stk::bracket_it("output directory") + " " + stk::bracket_it("args") + "\n"
                    + mpiCmd + stk::angle_it("outProcs") + " " + m_execName + " " + stk::dash_it(m_optionNames.infile) + " "
                             + stk::angle_it(m_optionNames.infile) + " " + stk::dash_it(m_optionNames.outputDirectory)
                             + " " + stk::bracket_it("output directory") + " " + stk::bracket_it("args") + "\n"
                    + mpiCmd + stk::angle_it("inProcs") + " " + m_execName + " " + stk::angle_it(m_optionNames.infile) + " "
                             + stk::bracket_it("output directory") + " --rebalance-to=<outProcs> " + stk::bracket_it("args") + "\n"
                    + "\n";
  return usage;
}

std::string Examples::get_long_examples()
{
  std::string examples = "Examples:\n\n";
  std::string tab = "  ";
  examples += tab + "To decompose for 16 processors:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " mesh.g\n";
  examples += "\n";
  examples += tab + "To decompose for 512 processors and put the decomposition into a directory named 'temp1':\n";
  examples += tab + tab + "> mpirun -n 512 " + m_execName + " mesh.g temp1\n";
  examples += "\n";
  examples += tab + "To decompose for 16 processors and use a relative contact search tolerance of 0.05:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " mesh.g " + stk::dash_it(m_optionNames.faceSearchRelTol) + "=0.05\n";
  examples += "\n";
  examples += tab + "To decompose for 16 processors with the RIB decomposition method:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " mesh.g " + stk::dash_it(m_optionNames.decompMethod) + "=rib\n";
  examples += "\n";
  examples += tab + "To rebalance a 16 processor mesh into 64 processors:\n";
  examples += tab + tab + "> mpirun -n 16 " + m_execName + " mesh.g " + stk::dash_it(m_optionNames.rebalanceTo) + "=64\n";
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

  stk::ErrorHandler orig = stk::set_assert_handler(stk::clean_error_handler);

  stk::parse_command_line(argc, argv, m_quickExample, m_longExamples, m_commandLineParser, m_comm);

  stk::set_assert_handler(orig);

  set_filenames(settings);
  set_processors(settings);
  set_app_type_defaults(settings);
  set_contact_search(settings);
  set_contact_search_tolerance(settings);
  set_fix_spiders(settings);
  set_fix_mechanisms(settings);
  set_decomp_method(settings);
  set_vertex_weight_block_multiplier(settings);
  set_cohesive_elements(settings);
  set_vertex_weight_method(settings);
  set_vertex_weight_field_name(settings);
  set_contact_search_edge_weight(settings);
  set_contact_search_vertex_weight_multiplier(settings);
  set_edge_weight_multiplier(settings);
  set_use_nested_decomp(settings);
  set_print_diagnostics(settings);
  set_logfile(settings);
}

std::string Parser::get_quick_error() const
{
  return m_quickError;
}

void Parser::add_options_to_parser()
{
  stk::CommandLineOption infile{m_optionNames.infile, "i",
                           "For balance: Undecomposed serial input mesh file\n"
                           "For rebalance: Mesh file decomposed to number MPI ranks"};
//                           "For rebalance: Mesh file that has been decomposed to the number of MPI ranks"};
  stk::CommandLineOption outputDirectory{m_optionNames.outputDirectory, "o",
                           "Output directory for decomposition"};
  stk::CommandLineOption logfile{m_optionNames.logfile, "l",
                           "Output log file path, one of: 'cout', 'cerr', or a file path."};
  stk::CommandLineOption printDiagnostics{m_optionNames.printDiagnostics, "d",
                           "Gather and print diagnostic output for each processor to help "
                           "evaluate the quality of the decomposition."};
  stk::CommandLineOption rebalanceTo{m_optionNames.rebalanceTo, "r",
                           "Rebalance the input mesh that is already decomposed for <inProcs> "
                           "processors into this number of processors.  The decomposition size "
                           "may be either increased or decreased."};

  std::ostringstream smStream;
  smStream << "Use settings suitable for solving Solid Mechanics problems. "
           << "This flag implies:" << std::endl
           << "    " << stk::dash_it(m_optionNames.faceSearchRelTol) << "=" << DefaultSettings::faceSearchRelTol << std::endl
           << "    " << stk::dash_it(m_optionNames.fixSpiders) << "=" << ((DefaultSettings::smFixSpiders) ? "on" : "off") << std::endl
           << "    " << stk::dash_it(m_optionNames.fixMechanisms) << "=" << ((DefaultSettings::fixMechanisms) ? "on" : "off") << std::endl
           << "    " << stk::dash_it(m_optionNames.vertexWeightMethod) << "=" << vertex_weight_method_name(DefaultSettings::smVertexWeightMethod) << std::endl
           << "    " << stk::dash_it(m_optionNames.edgeWeightMultiplier) << "=" << DefaultSettings::smGraphEdgeWeightMultiplier << std::endl
           << "    " << stk::dash_it(m_optionNames.contactSearchVertexWeightMultiplier) << "=" << DefaultSettings::smFaceSearchVertexMultiplier << std::endl
           << "    " << stk::dash_it(m_optionNames.contactSearchEdgeWeight) << "=" << DefaultSettings::smFaceSearchEdgeWeight << std::endl;
  stk::CommandLineOption smDefaults{m_optionNames.smDefaults, "", smStream.str()};

  std::ostringstream sdStream;
  sdStream << "Use settings suitable for solving Structural Dynamics problems. "
           << "This flag implies:" << std::endl
           << "    " << stk::dash_it(m_optionNames.faceSearchRelTol) << "=" << DefaultSettings::faceSearchRelTol << std::endl
           << "    " << stk::dash_it(m_optionNames.fixSpiders) << "=" << ((DefaultSettings::sdFixSpiders) ? "on" : "off") << std::endl
           << "    " << stk::dash_it(m_optionNames.fixMechanisms) << "=" << ((DefaultSettings::fixMechanisms) ? "on" : "off") << std::endl
           << "    " << stk::dash_it(m_optionNames.vertexWeightMethod) << "=" << vertex_weight_method_name(DefaultSettings::sdVertexWeightMethod) << std::endl
           << "    " << stk::dash_it(m_optionNames.edgeWeightMultiplier) << "=" << DefaultSettings::sdGraphEdgeWeightMultiplier << std::endl
           << "    " << stk::dash_it(m_optionNames.contactSearchVertexWeightMultiplier) << "=" << DefaultSettings::sdFaceSearchVertexMultiplier << std::endl
           << "    " << stk::dash_it(m_optionNames.contactSearchEdgeWeight) << "=" << DefaultSettings::sdFaceSearchEdgeWeight << std::endl;
  stk::CommandLineOption sdDefaults{m_optionNames.sdDefaults, "", sdStream.str()};

  stk::CommandLineOption faceSearchAbsTol{m_optionNames.faceSearchAbsTol, "",
                           "Use an absolute tolerance for face contact search. "
                           "Optionally provide a numeric tolerance value."};
  stk::CommandLineOption faceSearchRelTol{m_optionNames.faceSearchRelTol, "",
                           "Use a tolerance relative to the face size for face contact search. "
                           "Optionally provide a numeric tolerance value.  This is the global "
                           "default.  Values less than 0.5 are recommended."};
  stk::CommandLineOption contactSearch{m_optionNames.contactSearch, "",
                           "Use proximity search for contact [on|off]"};
  stk::CommandLineOption fixSpiders{m_optionNames.fixSpiders, "",
                           "Correct the decomp to group spider legs (large collection of beam "
                           "elements connected to a single node) onto fewer processors [on|off]"};
  stk::CommandLineOption fixMechanisms{m_optionNames.fixMechanisms, "",
                           "Remove mechanisms (partition components connected by a hinge) in "
                           "the decomp by reassigning element ownership [on|off]"};
  stk::CommandLineOption decompMethod{m_optionNames.decompMethod, "m",
                           "Use this geometric decomposition method [rcb|rib|multijagged] "
                           "or graph-based decomposition method [parmetis|scotch]. "
                           "Note that geometric methods do not use contact search and "
                           "ignore all search-related options, as well as ignoring spider elements."};
  stk::CommandLineOption vertexWeightBlockMultiplier{m_optionNames.vertexWeightBlockMultiplier, "",
                           "Specify a list of vertex weight multipliers through "
                           "comma-separated block_name:weight pairs to use for each element in "
                           "the block.  A multiplier of 1.0 is used for each unspecified block.\n"
                           " Syntax example: block_1:1.5,block_2:3"};
  stk::CommandLineOption cohesiveElements{m_optionNames.cohesiveElements, "",
                           "Specify a comma-separated list of blocks containing cohesive elements. "
                           "The elements in these blocks will not appear on a partition boundary. "
                           "This option requires parmetis as the decomp method\n"
                           " Syntax example: block_1,block_2"};
  stk::CommandLineOption useNested{m_optionNames.useNestedDecomp, "",
                           "If doing a rebalance, nest the new decomposition completely "
                           "within the boundaries of the input decomposition.  The new number "
                           "of processors must be an integer multiple of the input processors."};

  stk::CommandLineOption vertexWeightMethod{m_optionNames.vertexWeightMethod, "",
                           "Method used to calculate vertex weights given to the partitioner. "
                           "[constant|topology|connectivity|field]"};
  stk::CommandLineOption vertexWeightFieldName{m_optionNames.vertexWeightFieldName, "",
                           "Field name for field vertex weights given to the partitioner. "
                           "Using this option implies --vertex-weight-method=field."};
  stk::CommandLineOption contactSearchEdgeWeight{m_optionNames.contactSearchEdgeWeight, "",
                           "Graph edge weight to use between elements that are determined to be "
                           "in contact."};
  stk::CommandLineOption contactSearchVertexWeightMultiplier{m_optionNames.contactSearchVertexWeightMultiplier, "",
                           "Scale factor to be applied to graph vertex weights for elements that "
                           "are determined to be in contact."};
  stk::CommandLineOption edgeWeightMultiplier{m_optionNames.edgeWeightMultiplier, "",
                           "Scale factor to be applied to all graph edge weights.  This will be "
                           "automatically set to 1.0 for constant vertex weights, 1.0 for topology "
                           "vertex weights, and 10.0 for connectivity vertex weights."};


  m_commandLineParser.add_required_positional<std::string>(infile);
  m_commandLineParser.add_optional_positional<std::string>(outputDirectory, DefaultSettings::outputDirectory);
  m_commandLineParser.add_optional<std::string>(logfile, "<mesh>.<inProcs>_to_<outProcs>.log");
  m_commandLineParser.add_flag(printDiagnostics);
  m_commandLineParser.add_optional<unsigned>(rebalanceTo);
  m_commandLineParser.add_flag(smDefaults);
  m_commandLineParser.add_flag(sdDefaults);
  m_commandLineParser.add_optional_implicit(faceSearchAbsTol, DefaultSettings::faceSearchAbsTol);
  m_commandLineParser.add_optional_implicit(faceSearchRelTol, DefaultSettings::faceSearchRelTol);
  m_commandLineParser.add_optional(contactSearch, (DefaultSettings::useContactSearch) ? "on" : "off");
  m_commandLineParser.add_optional(fixSpiders, (DefaultSettings::fixSpiders) ? "on" : "off");
  m_commandLineParser.add_optional(fixMechanisms, (DefaultSettings::fixMechanisms) ? "on" : "off");
  m_commandLineParser.add_optional(decompMethod, DefaultSettings::decompMethod);
  m_commandLineParser.add_optional(vertexWeightBlockMultiplier, DefaultSettings::vertexWeightBlockMultiplier);
  m_commandLineParser.add_optional(cohesiveElements, DefaultSettings::cohesiveElements);
  m_commandLineParser.add_flag(useNested);

  m_commandLineParser.add_optional(vertexWeightMethod, vertex_weight_method_name(DefaultSettings::vertexWeightMethod));
  m_commandLineParser.add_optional<std::string>(vertexWeightFieldName);
  m_commandLineParser.add_optional<double>(contactSearchEdgeWeight, DefaultSettings::faceSearchEdgeWeight);
  m_commandLineParser.add_optional<double>(contactSearchVertexWeightMultiplier, DefaultSettings::faceSearchVertexMultiplier);
  m_commandLineParser.add_optional<double>(edgeWeightMultiplier, DefaultSettings::graphEdgeWeightMultiplier);

  m_commandLineParser.disallow_unrecognized();
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
  const std::string inputFilename = m_commandLineParser.get_option_value<std::string>(m_optionNames.infile);
  const std::string outputFilename = stk::io::construct_output_file_name(outputDirectory, inputFilename);

  settings.set_input_filename(inputFilename);
  settings.set_output_filename(outputFilename);
}

void Parser::set_processors(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.rebalanceTo)) {
    const unsigned numOutputProcs = m_commandLineParser.get_option_value<unsigned>(m_optionNames.rebalanceTo);
    STK_ThrowRequireMsg(numOutputProcs > 0, "Must specify a positive output processor count");
    settings.set_num_input_processors(stk::parallel_machine_size(m_comm));
    settings.set_num_output_processors(numOutputProcs);
    settings.set_is_rebalancing(true);
  }
  else {
    settings.set_num_input_processors(1);
    settings.set_num_output_processors(stk::parallel_machine_size(m_comm));
    settings.set_is_rebalancing(false);
  }
}

void Parser::set_logfile(BalanceSettings& settings) const
{

  if (m_commandLineParser.is_option_parsed(m_optionNames.logfile)) {
    settings.set_log_filename(m_commandLineParser.get_option_value<std::string>(m_optionNames.logfile));
  }
  else {
    const int inputRanks = settings.get_is_rebalancing() ? stk::parallel_machine_size(m_comm) : 1;
    const int outputRanks = settings.get_is_rebalancing() ? settings.get_num_output_processors()
                                                          : stk::parallel_machine_size(m_comm);
    settings.set_log_filename(stk::basename(stk::tailname(settings.get_input_filename())) + "." + std::to_string(inputRanks) +
                              "_to_" + std::to_string(outputRanks) + ".log");
  }
}

void Parser::set_app_type_defaults(BalanceSettings& settings) const
{
  bool useSM = m_commandLineParser.is_option_provided(m_optionNames.smDefaults);
  bool useSD = m_commandLineParser.is_option_provided(m_optionNames.sdDefaults);

  STK_ThrowRequireMsg( !(useSM && useSD), "Can't set default settings for multiple apps at the same time");

  if (useSM) {
    settings.setVertexWeightMethod(DefaultSettings::smVertexWeightMethod);
    settings.setGraphEdgeWeightMultiplier(DefaultSettings::smGraphEdgeWeightMultiplier);
    settings.setVertexWeightMultiplierForVertexInSearch(DefaultSettings::smFaceSearchVertexMultiplier);
    settings.setEdgeWeightForSearch(DefaultSettings::smFaceSearchEdgeWeight);
    settings.setShouldFixSpiders(DefaultSettings::smFixSpiders);
  }

  if (useSD) {
    settings.setVertexWeightMethod(DefaultSettings::sdVertexWeightMethod);
    settings.setGraphEdgeWeightMultiplier(DefaultSettings::sdGraphEdgeWeightMultiplier);
    settings.setVertexWeightMultiplierForVertexInSearch(DefaultSettings::sdFaceSearchVertexMultiplier);
    settings.setEdgeWeightForSearch(DefaultSettings::sdFaceSearchEdgeWeight);
    settings.setShouldFixSpiders(DefaultSettings::sdFixSpiders);
  }
}

void Parser::set_contact_search(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.contactSearch)) {
    std::string searchOption = m_commandLineParser.get_option_value<std::string>(m_optionNames.contactSearch);
    std::transform(searchOption.begin(), searchOption.end(), searchOption.begin(), ::tolower);

    STK_ThrowRequireMsg(searchOption == "on" || searchOption == "off",
        "Invalid contact search type (" + searchOption + ").  Must be one of: [on|off]");

    settings.setIncludeSearchResultsInGraph(searchOption == "on");
  }
}

void Parser::set_fix_spiders(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.fixSpiders)) {
    std::string fixSpiders = m_commandLineParser.get_option_value<std::string>(m_optionNames.fixSpiders);
    std::transform(fixSpiders.begin(), fixSpiders.end(), fixSpiders.begin(), ::tolower);

    STK_ThrowRequireMsg(fixSpiders == "on" || fixSpiders == "off",
        "Invalid spider fixing argument (" + fixSpiders + ").  Must be one of: [on|off]");

    settings.setShouldFixSpiders(fixSpiders == "on");
  }
}

void Parser::set_fix_mechanisms(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.fixMechanisms)) {
    std::string fixMechanisms = m_commandLineParser.get_option_value<std::string>(m_optionNames.fixMechanisms);
    std::transform(fixMechanisms.begin(), fixMechanisms.end(), fixMechanisms.begin(), ::tolower);

    STK_ThrowRequireMsg(fixMechanisms == "on" || fixMechanisms == "off",
        "Invalid mechanism fixing argument (" + fixMechanisms + ").  Must be one of: [on|off]");

    settings.setShouldFixMechanisms(fixMechanisms == "on");
  }
}

void Parser::set_contact_search_tolerance(BalanceSettings& settings) const
{
  bool useAbsTol = m_commandLineParser.is_option_provided(m_optionNames.faceSearchAbsTol);
  bool useRelTol = m_commandLineParser.is_option_provided(m_optionNames.faceSearchRelTol);

  STK_ThrowRequireMsg( !(useAbsTol && useRelTol), "Must not specify both an absolute and relative tolerance");

  if (useAbsTol) {
    settings.setToleranceForFaceSearch(m_commandLineParser.get_option_value<double>(m_optionNames.faceSearchAbsTol));
  }

  if (useRelTol) {
    settings.setToleranceFunctionForFaceSearch(
        std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>(
          m_commandLineParser.get_option_value<double>(m_optionNames.faceSearchRelTol)
        )
    );
  }
}

void Parser::set_decomp_method(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.decompMethod)) {
    settings.setDecompMethod(m_commandLineParser.get_option_value<std::string>(m_optionNames.decompMethod));
  }
}

void Parser::set_vertex_weight_block_multiplier(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.vertexWeightBlockMultiplier)) {
    const std::string blockMultiplierString =
        m_commandLineParser.get_option_value<std::string>(m_optionNames.vertexWeightBlockMultiplier);
    std::vector<std::string> blockSegments = stk::split_csv_string(blockMultiplierString);
    for (const std::string & blockSegment : blockSegments) {
      std::vector<std::string> multiplierSegments = stk::split_string(stk::trim_string(blockSegment), ':');
      STK_ThrowRequireMsg(multiplierSegments.size() == 2,
                      "Require block_name:value pairs for vertex weight block multiplier (" <<
                      stk::trim_string(blockSegment) << ")");
      const std::string blockName = stk::trim_string(multiplierSegments[0]);
      const double multiplier = std::stod(stk::trim_string(multiplierSegments[1]));
      settings.setVertexWeightBlockMultiplier(blockName, multiplier);
    }
  }
}

void Parser::set_cohesive_elements(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.cohesiveElements)) {
    const std::string blockString =
        m_commandLineParser.get_option_value<std::string>(m_optionNames.cohesiveElements);
    std::vector<std::string> blocks = stk::split_csv_string(blockString);
    for (const std::string & block : blocks) {
      const std::string blockName = stk::trim_string(block);
      settings.setCohesiveElements(blockName);
    }

    if (m_commandLineParser.is_option_parsed(m_optionNames.decompMethod)) {
      auto providedDecompMethod = m_commandLineParser.get_option_value<std::string>(m_optionNames.decompMethod);
      STK_ThrowRequireMsg(providedDecompMethod == "parmetis", "Cohesive elements can only be used with parmetis decomp method. Method " << providedDecompMethod << " was supplied instead.");
    }
  }
}

void Parser::set_use_nested_decomp(BalanceSettings& settings) const
{
  const bool useNestedDecomp = m_commandLineParser.is_option_provided(m_optionNames.useNestedDecomp);
  settings.set_use_nested_decomp(useNestedDecomp);

  if (useNestedDecomp) {
    const int inputNumProcs = stk::parallel_machine_size(m_comm);
    const int outputNumProcs = settings.get_num_output_processors();
    const bool isValidProcCount = (outputNumProcs % inputNumProcs) == 0;
    STK_ThrowRequireMsg(isValidProcCount, "Output number of processors (" << outputNumProcs << ") must be an integer "
                    << "multiple of input processors (" << inputNumProcs << ") to use a nested decomposition.");
   }
 }

void Parser::set_print_diagnostics(BalanceSettings &settings) const
{
  const bool printDiagnostics = m_commandLineParser.is_option_provided(m_optionNames.printDiagnostics);

  if (printDiagnostics) {
    settings.setShouldPrintDiagnostics(true);
  }
}

void Parser::set_vertex_weight_method(BalanceSettings &settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.vertexWeightMethod)) {
    const std::string vertexWeightMethodName = m_commandLineParser.get_option_value<std::string>(m_optionNames.vertexWeightMethod);
    // FIXME: case-insensitive comparison?  Need this for decomp method too?
    if (vertexWeightMethodName == vertex_weight_method_name(VertexWeightMethod::CONSTANT)) {
      settings.setVertexWeightMethod(VertexWeightMethod::CONSTANT);
      settings.setGraphEdgeWeightMultiplier(1.0);
    }
    else if (vertexWeightMethodName == vertex_weight_method_name(VertexWeightMethod::TOPOLOGY)) {
      settings.setVertexWeightMethod(VertexWeightMethod::TOPOLOGY);
      settings.setGraphEdgeWeightMultiplier(1.0);
    }
    else if (vertexWeightMethodName == vertex_weight_method_name(VertexWeightMethod::CONNECTIVITY)) {
      settings.setVertexWeightMethod(VertexWeightMethod::CONNECTIVITY);
      settings.setGraphEdgeWeightMultiplier(10.0);
    }
    else if (vertexWeightMethodName == vertex_weight_method_name(VertexWeightMethod::FIELD)) {
      settings.setVertexWeightMethod(VertexWeightMethod::FIELD);
    }
    else {
      STK_ThrowErrorMsg("Unrecognized vertex weight method: " << vertexWeightMethodName);
    }
  }
}

void Parser::set_vertex_weight_field_name(BalanceSettings& settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.vertexWeightFieldName)) {
    settings.setVertexWeightMethod(VertexWeightMethod::FIELD);
    settings.setGraphEdgeWeightMultiplier(1.0);
    const std::string vertexWeightFieldName = m_commandLineParser.get_option_value<std::string>(m_optionNames.vertexWeightFieldName);
    settings.setVertexWeightFieldName(vertexWeightFieldName);
  }
}

void Parser::set_contact_search_edge_weight(BalanceSettings &settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.contactSearchEdgeWeight)) {
    const double edgeWeight = m_commandLineParser.get_option_value<double>(m_optionNames.contactSearchEdgeWeight);
    settings.setEdgeWeightForSearch(edgeWeight);
  }
}

void Parser::set_contact_search_vertex_weight_multiplier(BalanceSettings &settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.contactSearchVertexWeightMultiplier)) {
    const double multiplier = m_commandLineParser.get_option_value<double>(m_optionNames.contactSearchVertexWeightMultiplier);
    settings.setVertexWeightMultiplierForVertexInSearch(multiplier);
  }
}

void Parser::set_edge_weight_multiplier(BalanceSettings &settings) const
{
  if (m_commandLineParser.is_option_parsed(m_optionNames.edgeWeightMultiplier)) {
    const double edgeWeightMultiplier = m_commandLineParser.get_option_value<double>(m_optionNames.edgeWeightMultiplier);
    settings.setGraphEdgeWeightMultiplier(edgeWeightMultiplier);
  }
}

}  // namespace balance
}  // namespace stk
