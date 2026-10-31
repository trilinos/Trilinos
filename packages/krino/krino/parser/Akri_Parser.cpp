// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Parser.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshInput_Parser.hpp>
#include <Akri_Simulation_Parser.hpp>
#include <Akri_YAML.hpp>

#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>
#include <stk_util/parallel/mpi_filebuf.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace krino {
namespace Parser {

#ifdef KRINO_HAVE_YAML

static const std::string node_type_name(YAML::NodeType::value type)
{
  switch (type)
  {
    case YAML::NodeType::Undefined:
      return "Undefined";

    case YAML::NodeType::Null:
      return "Null";

    case YAML::NodeType::Scalar:
      return "Scalar";

    case YAML::NodeType::Sequence:
      return "Sequence";

    case YAML::NodeType::Map:
      return "Map";

    default:
      return "Unknown";
  }
}

static void emit_yaml(YAML::Emitter& emout, const YAML::Node & node) {
  // recursive depth first
  YAML::NodeType::value type = node.Type();
  std::string out;
  switch (type)
    {
    case YAML::NodeType::Scalar:
      out = node.as<std::string>();
      emout << out;
      break;
    case YAML::NodeType::Sequence:
      emout << YAML::BeginSeq;
      for (unsigned int i = 0; i < node.size(); ++i) {
        const YAML::Node & subnode = node[i];
        emit_yaml(emout, subnode);
      }
      emout << YAML::EndSeq;
      break;
    case YAML::NodeType::Map:
      emout << YAML::BeginMap ;
      for (auto && entry : node) {
        const YAML::Node & key   = entry.first;
        const YAML::Node & value = entry.second;
        out = key.as<std::string>();
        emout << YAML::Key << out;
        emout << YAML::Value;
        emit_yaml(emout, value);
      }
      emout << YAML::EndMap ;
      break;
    case YAML::NodeType::Null:
      emout << " (empty) ";
      break;
    default:
      krinolog << "Warning: emit: unknown/unsupported node type" << stk::diag::dendl;
      break;
    }
}

static void emit_to_stream(std::ostream & os, const YAML::Node & node)
{
  YAML::Emitter out;
  emit_yaml(out, node);
  os << out.c_str() << std::endl;
}

static bool is_undefined_alias_error(const YAML::ParserException & e)
{
  return e.msg == YAML::ErrorMsg::UNKNOWN_ANCHOR || e.msg == YAML::ErrorMsg::ALIAS_NOT_FOUND;
}

std::string yaml_parse_error_message(const YAML::ParserException & e)
{
  std::string message = e.what();
  if (is_undefined_alias_error(e))
    message += "\nAn unquoted value that starts with '*' is a YAML alias.  If this is a name pattern, such as a "
        "part or subdomain name with a '*' wildcard, put it in quotes, e.g. - name: '*_dragon'.";
  return message;
}

static std::vector<std::string> split_lines(const std::string & text)
{
  std::vector<std::string> lines;
  std::istringstream textStream(text);
  std::string line;
  while (std::getline(textStream, line))
    lines.push_back(line);
  return lines;
}

std::string input_excerpt(const std::string & text, const int line, const int column, const int numContextLines)
{
  const std::vector<std::string> lines = split_lines(text);
  if (lines.empty())
    return std::string();

  const int numLines = lines.size();
  const int errorLine = std::clamp(line, 0, numLines-1);
  const int beginLine = std::max(errorLine-numContextLines, 0);
  const int endLine = std::min(errorLine+numContextLines, numLines-1);
  const int lineNumberWidth = std::to_string(endLine+1).size();

  std::ostringstream excerpt;
  for (int i=beginLine; i<=endLine; ++i)
  {
    const bool isErrorLine = (i == errorLine);
    excerpt << (isErrorLine ? "> " : "  ") << std::setw(lineNumberWidth) << i+1 << " | " << lines[i] << "\n";
    if (isErrorLine && line == errorLine && column >= 0)
      excerpt << std::string(2 + lineNumberWidth + 3 + column, ' ') << "^\n";
  }
  return excerpt.str();
}

static bool write_processed_input(const std::string & fileName, const std::string & text)
{
  std::ofstream file(fileName);
  file << text;
  return bool(file);
}

static YAML::Node parse_processed_input(const std::string & inputFileName, const std::string & processedInput)
{
  try
  {
    return YAML::Load(processedInput);
  }
  catch (YAML::ParserException &e)
  {
    std::ostringstream message;
    message << yaml_parse_error_message(e) << "\n"
        << "Input after aprepro processing, near line " << e.mark.line+1
        << " (line numbers refer to this processed text):\n"
        << input_excerpt(processedInput, e.mark.line, e.mark.column);

    // Every rank parses the same processed input and fails the same way, so only the root rank writes it.
    const std::string processedInputFileName = inputFileName + ".aprepro";
    if (0 == stk::EnvData::parallel_rank())
    {
      if (write_processed_input(processedInputFileName, processedInput))
        message << "The full processed input was written to " << processedInputFileName << ".";
      else
        message << "Failed to write the full processed input to " << processedInputFileName << ".";
    }

    throw std::runtime_error(message.str());
  }
}

// Logs the input file as written, before aprepro processing, and the aprepro settings, for
// provenance.  The root rank's log is the run's log, so only the root rank reads and logs it.
static void log_input_provenance(const std::string & inputFileName, const bool useAprepro, const std::string & apreproDefines)
{
  if (0 != stk::EnvData::parallel_rank())
    return;

  sierra::Env::outputP0() << "Aprepro processing: " << (useAprepro ? "on" : "off") << std::endl;
  if (useAprepro)
    sierra::Env::outputP0() << "Aprepro defines: " << (apreproDefines.empty() ? "(none)" : apreproDefines) << std::endl;

  std::ifstream inputFile(inputFileName);
  if (!inputFile)
    return;

  const std::string unprocessedInput{std::istreambuf_iterator<char>(inputFile), std::istreambuf_iterator<char>()};
  sierra::Env::outputP0()
      << "Input file " << inputFileName << " (before aprepro processing):" << std::endl
      << sierra::Env::subsection_separator() << std::endl
      << unprocessedInput;
  if (!unprocessedInput.empty() && unprocessedInput.back() != '\n')
    sierra::Env::outputP0() << std::endl;
  sierra::Env::outputP0()
      << sierra::Env::subsection_separator() << std::endl
      << "End of input file " << inputFileName << std::endl << std::endl;
}

static Node parse_yaml_file(const std::string & inputFileName)
{
  const std::string &aprepro = sierra::Env::get_param("aprepro");
  const bool use_aprepro = aprepro.empty() || aprepro != "off";
  const std::string &aprepro_defines = sierra::Env::get_param("define");

  log_input_provenance(inputFileName, use_aprepro, aprepro_defines);

  // parallel version of input file
  mpi_filebuf input_parallel_file(use_aprepro, aprepro_defines);
  if (!input_parallel_file.open(stk::EnvData::parallel_comm(), 0, std::ios::in, inputFileName.c_str())) {
    throw std::runtime_error("failed to open input file");
  }
  if(use_aprepro) {
    sierra::Env::outputP0()<<input_parallel_file.aprepro_parse_info();
    if(input_parallel_file.aprepro_parse_warning_count()) {
      sierra::Env::outputP0()<<input_parallel_file.aprepro_parse_warnings();
    }
    if(input_parallel_file.aprepro_parse_error_count()) {
      sierra::Env::outputP0()<<input_parallel_file.aprepro_parse_errors();
    }
  }

  std::istream input_stream(&input_parallel_file);
  const std::string processedInput{std::istreambuf_iterator<char>(input_stream), std::istreambuf_iterator<char>()};

  YAML::Node doc;

  try {
    doc = parse_processed_input(inputFileName, processedInput);
  }
  catch (std::runtime_error &e) {
    krinolog << "Parser Error: " << e.what() << stk::diag::dendl;
    throw std::runtime_error("YAML parse failed");
  }


  return Node(doc);
}

bool Node::is_scalar() const { return Type() == YAML::NodeType::Scalar; }
bool Node::is_map() const { return Type() == YAML::NodeType::Map; }
bool Node::is_sequence() const { return Type() == YAML::NodeType::Sequence; }
bool Node::is_null() const { return Type() == YAML::NodeType::Null; }

bool Node::is_empty_scalar() const
{
  if (!is_scalar())
    return false;

  return mNode.as<std::string>().empty();
}

Node Node::get_if_present(const std::string& key) const
{
  std::ostringstream err_msg;

  if (!*this)
  {
    return Node();
  }

  // A Null or Scalar node cannot be indexed by key.  (An Undefined node is
  // already falsy and returned by the check above, so it need not be listed
  // here: YAML::Node::operator bool() is defined as IsDefined().)
  if (Type() == YAML::NodeType::Null ||
      Type() == YAML::NodeType::Scalar)
  {
    emit_to_stream(err_msg, mNode);
    err_msg << "Check structure of input file.";

    stk::RuntimeDoomedAdHoc()
      << "parsing within non-searchable node of type = "
      << node_type_name(Type())
      << " for key= " << key
      << " at " << line_info()
      << " node= " << err_msg.str() << std::endl;

    return Node();
  }

  YAML::Node child = mNode[key];

  if (!child.IsDefined())
  {
    return Node();
  }

  return Node(child);
}

Node Node::get_null_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Null); }
Node Node::get_map_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Map); }
Node Node::get_sequence_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Sequence); }
Node Node::get_scalar_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Scalar); }

Node Node::get_options_if_present(const std::string& key) const
{
  const Node keywordNode = get_if_present(key);

  if (!keywordNode)
  {
    return Node();
  }

  // Nested syntax:
  //
  //   plane:
  //     normal: [...]
  //     offset: ...
  if (keywordNode.is_map())
  {
    return keywordNode;
  }

  // Flat syntax:
  //
  //   plane:
  //   normal: [...]
  //   offset: ...
  //
  // Depending on yaml-cpp / Parser behavior, the marker may appear
  // as Null or as an empty scalar.
  if (keywordNode.is_null() || keywordNode.is_empty_scalar())
  {
    return *this;
  }

  stk::RuntimeDoomedAdHoc()
      << "Invalid specification for '" << key << "'. "
      << "Expected either its options nested under it (a map) or no value with its options beside it. "
      << info();

  return Node();
}

template<>
bool Node::as<bool>() const
{
  std::string resultString = as<std::string>();
  std::transform(resultString.begin(), resultString.end(), resultString.begin(), ::toupper);

  if (resultString == "TRUE" || resultString == "YES" || resultString == "1")
    return true;
  else if (resultString == "FALSE" || resultString == "NO" || resultString == "0")
    return false;

  std::ostringstream err_msg;
  emit_to_stream(err_msg, mNode);
  stk::RuntimeDoomedAdHoc() << "Failed to parse boolean "
    << " at " << line_info()
    << " node= " << err_msg.str() << std::endl;
  return false;
}

std::string Node::info() const
{
  std::ostringstream sout;
  sout << "Node at " << line_info() << " => \n" ;
  emit_to_stream(sout, mNode);
  return sout.str();
}

void Node::report_conversion_failure(const std::string& key) const
{
  std::ostringstream err_msg;
  emit_to_stream(err_msg, mNode);

  stk::RuntimeDoomedAdHoc()
    << "Failed to parse value for key= " << key
    << " at " << line_info()
    << " node= " << err_msg.str() << std::endl;
}

std::string Node::line_info() const
{
  std::ostringstream sout;
  sout << "(line,column) = ("
       << mNode.Mark().line+1 << ", "
       << mNode.Mark().column+1 << ")";
  return sout.str();
}

Node Node::get_type_if_present(
    const std::string& key,
    YAML::NodeType::value expectedType) const
{
  std::ostringstream err_msg;

  const Node value = get_if_present(key);

  if (!value)
  {
    return Node();
  }

  if (value.Type() != expectedType)
  {
    // Emit the offending child value (not the whole parent map) and report its
    // location, so the error points at the mistyped entry rather than dumping
    // the enclosing subtree and pointing at the parent's line.
    emit_to_stream(err_msg, value.mNode);
    err_msg << "Check indentation of input file.";

    stk::RuntimeDoomedAdHoc()
      << "parsing expected type " << node_type_name(expectedType)
      << " got type = " << node_type_name(value.Type())
      << " for key= " << key
      << " at " << value.line_info()
      << " node= " << err_msg.str() << std::endl;

    return Node();
  }

  return value;
}

#else
static Node parse_yaml_file(const std::string & inputFileName) { throw std::runtime_error("YAML parser not enabled in krino build."); }
#endif

void parse(Simulation & simulation)
{/* %TRACE% */ Trace trace__("krino::YAML_Parser::parse()"); /* %TRACE% */
  const std::string inputFileName = stk::EnvData::getInputFileName();
  const Node doc = parse_yaml_file(inputFileName);

  MeshInput_Parser::parse(doc);
  Simulation_Parser::parse(simulation, doc);

  unsigned local_num_doomed = stk::get_doomed_count();
  unsigned num_doomed = 0;

  MPI_Allreduce(&local_num_doomed, &num_doomed, 1, MPI_UNSIGNED, MPI_SUM, stk::EnvData::parallel_comm());

  if (num_doomed)
  {
    throw std::runtime_error("Krino parse failed");
  }
}

}
}

