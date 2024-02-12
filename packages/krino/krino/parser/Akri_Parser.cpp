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

namespace krino {
namespace Parser {

#ifdef KRINO_HAVE_YAML

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
      std::cerr << "Warning: emit: unknown/unsupported node type" << std::endl;
      break;
    }
}

static void emit_to_stream(std::ostream & os, const YAML::Node & node)
{
  YAML::Emitter out;
  emit_yaml(out, node);
  os << out.c_str() << std::endl;
}

static YAML::Node parse_yaml_stream(std::istream & input_stream)
{
  YAML::Node doc;

  try {
    doc = YAML::Load(input_stream);
  }
  catch (YAML::ParserException &e) {
    throw std::runtime_error(e.what());
  }

  return doc;
}

static Node parse_yaml_file(const std::string & inputFileName)
{
  const std::string &aprepro = sierra::Env::get_param("aprepro");
  const bool use_aprepro = aprepro.empty() || aprepro != "off";
  const std::string &aprepro_defines = sierra::Env::get_param("define");

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

  YAML::Node doc;

  try {
    doc = parse_yaml_stream(input_stream);
    emit_to_stream(krinolog.getStream(), doc);
    krinolog << stk::diag::dendl;
  }
  catch (std::runtime_error &e) {
    krinolog << "Parser Error: " << e.what() << stk::diag::dendl;
    throw std::runtime_error("YAML parse failed");
  }


  return Node(doc);
}

bool Node::is_scalar() const { return Type() == YAML::NodeType::Scalar; }

Node Node::get_if_present(const std::string& key) const
{
  static std::string types[] = {"Null", "Scalar", "Sequence", "Map"};
  std::ostringstream err_msg;

  if (Type() == YAML::NodeType::Null || Type() == YAML::NodeType::Scalar)
  {
    emit_to_stream(err_msg, mNode);
    err_msg << "Check structure of input file.";
    stk::RuntimeDoomedAdHoc()
      << "parsing within non-searchable node of type = " << types[mNode.Type()]
      << " for key= " << key
      << " at " << line_info()
      << " node= " << err_msg.str() << std::endl;
    return Node();
  }
  return Node(mNode[key]);
}

Node Node::get_null_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Null); }
Node Node::get_map_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Map); }
Node Node::get_sequence_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Sequence); }
Node Node::get_scalar_if_present(const std::string& key) const { return get_type_if_present(key, YAML::NodeType::Scalar); }

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

std::string Node::line_info() const
{
  std::ostringstream sout;
  sout << "(line,column) = ("
       << mNode.Mark().line+1 << ", "
       << mNode.Mark().column+1 << ")";
  return sout.str();
}

Node Node::get_type_if_present(const std::string& key, YAML::NodeType::value type) const
{
  static std::string types[] = {"Null", "Scalar", "Sequence", "Map"};
  std::ostringstream err_msg;

  const Node value = get_if_present(key);
  if (value && (value.Type() != type))
  {
    emit_to_stream(err_msg, mNode);
    err_msg << "Check indentation of input file.";
    stk::RuntimeDoomedAdHoc()
      << "parsing expected type " << types[type] << " got type = " << types[value.Type()]
      << " for key= " << key
      << " at " << line_info()
      << " node= " << err_msg.str() << std::endl;
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

  int local_num_doomed = stk::get_doomed_count();
  int num_doomed = 0;

  MPI_Allreduce(&local_num_doomed, &num_doomed, 1, MPI_INT, MPI_SUM, stk::EnvData::parallel_comm());

  if (num_doomed)
  {
    throw std::runtime_error("Krino parse failed");
  }
}

}
}

