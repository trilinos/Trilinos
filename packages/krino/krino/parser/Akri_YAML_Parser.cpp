// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_YAML_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_MeshInput_Parser.hpp>
#include <Akri_Simulation_Parser.hpp>

#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>
#include <stk_util/parallel/mpi_filebuf.hpp>

#include <fstream>

#ifdef KRINO_HAVE_YAML

static void emit(YAML::Emitter& emout, const YAML::Node & node) {
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
        emit(emout, subnode);
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
        emit(emout, value);
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

/// uses Emitter to print node to stream
static void emit(std::ostream& sout, const YAML::Node & node) {
  YAML::Emitter out;
  emit(out,node);
  sout << out.c_str() << std::endl;
}

static YAML::Node parse_yaml(std::istream & input_stream)
{
  return YAML::Load(input_stream);
  YAML::Node doc;

  try {
    doc = YAML::Load(input_stream);
  }
  catch (YAML::ParserException &e) {
    throw std::runtime_error(e.what());
  }

  return doc;
}

static std::string line_info(const YAML::Node & node) {
  std::ostringstream sout;
  sout << "(line,column) = ("
       << node.Mark().line+1 << ", "
       << node.Mark().column+1 << ")";
  return sout.str();
}

#else

static YAML::Node parse_yaml(std::istream & input_stream) { throw std::runtime_error("YAML parser not enabled in krino build."); }
static void emit(std::ostream& sout, const YAML::Node & node) {}
static std::string line_info(const YAML::Node & node) { return std::string(); }

#endif

namespace krino{
namespace YAML_Parser {

std::string info(const YAML::Node & node) {
  std::ostringstream sout;
  sout << "Node at " << line_info(node) << " => \n" ;
  emit(sout, node);
  return sout.str();
}

const YAML::Node
get_if_present(const YAML::Node& node, const std::string& key)
{
  static std::string types[] = {"Null", "Scalar", "Sequence", "Map"};
  std::ostringstream err_msg;

  if (node.Type() == YAML::NodeType::Null || node.Type() == YAML::NodeType::Scalar)
  {
    emit(err_msg, node);
    err_msg << "Check structure of input file.";
    stk::RuntimeDoomedAdHoc()
      << "parsing within non-searchable node of type = " << types[node.Type()]
      << " for key= " << key
      << " at " << line_info(node)
      << " node= " << err_msg.str() << std::endl;
    YAML::Node empty;
    return empty;
  }
  return node[key];
}

const YAML::Node
get_type_if_present(const YAML::Node& node, const std::string& key, YAML::NodeType::value type)
{
  static std::string types[] = {"Null", "Scalar", "Sequence", "Map"};
  std::ostringstream err_msg;

  const YAML::Node value = get_if_present(node, key);
  if (value && (value.Type() != type))
  {
    emit(err_msg, node);
    err_msg << "Check indentation of input file.";
    stk::RuntimeDoomedAdHoc()
      << "parsing expected type " << types[type] << " got type = " << types[value.Type()]
      << " for key= " << key
      << " at " << line_info(node)
      << " node= " << err_msg.str() << std::endl;
  }
  return value;
}

// template specialization
template<>
bool get_if_present(const YAML::Node & node, const std::string& key, bool& result)
{
  std::string resultString;
  const bool haveNode = get_if_present(node, key, resultString);

  if (haveNode)
  {
    result = false;
    std::transform(resultString.begin(), resultString.end(), resultString.begin(), ::toupper);
    if (resultString == "TRUE" || resultString == "YES" || resultString == "1")
      result = true;
    else if (resultString == "FALSE" || resultString == "NO" || resultString == "0")
      result = false;
    else
    {
      std::ostringstream err_msg;
      emit(err_msg, node);
      stk::RuntimeDoomedAdHoc() << "Failed to parse boolean "
        << " for key= " << key
        << " at " << line_info(node)
        << " node= " << err_msg.str() << std::endl;
    }
    return true;
  }
  return false;
}


void parse()
{/* %TRACE% */ Trace trace__("krino::YAML_Parser::parse()"); /* %TRACE% */
    const std::string inputFileName = stk::EnvData::getInputFileName();
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
        //throw std::runtime_error(input_parallel_file.aprepro_parse_errors());
      sierra::Env::outputP0()<<input_parallel_file.aprepro_parse_errors();
    }
  }

  std::istream input_stream(&input_parallel_file);

  // proceed with reading input file "document" from YAML
  YAML::Node doc;

  try {
    doc = parse_yaml(input_stream);
    emit(krinolog.getStream(), doc);
    krinolog << stk::diag::dendl;
  }
  catch (std::runtime_error &e) {
    krinolog << "Parser Error: " << e.what() << stk::diag::dendl;
    throw std::runtime_error("YAML parse failed");
  }

  MeshInput_Parser::parse(doc);
  Simulation_Parser::parse(doc);

  int local_num_doomed = stk::get_doomed_count();
  int num_doomed = 0;

  MPI_Allreduce(&local_num_doomed, &num_doomed, 1, MPI_INT, MPI_SUM, stk::EnvData::parallel_comm());

  if (num_doomed)
  {
    throw std::runtime_error("Krino parse failed");
  }
}

} // namespace parser
} // namespace krino
