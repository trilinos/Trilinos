// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "TextMeshUtils.hpp"
#include <ctype.h>                                   // for toupper
#include <stddef.h>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <map>
#include <set>                                       // for set
#include <sstream>                                   // for operator<<, etc
#include <stk_io/IossBridge.hpp>                     // for is_part_io_part, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/FEMHelpers.hpp>              // for declare_element
#include <stk_mesh/base/Field.hpp>                   // for Field
#include <stk_mesh/base/GetEntities.hpp>             // for get_entities
#include <stk_mesh/base/MetaData.hpp>                // for MetaData, etc
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"       // for Cartesian
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for field_data
#include "stk_mesh/base/Types.hpp"                   // for EntityId, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/util/ReportHandler.hpp"           // for ThrowRequireMsg

namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{
TopologyMapping::TopologyMapping()
{
  initialize_topology_map();
}

void TopologyMapping::initialize_topology_map()
{
  m_nameToTopology = {
      {"NODE", stk::topology::NODE},
      {"LINE_2", stk::topology::LINE_2},
      {"LINE_3", stk::topology::LINE_3},
      {"TRI_3", stk::topology::TRI_3},
      {"TRI_4", stk::topology::TRI_4},
      {"TRI_6", stk::topology::TRI_6},
      {"QUAD_4", stk::topology::QUAD_4},
      {"QUAD_6", stk::topology::QUAD_6},
      {"QUAD_8", stk::topology::QUAD_8},
      {"QUAD_9", stk::topology::QUAD_9},
      {"PARTICLE", stk::topology::PARTICLE},
      {"LINE_2_1D", stk::topology::LINE_2_1D},
      {"LINE_3_1D", stk::topology::LINE_3_1D},
      {"BEAM_2", stk::topology::BEAM_2},
      {"BEAM_3", stk::topology::BEAM_3},
      {"SHELL_LINE_2", stk::topology::SHELL_LINE_2},
      {"SHELL_LINE_3", stk::topology::SHELL_LINE_3},
      {"SPRING_2", stk::topology::SPRING_2},
      {"SPRING_3", stk::topology::SPRING_3},
      {"TRI_3_2D", stk::topology::TRI_3_2D},
      {"TRI_4_2D", stk::topology::TRI_4_2D},
      {"TRI_6_2D", stk::topology::TRI_6_2D},
      {"QUAD_4_2D", stk::topology::QUAD_4_2D},
      {"QUAD_8_2D", stk::topology::QUAD_8_2D},
      {"QUAD_9_2D", stk::topology::QUAD_9_2D},
      {"SHELL_TRI_3", stk::topology::SHELL_TRI_3},
      {"SHELL_TRI_4", stk::topology::SHELL_TRI_4},
      {"SHELL_TRI_6", stk::topology::SHELL_TRI_6},
      {"SHELL_QUAD_4", stk::topology::SHELL_QUAD_4},
      {"SHELL_QUAD_8", stk::topology::SHELL_QUAD_8},
      {"SHELL_QUAD_9", stk::topology::SHELL_QUAD_9},
      {"TET_4", stk::topology::TET_4},
      {"TET_8", stk::topology::TET_8},
      {"TET_10", stk::topology::TET_10},
      {"TET_11", stk::topology::TET_11},
      {"PYRAMID_5", stk::topology::PYRAMID_5},
      {"PYRAMID_13", stk::topology::PYRAMID_13},
      {"PYRAMID_14", stk::topology::PYRAMID_14},
      {"WEDGE_6", stk::topology::WEDGE_6},
      {"WEDGE_12", stk::topology::WEDGE_12},
      {"WEDGE_15", stk::topology::WEDGE_15},
      {"WEDGE_18", stk::topology::WEDGE_18},
      {"HEX_8", stk::topology::HEX_8},
      {"HEX_20", stk::topology::HEX_20},
      {"HEX_27", stk::topology::HEX_27},
  };
}

stk::topology::topology_t TopologyMapping::topology(const std::string& name) const
{
  auto it = m_nameToTopology.find(name);
  return (it != m_nameToTopology.end() ? it->second : stk::topology::INVALID_TOPOLOGY);
}

PartIdMapping::PartIdMapping() : m_idsAssigned(false) {}

void PartIdMapping::register_part_name(const std::string& name)
{
  m_partNames.push_back(name);
  handle_block_part(name);
}

void PartIdMapping::register_part_name_with_id(const std::string& name, unsigned id)
{
  register_part_name(name);
  assign(name, id);
}

unsigned PartIdMapping::get_part_id(const std::string& name) const
{
  auto it = m_ids.find(name);
  ThrowRequireMsg(it != m_ids.end(), "PartIdMapping has no ID for invalid part name " << name);
  return it->second;
}

unsigned PartIdMapping::get(const std::string& name) const
{
  if (!m_idsAssigned) assign_ids();
  return get_part_id(name);
}

std::string PartIdMapping::get_part_name(unsigned id) const
{
  auto it = m_parts.find(id);
  ThrowRequireMsg(it != m_parts.end(), "PartIdMapping has no part name for invalid id " << id);
  return it->second;
}

std::string PartIdMapping::get(unsigned id) const
{
  if (!m_idsAssigned) assign_ids();
  return get_part_name(id);
}

unsigned PartIdMapping::size() const
{
  if (!m_idsAssigned) assign_ids();
  return m_ids.size();
}

std::vector<std::string> PartIdMapping::get_part_names_sorted_by_id() const
{
  if (!m_idsAssigned) assign_ids();

  std::vector<std::string> names;
  names.reserve(m_parts.size());

  for (auto iter : m_parts) {
    names.push_back(iter.second);
  }

  return names;
}

void PartIdMapping::handle_block_part(const std::string& name)
{
  const std::string blockPrefix = "BLOCK_";
  const unsigned prefixLength = blockPrefix.length();

  if (name.length() < prefixLength + 1) return;

  const std::string namePrefix = name.substr(0, prefixLength);
  const std::string nameSuffix = name.substr(prefixLength);

  if (namePrefix != blockPrefix) return;

  unsigned id;
  std::istringstream nameSuffixStream(nameSuffix);
  nameSuffixStream >> id;
  if (nameSuffixStream.fail()) {
    return;
  }
  assign(name, id);
}

void PartIdMapping::assign_ids() const
{
  unsigned nextPartId = 1;
  for (const std::string& name : m_partNames) {
    if (m_ids.find(name) == m_ids.end()) {
      while (is_assigned(nextPartId)) nextPartId++;

      assign(name, nextPartId);
    }
  }

  m_idsAssigned = true;
}

void PartIdMapping::assign(const std::string& name, unsigned id) const
{
  validate_name_and_id(name, id);
  m_ids[name] = id;
  m_parts[id] = name;
}

void PartIdMapping::validate_name_and_id(const std::string& name, unsigned id) const
{
  if (is_registered(name)) {
    ThrowRequireMsg(
        m_ids[name] == id, "Cannot assign part '" << name << "' two different ids: " << m_ids[name] << " and " << id);
  } else {
    ThrowRequireMsg(
        !is_assigned(id), "Part id " << id << " has already been assigned, cannot assign it to part '" << name << "'");
  }
}

bool PartIdMapping::is_registered(const std::string& name) const
{
  return m_ids.count(name) > 0;
}

bool PartIdMapping::is_assigned(unsigned id) const
{
  return m_parts.count(id) > 0;
}

void TextMeshData::add_element(const ElementData& elem)
{
  elementDataVec.push_back(elem);
  for (const stk::mesh::EntityId& nodeId : elem.nodeIds) {
    nodeIds.insert(nodeId);
    associate_node_with_proc(nodeId, elem.proc);
  }
}

const EntityIdSet TextMeshData::nodes_on_proc(int proc) const
{
  auto it = m_nodesOnProc.find(proc);
  return it != m_nodesOnProc.end() ? it->second : EntityIdSet();
}

unsigned TextMeshData::num_nodes_on_proc(int proc) const
{
  auto it = m_nodesOnProc.find(proc);
  return it != m_nodesOnProc.end() ? it->second.size() : 0;
}

const std::set<int> TextMeshData::procs_for_node(const stk::mesh::EntityId& nodeId) const
{
  auto it = m_procsForNode.find(nodeId);
  return it != m_procsForNode.end() ? it->second : std::set<int>();
}

void TextMeshData::associate_node_with_proc(const stk::mesh::EntityId& nodeId, int proc)
{
  m_procsForNode[nodeId].insert(proc);
  m_nodesOnProc[proc].insert(nodeId);
}

TextMeshLexer::TextMeshLexer() : m_currentIndex(0), m_token(""), m_isNumber(false) {}

void TextMeshLexer::set_input_string(const std::string& input)
{
  m_input = input;
  m_currentIndex = 0;
  read_next_token();
}

int TextMeshLexer::get_int()
{
  read_next_token();
  return std::stoi(m_oldToken);
}

unsigned TextMeshLexer::get_unsigned()
{
  read_next_token();
  return std::stoul(m_oldToken);
}

std::string TextMeshLexer::get_string()
{
  read_next_token();
  return make_upper_case(m_oldToken);
}

void TextMeshLexer::get_newline()
{
  read_next_token();
}

bool TextMeshLexer::has_token() const
{
  return m_token != "";
}
bool TextMeshLexer::has_newline() const
{
  return m_token == "\n";
}
bool TextMeshLexer::has_number() const
{
  return has_token() && m_isNumber;
}
bool TextMeshLexer::has_string() const
{
  return has_token() && !has_number() && !has_newline();
}

void TextMeshLexer::read_next_token()
{
  m_oldToken = m_token;

  if (char_is_newline()) {
    m_isNumber = false;
    m_token = "\n";
    m_currentIndex++;
    return;
  }

  m_token = "";
  m_isNumber = true;

  while (has_more_input()) {
    if (char_is_whitespace()) {
      m_currentIndex++;
      continue;
    }

    if (char_is_comma()) {
      m_currentIndex++;
      break;
    }

    if (char_is_newline()) {
      break;
    }

    m_isNumber &= char_is_digit();
    m_token += current_char();
    m_currentIndex++;
  }
}

bool TextMeshLexer::has_more_input()
{
  return m_currentIndex < m_input.size();
}

bool TextMeshLexer::char_is_whitespace()
{
  return current_char() == ' ';
}
bool TextMeshLexer::char_is_comma()
{
  return current_char() == ',';
}
bool TextMeshLexer::char_is_newline()
{
  return current_char() == '\n';
}
bool TextMeshLexer::char_is_digit()
{
  return std::isdigit(static_cast<unsigned char>(current_char()));
}

char TextMeshLexer::current_char()
{
  return m_input[m_currentIndex];
}

std::string TextMeshLexer::make_upper_case(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
  return str;
}

TextMeshParser::TextMeshParser(unsigned dim) : m_lineNumber(0)
{
  m_data.spatialDim = dim;
}

TextMeshData TextMeshParser::parse(const std::string& meshDescription)
{
  initialize_parse(meshDescription);
  parse_description();
  return m_data;
}

void TextMeshParser::initialize_parse(const std::string& meshDescription)
{
  m_lexer.set_input_string(meshDescription);
  m_lineNumber = 1;
  validate_required_field(m_lexer.has_token());
}

void TextMeshParser::parse_description()
{
  while (m_lexer.has_token()) {
    ElementData elemData = parse_element();
    m_data.add_element(elemData);

    validate_no_extra_fields();
    parse_newline();
  }
}

ElementData TextMeshParser::parse_element()
{
  ElementData elem;
  elem.proc = parse_proc_id();
  elem.identifier = parse_elem_id();
  elem.topology = parse_topology();
  elem.nodeIds = parse_node_ids(elem.topology);
  elem.partName = parse_part(elem.topology);
  return elem;
}

int TextMeshParser::parse_proc_id()
{
  validate_required_field(m_lexer.has_number());
  return parse_int();
}

stk::mesh::EntityId TextMeshParser::parse_elem_id()
{
  validate_required_field(m_lexer.has_number());
  return parse_unsigned();
}

stk::topology TextMeshParser::parse_topology()
{
  validate_required_field(m_lexer.has_string());
  std::string topologyName = parse_string();

  stk::topology topology = m_topologyMapping.topology(topologyName);
  validate_topology(topology, topologyName);

  return topology;
}

stk::mesh::EntityIdVector TextMeshParser::parse_node_ids(const stk::topology& topology)
{
  stk::mesh::EntityIdVector nodeIds;
  while (m_lexer.has_number()) {
    nodeIds.push_back(parse_unsigned());
  }
  validate_node_count(topology, nodeIds.size());
  return nodeIds;
}

std::string TextMeshParser::parse_part(const stk::topology& topology)
{
  std::string partName;

  if (m_lexer.has_string()) {
    partName = parse_string();
  } else {
    partName = "block_" + topology.name();
  }

  if (m_lexer.has_number()) {
    unsigned partId = parse_unsigned();
    m_data.partIds.register_part_name_with_id(partName, partId);
  } else {
    m_data.partIds.register_part_name(partName);
  }

  return partName;
}

int TextMeshParser::parse_int()
{
  return m_lexer.get_int();
}
unsigned TextMeshParser::parse_unsigned()
{
  return m_lexer.get_unsigned();
}
std::string TextMeshParser::parse_string()
{
  return m_lexer.get_string();
}

void TextMeshParser::parse_newline()
{
  m_lexer.get_newline();
  m_lineNumber++;
}

void TextMeshParser::validate_required_field(bool hasNextRequiredField)
{
  ThrowRequireMsg(hasNextRequiredField,
      "Error!  Each line must contain the following fields (with at least one node):  "
      "Processor, GlobalId, Element Topology, NodeIds.  Error on line "
          << m_lineNumber << ".");
}

void TextMeshParser::validate_no_extra_fields()
{
  ThrowRequireMsg(!m_lexer.has_token() || m_lexer.has_newline(),
      "Error!  Each line should not contain more than the following fields (with at least one node):  "
      "Processor, GlobalId, Element Topology, NodeIds, Part Name, PartId.  "
      "Error on line "
          << m_lineNumber << ".");
}

void TextMeshParser::validate_topology(const stk::topology& topology, const std::string& providedName)
{
  ThrowRequireMsg(topology != stk::topology::INVALID_TOPOLOGY,
      "Error!  Topology = >>" << providedName << "<< is invalid from line " << m_lineNumber << ".");
  ThrowRequireMsg(topology.defined_on_spatial_dimension(m_data.spatialDim),
      "Error on input line " << m_lineNumber << ".  Topology = " << topology
                             << " is not defined on spatial dimension = " << m_data.spatialDim << " set in MetaData.");
}

void TextMeshParser::validate_node_count(const stk::topology& topology, size_t numNodes)
{
  size_t numTopologyNodes = topology.num_nodes();
  ThrowRequireMsg(numNodes == numTopologyNodes, "Error!  The input line appears to contain "
                                                    << numNodes << " nodes, but the topology " << topology << " needs "
                                                    << numTopologyNodes << " nodes on line " << m_lineNumber << ".");
}

Coordinates::Coordinates() {}

Coordinates::Coordinates(const TextMeshData& data, const std::vector<double>& coordinates)
{
  set_coordinate_data(data, coordinates);
}

void Coordinates::set_coordinate_data(const TextMeshData& data, const std::vector<double>& coordinates)
{
  if (!coordinates.empty()) {
    validate_num_coordinates(data, coordinates);
    fill_coordinate_map(data, coordinates);
  }
}

const std::vector<double>& Coordinates::operator[](const stk::mesh::EntityId& nodeId) const
{
  auto it(m_nodalCoords.find(nodeId));
  return it->second;
}

void Coordinates::validate_num_coordinates(const TextMeshData& data, const std::vector<double>& coordinates)
{
  ThrowRequireMsg(coordinates.size() == data.nodeIds.size() * data.spatialDim,
      "Number of coordinates: " << coordinates.size() << ", Number of nodes: " << data.nodeIds.size()
                                << ", Spatial dimension: " << data.spatialDim);
}

void Coordinates::fill_coordinate_map(const TextMeshData& data, const std::vector<double>& coordinates)
{
  std::vector<double>::const_iterator coordIter = coordinates.begin();
  for (const stk::mesh::EntityId& nodeId : data.nodeIds) {
    m_nodalCoords[nodeId] = std::vector<double>(coordIter, coordIter + data.spatialDim);
    coordIter += data.spatialDim;
  }
}

}  // namespace unit_test_util
}  // namespace stk
