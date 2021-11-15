#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_TEXTMESHUTILS_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_TEXTMESHUTILS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "TextMesh.hpp"
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
using EntityIdSet = std::set<stk::mesh::EntityId>;

class TopologyMapping
{
 public:
  TopologyMapping();
  stk::topology::topology_t topology(const std::string& name) const;

 private:
  void initialize_topology_map();
  std::unordered_map<std::string, stk::topology::topology_t> m_nameToTopology;
};

class PartIdMapping
{
 public:
  PartIdMapping();

  void register_part_name(const std::string& name);

  void register_part_name_with_id(const std::string& name, unsigned id);

  unsigned get(const std::string& name) const;

  std::string get(unsigned id) const;

  unsigned size() const;

  std::vector<std::string> get_part_names_sorted_by_id() const;

 private:
  void handle_block_part(const std::string& name);

  void assign_ids() const;

  void assign(const std::string& name, unsigned id) const;

  void validate_name_and_id(const std::string& name, unsigned id) const;

  bool is_registered(const std::string& name) const;

  bool is_assigned(unsigned id) const;

  unsigned get_part_id(const std::string& name) const;
  std::string get_part_name(unsigned id) const;

  std::vector<std::string> m_partNames;
  mutable std::unordered_map<std::string, unsigned> m_ids;
  mutable std::map<unsigned, std::string> m_parts;
  mutable bool m_idsAssigned;
};

struct ElementData {
  int proc;
  stk::mesh::EntityId identifier;
  stk::topology topology;
  stk::mesh::EntityIdVector nodeIds;
  std::string partName = "";
};

struct TextMeshData {
  unsigned spatialDim;
  std::vector<ElementData> elementDataVec;
  PartIdMapping partIds;
  EntityIdSet nodeIds;

  void add_element(const ElementData& elem);

  const EntityIdSet nodes_on_proc(int proc) const;

  unsigned num_nodes_on_proc(int proc) const;

  const std::set<int> procs_for_node(const stk::mesh::EntityId& nodeId) const;

 private:
  void associate_node_with_proc(const stk::mesh::EntityId& nodeId, int proc);

  std::unordered_map<stk::mesh::EntityId, std::set<int>> m_procsForNode;
  std::unordered_map<int, EntityIdSet> m_nodesOnProc;
};

class TextMeshLexer
{
 public:
  TextMeshLexer();

  void set_input_string(const std::string& input);

  int get_int();

  unsigned get_unsigned();

  std::string get_string();

  void get_newline();

  bool has_token() const;
  bool has_newline() const;
  bool has_number() const;
  bool has_string() const;

 private:
  void read_next_token();

  bool has_more_input();
  bool char_is_whitespace();
  bool char_is_comma();
  bool char_is_newline();
  bool char_is_digit();

  char current_char();

  std::string make_upper_case(std::string str);

  std::string m_input;
  unsigned m_currentIndex;

  std::string m_oldToken;
  std::string m_token;

  bool m_isNumber;
};

class TextMeshParser
{
 public:
  TextMeshParser(unsigned dim);

  TextMeshData parse(const std::string& meshDescription);

 private:
  void initialize_parse(const std::string& meshDescription);

  void parse_description();

  ElementData parse_element();

  int parse_proc_id();

  stk::mesh::EntityId parse_elem_id();

  stk::topology parse_topology();

  stk::mesh::EntityIdVector parse_node_ids(const stk::topology& topology);

  std::string parse_part(const stk::topology& topology);

  int parse_int();
  unsigned parse_unsigned();
  std::string parse_string();

  void parse_newline();

  void validate_required_field(bool hasNextRequiredField);

  void validate_no_extra_fields();

  void validate_topology(const stk::topology& topology, const std::string& providedName);

  void validate_node_count(const stk::topology& topology, size_t numNodes);

  unsigned m_lineNumber;
  TextMeshData m_data;
  TextMeshLexer m_lexer;
  TopologyMapping m_topologyMapping;
};

class Coordinates
{
 public:
  Coordinates();
  Coordinates(const TextMeshData& data, const std::vector<double>& coordinates);

  const std::vector<double>& operator[](const stk::mesh::EntityId& nodeId) const;

  void set_coordinate_data(const TextMeshData& data, const std::vector<double>& coordinates);

 private:
  void validate_num_coordinates(const TextMeshData& data, const std::vector<double>& coordinates);

  void fill_coordinate_map(const TextMeshData& data, const std::vector<double>& coordinates);

  std::unordered_map<stk::mesh::EntityId, std::vector<double>> m_nodalCoords;
};

}  // namespace unit_test_util
}  // namespace stk

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_TEXTMESHUTILS_HPP_ */
