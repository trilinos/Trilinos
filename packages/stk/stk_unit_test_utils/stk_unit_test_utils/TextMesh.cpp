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
  stk::topology::topology_t topology(const std::string& name)
  {
    auto it = m_nameToTopology.find(name);
    return (it != m_nameToTopology.end() ? it->second : stk::topology::INVALID_TOPOLOGY);
  }

private:
  std::unordered_map<std::string, stk::topology::topology_t> m_nameToTopology =
  {
    {  "NODE"          , stk::topology::NODE         },
    {  "LINE_2"        , stk::topology::LINE_2       },
    {  "LINE_3"        , stk::topology::LINE_3       },
    {  "TRI_3"         , stk::topology::TRI_3        },
    {  "TRI_4"         , stk::topology::TRI_4        },
    {  "TRI_6"         , stk::topology::TRI_6        },
    {  "QUAD_4"        , stk::topology::QUAD_4       },
    {  "QUAD_6"        , stk::topology::QUAD_6       },
    {  "QUAD_8"        , stk::topology::QUAD_8       },
    {  "QUAD_9"        , stk::topology::QUAD_9       },
    {  "PARTICLE"      , stk::topology::PARTICLE     },
    {  "LINE_2_1D"     , stk::topology::LINE_2_1D    },
    {  "LINE_3_1D"     , stk::topology::LINE_3_1D    },
    {  "BEAM_2"        , stk::topology::BEAM_2       },
    {  "BEAM_3"        , stk::topology::BEAM_3       },
    {  "SHELL_LINE_2"  , stk::topology::SHELL_LINE_2 },
    {  "SHELL_LINE_3"  , stk::topology::SHELL_LINE_3 },
    {  "SPRING_2"      , stk::topology::SPRING_2     },
    {  "SPRING_3"      , stk::topology::SPRING_3     },
    {  "TRI_3_2D"      , stk::topology::TRI_3_2D     },
    {  "TRI_4_2D"      , stk::topology::TRI_4_2D     },
    {  "TRI_6_2D"      , stk::topology::TRI_6_2D     },
    {  "QUAD_4_2D"     , stk::topology::QUAD_4_2D    },
    {  "QUAD_8_2D"     , stk::topology::QUAD_8_2D    },
    {  "QUAD_9_2D"     , stk::topology::QUAD_9_2D    },
    {  "SHELL_TRI_3"   , stk::topology::SHELL_TRI_3  },
    {  "SHELL_TRI_4"   , stk::topology::SHELL_TRI_4  },
    {  "SHELL_TRI_6"   , stk::topology::SHELL_TRI_6  },
    {  "SHELL_QUAD_4"  , stk::topology::SHELL_QUAD_4 },
    {  "SHELL_QUAD_8"  , stk::topology::SHELL_QUAD_8 },
    {  "SHELL_QUAD_9"  , stk::topology::SHELL_QUAD_9 },
    {  "TET_4"         , stk::topology::TET_4        },
    {  "TET_8"         , stk::topology::TET_8        },
    {  "TET_10"        , stk::topology::TET_10       },
    {  "TET_11"        , stk::topology::TET_11       },
    {  "PYRAMID_5"     , stk::topology::PYRAMID_5    },
    {  "PYRAMID_13"    , stk::topology::PYRAMID_13   },
    {  "PYRAMID_14"    , stk::topology::PYRAMID_14   },
    {  "WEDGE_6"       , stk::topology::WEDGE_6      },
    {  "WEDGE_12"      , stk::topology::WEDGE_12     },
    {  "WEDGE_15"      , stk::topology::WEDGE_15     },
    {  "WEDGE_18"      , stk::topology::WEDGE_18     },
    {  "HEX_8"         , stk::topology::HEX_8        },
    {  "HEX_20"        , stk::topology::HEX_20       },
    {  "HEX_27"        , stk::topology::HEX_27       },
  };
};

class PartIdMapping
{
public:
  PartIdMapping()
    : m_idsAssigned(false)
  { }

  void register_part_name(const std::string& name)
  {
    m_partNames.push_back(name);
    handle_block_part(name);
  }

  void register_part_name_with_id(const std::string& name, unsigned id)
  {
    register_part_name(name);
    assign(name, id);
  }

  unsigned get(const std::string& name)
  {
    if (!m_idsAssigned) assign_ids();

    auto it = m_ids.find(name);
    ThrowRequireMsg(it != m_ids.end(), "PartIdMapping has no ID for invalid part name " << name);
    return it->second;
  }

private:
  void handle_block_part(const std::string& name)
  {
    const std::string blockPrefix = "BLOCK_";
    const unsigned prefixLength = blockPrefix.length();

    if (name.length() < prefixLength+1) return;

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

  void assign_ids()
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

  void assign(const std::string& name, unsigned id)
  {
    validate_name_and_id(name, id);
    m_ids[name] = id;
    m_assignedIds.insert(id);
  }

  void validate_name_and_id(const std::string& name, unsigned id)
  {
    if (is_registered(name)) {
      ThrowRequireMsg(m_ids[name] == id,
          "Cannot assign part '" << name << "' two different ids: " << m_ids[name] << " and " << id);
    }
    else {
      ThrowRequireMsg(!is_assigned(id),
          "Part id " << id << " has already been assigned, cannot assign it to part '" << name << "'");
    }
  }

  bool is_registered(const std::string& name)
  {
    return m_ids.count(name) > 0;
  }

  bool is_assigned(unsigned id) const
  {
    return m_assignedIds.count(id) > 0;
  }

  std::vector<std::string> m_partNames;
  std::set<unsigned> m_assignedIds;
  std::unordered_map<std::string, unsigned> m_ids;
  bool m_idsAssigned;
};

struct ElementData
{
  int proc;
  stk::mesh::EntityId identifier;
  stk::topology topology;
  stk::mesh::EntityIdVector nodeIds;
  std::string partName = "";
  unsigned partId;
};

struct TextMeshData
{
  unsigned spatialDim;
  std::vector<ElementData> elementDataVec;
  EntityIdSet nodeIds;

  void add_element(const ElementData& elem)
  {
    elementDataVec.push_back(elem);
    for (const stk::mesh::EntityId& nodeId : elem.nodeIds) {
      nodeIds.insert(nodeId);
      associate_node_with_proc(nodeId, elem.proc);
    }
  }

  const EntityIdSet nodes_on_proc(int proc) const
  {
    auto it = m_nodesOnProc.find(proc);
    return it != m_nodesOnProc.end() ? it->second : EntityIdSet();
  }

  const std::set<int> procs_for_node(const stk::mesh::EntityId& nodeId) const
  {
    auto it = m_procsForNode.find(nodeId);
    return it != m_procsForNode.end() ? it->second : std::set<int>();
  }

private:
  void associate_node_with_proc(const stk::mesh::EntityId& nodeId, int proc)
  {
    m_procsForNode[nodeId].insert(proc);
    m_nodesOnProc[proc].insert(nodeId);
  }

  std::unordered_map<stk::mesh::EntityId, std::set<int>> m_procsForNode;
  std::unordered_map<int, EntityIdSet> m_nodesOnProc;
};

class Lexer
{
public:
  Lexer()
    : m_token(""),
      m_isNumber(false)
  { }

  void set_input_string(const std::string& input)
  {
    m_input = input;
    m_currentIndex = 0;
    read_next_token();
  }

  int get_int()
  {
    read_next_token();
    return std::stoi(m_oldToken);
  }

  unsigned get_unsigned()
  {
    read_next_token();
    return std::stoul(m_oldToken);
  }

  std::string get_string()
  {
    read_next_token();
    return make_upper_case(m_oldToken);
  }

  void get_newline()
  {
    read_next_token();
  }

  bool has_token() const { return m_token != ""; }
  bool has_newline() const { return m_token == "\n"; }
  bool has_number() const { return has_token() && m_isNumber; }
  bool has_string() const { return has_token() && !has_number() && !has_newline(); }

private:
  void read_next_token()
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

  bool has_more_input()
  {
    return m_currentIndex < m_input.size();
  }

  bool char_is_whitespace() { return current_char() == ' '; }
  bool char_is_comma() { return current_char() == ','; }
  bool char_is_newline() { return current_char() == '\n'; }
  bool char_is_digit() { return std::isdigit(static_cast<unsigned char>(current_char())); }

  char current_char() { return m_input[m_currentIndex]; }

  std::string make_upper_case(std::string str)
  {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
  }

  std::string m_input;
  unsigned m_currentIndex;

  std::string m_oldToken;
  std::string m_token;

  bool m_isNumber;
};

class TextMeshParser
{
public:
  TextMeshParser(unsigned dim)
  {
    m_data.spatialDim = dim;
  }

  TextMeshData parse(const std::string& meshDescription)
  {
    initialize_parse(meshDescription);
    parse_description();
    set_part_ids();
    return m_data;
  }

private:
  void initialize_parse(const std::string& meshDescription)
  {
    m_lexer.set_input_string(meshDescription);
    m_lineNumber = 1;
    validate_required_field(m_lexer.has_token());
  }

  void parse_description()
  {
    while (m_lexer.has_token()) {
      ElementData elemData = parse_element();
      m_data.add_element(elemData);

      validate_no_extra_fields();
      parse_newline();
    }
  }

  ElementData parse_element()
  {
    ElementData elem;
    elem.proc = parse_proc_id();
    elem.identifier = parse_elem_id();
    elem.topology = parse_topology();
    elem.nodeIds = parse_node_ids(elem.topology);
    elem.partName = parse_part(elem.topology);
    return elem;
  }

  int parse_proc_id()
  {
    validate_required_field(m_lexer.has_number());
    return parse_int();
  }

  stk::mesh::EntityId parse_elem_id()
  {
    validate_required_field(m_lexer.has_number());
    return parse_unsigned();
  }

  stk::topology parse_topology()
  {
    validate_required_field(m_lexer.has_string());
    std::string topologyName = parse_string();

    stk::topology topology = m_topologyMapping.topology(topologyName);
    validate_topology(topology, topologyName);

    return topology;
  }

  stk::mesh::EntityIdVector parse_node_ids(const stk::topology& topology)
  {
    stk::mesh::EntityIdVector nodeIds;
    while (m_lexer.has_number()) {
      nodeIds.push_back(parse_unsigned());
    }
    validate_node_count(topology, nodeIds.size());
    return nodeIds;
  }

  std::string parse_part(const stk::topology& topology)
  {
    std::string partName;

    if (m_lexer.has_string()) {
      partName = parse_string();
    }
    else {
      partName = "block_" + topology.name();
    }

    if (m_lexer.has_number()) {
      unsigned partId = parse_unsigned();
      m_partIds.register_part_name_with_id(partName, partId);
    }
    else {
      m_partIds.register_part_name(partName);
    }

    return partName;
  }

  int parse_int() { return m_lexer.get_int(); }
  unsigned parse_unsigned() { return m_lexer.get_unsigned(); }
  std::string parse_string() { return m_lexer.get_string(); }

  void parse_newline()
  {
    m_lexer.get_newline();
    m_lineNumber++;
  }

  void validate_required_field(bool hasNextRequiredField)
  {
    ThrowRequireMsg(hasNextRequiredField,
                    "Error!  Each line must contain the following fields (with at least one node):  "
                    "Processor, GlobalId, Element Topology, NodeIds.  Error on line " << m_lineNumber << ".");
  }

  void validate_no_extra_fields()
  {
    ThrowRequireMsg(!m_lexer.has_token() || m_lexer.has_newline(),
                    "Error!  Each line should not contain more than the following fields (with at least one node):  "
                    "Processor, GlobalId, Element Topology, NodeIds, Part Name, PartId.  "
                    "Error on line " << m_lineNumber << ".");
  }

  void validate_topology(const stk::topology& topology, const std::string& providedName)
  {
    ThrowRequireMsg(topology != stk::topology::INVALID_TOPOLOGY,
                    "Error!  Topology = >>" << providedName << "<< is invalid from line " << m_lineNumber << ".");
    ThrowRequireMsg(topology.defined_on_spatial_dimension(m_data.spatialDim),
                    "Error on input line " << m_lineNumber << ".  Topology = " << topology
                    << " is not defined on spatial dimension = " << m_data.spatialDim
                    << " set in MetaData.");
  }

  void validate_node_count(const stk::topology& topology, size_t numNodes)
  {
    size_t numTopologyNodes = topology.num_nodes();
    ThrowRequireMsg(numNodes == numTopologyNodes,
                    "Error!  The input line appears to contain " << numNodes << " nodes, but the topology "
                    << topology << " needs " << numTopologyNodes << " nodes on line " << m_lineNumber << ".");
  }

  void set_part_ids()
  {
    for (ElementData& elem : m_data.elementDataVec) {
      elem.partId = m_partIds.get(elem.partName);
    }
  }

  unsigned m_lineNumber;
  TextMeshData m_data;
  Lexer m_lexer;
  TopologyMapping m_topologyMapping;
  PartIdMapping m_partIds;
};

class MetaDataInitializer
{
public:
  MetaDataInitializer(const TextMeshData& d, stk::mesh::MetaData& m)
    : m_data(d), m_meta(m)
  { }

  void setup()
  {
    declare_parts();
    declare_coordinate_field();
  }

private:
  void declare_parts()
  {
    for (const ElementData& elementData : m_data.elementDataVec) {
      if (m_meta.get_part(elementData.partName) == nullptr) {
        stk::mesh::Part& part = m_meta.declare_part_with_topology(elementData.partName, elementData.topology);

        stk::io::put_io_part_attribute(part);
        m_meta.set_part_id(part, elementData.partId);
      }
    }
  }

  void declare_coordinate_field()
  {
    if (m_data.spatialDim == 3) {
      declare_coordinate_field_with_type<stk::mesh::CoordinatesField>();
    }
    else if (m_data.spatialDim == 2) {
      declare_coordinate_field_with_type<stk::mesh::Field<double, stk::mesh::Cartesian2d>>();
    }
  }

  template<typename F>
  void declare_coordinate_field_with_type()
  {
    F& coordsField = m_meta.declare_field<F>(stk::topology::NODE_RANK, m_meta.coordinate_field_name());
    stk::mesh::put_field_on_mesh(coordsField, m_meta.universal_part(), m_data.spatialDim, static_cast<double*>(nullptr));
  }

  const TextMeshData& m_data;
  stk::mesh::MetaData& m_meta;
};

class BulkDataInitializer
{
public:
  BulkDataInitializer(const TextMeshData& d, stk::mesh::BulkData& b)
    : m_data(d),
      m_bulk(b),
      m_meta(m_bulk.mesh_meta_data())
  { }

  void setup()
  {
    m_bulk.modification_begin();
    for (const ElementData& elementData : m_data.elementDataVec) {
      if (is_locally_owned(elementData)) {
        add_element(elementData);
      }
    }
    setup_node_sharing();
    m_bulk.modification_end();
  }

private:
  bool is_locally_owned(const ElementData& elem)
  {
    return elem.proc == m_bulk.parallel_rank();
  }

  void add_element(const ElementData& elem)
  {
    stk::mesh::PartVector parts;
    parts.push_back(m_meta.get_part(elem.partName));
    parts.push_back(&m_meta.get_topology_root_part(elem.topology));
    stk::mesh::declare_element(m_bulk, parts, elem.identifier, elem.nodeIds);
  }

  void setup_node_sharing()
  {
    for (const stk::mesh::EntityId& nodeId : m_data.nodes_on_proc(m_bulk.parallel_rank())) {
      for (int proc : m_data.procs_for_node(nodeId)) {
        if (proc == m_bulk.parallel_rank()) continue;

        const stk::mesh::Entity& node = m_bulk.get_entity(stk::topology::NODE_RANK, nodeId);
        m_bulk.add_node_sharing(node, proc);
      }
    }
  }

  const TextMeshData& m_data;

  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;
};

class Coordinates
{
public:
  Coordinates(const TextMeshData& data, const std::vector<double> coordinates)
  {
    validate_num_coordinates(data, coordinates);
    fill_coordinate_map(data, coordinates);
  }

  const std::vector<double>& operator[](const stk::mesh::EntityId& nodeId) const
  {
    auto it(m_nodalCoords.find(nodeId));
    return it->second;
  }

private:
  void validate_num_coordinates(const TextMeshData& data, const std::vector<double>& coordinates)
  {
    ThrowRequireMsg(coordinates.size() == data.nodeIds.size()*data.spatialDim,
                    "Number of coordinates: " << coordinates.size()
                    << ", Number of nodes: " << data.nodeIds.size()
                    << ", Spatial dimension: " << data.spatialDim);
  }

  void fill_coordinate_map(const TextMeshData& data, const std::vector<double>& coordinates)
  {
    std::vector<double>::const_iterator coordIter = coordinates.begin();
    for (const stk::mesh::EntityId& nodeId : data.nodeIds) {
      m_nodalCoords[nodeId] = std::vector<double>(coordIter, coordIter+data.spatialDim);
      coordIter += data.spatialDim;
    }
  }

  std::unordered_map<stk::mesh::EntityId, std::vector<double>> m_nodalCoords;
};

class CoordinateInitializer
{
public:
  CoordinateInitializer(const TextMeshData& d, stk::mesh::BulkData& b)
    : m_data(d),
      m_bulk(b),
      m_meta(m_bulk.mesh_meta_data()),
      m_coordsField(*m_meta.coordinate_field())
  { }

  void setup(const std::vector<double>& coordinates)
  {
    Coordinates coords(m_data, coordinates);
    fill_node_list();
    fill_coordinate_field(coords);
    communicate_coordinate_field();
  }

private:
  void fill_node_list()
  {
    stk::mesh::get_selected_entities(m_meta.universal_part(), m_bulk.buckets(stk::topology::NODE_RANK), m_nodes, true);
  }

  void fill_coordinate_field(const Coordinates& coordinates)
  {
    for (stk::mesh::Entity& node : m_nodes) {
      stk::mesh::EntityId nodeId = m_bulk.identifier(node);
      double* nodalCoordsLocation = static_cast<double*>(stk::mesh::field_data(m_coordsField, node));
      copy_coordinates(coordinates[nodeId], nodalCoordsLocation);
    }
  }

  void copy_coordinates(const std::vector<double>& coords, double* dest)
  {
    for (unsigned i=0; i<coords.size(); i++) {
      dest[i] = coords[i];
    }
  }

  void communicate_coordinate_field()
  {
    if (m_bulk.is_automatic_aura_on()) {
      stk::mesh::communicate_field_data(m_bulk, {&m_coordsField});
    }
  }

  const TextMeshData& m_data;

  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;

  const stk::mesh::FieldBase& m_coordsField;
  stk::mesh::EntityVector m_nodes;
};

class TextMesh
{
public:
  TextMesh(stk::mesh::BulkData& b, const std::string& meshDesc)
    : m_bulk(b),
      m_meta(m_bulk.mesh_meta_data()),
      m_parser(m_meta.spatial_dimension())
  {
    validate_spatial_dim(m_meta.spatial_dimension());
    m_data = m_parser.parse(meshDesc);
  }

  void setup_mesh()
  {
    MetaDataInitializer metaInit(m_data, m_meta);
    metaInit.setup();

    BulkDataInitializer bulkInit(m_data, m_bulk);
    bulkInit.setup();
  }

  void setup_coordinates(const std::vector<double>& coordinates)
  {
    CoordinateInitializer coordInit(m_data, m_bulk);
    coordInit.setup(coordinates);
  }

private:
  void validate_spatial_dim(unsigned spatialDim)
  {
    ThrowRequireMsg(spatialDim == 2 || spatialDim == 3, "Error!  Spatial dimension not defined to be 2 or 3!");
  }

  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;

  TextMeshParser m_parser;
  TextMeshData m_data;
};

void setup_text_mesh(stk::mesh::BulkData& bulk, const std::string& meshDesc)
{
  TextMesh mesh(bulk, meshDesc);
  mesh.setup_mesh();
}

void setup_text_mesh(stk::mesh::BulkData& bulk, const std::string& meshDesc, const std::vector<double>& coordinates)
{
  TextMesh mesh(bulk, meshDesc);
  mesh.setup_mesh();
  mesh.setup_coordinates(coordinates);
}

} // namespace unit_test_util
} // namespace stk
