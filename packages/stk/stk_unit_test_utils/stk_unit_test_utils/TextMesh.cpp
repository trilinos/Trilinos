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
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
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

typedef std::set<stk::mesh::EntityId> EntityIdSet;

class TopologyMapping
{
public:
  stk::topology::topology_t topology(const std::string& name)
  {
    auto it = nameToTopology.find(name);
    return (it != nameToTopology.end() ? it->second : stk::topology::INVALID_TOPOLOGY);
  }

private:
  std::unordered_map<std::string, stk::topology::topology_t> nameToTopology =
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

struct ElementData
{
  int proc;
  stk::topology topology;
  stk::mesh::EntityId identifier;
  stk::mesh::EntityIdVector nodeIds;
  std::string partName = "";
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
    auto it = nodesOnProc.find(proc);
    return it != nodesOnProc.end() ? it->second : EntityIdSet();
  }

  const std::set<int> procs_for_node(const stk::mesh::EntityId& nodeId) const
  {
    auto it = procsForNode.find(nodeId);
    return it != procsForNode.end() ? it->second : std::set<int>();
  }

private:
  void associate_node_with_proc(const stk::mesh::EntityId& nodeId, int proc)
  {
    procsForNode[nodeId].insert(proc);
    nodesOnProc[proc].insert(nodeId);
  }

  std::unordered_map<stk::mesh::EntityId, std::set<int>> procsForNode;
  std::unordered_map<int, EntityIdSet> nodesOnProc;
};

class TextMeshParser
{
public:
  TextMeshParser(unsigned dim)
    : spatialDim(dim)
  {
    validate_spatial_dim();
  }

  TextMeshData parse(const std::string& meshDescription)
  {
    TextMeshData data;
    data.spatialDim = spatialDim;

    const std::vector<std::string> lines = split(meshDescription, '\n');

    lineNumber = 1;
    for (const std::string& line : lines) {
      ElementData elem = parse_element(line);
      data.add_element(elem);
      lineNumber++;
    }

    return data;
  }

private:
  void validate_spatial_dim()
  {
    ThrowRequireMsg(spatialDim == 2 || spatialDim == 3, "Error!  Spatial dimension not defined to be 2 or 3!");
  }

  std::vector<std::string> split(const std::string &s, char delim)
  {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      tokens.push_back(item);
    }
    return tokens;
  }

  ElementData parse_element(const std::string& rawLine)
  {
    std::string line = make_upper_case(remove_spaces(rawLine));

    const std::vector<std::string> tokens = split(line, ',');
    validate_min_token_count(tokens.size());

    ElementData elementData;
    elementData.proc = parse_proc(tokens[0]);
    elementData.identifier = parse_identifier(tokens[1]);
    elementData.topology = parse_topology(tokens[2]);

    unsigned numNodes = elementData.topology.num_nodes();

    validate_token_count(tokens.size(), numNodes, elementData.topology.name());

    for (unsigned i=0; i<numNodes; ++i) {
      stk::mesh::EntityId nodeId = parse_identifier(tokens[3+i]);
      elementData.nodeIds.push_back(nodeId);
    }

    if(tokens.size() == numNodes+4) {
      elementData.partName = tokens[3+numNodes];
    }
    else {
      elementData.partName = "block_" + elementData.topology.name();
    }

    return elementData;
  }

  std::string make_upper_case(std::string str)
  {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
  }

  std::string remove_spaces(std::string str)
  {
    std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
    str.erase(end_pos, str.end());
    return str;
  }

  void validate_min_token_count(size_t numTokens)
  {
    ThrowRequireMsg(numTokens >= 4,
                    "Error!  Each line must contain the following fields (with at least one node):  "
                    "Processor, GlobalId, Element Topology, NodeIds.  Error on line " << lineNumber << ".");
  }

  int parse_proc(const std::string& token)
  {
    return std::stoi(token);
  }

  stk::mesh::EntityId parse_identifier(const std::string& token)
  {
    return static_cast<stk::mesh::EntityId>(std::stoul(token));
  }

  stk::topology parse_topology(const std::string& token)
  {
    stk::topology topology = topologyMapping.topology(token);
    validate_topology(topology, token);
    return topology;
  }

  void validate_topology(const stk::topology& topology, const std::string& token)
  {
    ThrowRequireMsg(topology != stk::topology::INVALID_TOPOLOGY,
                    "Error!  Topology = >>" << token << "<< is invalid from line " << lineNumber << ".");
    ThrowRequireMsg(topology.defined_on_spatial_dimension(spatialDim),
                    "Error on input line " << lineNumber << ".  Topology = " << topology
                    << " is not defined on spatial dimension = " << spatialDim
                    << " set in MetaData.");
  }

  void validate_token_count(size_t numTokens, unsigned numNodes, const std::string& topologyName)
  {
    ThrowRequireMsg(numTokens == numNodes+3 || numTokens == numNodes+4,
                    "Error!  The input line appears to contain " << numTokens-3 << " nodes, but the topology "
                    << topologyName << " needs " << numNodes << " nodes on line " << lineNumber << ".");
  }

  unsigned spatialDim;
  TopologyMapping topologyMapping;
  size_t lineNumber;
};

class MetaDataInitializer
{
public:
  MetaDataInitializer(const TextMeshData& d, stk::mesh::MetaData& m)
    : data(d), meta(m)
  { }

  void setup()
  {
    declare_parts();
    declare_coordinate_field();
  }

private:
  void declare_parts()
  {
    for (const ElementData& elementData : data.elementDataVec) {
      stk::mesh::Part& part = meta.declare_part_with_topology(elementData.partName, elementData.topology);

      if (!stk::io::is_part_io_part(part)) {
        stk::io::put_io_part_attribute(part);
      }
    }
  }

  void declare_coordinate_field()
  {
    if (data.spatialDim == 3) {
      declare_coordinate_field_with_type<stk::mesh::CoordinatesField>();
    }
    else if (data.spatialDim == 2) {
      declare_coordinate_field_with_type<stk::mesh::Field<double, stk::mesh::Cartesian2d>>();
    }
  }

  template<typename F>
  void declare_coordinate_field_with_type()
  {
    F& coordsField = meta.declare_field<F>(stk::topology::NODE_RANK, meta.coordinate_field_name());
    stk::mesh::put_field_on_mesh(coordsField, meta.universal_part(), data.spatialDim, static_cast<double*>(nullptr));
  }

  const TextMeshData& data;
  stk::mesh::MetaData& meta;
};

class BulkDataInitializer
{
public:
  BulkDataInitializer(const TextMeshData& d, stk::mesh::BulkData& b)
    : data(d),
      bulk(b),
      meta(bulk.mesh_meta_data())
  { }

  void setup()
  {
    bulk.modification_begin();
    for (const ElementData& elementData : data.elementDataVec) {
      if (is_locally_owned(elementData)) {
        add_element(elementData);
      }
    }
    setup_node_sharing();
    bulk.modification_end();
  }

private:
  bool is_locally_owned(const ElementData& elem)
  {
    return elem.proc == bulk.parallel_rank();
  }

  void add_element(const ElementData& elem)
  {
    stk::mesh::Part* part = meta.get_part(elem.partName);
    stk::mesh::declare_element(bulk, *part, elem.identifier, elem.nodeIds);
  }

  void setup_node_sharing()
  {
    for (const stk::mesh::EntityId& nodeId : data.nodes_on_proc(bulk.parallel_rank())) {
      for (int proc : data.procs_for_node(nodeId)) {
        if (proc == bulk.parallel_rank()) continue;

        const stk::mesh::Entity& node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
        bulk.add_node_sharing(node, proc);
      }
    }
  }

  const TextMeshData& data;

  stk::mesh::BulkData& bulk;
  stk::mesh::MetaData& meta;
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
    auto it(nodalCoords.find(nodeId));
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
      nodalCoords[nodeId] = std::vector<double>(coordIter, coordIter+data.spatialDim);
      coordIter += data.spatialDim;
    }
  }

  std::unordered_map<stk::mesh::EntityId, std::vector<double>> nodalCoords;
};

class CoordinateInitializer
{
public:
  CoordinateInitializer(const TextMeshData& d, stk::mesh::BulkData& b)
    : data(d),
      bulk(b),
      meta(bulk.mesh_meta_data()),
      coordsField(*meta.coordinate_field())
  { }

  void setup(const std::vector<double>& coordinates)
  {
    Coordinates coords(data, coordinates);
    fill_node_list();
    fill_coordinate_field(coords);
    communicate_coordinate_field();
  }

private:
  void fill_node_list()
  {
    stk::mesh::get_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::NODE_RANK), nodes, true);
  }

  void fill_coordinate_field(const Coordinates& coordinates)
  {
    for (stk::mesh::Entity& node : nodes) {
      stk::mesh::EntityId nodeId = bulk.identifier(node);
      double* nodalCoordsLocation = static_cast<double*>(stk::mesh::field_data(coordsField, node));
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
    if (bulk.is_automatic_aura_on()) {
      stk::mesh::communicate_field_data(bulk, {&coordsField});
    }
  }

  const TextMeshData& data;

  stk::mesh::BulkData& bulk;
  stk::mesh::MetaData& meta;

  const stk::mesh::FieldBase& coordsField;
  stk::mesh::EntityVector nodes;
};

class TextMesh
{
public:
  TextMesh(stk::mesh::BulkData& b, const std::string& meshDesc)
    : bulk(b),
      meta(bulk.mesh_meta_data()),
      parser(meta.spatial_dimension()),
      data(parser.parse(meshDesc))
  { }

  void setup_mesh()
  {
    MetaDataInitializer metaInit(data, meta);
    metaInit.setup();

    BulkDataInitializer bulkInit(data, bulk);
    bulkInit.setup();
  }

  void setup_coordinates(const std::vector<double>& coordinates)
  {
    CoordinateInitializer coordInit(data, bulk);
    coordInit.setup(coordinates);
  }

private:
  stk::mesh::BulkData& bulk;
  stk::mesh::MetaData& meta;

  TextMeshParser parser;
  TextMeshData data;
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
