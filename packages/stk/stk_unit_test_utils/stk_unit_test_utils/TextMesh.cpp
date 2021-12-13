// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_util/util/ReportHandler.hpp"           // for ThrowRequireMsg

#include "TextMeshUtils.hpp"
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

namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{
class StkTopologyMapping : public text_mesh::TopologyMapping<stk::topology>
{
 public:
  stk::topology invalid_topology() const override { return stk::topology::INVALID_TOPOLOGY; }

  void initialize_topology_map() override
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
};

using TextMeshData = text_mesh::TextMeshData<stk::mesh::EntityId, stk::topology>;
using ElementData = text_mesh::ElementData<stk::mesh::EntityId, stk::topology>;
using Coordinates = text_mesh::Coordinates<stk::mesh::EntityId, stk::topology>;
using TextMeshParser = text_mesh::TextMeshParser<stk::mesh::EntityId, StkTopologyMapping>;

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
        m_meta.set_part_id(part, m_data.partIds.get(elementData.partName));
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
    Coordinates coords;
    coords.set_coordinate_data(m_data, coordinates);
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
