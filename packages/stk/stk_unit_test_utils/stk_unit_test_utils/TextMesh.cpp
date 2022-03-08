// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_util/util/ReportHandler.hpp"           // for ThrowRequireMsg

#include "TextMeshUtils.hpp"
#include "TextMesh.hpp"
#include "TextMeshStkTopologyMapping.hpp"

#include <ctype.h>                                   // for toupper
#include <stddef.h>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <map>
#include <set>                                       // for set
#include <sstream>                                   // for operator<<, etc
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector

#include <stk_io/IossBridge.hpp>                     // for is_part_io_part, etc
#include <stk_mesh/base/BulkData.hpp>                // for BulkData
#include <stk_mesh/base/FEMHelpers.hpp>              // for declare_element
#include <stk_mesh/base/Field.hpp>                   // for Field
#include <stk_mesh/base/GetEntities.hpp>             // for get_entities
#include <stk_mesh/base/MetaData.hpp>                // for MetaData, etc
#include "stk_mesh/base/TopologyDimensions.hpp"      // for ElementNode
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
using Topology = StkTopologyMapEntry;
using TextMeshData = text_mesh::TextMeshData<stk::mesh::EntityId, StkTopologyMapEntry>;
using ElementData = text_mesh::ElementData<stk::mesh::EntityId, StkTopologyMapEntry>;
using SidesetData = text_mesh::SidesetData<stk::mesh::EntityId, StkTopologyMapEntry>;
using NodesetData = text_mesh::NodesetData<stk::mesh::EntityId>;
using Coordinates = text_mesh::Coordinates<stk::mesh::EntityId>;
using TextMeshParser = text_mesh::TextMeshParser<stk::mesh::EntityId, StkTopologyMapping>;
using SideBlockInfo = text_mesh::SideBlockInfo;
using SplitType = text_mesh::SplitType;

inline std::ostream& operator<<(std::ostream& out, const StkTopologyMapEntry& t)
{
  return out << t.name();
}

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
    declare_nodeset_distribution_factor_fields();
    declare_sideset_distribution_factor_fields();
  }

private:
 void declare_nodeset_distribution_factor_fields()
 {
   for (const NodesetData& nodesetData : m_data.nodesets.get_group_data()) {
     stk::mesh::Part* part = m_meta.get_part(nodesetData.name);

     std::string nodesetDistFieldName = "distribution_factors_" + nodesetData.name;

     stk::mesh::Field<double>& distributionFactorsFieldPerNodeset =
         m_meta.declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, nodesetDistFieldName);

     stk::io::set_field_role(distributionFactorsFieldPerNodeset, Ioss::Field::MESH);
     stk::mesh::put_field_on_mesh(distributionFactorsFieldPerNodeset, *part,
         (stk::mesh::FieldTraits<stk::mesh::Field<double>>::data_type*) nullptr);
   }
 }

 void declare_sideblock_distribution_factor_field(
     const SideBlockInfo& sideBlock, stk::mesh::Field<double, stk::mesh::ElementNode>* distributionFactorsField)
 {
   stk::mesh::Part* sideBlockPart = m_meta.get_part(sideBlock.name);

   if (nullptr != distributionFactorsField) {
     stk::io::set_distribution_factor_field(*sideBlockPart, *distributionFactorsField);
     stk::mesh::put_field_on_mesh(*distributionFactorsField, *sideBlockPart, sideBlock.numNodesPerSide,
         (stk::mesh::FieldTraits<stk::mesh::Field<double, stk::mesh::ElementNode>>::data_type*) nullptr);
   }
 }

 stk::mesh::Field<double, stk::mesh::ElementNode>* declare_sideset_distribution_factor_field(
     const SidesetData& sidesetData)
 {
   stk::mesh::Field<double, stk::mesh::ElementNode>* distributionFactorsField = nullptr;
   stk::mesh::Part* sidesetPart = m_meta.get_part(sidesetData.name);

   SplitType splitType = sidesetData.get_split_type();
   if (splitType != SplitType::NO_SPLIT) {
     std::string fieldName = sidesetData.name + "_df";

     distributionFactorsField =
         &m_meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode>>(m_meta.side_rank(), fieldName);

     stk::io::set_field_role(*distributionFactorsField, Ioss::Field::MESH);
     stk::io::set_distribution_factor_field(*sidesetPart, *distributionFactorsField);
   }

   return distributionFactorsField;
 }

 void declare_sideset_distribution_factor_fields()
 {
   for (const SidesetData& sidesetData : m_data.sidesets.get_group_data()) {
     stk::mesh::Field<double, stk::mesh::ElementNode>* distributionFactorsField =
         declare_sideset_distribution_factor_field(sidesetData);
     std::vector<SideBlockInfo> sideBlocks = sidesetData.get_side_block_info();

     for (const auto& sideBlock : sideBlocks) {
       declare_sideblock_distribution_factor_field(sideBlock, distributionFactorsField);
     }
   }
 }

 void declare_block_parts()
 {
   for (const ElementData& elementData : m_data.elementDataVec) {
     if (m_meta.get_part(elementData.partName) == nullptr) {
       stk::mesh::Part& part = m_meta.declare_part_with_topology(elementData.partName, elementData.topology.topology);

       stk::io::put_io_part_attribute(part);
       m_meta.set_part_id(part, m_data.partIds.get(elementData.partName));
     }
   }
 }

 void declare_sideblock_part(stk::mesh::Part& sidesetPart, const SideBlockInfo& sideBlock)
 {
   stk::mesh::Part& sideBlockPart = m_meta.declare_part(sideBlock.name, m_meta.side_rank());

   stk::io::put_io_part_attribute(sideBlockPart);
   m_meta.set_part_id(sideBlockPart, sidesetPart.id());

   if (sidesetPart.mesh_meta_data_ordinal() != sideBlockPart.mesh_meta_data_ordinal()) {
     m_meta.declare_part_subset(sidesetPart, sideBlockPart);
   }
 }

 void declare_sideset_part(const SidesetData& sidesetData)
 {
   if (m_meta.get_part(sidesetData.name) == nullptr) {
     stk::mesh::Part& sidesetPart = m_meta.declare_part(sidesetData.name, m_meta.side_rank());

     stk::io::put_io_part_attribute(sidesetPart);
     m_meta.set_part_id(sidesetPart, sidesetData.id);

     std::vector<SideBlockInfo> sideBlocks = sidesetData.get_side_block_info();

     for (const auto& sideBlock : sideBlocks) {
       declare_sideblock_part(sidesetPart, sideBlock);
     }
   }
 }

 void declare_sideset_parts()
 {
   for (const SidesetData& sidesetData : m_data.sidesets.get_group_data()) {
     declare_sideset_part(sidesetData);
   }
 }

 void declare_nodeset_parts()
 {
   for (const NodesetData& nodesetData : m_data.nodesets.get_group_data()) {
     if (m_meta.get_part(nodesetData.name) == nullptr) {
       stk::mesh::Part& part = m_meta.declare_part(nodesetData.name, stk::topology::NODE_RANK);

       stk::io::put_io_part_attribute(part);
       m_meta.set_part_id(part, nodesetData.id);
     }
   }
 }

 void declare_parts()
 {
   declare_block_parts();
   declare_sideset_parts();
   declare_nodeset_parts();
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
    setup_blocks();
    setup_node_sharing();
    setup_sidesets();
    setup_nodesets();
    m_bulk.modification_end();
  }

 private:
  bool is_locally_owned(const ElementData& elem) { return elem.proc == m_bulk.parallel_rank(); }

  void add_element(const ElementData& elem)
  {
    stk::mesh::PartVector parts;
    parts.push_back(m_meta.get_part(elem.partName));
    parts.push_back(&m_meta.get_topology_root_part(elem.topology.topology));
    stk::mesh::declare_element(m_bulk, parts, elem.identifier, elem.nodeIds);
  }

  void setup_sideblock(const SideBlockInfo& sideBlock, const SidesetData& sidesetData, stk::mesh::SideSet& stkSideSet)
  {
    stk::mesh::Part* sideBlockPart = m_meta.get_part(sideBlock.name);
    std::vector<size_t> sideIndices =
        sidesetData.get_sideblock_indices_local_to_proc(sideBlock, m_bulk.parallel_rank());
    for (size_t sideIndex : sideIndices) {
      const SidesetData::DataType& elemSidePair = sidesetData.data[sideIndex];
      const stk::mesh::Entity elem = m_bulk.get_entity(stk::topology::ELEM_RANK, elemSidePair.first);
      if (m_bulk.is_valid(elem)) {
        int sideOrdinal = elemSidePair.second - 1;
        stkSideSet.add({elem, sideOrdinal});
        if (m_bulk.bucket(elem).owned()) {
          m_bulk.declare_element_side(elem, sideOrdinal, stk::mesh::PartVector{sideBlockPart});
        }
      }
    }
  }

  void setup_sidesets()
  {
    const bool fromInput = true;
    for (const SidesetData& sidesetData : m_data.sidesets.get_group_data()) {
      stk::mesh::Part* part = m_meta.get_part(sidesetData.name);
      stk::mesh::SideSet& stkSideSet = m_bulk.create_sideset(*part, fromInput);

      std::vector<SideBlockInfo> sideBlocks = sidesetData.get_side_block_info();

      for (const auto& sideBlock : sideBlocks) {
        setup_sideblock(sideBlock, sidesetData, stkSideSet);
      }
    }
  }

  void setup_nodesets()
  {
    for (const NodesetData& nodesetData : m_data.nodesets.get_group_data()) {
      stk::mesh::Part* part = m_meta.get_part(nodesetData.name);
      stk::mesh::PartVector addParts(1, part);

      for (const NodesetData::DataType& nodeId : nodesetData.data) {
        stk::mesh::Entity const node = m_bulk.get_entity(stk::topology::NODE_RANK, nodeId);
        if (m_bulk.is_valid(node)) {
          m_bulk.change_entity_parts(node, addParts);
        }
      }
    }
  }

  void setup_blocks()
  {
    for (const ElementData& elementData : m_data.elementDataVec) {
      if (is_locally_owned(elementData)) {
        add_element(elementData);
      }
    }
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

  void setup()
  {
    if (m_data.coords.has_coordinate_data()) {
      fill_node_list();
      fill_coordinate_field(m_data.coords);
      communicate_coordinate_field();
    }
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

    CoordinateInitializer coordInit(m_data, m_bulk);
    coordInit.setup();
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

std::string get_full_text_mesh_desc(const std::string& textMeshDesc, const std::vector<double>& coordVec)
{
  std::stringstream coords;
  coords << "|coordinates:";

  for (double coord : coordVec) {
    coords << coord << ",";
  }

  std::string meshDesc = textMeshDesc + coords.str();
  return meshDesc;
}

void setup_text_mesh(stk::mesh::BulkData& bulk, const std::string& meshDesc)
{
  TextMesh mesh(bulk, meshDesc);
  mesh.setup_mesh();
}

} // namespace unit_test_util
} // namespace stk
