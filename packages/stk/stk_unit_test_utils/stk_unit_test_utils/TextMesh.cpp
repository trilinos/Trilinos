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
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/CompositeRank.hpp"
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for field_data
#include "stk_mesh/base/Types.hpp"                   // for EntityId, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"

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
using ElementDataLess = text_mesh::ElementDataLess<stk::mesh::EntityId, StkTopologyMapEntry>;
using SidesetData = text_mesh::SidesetData<stk::mesh::EntityId, StkTopologyMapEntry>;
using NodesetData = text_mesh::NodesetData<stk::mesh::EntityId>;
using AssemblyData = text_mesh::AssemblyData;
using Coordinates = text_mesh::Coordinates<stk::mesh::EntityId>;
using TextMeshParser = text_mesh::TextMeshParser<stk::mesh::EntityId, StkTopologyMapping>;
using ErrorHandler   = text_mesh::ErrorHandler;
using SideBlockInfo = text_mesh::SideBlockInfo;
using SplitType = text_mesh::SplitType;
using AssemblyType = text_mesh::AssemblyType;

class MetaDataInitializer
{
public:
  MetaDataInitializer(const StkTopologyMapping& t, const TextMeshData& d, stk::mesh::MetaData& m)
    : m_topologyMapping(t), m_data(d), m_meta(m)
  { }

  void setup()
  {
    declare_parts();
    declare_coordinate_field<double>();
    declare_nodeset_distribution_factor_fields();
    declare_sideset_distribution_factor_fields();
  }

private:
 void declare_nodeset_distribution_factor_fields()
 {
   for (const NodesetData& nodesetData : m_data.nodesets.get_group_data()) {
     stk::mesh::Part* part = m_meta.get_part(nodesetData.name);

     std::string nodesetDistFieldName = "distribution_factors_" + nodesetData.name;

     stk::mesh::Field<double> * distributionFactorsFieldPerNodeset =
         &m_meta.declare_field<double>(stk::topology::NODE_RANK, nodesetDistFieldName);

     stk::io::set_field_role(*distributionFactorsFieldPerNodeset, Ioss::Field::MESH);
     stk::mesh::put_field_on_mesh(*distributionFactorsFieldPerNodeset, *part, nullptr);
   }
 }

 void declare_sideblock_distribution_factor_field(const SideBlockInfo& sideBlock,
                                                  stk::mesh::Field<double>* distributionFactorsField)
 {
   stk::mesh::Part* sideBlockPart = m_meta.get_part(sideBlock.name);

   if (nullptr != distributionFactorsField) {
     stk::io::set_distribution_factor_field(*sideBlockPart, *distributionFactorsField);
     stk::mesh::put_field_on_mesh(*distributionFactorsField, *sideBlockPart, sideBlock.numNodesPerSide, nullptr);
   }
 }

 stk::mesh::Field<double>* declare_sideset_distribution_factor_field(const SidesetData& sidesetData)
 {
   stk::mesh::Field<double>* distributionFactorsField = nullptr;
   stk::mesh::Part* sidesetPart = m_meta.get_part(sidesetData.name);

   SplitType splitType = sidesetData.get_split_type();
   if (splitType != SplitType::NO_SPLIT) {
     std::string fieldName = sidesetData.name + "_df";

     distributionFactorsField = &m_meta.declare_field<double>(m_meta.side_rank(), fieldName);

     stk::io::set_field_role(*distributionFactorsField, Ioss::Field::MESH);
     stk::io::set_distribution_factor_field(*sidesetPart, *distributionFactorsField);
   }

   return distributionFactorsField;
 }

 void declare_sideset_distribution_factor_fields()
 {
   for (const SidesetData& sidesetData : m_data.sidesets.get_group_data()) {
     stk::mesh::Field<double>* distributionFactorsField = declare_sideset_distribution_factor_field(sidesetData);
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
   stk::topology sideTopology = m_topologyMapping.topology(sideBlock.sideTopology).topology;

   stk::mesh::Part* sideBlockPart = nullptr;

   if(stk::topology::INVALID_TOPOLOGY == sideTopology) {
     sideBlockPart = &m_meta.declare_part(sideBlock.name, m_meta.side_rank());
   } else {
     sideBlockPart = &m_meta.declare_part_with_topology(sideBlock.name, sideTopology);
   }

   STK_ThrowRequire(nullptr != sideBlockPart);
   stk::io::put_io_part_attribute(*sideBlockPart);
   m_meta.set_part_id(*sideBlockPart, sidesetPart.id());

   if (sidesetPart.mesh_meta_data_ordinal() != sideBlockPart->mesh_meta_data_ordinal()) {
     m_meta.declare_part_subset(sidesetPart, *sideBlockPart);
   }
 }

 void compute_block_membership(const SidesetData& sidesetData, const SideBlockInfo& sideBlock,
                               std::vector<const stk::mesh::Part*>& sideBlockTouchingBlockParts)
 {
   std::vector<int> l_blockIndex(m_data.partIds.size(), 0);
   std::vector<int> g_blockIndex(m_data.partIds.size(), 0);

   if (l_blockIndex.size() == 1) {
     l_blockIndex[0] = 1;
   } else {
     std::vector<std::string> partNames = m_data.partIds.get_part_names_sorted_by_id();
     std::sort(partNames.begin(), partNames.end());
     stk::mesh::BulkData& bulk = m_meta.mesh_bulk_data();

     std::vector<size_t> sideIndices = sidesetData.get_sideblock_indices_local_to_proc(sideBlock, bulk.parallel_rank());
     for (size_t sideIndex : sideIndices) {
       stk::mesh::EntityId elemId = sidesetData.data[sideIndex].first;
       auto elemIter = text_mesh::bound_search(m_data.elementDataVec.begin(), m_data.elementDataVec.end(), elemId, ElementDataLess());
       STK_ThrowRequire(elemIter != m_data.elementDataVec.end());

       std::string partName = elemIter->partName;

       auto partNameIter = text_mesh::bound_search(partNames.begin(), partNames.end(), partName);
       STK_ThrowRequire(partNameIter != partNames.end());

       unsigned index = std::distance(partNames.begin(), partNameIter);
       STK_ThrowRequire(index < m_data.partIds.size());

       l_blockIndex[index] = 1;
     }

     stk::all_reduce_max(bulk.parallel(), l_blockIndex.data(), g_blockIndex.data(), l_blockIndex.size());

     for(unsigned i=0; i<g_blockIndex.size(); i++) {
       if(g_blockIndex[i] == 1) {
         stk::mesh::Part* part = m_meta.get_part(partNames[i]);
         STK_ThrowRequire(nullptr != part);
         sideBlockTouchingBlockParts.push_back(part);
       }
     }
   }
 }

 void setup_surface_to_block_mapping(const SidesetData& sidesetData)
 {
   stk::mesh::Part* sidesetPart = m_meta.get_part(sidesetData.name);
   STK_ThrowRequire(nullptr != sidesetPart);

   std::vector<SideBlockInfo> sideBlocks = sidesetData.get_side_block_info();

   std::vector<const stk::mesh::Part*> touchingBlockParts;

   for (const auto& sideBlock : sideBlocks) {
     stk::mesh::Part* sideBlockPart = m_meta.get_part(sideBlock.name);
     std::vector<const stk::mesh::Part*> sideBlockTouchingBlockParts;

     if(!sideBlock.touchingBlock.empty()) {
       const stk::mesh::Part* touchingBlockPart = m_meta.get_part(sideBlock.touchingBlock);
       STK_ThrowRequire(nullptr != touchingBlockPart);
       sideBlockTouchingBlockParts.push_back(touchingBlockPart);
     } else {
       compute_block_membership(sidesetData, sideBlock, sideBlockTouchingBlockParts);
     }

     std::sort(sideBlockTouchingBlockParts.begin(), sideBlockTouchingBlockParts.end(), stk::mesh::PartLess());
     stk::util::insert_keep_sorted_and_unique(sideBlockTouchingBlockParts, touchingBlockParts, stk::mesh::PartLess());
     m_meta.set_surface_to_block_mapping(sideBlockPart, sideBlockTouchingBlockParts);
   }

   m_meta.set_surface_to_block_mapping(sidesetPart, touchingBlockParts);
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
     setup_surface_to_block_mapping(sidesetData);
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

 void declare_assembly_parts()
 {
   for (const AssemblyData& assemblyData : m_data.assemblies.get_group_data()) {
     if (m_meta.get_part(assemblyData.name) == nullptr) {
       stk::mesh::Part& part = m_meta.declare_part(assemblyData.name);
       stk::io::put_assembly_io_part_attribute(part);
       m_meta.set_part_id(part, assemblyData.id);
     }
   }
 }

 void declare_parts()
 {
   declare_block_parts();
   declare_sideset_parts();
   declare_nodeset_parts();
   declare_assembly_parts();
 }

  template<typename T>
  void declare_coordinate_field()
  {
    stk::mesh::Field<T>& coordsField = m_meta.declare_field<T>(stk::topology::NODE_RANK, m_meta.coordinate_field_name());
    stk::mesh::put_field_on_mesh(coordsField, m_meta.universal_part(), m_data.spatialDim, nullptr);
    stk::io::set_field_output_type(coordsField, stk::io::FieldOutputType::VECTOR_3D);
  }

  const StkTopologyMapping &m_topologyMapping;
  const TextMeshData& m_data;
  stk::mesh::MetaData& m_meta;
};

class BulkDataInitializer
{
public:
  BulkDataInitializer(const StkTopologyMapping& t, const TextMeshData& d, stk::mesh::BulkData& b)
    : m_topologyMapping(t),
      m_data(d),
      m_bulk(b),
      m_meta(m_bulk.mesh_meta_data())
  { }

  void setup()
  {
    m_bulk.modification_begin();
    setup_blocks();
    setup_node_sharing();
    m_bulk.modification_end();

    m_bulk.initialize_face_adjacent_element_graph();

    m_bulk.modification_begin();
    setup_sidesets();
    setup_nodesets();
    setup_assemblies();
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

  stk::topology get_side_block_topology_from_entries(const SidesetData& ss, const SideBlockInfo& sideBlock)
  {
    std::vector<size_t> sideIndices = ss.get_sideblock_indices_local_to_proc(sideBlock, m_bulk.parallel_rank());
    stk::topology invalidTopo = stk::topology::INVALID_TOPOLOGY;
    stk::topology blockSideTopo = stk::topology::INVALID_TOPOLOGY;
    int heterogenousTopo = 0;
    
    for (size_t sideIndex : sideIndices) {
      stk::mesh::EntityId elemId = ss.data[sideIndex].first;
      int elemSide = ss.data[sideIndex].second - 1;

      auto elemIter = text_mesh::bound_search(m_data.elementDataVec.begin(), m_data.elementDataVec.end(), elemId, ElementDataLess());
      STK_ThrowRequire(elemIter != m_data.elementDataVec.end());

      Topology blockTopo = elemIter->topology;
      Topology sideTopo = blockTopo.side_topology(elemSide);

      if(stk::topology::INVALID_TOPOLOGY != blockSideTopo && sideTopo.topology != blockSideTopo) {
	blockSideTopo = stk::topology::INVALID_TOPOLOGY;
	heterogenousTopo = 1;
	break;      
      }

      blockSideTopo = sideTopo.topology;	
    }

    int topoId = (stk::topology::topology_t) blockSideTopo;

    if (m_bulk.parallel_size()) {
      std::vector<int> l_topoData{topoId, heterogenousTopo};
      std::vector<int> g_topoData(2);
      stk::all_reduce_max(m_bulk.parallel(), l_topoData.data(), g_topoData.data(), 2);

      topoId = g_topoData[0];
      heterogenousTopo = g_topoData[1];
    }

    blockSideTopo = heterogenousTopo ? invalidTopo : stk::topology( (stk::topology::topology_t)topoId );

    return blockSideTopo;
  }

  bool set_sideset_topology(const SidesetData& ss, stk::mesh::Part* part,
			    const std::string& sideBlockTopo, stk::topology stkTopology,
			    bool printWarning = false)
  {
    if (stkTopology == stk::topology::INVALID_TOPOLOGY) {
      if (printWarning) {
	std::ostringstream os;
	os<<"TextMesh WARNING: failed to obtain sensible topology for sideset: " << ss.name<<", iossTopology: "<<sideBlockTopo<<", stk-part: "<<part->name()<<", rank: "<<part->primary_entity_rank()<<", stk-topology: "<<stkTopology<<". Probably because this SideSet is empty on this MPI rank. Unable to set correct stk topology hierarchy. Proceeding, but beware of unexpected behavior."<<std::endl;
	std::cerr<<os.str();
      }
    }
    else {
      for(auto subsetPart : part->subsets()) {
	if(subsetPart->topology() != stkTopology) {
	  return false;
	}
      }
	
      stk::mesh::set_topology(*part, stkTopology);
      return true;      
    }

    return false;    
  }

  
  void setup_sideset_topology(const SidesetData& ss, stk::mesh::Part* part)
  {
    if(nullptr == part) return;
    
    std::vector<SideBlockInfo> sideBlocks = ss.get_side_block_info();
    if(sideBlocks.size() != 1u) return;
    
    const SideBlockInfo& sideBlock = sideBlocks[0];

    if(sideBlock.name != ss.name) {
      stk::topology stkTopology = m_topologyMapping.topology(sideBlock.sideTopology).topology;
      if(set_sideset_topology(ss, part, sideBlock.sideTopology, stkTopology, true)) return;
    }

    stk::topology blockSideTopo = get_side_block_topology_from_entries(ss, sideBlock);
    set_sideset_topology(ss, part, sideBlock.sideTopology, blockSideTopo);
  }
  
  void setup_sidesets()
  {
    const bool fromInput = true;
    for (const SidesetData& sidesetData : m_data.sidesets.get_group_data()) {
      stk::mesh::Part* part = m_meta.get_part(sidesetData.name);
      stk::mesh::SideSet& stkSideSet = m_bulk.create_sideset(*part, fromInput);

      setup_sideset_topology(sidesetData, part);
      
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
        if (m_bulk.is_valid(node) && m_bulk.parallel_owner_rank(node) == m_bulk.parallel_rank()) {
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

  stk::mesh::EntityRank get_assembly_rank_from_type(const AssemblyType type)
  {
    if(type == AssemblyType::BLOCK) {
      return stk::topology::ELEM_RANK;
    } else if(type == AssemblyType::NODESET) {
      return stk::topology::NODE_RANK;
    } else if(type == AssemblyType::SIDESET) {
      return m_meta.side_rank();
    }

    return stk::topology::INVALID_RANK;
  }

  void setup_assemblies()
  {
    for (const AssemblyData& assemblyData : m_data.assemblies.get_group_data()) {
      stk::mesh::Part* part = m_meta.get_part(assemblyData.name);
      if(nullptr != part) {
        for(const std::string& member : assemblyData.data) {
          stk::mesh::Part* child = m_meta.get_part(member);

          if(nullptr != child) {
            m_meta.declare_part_subset(*part, *child);
          }
        }

        stk::mesh::EntityRank assemblyCompositeRank = stk::mesh::CompositeRank::get_rank(part);
        stk::mesh::EntityRank assemblyParsedRank = get_assembly_rank_from_type(assemblyData.get_assembly_type());

        if(stk::topology::INVALID_RANK != assemblyCompositeRank) {
          m_meta.declare_part(assemblyData.name, assemblyCompositeRank);

          STK_ThrowRequire((assemblyParsedRank == assemblyCompositeRank) ||
                       (assemblyParsedRank == stk::topology::INVALID_RANK));
        }
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

  const StkTopologyMapping &m_topologyMapping;
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
    m_topologyMapping.initialize_topology_map();
    m_data = m_parser.parse(meshDesc);
  }

  void setup_mesh()
  {
    MetaDataInitializer metaInit(m_topologyMapping, m_data, m_meta);
    metaInit.setup();

    BulkDataInitializer bulkInit(m_topologyMapping, m_data, m_bulk);
    bulkInit.setup();

    CoordinateInitializer coordInit(m_data, m_bulk);
    coordInit.setup();
  }

private:
  void validate_spatial_dim(unsigned spatialDim)
  {
    STK_ThrowRequireMsg(spatialDim == 2 || spatialDim == 3, "Error!  Spatial dimension not defined to be 2 or 3!");
  }

  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;

  TextMeshParser m_parser;
  TextMeshData m_data;
  StkTopologyMapping m_topologyMapping;
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

namespace simple_fields {

std::string get_full_text_mesh_desc(const std::string& textMeshConnectivityDesc,
                                    const std::vector<double>& coordVec)
{
  return stk::unit_test_util::get_full_text_mesh_desc(textMeshConnectivityDesc, coordVec);
}

void setup_text_mesh(stk::mesh::BulkData& bulkData, const std::string& meshDesc)
{
  return stk::unit_test_util::setup_text_mesh(bulkData, meshDesc);
}

}

} // namespace unit_test_util
} // namespace stk
