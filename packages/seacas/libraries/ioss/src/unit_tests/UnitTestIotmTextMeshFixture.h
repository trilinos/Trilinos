// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <Ionit_Initializer.h>
#include <Ioss_DBUsage.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_Region.h>

#include "Ioss_DatabaseIO.h"     // for DatabaseIO
#include "Ioss_EntityType.h"     // for EntityType, etc
#include "Ioss_Field.h"          // for Field, etc
#include "Ioss_GroupingEntity.h" // for GroupingEntity
#include "Ioss_IOFactory.h"      // for IOFactory
#include "Ioss_MeshType.h"       // for MeshType, etc

#include <gtest/gtest.h>
#include <mpi.h>
#include <string>
#include <unordered_map>

#include <fmt/ostream.h>

#include <string>
#include <vector>
#include <strings.h>

#include "Ioss_CommSet.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_StandardElementTypes.h"
#include "Ioss_ParallelUtils.h"

#include "text_mesh/Iotm_TextMeshTopologyMapping.h"

#define ThrowRequireWithMsg(expr, message)                                   \
  do {                                                                       \
    if ( !(expr) ) {                                                         \
      std::ostringstream internal_throw_require_oss;                         \
      internal_throw_require_oss << message;                                 \
      IOSS_ERROR(internal_throw_require_oss );                               \
    }                                                                        \
  } while (false)

namespace Iotm{ namespace unit_test {

inline
void filename_substitution(std::string &filename)
{
  // See if filename contains "%P" which is replaced by the number of processors...
  // Assumes that %P only occurs once...
  // filename is changed.
  size_t pos = filename.find("%P");
  if (pos != std::string::npos) {
    // Found the characters...  Replace with the processor count...
    int num_proc = std::max(1, get_parallel_size());
    std::string tmp(filename, 0, pos);
    tmp += std::to_string(num_proc);
    tmp += filename.substr(pos+2);
    filename = tmp;
  }
}

class TextMeshFixture : public ::testing::Test
{
 protected:
  using EntityId = int64_t;
  using EntityIdVector = std::vector<EntityId>;
  using EntityIdSet = std::set<EntityId>;

  using PartNameId = std::pair<std::string, unsigned>;

  struct ElementInfo {
    const Ioss::ElementTopology* topology;
    EntityIdVector connectivity;

    ElementInfo() : topology(nullptr) { }
    ElementInfo(const Ioss::ElementTopology* topology_, const EntityIdVector& connectivity_)
    : topology(topology_)
    , connectivity(connectivity_) {}
  };

  struct PartInfo {
    std::string blockName;
    EntityIdSet ids;
  };

  TextMeshFixture(unsigned spatialDimension)
  : m_spatialDimension(spatialDimension)
  {
    Ioss::Init::Initializer io;
  }

  ~TextMeshFixture()
  {
    delete m_region;
  }

  int get_parallel_size()
  {
    return Ioss::ParallelUtils(get_comm()).parallel_size();
  }

  int get_parallel_rank()
  {
    return Ioss::ParallelUtils(get_comm()).parallel_rank();
  }

  void fill_mesh(const std::string& meshDesc)
  {
    std::pair<std::string, std::string> result = get_database_type_and_filename(meshDesc);

    std::string type = result.first;
    std::string filename = result.second;

    EXPECT_EQ("textmesh", type);

    create_database(filename, type);
    create_ioss_region();
  }

  MPI_Comm get_comm() const { return communicator;}
  void set_comm(MPI_Comm comm) {communicator = comm;}

  std::string get_mesh_desc(const std::string& textMeshDesc)
  {
    std::string header = "textmesh:";
    std::string meshDesc = header + textMeshDesc;
    return meshDesc;
  }

  std::string get_mesh_desc(const std::string& textMeshDesc, const std::vector<double>& coordVec)
  {
    std::stringstream coords;
    coords << "|coordinates:";

    for (double coord : coordVec) {
      coords << coord << ",";
    }

    std::string meshDesc = get_mesh_desc(textMeshDesc) + coords.str();
    return meshDesc;
  }

  std::string get_mesh_desc(const std::string& textMeshDesc, unsigned dimension)
  {
    std::stringstream dim;
    dim << "|dimension:" << dimension;

    std::string meshDesc = get_mesh_desc(textMeshDesc) + dim.str();
    return meshDesc;
  }

  std::string get_mesh_desc(const std::string& textMeshDesc, const std::vector<double>& coordVec, unsigned dimension)
  {
    std::stringstream dim;
    dim << "|dimension:" << dimension;

    std::string meshDesc = get_mesh_desc(textMeshDesc, coordVec) + dim.str();
    return meshDesc;
  }

  void verify_shared_nodes(const EntityIdVector& nodeIds, int sharingProc)
  {
    EXPECT_EQ(nodeIds.size(), get_node_sharing_count());

    for (EntityId nodeId : nodeIds) {
      EXPECT_TRUE(node_is_shared_with_proc(nodeId, sharingProc));
    }
  }

  void verify_num_elements(size_t goldCount)
  {
    size_t count = get_element_count();
    EXPECT_EQ(goldCount, count);
  }

  void verify_single_element(EntityId elemId, const std::string& textMeshTopologyName, const EntityIdVector& nodeIds)
  {
    Ioss::ElementTopology* topology = m_topologyMapping.topology(textMeshTopologyName).topology;
    ElementInfo info = get_element_info(elemId);
    EXPECT_TRUE(is_valid_element(info));
    EXPECT_EQ(topology, info.topology);
    verify_nodes_on_element(info, nodeIds);
  }

  void verify_part_membership(const std::vector<PartInfo> golds)
  {
    for (const PartInfo& gold : golds) {
      Ioss::ElementBlock* block = get_element_block(gold.blockName);

      verify_block(block);
      verify_elements_on_block(block, gold.ids);
    }
  }

  void verify_part_ids(const std::vector<PartNameId>& golds)
  {
    for (const PartNameId& gold : golds) {
      Ioss::ElementBlock* block = get_element_block(gold.first);

      verify_block(block);
      unsigned id = block->get_property("id").get_int();
      EXPECT_EQ(id, gold.second);
    }
  }

  void verify_nodes_on_element(const ElementInfo& info, const EntityIdVector& goldNodeIds)
  {
    EXPECT_EQ(goldNodeIds, info.connectivity);
  }

  void verify_coordinates(const EntityIdVector& goldNodeIds, const std::vector<double>& goldCoordinates)
  {
    CoordinateVerifier cv(*m_region, goldNodeIds, goldCoordinates);
    cv.verify();
  }

  void setup_text_mesh(const std::string& textMeshDesc)
  {
    fill_mesh(get_mesh_desc(textMeshDesc, m_spatialDimension));
  }

  void setup_text_mesh(const std::string& textMeshDesc, const std::vector<double>& coordinates)
  {
    fill_mesh(get_mesh_desc(textMeshDesc, coordinates, m_spatialDimension));
  }

  std::string get_topology_name(const std::string& textMeshTopologyName)
  {
    return m_topologyMapping.topology(textMeshTopologyName).name();
  }

  size_t db_api_int_size() const
  {
    assert(m_database != nullptr);
    return m_database->int_byte_size_api();
  }

  size_t get_node_sharing_count() const
  {
    ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

    Ioss::CommSet* io_cs = m_region->get_commset("commset_node");
    size_t numSharings = io_cs->get_field("entity_processor").raw_count();

    return numSharings;
  }

  template<typename INT>
  EntityIdVector get_element_ids_from_block_impl(const Ioss::ElementBlock* block) const
  {
    EntityIdVector elemIds;

    std::vector<INT> ids ;

    block->get_field_data("ids", ids);

    for(INT id : ids) {
      elemIds.push_back(static_cast<EntityId>(id));
    }

    return elemIds;
  }

  EntityIdVector get_element_ids_from_block(const Ioss::ElementBlock* block) const
  {
    if(db_api_int_size() == 4) {
      return get_element_ids_from_block_impl<int>(block);
    } else {
      return get_element_ids_from_block_impl<int64_t>(block);
    }
  }

  template<typename INT>
  ElementInfo get_element_info_from_block_impl(EntityId elemId, const Ioss::ElementBlock* block) const
  {
    const Ioss::ElementTopology* topo = nullptr;
    EntityIdVector elemConn;

    std::vector<INT> connectivity ;
    std::vector<INT> elemIds ;

    block->get_field_data("ids", elemIds);
    block->get_field_data("connectivity", connectivity);

    topo = block->topology();

    size_t elementCount = elemIds.size();
    int nodesPerElem = topo->number_nodes();

    for(size_t i=0; i<elementCount; ++i) {
      INT *conn = &connectivity[i*nodesPerElem];
      EntityId id = static_cast<EntityId>(elemIds[i]);

      if(id == elemId) {
        for(int j=0; j<nodesPerElem; j++) {
          elemConn.push_back(conn[j]);
        }

        return ElementInfo(topo, elemConn);
      }
    }

    return ElementInfo(nullptr, EntityIdVector());
  }

  template<typename INT>
  ElementInfo get_element_info_impl(EntityId elemId) const
  {
    ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");
    ElementInfo elemInfo;

    EntityIdVector elemConn;

    const Ioss::ElementBlockContainer& elemBlocks = m_region->get_element_blocks();
    bool found = false;

    for(const Ioss::ElementBlock* block : elemBlocks) {
      ElementInfo info = get_element_info_from_block_impl<INT>(elemId, block);

      if(is_valid_element(info)) {
        ThrowRequireWithMsg(!found, "Element with id " << elemId << " exists in more than one block!");
        found = true;
        elemInfo = info;
      }
    }

    return elemInfo;
  }

  ElementInfo get_element_info(EntityId elemId) const
  {
    if(db_api_int_size() == 4) {
      return get_element_info_impl<int>(elemId);
    } else {
      return get_element_info_impl<int64_t>(elemId);
    }
  }

  template<typename INT>
  size_t get_element_count_impl() const
  {
    ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

    const Ioss::ElementBlockContainer& elemBlocks = m_region->get_element_blocks();
    size_t count = 0;

    for(const Ioss::ElementBlock* block : elemBlocks) {
      std::vector<INT> elemIds ;
      block->get_field_data("ids", elemIds);
      count += elemIds.size();
    }

    return count;
  }

  size_t get_element_count() const
  {
    if(db_api_int_size() == 4) {
      return get_element_count_impl<int>();
    } else {
      return get_element_count_impl<int64_t>();
    }
  }

  template<typename INT>
  bool node_is_shared_with_proc_impl(EntityId nodeId, int sharingProc) const
  {
    ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

    Ioss::CommSet* io_cs = m_region->get_commset("commset_node");
    size_t numSharings = io_cs->get_field("entity_processor").raw_count();

    std::vector<INT> entityProc;
    io_cs->get_field_data("entity_processor", entityProc);

    for (size_t i = 0; i < numSharings; ++i) {
        EntityId iossNodeId = entityProc[i*2];
        int iossSharingProc = entityProc[i*2+1];

        if(iossNodeId == nodeId && iossSharingProc == sharingProc){
          return true;
        }
    }

    return false;
  }

  bool node_is_shared_with_proc(EntityId nodeId, int sharingProc) const
  {
    if(db_api_int_size() == 4) {
      return node_is_shared_with_proc_impl<int>(nodeId, sharingProc);
    } else {
      return node_is_shared_with_proc_impl<int64_t>(nodeId, sharingProc);
    }
  }

  bool is_valid_element(const ElementInfo& info) const
  {
    bool validTopology = info.topology != nullptr && info.topology != Ioss::ElementTopology::factory(Ioss::Unknown::name);
    bool validConnectivitySize = info.connectivity.size() != 0;
    bool validNumNodes = validTopology ? info.topology->number_nodes() == static_cast<int>(info.connectivity.size()) : false;

    return validConnectivitySize && validNumNodes;
  }

  Ioss::ElementBlock* get_element_block(const std::string& blockName) const
  {
    ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

    const Ioss::ElementBlockContainer& elemBlocks = m_region->get_element_blocks();
    Ioss::ElementBlock* elemBlock = nullptr;

    for(Ioss::ElementBlock* block : elemBlocks) {
      if(strcasecmp(block->name().c_str(), blockName.c_str()) == 0) {
        elemBlock = block;
      }
    }

    return elemBlock;
  }

  void verify_block(Ioss::ElementBlock* block)
  {
    ASSERT_TRUE(block != nullptr);
  }

  void verify_elements_on_block(const Ioss::ElementBlock* block, const std::set<EntityId>& goldIds)
  {
    EntityIdVector elemIds = get_element_ids_from_block(block);

    ASSERT_EQ(goldIds.size(), elemIds.size());
    for (EntityId elemId : elemIds) {
      EXPECT_EQ(1u, goldIds.count(elemId));
    }
  }

  void create_ioss_region()
  {
    if (m_region == nullptr) {
      EXPECT_TRUE(m_database != nullptr);
      m_region = new Ioss::Region(m_database, "input_model");
      EXPECT_TRUE(m_region != nullptr);
    }
  }

  void create_database(const std::string& fileName, const std::string& meshType)
  {
    if(m_database == nullptr) {
      Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL;

      std::string meshFileName = fileName;
      filename_substitution(meshFileName);
      m_database = Ioss::IOFactory::create(meshType, meshFileName, db_usage, get_comm(), m_propertyManager);
      EXPECT_TRUE(m_database != nullptr);
      EXPECT_TRUE(m_database->ok(true));
      EXPECT_EQ(m_database->get_format(), "TextMesh");
    }
  }

  std::pair<std::string, std::string> get_database_type_and_filename(const std::string& meshDesc)
  {
    std::string type;
    std::string filename;

    size_t colon = meshDesc.find(':');
    if (colon != std::string::npos && colon > 0) {
      type = meshDesc.substr(0, colon);
      filename = meshDesc.substr(colon + 1);
    } else {
      type = "textmesh";
      filename = meshDesc;
    }

    return std::make_pair(type, filename);
  }

  class CoordinateVerifier
  {
   public:
    CoordinateVerifier(const Ioss::Region& r, const EntityIdVector& ids, const std::vector<double>& coords)
    : region(r),
      spatialDim(region.get_property("spatial_dimension").get_int()),
      goldNodeIds(ids),
      goldCoordinates(coords)
    {
      fill_coordinates_from_ioss();
    }

    void verify()
    {
      verify_num_nodes();

      for (size_t nodeIndex = 0; nodeIndex < goldNodeIds.size(); nodeIndex++) {
        EntityId nodeId = goldNodeIds[nodeIndex];
        EXPECT_TRUE(is_valid_node(nodeId));

        const double* nodalCoords = get_nodal_coordinates(nodeId);
        const double* goldCoords = &goldCoordinates[nodeIndex * spatialDim];

        verify_nodal_coordinates(nodeId, goldCoords, nodalCoords);
      }
    }

   private:
    template <typename T>
    size_t field_data_from_ioss(Ioss::GroupingEntity *io_entity,
                              const std::string &io_fld_name,
                              std::vector<T>& io_field_data)
    {
      size_t io_entity_count = 0;
      if (io_entity->field_exists(io_fld_name)) {
        const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
        io_entity_count = io_entity->get_field_data(io_field.get_name(), io_field_data);
      }
      return io_entity_count;
    }

    size_t db_api_int_size() const
    {
      return region.get_database()->int_byte_size_api();
    }

    template <typename INT>
    EntityIdVector get_node_ids_impl() const
    {
      const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
      assert(node_blocks.size() == 1);

      Ioss::NodeBlock *nb = node_blocks[0];

      std::vector<INT> ids;
      nb->get_field_data("ids", ids);

      EntityIdVector nodeIds;
      for(INT id : ids) {
        nodeIds.push_back(static_cast<EntityId>(id));
      }

      return nodeIds;
    }

    EntityIdVector get_node_ids() const
    {
      if(db_api_int_size() == 4) {
        return get_node_ids_impl<int>();
      } else {
        return get_node_ids_impl<int64_t>();
      }
    }

    template <typename INT>
    bool is_valid_node_impl(EntityId nodeId) const
    {
      EntityIdVector ids = get_node_ids_impl<INT>();
      auto iter = std::find(ids.begin(), ids.end(), INT(nodeId));
      return iter != ids.end();
    }

    bool is_valid_node(EntityId nodeId) const
    {
      if(db_api_int_size() == 4) {
        return is_valid_node_impl<int>(nodeId);
      } else {
        return is_valid_node_impl<int64_t>(nodeId);
      }
    }

    void verify_num_nodes()
    {
      EntityIdVector ids = get_node_ids();
      EXPECT_EQ(goldNodeIds.size(), ids.size());
    }

    void fill_coordinate_map(const EntityIdVector& nodeIds, const std::vector<double>& coordinates)
    {
      std::vector<double>::const_iterator coordIter = coordinates.begin();
      for (const EntityId& nodeId : nodeIds) {
        m_nodalCoords[nodeId] = std::vector<double>(coordIter, coordIter + spatialDim);
        coordIter += spatialDim;
      }
    }

    void fill_coordinates_from_ioss()
    {
      std::vector<double> iossCoordinates;

      const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
      assert(node_blocks.size() == 1);

      Ioss::NodeBlock *nb = node_blocks[0];

      size_t node_count = nb->get_property("entity_count").get_int();

      EntityIdVector nodeIds = get_node_ids();

      size_t numIossNodes = field_data_from_ioss<double>(nb, "mesh_model_coordinates", iossCoordinates);
      ThrowRequireWithMsg(node_count == numIossNodes, "Node count mismatch");
      ThrowRequireWithMsg(iossCoordinates.size() == numIossNodes*spatialDim, "Invalid coordinate data size");

      fill_coordinate_map(nodeIds, iossCoordinates);
    }

    const std::vector<double>& operator[](const EntityId nodeId) const
    {
      auto it(m_nodalCoords.find(nodeId));
      return it->second;
    }

    const double* get_nodal_coordinates(const EntityId& nodeId) const
    {
      return (*this)[nodeId].data();
    }

    void verify_nodal_coordinates(const EntityId& nodeId, const double* goldCoords, const double* nodalCoords)
    {
      for (unsigned i = 0; i < spatialDim; i++) {
        EXPECT_NEAR(goldCoords[i], nodalCoords[i], 1.0e-9) << error_message(nodeId, i);
      }
    }

    std::string error_message(const EntityId& nodeId, unsigned coordIndex)
    {
      std::stringstream message;
      message << "Proc " << region.get_database()->util().parallel_rank() << ", Node ID " << nodeId << ", coord index " << coordIndex;
      return message.str();
    }

    const Ioss::Region& region;

    const unsigned spatialDim;

    const EntityIdVector& goldNodeIds;
    const std::vector<double>& goldCoordinates;
    std::unordered_map<EntityId, std::vector<double>> m_nodalCoords;
  };

  unsigned m_spatialDimension = 3;
  Ioss::PropertyManager m_propertyManager;
  Ioss::DatabaseIO* m_database = nullptr;
  Ioss::Region* m_region = nullptr;
  MPI_Comm communicator = MPI_COMM_WORLD;
  IossTopologyMapping m_topologyMapping;
};

}}

using TextMeshFixture = Iotm::unit_test::TextMeshFixture;
using EntityIdVector = TextMeshFixture::EntityIdVector;

namespace
{
class TestTextMesh : public Iotm::unit_test::TextMeshFixture
{
protected:
  TestTextMesh() : TextMeshFixture(3) { }
};

class TestTextMesh2d : public Iotm::unit_test::TextMeshFixture
{
protected:
  TestTextMesh2d() : TextMeshFixture(2) {}
};

class TestTextMesh1d : public Iotm::unit_test::TextMeshFixture
{
protected:
  TestTextMesh1d() : TextMeshFixture(1) {}
};
}

#endif /* SEACAS_LIBRARIES_IOSS_SRC_UNIT_TESTS_UNITTESTIOTMTEXTMESHFIXTURE_H_ */
