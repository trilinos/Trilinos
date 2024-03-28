// Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "Ioss_CodeTypes.h"

#include "Ionit_Initializer.h"
#include "Ioss_Assembly.h"
#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_EntityType.h"     // for EntityType, etc
#include "Ioss_Field.h"          // for Field, etc
#include "Ioss_GroupingEntity.h" // for GroupingEntity
#include "Ioss_IOFactory.h"      // for IOFactory
#include "Ioss_MeshType.h"       // for MeshType, etc
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_StandardElementTypes.h"

#include <gtest/gtest.h>
#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif
#include <string>
#include <unordered_map>

#include <memory>
#include <string>
#include <strings.h>
#include <vector>

#include "text_mesh/Iotm_TextMeshNodeset.h"
#include "text_mesh/Iotm_TextMeshSideset.h"
#include "text_mesh/Iotm_TextMeshTopologyMapping.h"
#include "text_mesh/Iotm_TextMeshUtils.h"

#define ThrowRequireWithMsg(expr, message)                                                         \
  do {                                                                                             \
    if (!(expr)) {                                                                                 \
      std::ostringstream internal_throw_require_oss;                                               \
      internal_throw_require_oss << message;                                                       \
      throw std::logic_error(internal_throw_require_oss.str());                                    \
    }                                                                                              \
  } while (false)

using Topology           = Iotm::TopologyMapEntry;
using TopologyMapping    = Iotm::IossTopologyMapping;
using EntityId           = int64_t;
using EntityIdVector     = std::vector<EntityId>;
using EntityIdSet        = std::set<EntityId>;
using TextMeshData       = Iotm::text_mesh::TextMeshData<EntityId, Topology>;
using ElementData        = Iotm::text_mesh::ElementData<EntityId, Topology>;
using SidesetData        = Iotm::text_mesh::SidesetData<EntityId, Topology>;
using NodesetData        = Iotm::text_mesh::NodesetData<EntityId>;
using Coordinates        = Iotm::text_mesh::Coordinates<EntityId>;
using TextMeshParser     = Iotm::text_mesh::TextMeshParser<EntityId, TopologyMapping>;
using SideAdjacencyGraph = Iotm::text_mesh::SideAdjacencyGraph<EntityId, Topology>;
using SideBlockInfo      = Iotm::text_mesh::SideBlockInfo;
using SideEntry          = std::pair<EntityId, int>;
using SideVector         = std::vector<SideEntry>;
using SplitType          = Iotm::text_mesh::SplitType;

struct Adjacency
{
  using NeighborVector       = std::vector<std::pair<int, SideAdjacencyGraph::IndexType>>;
  using SimpleNeighborVector = std::vector<SideAdjacencyGraph::IndexType>;

  size_t         elemIndex;
  NeighborVector neighborIndices;

  Adjacency(size_t elemIndex_, const NeighborVector &neighborIndices_)
      : elemIndex(elemIndex_), neighborIndices(neighborIndices_)
  {
  }

  Adjacency(size_t elemIndex_, const SimpleNeighborVector &neighborIndices_)
      : elemIndex(elemIndex_), neighborIndices(get_full_neighbor_vector(neighborIndices_))
  {
  }

  NeighborVector get_full_neighbor_vector(const SimpleNeighborVector &simpleNeighborVector)
  {
    NeighborVector fullNeighborVector;
    fullNeighborVector.reserve(simpleNeighborVector.size());

    for (unsigned i = 0; i < simpleNeighborVector.size(); i++) {
      fullNeighborVector.push_back(std::make_pair(static_cast<int>(i), simpleNeighborVector[i]));
    }
    return fullNeighborVector;
  }
};

struct SideEntryLess
{
  inline bool operator()(const SideEntry &lhs, const SideEntry &rhs) const
  {
    if (lhs.first < rhs.first)
      return true;
    else if (lhs.first == rhs.first && lhs.second < rhs.second)
      return true;

    return false;
  }
};

namespace Iotm {
  namespace unit_test {

    class AssemblyTreeGraph
    {
    public:
      AssemblyTreeGraph()                          = delete;
      AssemblyTreeGraph(const AssemblyTreeGraph &) = delete;

      AssemblyTreeGraph(Ioss::Region *region) : m_region(region) {}

    public:
      std::vector<std::string> get_unique_leaf_members(const std::string &name)
      {
        m_leafMembers.clear();
        m_visitedAssemblies.clear();

        for (Ioss::Assembly *assembly : m_region->get_assemblies()) {
          m_visitedAssemblies[assembly] = false;
        }

        traverse_tree(m_region->get_assembly(name));

        std::sort(m_leafMembers.begin(), m_leafMembers.end(), std::less<std::string>());
        auto endIter = std::unique(m_leafMembers.begin(), m_leafMembers.end());
        m_leafMembers.resize(endIter - m_leafMembers.begin());

        return m_leafMembers;
      }

    private:
      void traverse_tree(const Ioss::Assembly *assembly)
      {
        // Walk the tree without cyclic dependency
        if (assembly != nullptr) {
          if (m_visitedAssemblies[assembly] == false) {
            m_visitedAssemblies[assembly] = true;

            const Ioss::EntityType assemblyType = assembly->get_member_type();
            if (Ioss::ASSEMBLY != assemblyType) {
              for (const Ioss::GroupingEntity *ge : assembly->get_members()) {
                m_leafMembers.push_back(ge->name());
              }
            }
            else {
              for (const Ioss::GroupingEntity *ge : assembly->get_members()) {
                const Ioss::Assembly *assemblyMember = dynamic_cast<const Ioss::Assembly *>(ge);
                ThrowRequireWithMsg(nullptr != assemblyMember,
                                    "Non-assembly member: " << ge->name()
                                                            << " in ASSEMBLY rank assembly: "
                                                            << assembly->name());
                traverse_tree(assemblyMember);
              }
            }
          }
        }
      }

      Ioss::Region                                            *m_region = nullptr;
      mutable std::unordered_map<const Ioss::Assembly *, bool> m_visitedAssemblies;
      mutable std::vector<std::string>                         m_leafMembers;
    };

    class TextMeshFixture : public ::testing::Test
    {
    protected:
      using PartNameId = std::pair<std::string, unsigned>;

      struct ElementInfo
      {
        const Ioss::ElementTopology *topology;
        EntityIdVector               connectivity;

        ElementInfo() : topology(nullptr) {}
        ElementInfo(const Ioss::ElementTopology *topology_, const EntityIdVector &connectivity_)
            : topology(topology_), connectivity(connectivity_)
        {
        }
      };

      struct PartInfo
      {
        std::string blockName;
        EntityIdSet ids;
      };

      TextMeshFixture(unsigned spatialDimension) : m_spatialDimension(spatialDimension)
      {
        Ioss::Init::Initializer io;
        m_topologyMapping.initialize_topology_map();
      }

      ~TextMeshFixture() { delete m_region; }

      int get_parallel_size() { return Ioss::ParallelUtils(get_comm()).parallel_size(); }

      int get_parallel_rank() { return Ioss::ParallelUtils(get_comm()).parallel_rank(); }

      void fill_mesh(const std::string &meshDesc)
      {
        std::pair<std::string, std::string> result = get_database_type_and_filename(meshDesc);

        std::string type     = result.first;
        std::string filename = result.second;

        EXPECT_EQ("textmesh", type);

        create_database(filename, type);
        create_ioss_region();
      }

      Ioss_MPI_Comm get_comm() const { return Ioss::ParallelUtils::comm_world(); }

      std::string get_mesh_desc(const std::string &textMeshDesc)
      {
        std::string header   = "textmesh:";
        std::string meshDesc = header + textMeshDesc;
        return meshDesc;
      }

      std::string get_mesh_desc(const std::string &textMeshDesc, unsigned dimension)
      {
        std::stringstream dim;
        dim << "|dimension:" << dimension;

        std::string meshDesc = get_mesh_desc(textMeshDesc) + dim.str();
        return meshDesc;
      }

      void verify_shared_nodes(const EntityIdVector &nodeIds, int sharingProc)
      {
        EXPECT_EQ(nodeIds.size(), get_node_sharing_count(sharingProc));

        for (EntityId nodeId : nodeIds) {
          EXPECT_TRUE(node_is_shared_with_proc(nodeId, sharingProc));
        }
      }

      void verify_num_elements(size_t goldCount)
      {
        size_t count = get_element_count();
        EXPECT_EQ(goldCount, count);
      }

      void verify_single_element(EntityId elemId, const std::string &textMeshTopologyName,
                                 const EntityIdVector &nodeIds)
      {
        Topology    topology = m_topologyMapping.topology(textMeshTopologyName);
        ElementInfo info     = get_element_info(elemId);
        EXPECT_TRUE(is_valid_element(info));
        EXPECT_EQ(topology, info.topology);
        verify_nodes_on_element(info, nodeIds);
      }

      void verify_num_sidesets(size_t goldCount)
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");
        size_t count = m_region->get_sidesets().size();
        EXPECT_EQ(goldCount, count);
      }

      void verify_sideset_subset(const Ioss::SideSet *sideset, const unsigned id,
                                 const std::vector<std::string> &subsetNames)
      {
        EXPECT_TRUE(nullptr != sideset);
        EXPECT_EQ(id, sideset->get_property("id").get_int());

        if (subsetNames.empty()) {
          EXPECT_EQ(1u, sideset->get_side_blocks().size());
        }
        else {
          EXPECT_EQ(subsetNames.size(), sideset->get_side_blocks().size());
        }

        for (std::string subsetName : subsetNames) {
          std::transform(subsetName.begin(), subsetName.end(), subsetName.begin(), ::toupper);
          Ioss::SideBlock *sideBlock = sideset->get_side_block(subsetName);
          EXPECT_TRUE(nullptr != sideBlock);
          EXPECT_EQ(id, sideBlock->get_property("id").get_int());
        }
      }

      void verify_single_sideset(const std::string &name, const unsigned id,
                                 const SideVector &goldElemSidePairs)
      {
        verify_single_sideset(name, id, std::vector<std::string>{}, goldElemSidePairs);
      }

      void verify_single_sideset(const std::string &name, const unsigned id,
                                 const std::vector<std::string> &subsets,
                                 const SideVector               &goldElemSidePairs)
      {
        Ioss::SideSet *sideset = get_sideset(name);
        verify_sideset_subset(sideset, id, subsets);

        EXPECT_TRUE(nullptr != sideset);

        SideVector elemSidePairs = get_element_side_pairs_from_sideset(sideset);
        std::sort(elemSidePairs.begin(), elemSidePairs.end(), SideEntryLess());

        for (const SideEntry &sideEntry : goldElemSidePairs) {
          EntityId elemId = sideEntry.first;
          int      side   = sideEntry.second;

          ElementInfo info = get_element_info(elemId);
          EXPECT_TRUE(is_valid_element(info));

          EXPECT_TRUE(side > 0);
          EXPECT_TRUE(side <= (int)info.topology->number_boundaries());

          EXPECT_TRUE(std::binary_search(elemSidePairs.begin(), elemSidePairs.end(), sideEntry,
                                         SideEntryLess()));
        }
      }

      void verify_num_nodesets(size_t goldCount)
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");
        size_t count = m_region->get_nodesets().size();
        EXPECT_EQ(goldCount, count);
      }

      void verify_single_nodeset(const std::string &name, const unsigned id,
                                 const EntityIdVector &goldNodeIds)
      {
        Ioss::NodeSet *nodeset = get_nodeset(name);
        EXPECT_TRUE(nullptr != nodeset);
        EXPECT_EQ(id, nodeset->get_property("id").get_int());

        EntityIdVector nodeIds = get_node_ids_from_nodeset(nodeset);
        std::sort(nodeIds.begin(), nodeIds.end());

        for (EntityId node : goldNodeIds) {
          EXPECT_TRUE(std::binary_search(nodeIds.begin(), nodeIds.end(), node));
        }
      }

      void verify_num_assemblies(size_t goldCount)
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");
        size_t count = m_region->get_assemblies().size();
        EXPECT_EQ(goldCount, count);
      }

      void verify_single_assembly(const std::string &name, const unsigned id,
                                  const std::vector<std::string> &goldMembers)
      {
        Ioss::Assembly *assembly = get_assembly(name);
        EXPECT_TRUE(nullptr != assembly);
        EXPECT_EQ(id, assembly->get_property("id").get_int());

        AssemblyTreeGraph        graph(m_region);
        std::vector<std::string> leafMembers = graph.get_unique_leaf_members(name);
        EXPECT_EQ(goldMembers.size(), leafMembers.size());

        for (size_t i = 0; i < goldMembers.size(); i++) {
          const std::string &goldMember = goldMembers[i];
          const std::string &leafMember = leafMembers[i];
          EXPECT_EQ(0, strcasecmp(goldMember.c_str(), leafMember.c_str()))
              << "Comparison failure for " << name << ": " << goldMember << " <-> " << leafMember;
        }
      }

      void verify_part_membership(const std::vector<PartInfo> golds)
      {
        for (const PartInfo &gold : golds) {
          Ioss::ElementBlock *block = get_element_block(gold.blockName);

          verify_block(block);
          verify_elements_on_block(block, gold.ids);
        }
      }

      void verify_part_ids(const std::vector<PartNameId> &golds)
      {
        for (const PartNameId &gold : golds) {
          Ioss::ElementBlock *block = get_element_block(gold.first);

          verify_block(block);
          unsigned id = block->get_property("id").get_int();
          EXPECT_EQ(id, gold.second);
        }
      }

      void verify_nodes_on_element(const ElementInfo &info, const EntityIdVector &goldNodeIds)
      {
        EXPECT_EQ(goldNodeIds, info.connectivity);
      }

      void verify_coordinates(const EntityIdVector      &goldNodeIds,
                              const std::vector<double> &goldCoordinates)
      {
        CoordinateVerifier cv(*m_region, goldNodeIds, goldCoordinates);
        cv.verify();
      }

      void setup_text_mesh(const std::string &textMeshDesc)
      {
        fill_mesh(get_mesh_desc(textMeshDesc, m_spatialDimension));
      }

      std::string get_topology_name(const std::string &textMeshTopologyName)
      {
        return m_topologyMapping.topology(textMeshTopologyName).name();
      }

      size_t db_api_int_size() const
      {
        assert(m_database != nullptr);
        return m_database->int_byte_size_api();
      }

      size_t get_node_sharing_count(int sharingProc) const
      {
        if (db_api_int_size() == 4) {
          return get_node_sharing_count_impl<int>(sharingProc);
        }
        else {
          return get_node_sharing_count_impl<int64_t>(sharingProc);
        }
      }

      template <typename INT> size_t get_node_sharing_count_impl(int sharingProc) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        Ioss::CommSet *io_cs       = m_region->get_commset("commset_node");
        size_t         numSharings = io_cs->get_field("entity_processor").raw_count();

        std::vector<INT> entityProc;
        io_cs->get_field_data("entity_processor", entityProc);

        size_t count = 0;

        for (size_t i = 0; i < numSharings; ++i) {
          int iossSharingProc = entityProc[i * 2 + 1];

          if (iossSharingProc == sharingProc) {
            count++;
          }
        }

        return count;
      }

      template <typename INT>
      EntityIdVector get_element_ids_from_block_impl(const Ioss::ElementBlock *block) const
      {
        EntityIdVector elemIds;

        std::vector<INT> ids;

        block->get_field_data("ids", ids);

        for (INT id : ids) {
          elemIds.push_back(static_cast<EntityId>(id));
        }

        return elemIds;
      }

      EntityIdVector get_element_ids_from_block(const Ioss::ElementBlock *block) const
      {
        if (db_api_int_size() == 4) {
          return get_element_ids_from_block_impl<int>(block);
        }
        else {
          return get_element_ids_from_block_impl<int64_t>(block);
        }
      }

      template <typename INT>
      EntityIdVector get_node_ids_from_nodeset_impl(const Ioss::NodeSet *ns) const
      {
        EntityIdVector nodeIds;

        std::vector<INT> ids;

        ns->get_field_data("ids", ids);

        for (INT id : ids) {
          nodeIds.push_back(static_cast<EntityId>(id));
        }

        return nodeIds;
      }

      EntityIdVector get_node_ids_from_nodeset(const Ioss::NodeSet *ns) const
      {
        if (db_api_int_size() == 4) {
          return get_node_ids_from_nodeset_impl<int>(ns);
        }
        else {
          return get_node_ids_from_nodeset_impl<int64_t>(ns);
        }
      }

      template <typename INT>
      SideVector get_element_side_pairs_from_sideset_impl(const Ioss::SideSet *ss) const
      {
        SideVector elemSides;

        for (const Ioss::SideBlock *sb : ss->get_side_blocks()) {
          std::vector<INT> elemSideVec;
          sb->get_field_data("element_side", elemSideVec);

          for (unsigned i = 0; i < sb->entity_count(); i++) {
            EntityId elem = elemSideVec[2 * i + 0];
            int      side = elemSideVec[2 * i + 1];
            elemSides.push_back({elem, side});
          }
        }

        return elemSides;
      }

      SideVector get_element_side_pairs_from_sideset(const Ioss::SideSet *ss) const
      {
        if (db_api_int_size() == 4) {
          return get_element_side_pairs_from_sideset_impl<int>(ss);
        }
        else {
          return get_element_side_pairs_from_sideset_impl<int64_t>(ss);
        }
      }

      template <typename INT>
      ElementInfo get_element_info_from_block_impl(EntityId                  elemId,
                                                   const Ioss::ElementBlock *block) const
      {
        const Ioss::ElementTopology *topo = nullptr;
        EntityIdVector               elemConn;

        std::vector<INT> connectivity;
        std::vector<INT> elemIds;

        block->get_field_data("ids", elemIds);
        block->get_field_data("connectivity", connectivity);

        topo = block->topology();

        size_t elementCount = elemIds.size();
        int    nodesPerElem = topo->number_nodes();

        for (size_t i = 0; i < elementCount; ++i) {
          INT *conn = &connectivity[i * nodesPerElem];
          auto id   = static_cast<EntityId>(elemIds[i]);

          if (id == elemId) {
            for (int j = 0; j < nodesPerElem; j++) {
              elemConn.push_back(conn[j]);
            }

            return ElementInfo(topo, elemConn);
          }
        }

        return ElementInfo(nullptr, EntityIdVector());
      }

      template <typename INT> ElementInfo get_element_info_impl(EntityId elemId) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");
        ElementInfo elemInfo;

        EntityIdVector elemConn;

        const Ioss::ElementBlockContainer &elemBlocks = m_region->get_element_blocks();
        bool                               found      = false;

        for (const Ioss::ElementBlock *block : elemBlocks) {
          ElementInfo info = get_element_info_from_block_impl<INT>(elemId, block);

          if (is_valid_element(info)) {
            ThrowRequireWithMsg(!found,
                                "Element with id " << elemId << " exists in more than one block!");
            found    = true;
            elemInfo = info;
          }
        }

        return elemInfo;
      }

      ElementInfo get_element_info(EntityId elemId) const
      {
        if (db_api_int_size() == 4) {
          return get_element_info_impl<int>(elemId);
        }
        else {
          return get_element_info_impl<int64_t>(elemId);
        }
      }

      template <typename INT> size_t get_element_count_impl() const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        const Ioss::ElementBlockContainer &elemBlocks = m_region->get_element_blocks();
        size_t                             count      = 0;

        for (const Ioss::ElementBlock *block : elemBlocks) {
          std::vector<INT> elemIds;
          block->get_field_data("ids", elemIds);
          count += elemIds.size();
        }

        return count;
      }

      size_t get_element_count() const
      {
        if (db_api_int_size() == 4) {
          return get_element_count_impl<int>();
        }
        else {
          return get_element_count_impl<int64_t>();
        }
      }

      template <typename INT>
      bool node_is_shared_with_proc_impl(EntityId nodeId, int sharingProc) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        Ioss::CommSet *io_cs       = m_region->get_commset("commset_node");
        size_t         numSharings = io_cs->get_field("entity_processor").raw_count();

        std::vector<INT> entityProc;
        io_cs->get_field_data("entity_processor", entityProc);

        for (size_t i = 0; i < numSharings; ++i) {
          EntityId iossNodeId      = entityProc[i * 2];
          int      iossSharingProc = entityProc[i * 2 + 1];

          if (iossNodeId == nodeId && iossSharingProc == sharingProc) {
            return true;
          }
        }

        return false;
      }

      bool node_is_shared_with_proc(EntityId nodeId, int sharingProc) const
      {
        if (db_api_int_size() == 4) {
          return node_is_shared_with_proc_impl<int>(nodeId, sharingProc);
        }
        else {
          return node_is_shared_with_proc_impl<int64_t>(nodeId, sharingProc);
        }
      }

      bool is_valid_element(const ElementInfo &info) const
      {
        bool validTopology = info.topology != nullptr &&
                             info.topology != Ioss::ElementTopology::factory(Ioss::Unknown::name);
        bool validConnectivitySize = info.connectivity.size() != 0;
        bool validNumNodes         = validTopology ? info.topology->number_nodes() ==
                                                 static_cast<int>(info.connectivity.size())
                                                   : false;

        return validConnectivitySize && validNumNodes;
      }

      Ioss::ElementBlock *get_element_block(const std::string &blockName) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        const Ioss::ElementBlockContainer &elemBlocks = m_region->get_element_blocks();
        Ioss::ElementBlock                *elemBlock  = nullptr;

        for (Ioss::ElementBlock *block : elemBlocks) {
          if (strcasecmp(block->name().c_str(), blockName.c_str()) == 0) {
            elemBlock = block;
          }
        }

        return elemBlock;
      }

      Ioss::NodeSet *get_nodeset(const std::string &name) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        const Ioss::NodeSetContainer &nodesets = m_region->get_nodesets();
        Ioss::NodeSet                *nodeset  = nullptr;

        for (Ioss::NodeSet *ns : nodesets) {
          if (strcasecmp(ns->name().c_str(), name.c_str()) == 0) {
            nodeset = ns;
          }
        }

        return nodeset;
      }

      Ioss::SideSet *get_sideset(const std::string &name) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        const Ioss::SideSetContainer &sidesets = m_region->get_sidesets();
        Ioss::SideSet                *sideset  = nullptr;

        for (Ioss::SideSet *ss : sidesets) {
          if (strcasecmp(ss->name().c_str(), name.c_str()) == 0) {
            sideset = ss;
          }
        }

        return sideset;
      }

      Ioss::Assembly *get_assembly(const std::string &name) const
      {
        ThrowRequireWithMsg(m_region != nullptr, "Ioss region has not been created");

        const Ioss::AssemblyContainer &assemblies = m_region->get_assemblies();
        Ioss::Assembly                *assembly   = nullptr;

        for (Ioss::Assembly *ass : assemblies) {
          if (strcasecmp(ass->name().c_str(), name.c_str()) == 0) {
            assembly = ass;
          }
        }

        return assembly;
      }

      void verify_block(Ioss::ElementBlock *block) { ASSERT_TRUE(block != nullptr); }

      void verify_elements_on_block(const Ioss::ElementBlock *block,
                                    const std::set<EntityId> &goldIds)
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

      void create_database(const std::string &fileName, const std::string &meshType)
      {
        if (m_database == nullptr) {
          Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL;

          std::string meshFileName = fileName;
          filename_substitution(meshFileName);
          m_database = Ioss::IOFactory::create(meshType, meshFileName, db_usage, get_comm(),
                                               m_propertyManager);
          EXPECT_TRUE(m_database != nullptr);
          EXPECT_TRUE(m_database->ok(true));
          EXPECT_EQ(m_database->get_format(), "TextMesh");
        }
      }

      void filename_substitution(std::string &filename)
      {
        // See if filename contains "%P" which is replaced by the number of processors...
        // Assumes that %P only occurs once...
        // filename is changed.
        size_t pos = filename.find("%P");
        if (pos != std::string::npos) {
          // Found the characters...  Replace with the processor count...
          int         num_proc = std::max(1, get_parallel_size());
          std::string tmp(filename, 0, pos);
          tmp += std::to_string(num_proc);
          tmp += filename.substr(pos + 2);
          filename = tmp;
        }
      }

      std::pair<std::string, std::string>
      get_database_type_and_filename(const std::string &meshDesc)
      {
        std::string type;
        std::string filename;

        size_t colon = meshDesc.find(':');
        if (colon != std::string::npos && colon > 0) {
          type     = meshDesc.substr(0, colon);
          filename = meshDesc.substr(colon + 1);
        }
        else {
          type     = "textmesh";
          filename = meshDesc;
        }

        return std::make_pair(type, filename);
      }

      class CoordinateVerifier
      {
      public:
        CoordinateVerifier(const Ioss::Region &r, const EntityIdVector &ids,
                           const std::vector<double> &coords)
            : region(r), spatialDim(region.get_property("spatial_dimension").get_int()),
              goldNodeIds(ids), goldCoordinates(coords)
        {
          fill_coordinates_from_ioss();
        }

        void verify()
        {
          verify_num_nodes();

          for (size_t nodeIndex = 0; nodeIndex < goldNodeIds.size(); nodeIndex++) {
            EntityId nodeId = goldNodeIds[nodeIndex];
            EXPECT_TRUE(is_valid_node(nodeId));

            const double *nodalCoords = get_nodal_coordinates(nodeId);
            const double *goldCoords  = &goldCoordinates[nodeIndex * spatialDim];

            verify_nodal_coordinates(nodeId, goldCoords, nodalCoords);
          }
        }

      private:
        template <typename T>
        size_t field_data_from_ioss(Ioss::GroupingEntity *io_entity, const std::string &io_fld_name,
                                    std::vector<T> &io_field_data)
        {
          size_t io_entity_count = 0;
          if (io_entity->field_exists(io_fld_name)) {
            const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
            io_entity_count = io_entity->get_field_data(io_field.get_name(), io_field_data);
          }
          return io_entity_count;
        }

        size_t db_api_int_size() const { return region.get_database()->int_byte_size_api(); }

        template <typename INT> EntityIdVector get_node_ids_impl() const
        {
          const Ioss::NodeBlockContainer &node_blocks = region.get_node_blocks();
          assert(node_blocks.size() == 1);

          Ioss::NodeBlock *nb = node_blocks[0];

          std::vector<INT> ids;
          nb->get_field_data("ids", ids);

          EntityIdVector nodeIds;
          for (INT id : ids) {
            nodeIds.push_back(static_cast<EntityId>(id));
          }

          return nodeIds;
        }

        EntityIdVector get_node_ids() const
        {
          if (db_api_int_size() == 4) {
            return get_node_ids_impl<int>();
          }
          else {
            return get_node_ids_impl<int64_t>();
          }
        }

        template <typename INT> bool is_valid_node_impl(EntityId nodeId) const
        {
          EntityIdVector ids  = get_node_ids_impl<INT>();
          auto           iter = std::find(ids.begin(), ids.end(), INT(nodeId));
          return iter != ids.end();
        }

        bool is_valid_node(EntityId nodeId) const
        {
          if (db_api_int_size() == 4) {
            return is_valid_node_impl<int>(nodeId);
          }
          else {
            return is_valid_node_impl<int64_t>(nodeId);
          }
        }

        void verify_num_nodes()
        {
          EntityIdVector ids = get_node_ids();
          EXPECT_EQ(goldNodeIds.size(), ids.size());
        }

        void fill_coordinate_map(const EntityIdVector      &nodeIds,
                                 const std::vector<double> &coordinates)
        {
          std::vector<double>::const_iterator coordIter = coordinates.begin();
          for (const EntityId &nodeId : nodeIds) {
            m_nodalCoords[nodeId] = std::vector<double>(coordIter, coordIter + spatialDim);
            coordIter += spatialDim;
          }
        }

        void fill_coordinates_from_ioss()
        {
          std::vector<double> iossCoordinates;

          const Ioss::NodeBlockContainer &node_blocks = region.get_node_blocks();
          assert(node_blocks.size() == 1);

          Ioss::NodeBlock *nb = node_blocks[0];

          size_t node_count = nb->get_property("entity_count").get_int();

          EntityIdVector nodeIds = get_node_ids();

          size_t numIossNodes =
              field_data_from_ioss<double>(nb, "mesh_model_coordinates", iossCoordinates);
          ThrowRequireWithMsg(node_count == numIossNodes, "Node count mismatch");
          ThrowRequireWithMsg(iossCoordinates.size() == numIossNodes * spatialDim,
                              "Invalid coordinate data size");

          fill_coordinate_map(nodeIds, iossCoordinates);
        }

        const std::vector<double> &operator[](const EntityId nodeId) const
        {
          auto it = m_nodalCoords.find(nodeId);
          return it->second;
        }

        const double *get_nodal_coordinates(const EntityId &nodeId) const
        {
          return Data((*this)[nodeId]);
        }

        void verify_nodal_coordinates(const EntityId &nodeId, const double *goldCoords,
                                      const double *nodalCoords)
        {
          for (unsigned i = 0; i < spatialDim; i++) {
            EXPECT_NEAR(goldCoords[i], nodalCoords[i], 1.0e-9) << error_message(nodeId, i);
          }
        }

        std::string error_message(const EntityId &nodeId, unsigned coordIndex)
        {
          std::stringstream message;
          message << "Proc " << region.get_database()->util().parallel_rank() << ", Node ID "
                  << nodeId << ", coord index " << coordIndex;
          return message.str();
        }

        const Ioss::Region &region;

        const unsigned spatialDim;

        const EntityIdVector                             &goldNodeIds;
        const std::vector<double>                        &goldCoordinates;
        std::unordered_map<EntityId, std::vector<double>> m_nodalCoords;
      };

      unsigned              m_spatialDimension = 3;
      Ioss::PropertyManager m_propertyManager;
      Ioss::DatabaseIO     *m_database = nullptr;
      Ioss::Region         *m_region   = nullptr;
      IossTopologyMapping   m_topologyMapping;
    };

  } // namespace unit_test
} // namespace Iotm

namespace {
  class TestTextMesh : public Iotm::unit_test::TextMeshFixture
  {
  protected:
    TestTextMesh() : TextMeshFixture(3) {}
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

  class TestTextMeshGraph : public Iotm::unit_test::TextMeshFixture
  {
  protected:
    TestTextMeshGraph() : TextMeshFixture(3) {}

    class TextMeshGraph : public SideAdjacencyGraph
    {
    public:
      TextMeshGraph(const TextMeshData &data) : m_data(data) {}

      size_t get_num_elements() const override { return m_data.elementDataVec.size(); }

      int get_element_proc(const size_t elemIndex) const override
      {
        const ElementData &elemData = m_data.elementDataVec[elemIndex];
        return elemData.proc;
      }

      bool element_has_any_node_on_proc(const size_t elemIndex, int proc) const override
      {
        const ElementData &elemData = m_data.elementDataVec[elemIndex];

        for (const EntityId &nodeId : elemData.nodeIds) {
          const std::set<int> &procsForNode = m_data.procs_for_node(nodeId);
          if (procsForNode.count(proc) > 0) {
            return true;
          }
        }

        return false;
      }

      const std::string &get_element_block_name(const size_t elemIndex) const override
      {
        const ElementData &elemData = m_data.elementDataVec[elemIndex];
        return elemData.partName;
      }

      const std::vector<EntityId> &get_element_node_ids(const size_t elemIndex) const override
      {
        const ElementData &elemData = m_data.elementDataVec[elemIndex];
        return elemData.nodeIds;
      }

      const Topology &get_element_topology(const size_t elemIndex) const override
      {
        const ElementData &elemData = m_data.elementDataVec[elemIndex];
        return elemData.topology;
      }

      EntityId get_element_id(const size_t elemIndex) const override
      {
        const ElementData &elemData = m_data.elementDataVec[elemIndex];
        return elemData.identifier;
      }

    private:
      const TextMeshData &m_data;
    };

    void dump_graph(std::ostream &out = std::cout) { m_graph->dump(m_data.elementDataVec, out); }

    void setup_text_mesh_graph(const std::string              &meshDesc,
                               const std::vector<std::string> &selectedBlocks = {},
                               int                             proc = SideAdjacencyGraph::ANY_PROC)
    {
      TextMeshParser parser;
      m_data  = parser.parse(meshDesc);
      m_graph = std::make_shared<TextMeshGraph>(m_data);
      m_graph->create_graph(selectedBlocks, proc);
    }

    void verify_side_adjacency(const std::vector<Adjacency> &goldNeighbors)
    {
      EXPECT_EQ(m_graph->size(), goldNeighbors.size());
      for (size_t i = 0; i < goldNeighbors.size(); ++i) {
        const auto &graphNeighborIndices = (*m_graph)[goldNeighbors[i].elemIndex];
        const auto &goldNeighborIndices  = goldNeighbors[i].neighborIndices;

        unsigned numActualGoldConnections = 0;
        for (const auto &entry : goldNeighborIndices) {
          if (entry.second >= 0) {
            numActualGoldConnections++;
          }
        }

        EXPECT_EQ(numActualGoldConnections, graphNeighborIndices.connections.size());

        for (const auto &entry : goldNeighborIndices) {
          int                           side              = entry.first + 1;
          SideAdjacencyGraph::IndexType neighborElemIndex = entry.second;

          if (neighborElemIndex >= 0) {
            EXPECT_LT(0, graphNeighborIndices.sideReference[side - 1]);
            EXPECT_TRUE(graphNeighborIndices.has_any_connection(side, neighborElemIndex));
          }
          else {
            EXPECT_EQ(0, graphNeighborIndices.sideReference[side - 1]);
            EXPECT_FALSE(graphNeighborIndices.has_any_connection(side));
          }
        }
      }
    }

    TextMeshData                   m_data;
    std::shared_ptr<TextMeshGraph> m_graph;
  };

  using TestTextMeshSkin = TestTextMesh;
} // namespace
