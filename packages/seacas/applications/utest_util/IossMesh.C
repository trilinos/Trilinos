// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <algorithm>
#include <string>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "IossMesh.h"

#include "Ioss_Assembly.h"
#include "Ioss_CommSet.h"
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

namespace {
  size_t db_api_int_size(const Ioss::GroupingEntity *entity)
  {
    return entity->get_database()->int_byte_size_api();
  }

  template <typename INT>
  void populate_element_data(Ioss::Region                                  *region,
                             std::vector<utest_util::IossElementData>      &elemData,
                             std::vector<utest_util::IossElementBlockData> &elemBlockData)
  {
    elemData.clear();
    size_t localId = 0;

    for (Ioss::ElementBlock *entity : region->get_element_blocks()) {
      std::vector<INT> elem_ids;
      std::vector<INT> connectivity;
      std::vector<INT> connectivity_raw;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);
      entity->get_field_data("connectivity_raw", connectivity_raw);

      EXPECT_TRUE(connectivity.size() == connectivity_raw.size());

      utest_util::IossElementBlockData elemBlock;
      elemBlock.block = entity;

      if (elem_ids.empty()) {
        elemBlockData.push_back(elemBlock);
        continue;
      }

      const Ioss::ElementTopology *topology = entity->topology();

      size_t element_count  = elem_ids.size();
      int    nodes_per_elem = topology->number_nodes();

      for (size_t i = 0; i < element_count; ++i, ++localId) {
        INT *conn     = &connectivity[i * nodes_per_elem];
        INT *conn_raw = &connectivity_raw[i * nodes_per_elem];

        INT id = elem_ids[i];

        utest_util::IossElementData elem;
        elem.id       = id;
        elem.topology = topology;
        elem.block    = entity;
        elem.localId  = localId;

        elemBlock.localIds.push_back(localId);

        for (int j = 0; j < nodes_per_elem; j++) {
          elem.nodeIds.push_back(conn[j]);

          INT zeroBasedLocalConn = conn_raw[j] - 1;
          assert(zeroBasedLocalConn >= 0);
          elem.localConnectivity.push_back(zeroBasedLocalConn);
        }

        elemData.push_back(elem);
      }

      elemBlockData.push_back(elemBlock);
    }
  }

  template <typename INT>
  void populate_node_data(Ioss::Region *region, std::vector<utest_util::IossNodeData> &nodeData)
  {
    nodeData.clear();
    size_t localId = 0;

    Ioss::NodeBlock *nb = region->get_node_block("nodeblock_1");
    EXPECT_TRUE(nb != nullptr);

    std::vector<INT> node_ids;

    nb->get_field_data("ids", node_ids);

    size_t node_count = node_ids.size();

    for (size_t i = 0; i < node_count; ++i, ++localId) {
      INT id = node_ids[i];

      utest_util::IossNodeData node;
      node.id      = id;
      node.localId = localId;

      nodeData.push_back(node);
    }
  }
} // namespace

namespace utest_util {

  IossMesh::IossMesh() : m_communicator(Ioss::ParallelUtils::comm_world()) {}

  IossMesh::IossMesh(Ioss_MPI_Comm comm) : m_communicator(comm) {}

  int IossMesh::get_parallel_size() { return Ioss::ParallelUtils(get_comm()).parallel_size(); }

  int IossMesh::get_parallel_rank() { return Ioss::ParallelUtils(get_comm()).parallel_rank(); }

  void IossMesh::fill_mesh(const std::string &meshDesc, const std::string &regionName)
  {
    std::pair<std::string, std::string> result = get_database_type_and_filename(meshDesc);

    std::string type     = result.first;
    std::string filename = result.second;

    create_database(filename, type);
    create_ioss_region(regionName);

    fill_element_data();
    fill_node_data();

    m_numGlobalNodes         = get_region()->get_optional_property("global_node_count", 1);
    m_numGlobalElements      = get_region()->get_optional_property("global_element_count", 1);
    m_numGlobalElementBlocks = get_region()->get_element_blocks().size();

    Ioss::NodeBlock *nb = m_region->get_node_block("nodeblock_1");
    EXPECT_TRUE(nb != nullptr);
    m_spatialDimension = nb->get_property("component_degree").get_int();
  }

  void IossMesh::fill_mesh(const std::string &meshDesc) { fill_mesh(meshDesc, "input_region"); }

  void IossMesh::create_ioss_region(const std::string &regionName)
  {
    if (m_region == nullptr) {
      EXPECT_TRUE(m_database != nullptr);
      m_region = new Ioss::Region(m_database, regionName);
      EXPECT_TRUE(m_region != nullptr);
    }
  }

  void IossMesh::filename_substitution(std::string &filename)
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

  void IossMesh::create_database(const std::string &fileName, const std::string &meshType,
                                 Ioss::DatabaseUsage dbUsage)
  {
    if (m_database == nullptr) {
      std::string meshFileName = fileName;
      filename_substitution(meshFileName);
      m_database =
          Ioss::IOFactory::create(meshType, meshFileName, dbUsage, get_comm(), m_propertyManager);
      EXPECT_TRUE(m_database != nullptr);
      EXPECT_TRUE(m_database->ok(true));
    }
  }

  std::pair<std::string, std::string>
  IossMesh::get_database_type_and_filename(const std::string &meshDesc)
  {
    std::string type;
    std::string filename;

    size_t colon = meshDesc.find(':');
    if (colon != std::string::npos && colon > 0) {
      type     = meshDesc.substr(0, colon);
      filename = meshDesc.substr(colon + 1);
    }
    else {
      type     = "exodusII";
      filename = meshDesc;
    }

    return std::make_pair(type, filename);
  }

  void IossMesh::fill_element_data()
  {
    EXPECT_TRUE(m_database != nullptr) << "Database has not been read";
    EXPECT_TRUE(m_region != nullptr) << "Region has not been constructed";

    bool ints64bit = db_api_int_size(m_region) == 8;

    if (ints64bit) {
      populate_element_data<int64_t>(m_region, m_elementData, m_elementBlockData);
    }
    else {
      populate_element_data<int>(m_region, m_elementData, m_elementBlockData);
    }
  }

  IossElementData IossMesh::get_local_element(size_t index) const
  {
    if (index >= m_elementData.size())
      return IossElementData();

    return m_elementData[index];
  }

  IossElementData IossMesh::get_global_element(EntityId id) const
  {
    auto lowerBound = std::lower_bound(m_elementData.begin(), m_elementData.end(), id);

    if (lowerBound != m_elementData.end() && lowerBound->id == id) {
      return *lowerBound;
    }

    return IossElementData();
  }

  void IossMesh::fill_node_data()
  {
    EXPECT_TRUE(m_database != nullptr) << "Database has not been read";
    EXPECT_TRUE(m_region != nullptr) << "Region has not been constructed";

    bool ints64bit = db_api_int_size(m_region) == 8;

    if (ints64bit) {
      populate_node_data<int64_t>(m_region, m_nodeData);
    }
    else {
      populate_node_data<int>(m_region, m_nodeData);
    }
  }

  IossNodeData IossMesh::get_local_node(size_t index) const
  {
    if (index >= m_nodeData.size())
      return IossNodeData();

    return m_nodeData[index];
  }

  IossNodeData IossMesh::get_global_node(EntityId id) const
  {
    auto lowerBound = std::lower_bound(m_nodeData.begin(), m_nodeData.end(), id);

    if (lowerBound != m_nodeData.end() && lowerBound->id == id) {
      return *lowerBound;
    }

    return IossNodeData();
  }
} // namespace utest_util
