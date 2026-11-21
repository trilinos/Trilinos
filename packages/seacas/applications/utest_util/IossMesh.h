// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif
#include <string>
#include <unordered_map>

#include <memory>
#include <string>
#if defined(_WIN32) && !defined(__MINGW32__)
#include <string.h>
#define strcasecmp  _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif
#include <limits>
#include <vector>

#include "Ioss_CodeTypes.h"

#include "Ionit_Initializer.h"

#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_StandardElementTypes.h"

#include "IossMeshTypes.h"

namespace utest_util {

  class IossMesh
  {
  public:
    IossMesh();
    IossMesh(Ioss_MPI_Comm comm);

    virtual ~IossMesh() { delete m_region; }

    int get_parallel_size();

    int get_parallel_rank();

    Ioss_MPI_Comm get_comm() const { return m_communicator; }

    void fill_mesh(const std::string &meshDesc);
    void fill_mesh(const std::string &meshDesc, const std::string &regionName);

    Ioss::DatabaseIO *get_database() { return m_database; }

    Ioss::Region *get_region() { return m_region; }

    int64_t get_num_global_nodes() const { return m_numGlobalNodes; }
    int64_t get_num_global_elements() const { return m_numGlobalElements; }
    int64_t get_num_global_element_blocks() const { return m_numGlobalElementBlocks; }

    size_t get_num_local_elements() const { return m_elementData.size(); }
    size_t get_num_local_nodes() const { return m_nodeData.size(); };

    IossElementData get_local_element(size_t index) const;
    IossElementData get_global_element(EntityId id) const;

    IossNodeData get_local_node(size_t index) const;
    IossNodeData get_global_node(EntityId id) const;

    unsigned get_spatial_dimension() const { return m_spatialDimension; }

  protected:
    Ioss_MPI_Comm         m_communicator;
    Ioss::PropertyManager m_propertyManager;
    Ioss::DatabaseIO     *m_database = nullptr;
    Ioss::Region         *m_region   = nullptr;

    int64_t m_numGlobalNodes{0};
    int64_t m_numGlobalElements{0};
    int64_t m_numGlobalElementBlocks{0};

    unsigned m_spatialDimension{0};

    std::vector<IossNodeData>         m_nodeData;
    std::vector<IossElementData>      m_elementData;
    std::vector<IossElementBlockData> m_elementBlockData;

  protected:
    void filename_substitution(std::string &filename);

    std::pair<std::string, std::string> get_database_type_and_filename(const std::string &meshDesc);

    void create_ioss_region(const std::string &regionName = "input_region");

    void create_database(const std::string &fileName, const std::string &meshType,
                         Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL);

    void fill_element_data();

    void fill_node_data();
  };

} // namespace utest_util
