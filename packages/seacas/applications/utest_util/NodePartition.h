// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <gtest/gtest.h>
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
#include <ostream>
#include <vector>

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "Ioss_CodeTypes.h" // for Ioss_MPI_Comm
#include "Ioss_DBUsage.h"   // for DatabaseUsage
#include "Ioss_ElementTopology.h"
#include "Ioss_IOFactory.h" // for IOFactory
#include "Ioss_MeshCopyOptions.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h" // for IOFactory

#include "IossMesh.h"
#include "Partition.h"
#include "PartitionTypes.h"

namespace utest_util {

  template <typename INT> class NodePartition : public Partition
  {
  public:
    NodePartition(IossMesh *mesh, const std::vector<EntityProc> &procAssign, const int nProc);

    void write_nemesis_data(int exoid) const override;

    size_t get_num_local_entities() const override;

    int get_processor(size_t localId) const override;

    size_t api_size() const override { return sizeof(INT); }

    std::string type() const { return "nodal"; }

    virtual ~NodePartition() = default;

  private:
    /* Nodal */
    std::vector<std::vector<INT>> int_nodes{};
    std::vector<std::vector<INT>> bor_nodes{};
    std::vector<std::vector<INT>> ext_nodes{};
    std::vector<std::vector<INT>> ext_procs{};

    /* Elemental */
    std::vector<std::vector<INT>> int_elems{};
    std::vector<std::vector<INT>> bor_elems{};

  private:
    void categorize_elements(const std::vector<std::vector<INT>> &sur_elem, const int max_nsur);
    void categorize_nodes(const std::vector<std::vector<INT>> &sur_elem, const int max_nsur);

    void generate_partition_maps();
  };
} // namespace utest_util
