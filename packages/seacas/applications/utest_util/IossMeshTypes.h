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
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_StandardElementTypes.h"

namespace utest_util {

  using EntityId = int64_t;

  struct IossNodeData
  {
    EntityId id{0};

    // Local index into vector held by IossMesh
    size_t localId{std::numeric_limits<size_t>::max()};

    operator EntityId() const { return id; }

    bool is_valid() const { return (id > 0) && (localId != std::numeric_limits<size_t>::max()); }
  };

  struct IossNodeDataLess
  {
    bool operator()(const IossNodeData &lhs, const IossNodeData &rhs) { return lhs.id < rhs.id; };

    bool operator()(const IossNodeData &lhs, const EntityId rhs) { return lhs.id < rhs; };

    bool operator()(const EntityId lhs, const IossNodeData &rhs) { return lhs < rhs.id; };

    bool operator()(const EntityId lhs, const EntityId rhs) { return lhs < rhs; };
  };

  struct IossElementData
  {
    EntityId id{0};

    // Local index into vector held by IossMesh
    size_t localId{std::numeric_limits<size_t>::max()};

    const Ioss::ElementTopology *topology{nullptr};
    const Ioss::ElementBlock    *block{nullptr};

    // Global node id connectivity
    std::vector<EntityId> nodeIds{};

    // Zero based local connectivity
    std::vector<EntityId> localConnectivity{};

    operator EntityId() const { return id; }

    bool is_valid() const
    {
      return (id > 0) && (localId != std::numeric_limits<size_t>::max()) && (topology != nullptr) &&
             (block != nullptr) && !nodeIds.empty();
    }
  };

  struct IossElementDataLess
  {
    bool operator()(const IossElementData &lhs, const IossElementData &rhs)
    {
      return lhs.id < rhs.id;
    };

    bool operator()(const IossElementData &lhs, const EntityId rhs) { return lhs.id < rhs; };

    bool operator()(const EntityId lhs, const IossElementData &rhs) { return lhs < rhs.id; };

    bool operator()(const EntityId lhs, const EntityId rhs) { return lhs < rhs; };
  };

  struct IossElementBlockData
  {
    const Ioss::ElementBlock *block{nullptr};
    std::vector<size_t>       localIds;
  };

} // namespace utest_util
