// Copyright(C) 1999-2020, 2022, 2023, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <CJ_ObjectType.h>
#include <cstring>
#include <smart_assert.h>
#include <string>
#include <vector>

namespace Excn {
  enum class InOut { IN_ = 1, OUT_ = 2 };

  using IntVector = std::vector<int>;

  struct Variables
  {
    Variables(ObjectType otype, int part_count, bool arg_add_status = false)
        : objectType(otype), addStatus(arg_add_status)
    {
      inputIndex_.resize(part_count);
      SMART_ASSERT(otype == ObjectType::EBLK || otype == ObjectType::NSET ||
                   otype == ObjectType::SSET || otype == ObjectType::NODE ||
                   otype == ObjectType::GLOBAL);
    }

    size_t count() const { return names_.size() + (addStatus ? 1 : 0); }

    size_t in_count(int part) const { return inputIndex_[part].size(); }

    const char *label() const
    {
      switch (objectType) {
      case ObjectType::EBLK: return "element";
      case ObjectType::NSET: return "nodeset";
      case ObjectType::GLOBAL: return "global";
      case ObjectType::NODE: return "nodal";
      case ObjectType::SSET: return "sideset";
      default: return "UNKNOWN";
      }
    }

    ex_entity_type type() const
    {
      switch (objectType) {
      case ObjectType::EBLK: return EX_ELEM_BLOCK;
      case ObjectType::NSET: return EX_NODE_SET;
      case ObjectType::SSET: return EX_SIDE_SET;
      case ObjectType::NODE: return EX_NODAL;
      case ObjectType::GLOBAL: return EX_GLOBAL;
      default: return EX_INVALID;
      }
    }

    IntVector &input_index(int part) { return inputIndex_[part]; }

    bool add_status() const { return addStatus; }

    ObjectType             objectType;
    bool                   addStatus;
    std::vector<IntVector> inputIndex_{};
    StringVector           names_{};
  };
} // namespace Excn
