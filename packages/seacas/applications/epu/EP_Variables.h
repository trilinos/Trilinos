/*
 * Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "EP_ObjectType.h"
#include <cstring>
#include <smart_assert.h>
#include <string>
#include <vector>

namespace Excn {
  enum class InOut { IN = 1, OUT = 2 };

  using IntVector = std::vector<int>;

  struct Variables
  {
    explicit Variables(ObjectType otype) : objectType(otype)
    {
      SMART_ASSERT(otype == Excn::ObjectType::EBLK || otype == Excn::ObjectType::NSET ||
                   otype == Excn::ObjectType::SSET || otype == Excn::ObjectType::NODE ||
                   otype == Excn::ObjectType::EDBLK || otype == Excn::ObjectType::FABLK ||
                   otype == Excn::ObjectType::GLOBAL);
    }

    int count(InOut in_out = InOut::IN) const
    {
      int ret_val = 0;
      switch (in_out) {
      case InOut::IN: ret_val = index_.size() - (addProcessorId ? 1 : 0); break;
      case InOut::OUT: ret_val = outputCount; break;
      }
      return ret_val;
    }

    const char *label() const
    {
      switch (objectType) {
      case Excn::ObjectType::EBLK: return "element";
      case Excn::ObjectType::NSET: return "nodeset";
      case Excn::ObjectType::GLOBAL: return "global";
      case Excn::ObjectType::NODE: return "nodal";
      case Excn::ObjectType::SSET: return "sideset";
      case Excn::ObjectType::EDBLK: return "edgeblock";
      case Excn::ObjectType::FABLK: return "faceblock";
      default: return "UNKNOWN";
      }
    }

    ex_entity_type type() const
    {
      switch (objectType) {
      case Excn::ObjectType::EBLK: return EX_ELEM_BLOCK;
      case Excn::ObjectType::NSET: return EX_NODE_SET;
      case Excn::ObjectType::SSET: return EX_SIDE_SET;
      case Excn::ObjectType::NODE: return EX_NODAL;
      case Excn::ObjectType::GLOBAL: return EX_GLOBAL;
      case Excn::ObjectType::EDBLK: return EX_EDGE_BLOCK;
      case Excn::ObjectType::FABLK: return EX_FACE_BLOCK;
      default: return EX_INVALID;
      }
    }

    bool add_processor_id() const { return addProcessorId; }

    ObjectType  objectType;
    int         outputCount{0};
    bool        addProcessorId{false};
    IntVector   index_{};
    std::string type_{};
  };
} // namespace Excn
