// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef SEACAS_Variables_H
#define SEACAS_Variables_H
#include <CJ_ObjectType.h>
#include <cstring>
#include <smart_assert.h>
#include <string>
#include <vector>

namespace Excn {
  enum InOut { IN_ = 1, OUT_ = 2 };

  using IntVector = std::vector<int>;

  struct Variables
  {
    explicit Variables(ObjectType otype, bool arg_add_status = false)
        : objectType(otype), outputCount(0), addStatus(arg_add_status)
    {
      SMART_ASSERT(otype == EBLK || otype == NSET || otype == SSET || otype == NODE ||
                   otype == GLOBAL);
    }

    int count(InOut in_out = IN_) const
    {
      int ret_val = 0;
      switch (in_out) {
      case IN_: ret_val = index_.size() - (addStatus ? 1 : 0); break;
      case OUT_: ret_val = outputCount; break;
      }
      return ret_val;
    }

    const char *label() const
    {
      switch (objectType) {
      case EBLK: return "element";
      case NSET: return "nodeset";
      case GLOBAL: return "global";
      case NODE: return "nodal";
      case SSET: return "sideset";
      default: return "UNKNOWN";
      }
    }

    ex_entity_type type() const
    {
      switch (objectType) {
      case EBLK: return EX_ELEM_BLOCK;
      case NSET: return EX_NODE_SET;
      case SSET: return EX_SIDE_SET;
      case NODE: return EX_NODAL;
      case GLOBAL: return EX_GLOBAL;
      default: return EX_INVALID;
      }
    }

    bool add_status() const { return addStatus; }

    ObjectType  objectType;
    int         outputCount;
    bool        addStatus;
    IntVector   index_{};
    std::string type_{};
  };
} // namespace Excn

#endif
