// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "IossMeshTypes.h"

namespace utest_util {

  struct EntityProc
  {
    EntityProc(EntityId id_, int proc_) : id(id_), proc(proc_) {}
    EntityProc() : id(0), proc(-1) {}

    EntityId id{0};
    int      proc{-1};
  };

  class EntityProcLess
  {
  public:
    EntityProcLess() = default;
    bool operator()(const EntityProc &lhs, const EntityProc &rhs) const { return lhs.id < rhs.id; }

    bool operator()(const EntityProc &lhs, const EntityId &rhs) const { return lhs.id < rhs; }

    bool operator()(const EntityId &lhs, const EntityProc &rhs) const { return lhs < rhs.id; }

    bool operator()(const EntityId &lhs, const EntityId &rhs) const { return lhs < rhs; }
  };

  class EntityProcEqual
  {
  public:
    EntityProcEqual() = default;
    bool operator()(const EntityProc &lhs, const EntityProc &rhs) const
    {
      return lhs.id == rhs.id && lhs.proc == rhs.proc;
    }
  };
} // namespace utest_util
