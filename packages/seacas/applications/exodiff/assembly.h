// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "exo_entity.h"
#include <iostream>

#include <string>

template <typename INT> class ExoII_Read;

template <typename INT> class Assembly : public Exo_Entity
{
public:
  Assembly();
  Assembly(int file_id, size_t assembly_id);
  Assembly(const Assembly &)                  = delete;
  const Assembly &operator=(const Assembly &) = delete;

  ex_entity_type                   Type() const { return assembly_type; }
  const std::vector<ex_entity_id> &Entities() const { return entities; }
  size_t                           Size() const override { return entities.size(); }

private:
  int  Check_State() const override;
  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Assembly"; }
  const char *short_label() const override { return "assembly"; }

  int                       entity_count{0};
  ex_entity_type            assembly_type{EX_INVALID};
  std::vector<ex_entity_id> entities{};

  friend class ExoII_Read<INT>;
};
