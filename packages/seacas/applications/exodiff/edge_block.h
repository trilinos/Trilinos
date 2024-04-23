// Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "exo_entity.h"
#include <iostream>

template <typename INT> class ExoII_Read;

template <typename INT> class Edge_Block : public Exo_Entity
{
public:
  Edge_Block();
  Edge_Block(int file_id, size_t id);
  Edge_Block(int file_id, size_t id, size_t ne);
  Edge_Block(const Edge_Block &)                  = delete; // Not written.
  const Edge_Block &operator=(const Edge_Block &) = delete; // Not written.

  ~Edge_Block() override;

  size_t Edge_Index(size_t position) const;

private:
  int  Check_State() const override;
  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Edgeblock"; }
  const char *short_label() const override { return "edgeblock"; }

  std::string elmt_type{};
  int         num_edges_per_elmt{-1};

  friend class ExoII_Read<INT>;
};
