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

template <typename INT> class Exo_Block : public Exo_Entity
{
public:
  Exo_Block();
  Exo_Block(int file_id, size_t exo_block_id);
  Exo_Block(int file_id, size_t id, const char *type, size_t num_e, size_t num_npe);
  Exo_Block(const Exo_Block &)                  = delete;
  const Exo_Block &operator=(const Exo_Block &) = delete;

  std::string Load_Connectivity();
  std::string Free_Connectivity();

  // Access functions:
  const std::string &Element_Type() const { return elmt_type; }
  size_t             Num_Nodes_per_Element() const { return num_nodes_per_elmt; }

  // Block description access functions:
  const std::vector<INT> &Connectivity() const { return conn; }  // 1-offset connectivity
  const INT              *Connectivity(size_t elmt_index) const; // 1-offset connectivity

  // Misc:
  int64_t offset() const { return offset_; }
  void    offset(int64_t off) { offset_ = off; }

private:
  int  Check_State() const override;
  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Element Block"; }
  const char *short_label() const override { return "block"; }

  std::string      elmt_type{};
  int              num_nodes_per_elmt{-1};
  int64_t          offset_{0};
  std::vector<INT> conn{}; // Array; holds a matrix, num_elmts by num_nodes_per_elmt.

  friend class ExoII_Read<INT>;
};
