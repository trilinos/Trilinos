// Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "exo_entity.h"
#include "side_set.h" // for Side_Set

template <typename INT> class ExoII_Read;

template <typename INT> class Node_Set : public Exo_Entity
{
public:
  Node_Set();
  Node_Set(int file_id, size_t id);
  Node_Set(int file_id, size_t id, size_t nnodes, size_t ndfs = 0);
  Node_Set(const Node_Set &)                  = delete;
  const Node_Set &operator=(const Node_Set &) = delete;

  ~Node_Set() override;

  void       apply_map(const std::vector<INT> &node_map);
  const INT *Nodes() const;
  size_t     Node_Id(size_t position) const;
  size_t     Node_Index(size_t position) const;

  const double *Distribution_Factors() const;
  void          Free_Distribution_Factors() const;

private:
  int  Check_State() const override;
  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Nodeset"; }
  const char *short_label() const override { return "nodeset"; }

  void load_nodes(const std::vector<INT> &node_map) const;

  size_t num_dist_factors{0};

  mutable INT    *nodes{nullptr};     // Array.
  mutable INT    *nodeIndex{nullptr}; // An index array which orders the nodelist in sorted order.
  mutable double *dist_factors{nullptr}; // Array.

  friend class ExoII_Read<INT>;
};
