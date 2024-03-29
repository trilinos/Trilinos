// Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, etc
#include "exodusII.h"           // for ex_set, etc
#include "fmt/format.h"
#include "iqsort.h" // for index_qsort
#include "node_set.h"
#include "smart_assert.h" // for SMART_ASSERT
#include <cstdlib>        // for exit
#include <vector>         // for vector

template <typename INT> Node_Set<INT>::Node_Set() : Exo_Entity() {}

template <typename INT> Node_Set<INT>::Node_Set(int file_id, size_t id) : Exo_Entity(file_id, id) {}

template <typename INT>
Node_Set<INT>::Node_Set(int file_id, size_t id, size_t nnodes, size_t ndfs)
    : Exo_Entity(file_id, id, nnodes), num_dist_factors(ndfs)
{
}

template <typename INT> Node_Set<INT>::~Node_Set()
{
  delete[] nodes;
  delete[] nodeIndex;
  delete[] dist_factors;
}

template <typename INT> EXOTYPE Node_Set<INT>::exodus_type() const { return EX_NODE_SET; }

template <typename INT> const INT *Node_Set<INT>::Nodes() const
{
  // See if already loaded...
  if (!nodes) {
    std::vector<INT> tmp;
    load_nodes(tmp);
  }
  return nodes;
}

template <typename INT> size_t Node_Set<INT>::Node_Id(size_t position) const
{
  if (numEntity == 0) {
    return 0;
  }

  // See if already loaded...
  if (!nodes) {
    std::vector<INT> tmp;
    load_nodes(tmp);
  }
  SMART_ASSERT(position < numEntity);
  return nodes[nodeIndex[position]];
}

template <typename INT> size_t Node_Set<INT>::Node_Index(size_t position) const
{
  if (numEntity == 0) {
    return 0;
  }

  // See if already loaded...
  if (!nodes) {
    std::vector<INT> tmp;
    load_nodes(tmp);
  }
  SMART_ASSERT(position < numEntity);
  SMART_ASSERT(nodeIndex != nullptr);
  return nodeIndex[position];
}

template <typename INT> void Node_Set<INT>::apply_map(const std::vector<INT> &node_map)
{
  SMART_ASSERT(!node_map.empty());
  if (nodes != nullptr) {
    delete[] nodes;
    nodes = nullptr;
    delete[] nodeIndex;
    nodeIndex = nullptr;
  }
  load_nodes(node_map);
}

template <typename INT> void Node_Set<INT>::load_nodes(const std::vector<INT> &node_map) const
{
  if (numEntity > 0) {
    nodes = new INT[numEntity];
    SMART_ASSERT(nodes != nullptr);
    nodeIndex = new INT[numEntity];
    SMART_ASSERT(nodeIndex != nullptr);
    ex_get_set(fileId, EX_NODE_SET, id_, nodes, nullptr);

    if (!node_map.empty()) {
      for (size_t i = 0; i < numEntity; i++) {
        nodes[i] = 1 + node_map[nodes[i] - 1];
      }
    }

    for (size_t i = 0; i < numEntity; i++) {
      nodeIndex[i] = i;
    }
    if (interFace.nsmap_flag) {
      index_qsort(nodes, nodeIndex, numEntity);
    }
  }
}

template <typename INT> const double *Node_Set<INT>::Distribution_Factors() const
{
  if ((dist_factors == nullptr) && num_dist_factors > 0) {
    dist_factors = new double[num_dist_factors];
    SMART_ASSERT(dist_factors != nullptr);
    ex_get_set_dist_fact(fileId, EX_NODE_SET, id_, dist_factors);
  }
  return dist_factors;
}

template <typename INT> void Node_Set<INT>::Free_Distribution_Factors() const
{
  if (dist_factors) {
    delete[] dist_factors;
    dist_factors = nullptr;
  }
}

template <typename INT> int Node_Set<INT>::Check_State() const
{
  SMART_ASSERT(id_ >= EX_INVALID_ID);
  SMART_ASSERT(!(id_ == EX_INVALID_ID && numEntity > 0));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && num_dist_factors > 0));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && nodes));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && dist_factors));

  return 1;
}

template <typename INT> void Node_Set<INT>::entity_load_params()
{
  std::vector<ex_set> sets(1);
  sets[0].id                       = id_;
  sets[0].type                     = EX_NODE_SET;
  sets[0].entry_list               = nullptr;
  sets[0].extra_list               = nullptr;
  sets[0].distribution_factor_list = nullptr;

  int err = ex_get_sets(fileId, 1, Data(sets));

  if (err < 0) {
    Error(fmt::format("Failed to get nodeset parameters for nodeset {}. !  Aborting...\n", id_));
  }

  numEntity        = sets[0].num_entry;
  num_dist_factors = sets[0].num_distribution_factor;
}

template class Node_Set<int>;
template class Node_Set<int64_t>;
