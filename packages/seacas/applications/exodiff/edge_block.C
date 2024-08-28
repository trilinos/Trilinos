// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, etc
#include "edge_block.h"
#include "exodusII.h"     // for ex_set, etc
#include "iqsort.h"       // for index_qsort
#include "smart_assert.h" // for SMART_ASSERT
#include <cstdlib>        // for exit
#include <vector>         // for vector

template <typename INT> Edge_Block<INT>::Edge_Block() : Exo_Entity() {}

template <typename INT>
Edge_Block<INT>::Edge_Block(int file_id, size_t id) : Exo_Entity(file_id, id)
{
  SMART_ASSERT((int)id != EX_INVALID_ID);
}

template <typename INT>
Edge_Block<INT>::Edge_Block(int file_id, size_t id, size_t ne) : Exo_Entity(file_id, id, ne)
{
  SMART_ASSERT(id > 0);
}

template <typename INT> Edge_Block<INT>::~Edge_Block() { SMART_ASSERT(Check_State()); }

template <typename INT> EXOTYPE Edge_Block<INT>::exodus_type() const { return EX_EDGE_BLOCK; }

template <typename INT> void Edge_Block<INT>::entity_load_params()
{
  int      num_attr;
  ex_block block{};
  block.id   = id_;
  block.type = EX_EDGE_BLOCK;
  int err    = ex_get_block_param(fileId, &block);

  if (err < 0) {
    Error("Edge_Block<INT>::entity_load_params(): Failed to get edge"
          " block parameters!  Aborting...\n");
  }

  numEntity          = block.num_entry;
  num_edges_per_elmt = block.num_edges_per_entry;
  num_attr           = block.num_attribute;
  elmt_type          = block.topology;

  if (num_edges_per_elmt < 0 || num_attr < 0) {
    Error(fmt::format(
		      fmt::runtime("Edge_Block<INT>::entity_load_params(): Data appears corrupt for edge block {}!\n"
				   "\tnum elmts          = {}\n"
				   "\tnum edges per elmt = {}\n"
				   "\tnum attributes     = {}\n"
				   " ... Aborting...\n"),
		      fmt::group_digits(numEntity), num_edges_per_elmt, num_attr));
  }
}

template <typename INT> size_t Edge_Block<INT>::Edge_Index(size_t position) const
{
  SMART_ASSERT(position < numEntity);
  return position;
}

template <typename INT> int Edge_Block<INT>::Check_State() const
{
  SMART_ASSERT(id_ >= EX_INVALID_ID);
  SMART_ASSERT(!(id_ == EX_INVALID_ID && numEntity > 0));

  return 1;
}

template class Edge_Block<int>;
template class Edge_Block<int64_t>;
