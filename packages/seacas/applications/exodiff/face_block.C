// Copyright(C) 1999-, 20212021,  National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, etc
#include "exodusII.h"           // for ex_set, etc
#include "face_block.h"
#include "iqsort.h"       // for index_qsort
#include "smart_assert.h" // for SMART_ASSERT
#include <cstdlib>        // for exit
#include <vector>         // for vector

template <typename INT> Face_Block<INT>::Face_Block() : Exo_Entity() {}

template <typename INT>
Face_Block<INT>::Face_Block(int file_id, size_t id) : Exo_Entity(file_id, id)
{
  SMART_ASSERT((int)id != EX_INVALID_ID);
}

template <typename INT>
Face_Block<INT>::Face_Block(int file_id, size_t id, size_t ne) : Exo_Entity(file_id, id, ne)
{
  SMART_ASSERT(id > 0);
}

template <typename INT> Face_Block<INT>::~Face_Block() { SMART_ASSERT(Check_State()); }

template <typename INT> EXOTYPE Face_Block<INT>::exodus_type() const { return EX_FACE_BLOCK; }

template <typename INT> void Face_Block<INT>::entity_load_params()
{
  int      num_attr;
  ex_block block{};
  block.id   = id_;
  block.type = EX_FACE_BLOCK;
  int err    = ex_get_block_param(fileId, &block);

  if (err < 0) {
    Error("Face_Block<INT>::entity_load_params(): Failed to get face"
          " block parameters!  Aborting...\n");
  }

  numEntity          = block.num_entry;
  num_faces_per_elmt = block.num_faces_per_entry;
  num_attr           = block.num_attribute;
  elmt_type          = block.topology;

  if (num_faces_per_elmt < 0 || num_attr < 0) {
    Error(fmt::format(
		      fmt::runtime("Face_Block<INT>::entity_load_params(): Data appears corrupt for face block {}!\n"
				   "\tnum elmts          = {}\n"
				   "\tnum faces per elmt = {}\n"
				   "\tnum attributes     = {}\n"
				   " ... Aborting...\n"),
        fmt::group_digits(numEntity), num_faces_per_elmt, num_attr));
  }
}

template <typename INT> size_t Face_Block<INT>::Face_Index(size_t position) const
{
  SMART_ASSERT(position < numEntity);
  return position;
}

template <typename INT> int Face_Block<INT>::Check_State() const
{
  SMART_ASSERT(id_ >= EX_INVALID_ID);
  SMART_ASSERT(!(id_ == EX_INVALID_ID && numEntity > 0));

  return 1;
}

template class Face_Block<int>;
template class Face_Block<int64_t>;
