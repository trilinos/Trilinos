// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, interFace
#include "assembly.h"
#include "exodusII.h" // for ex_block, etc
#include "fmt/ostream.h"
#include "smart_assert.h" // for SMART_ASSERT
#include <cstdlib>        // for exit, nullptr
#include <string>         // for string, char_traits

template <typename INT> Assembly<INT>::Assembly() : Exo_Entity() {}

template <typename INT>
Assembly<INT>::Assembly(int file_id, size_t assembly_id) : Exo_Entity(file_id, assembly_id)
{
  SMART_ASSERT(file_id >= 0);
  SMART_ASSERT((int)assembly_id > EX_INVALID_ID);

  initialize(file_id, assembly_id);
}

template <typename INT> EXOTYPE Assembly<INT>::exodus_type() const { return EX_ASSEMBLY; }

template <typename INT> void Assembly<INT>::entity_load_params()
{
  ex_assembly assembly{};
  assembly.id = id_;
  int err     = ex_get_assembly(fileId, &assembly);

  if (err < 0) {
    Error("Assembly<INT>::entity_load_params(): Failed to get assembly"
          " parameters!  Aborting...\n");
  }

  entity_count  = assembly.entity_count;
  assembly_type = assembly.type;
  entities.resize(entity_count);

  // Now fill in the entities list.
  assembly.entity_list = Data(entities);
  err                  = ex_get_assembly(fileId, &assembly);

  if (err < 0) {
    Error("Assembly<INT>::entity_load_params(): Failed to get assembly"
          " parameters!  Aborting...\n");
  }
}

template <typename INT> int Assembly<INT>::Check_State() const
{
  SMART_ASSERT(id_ >= EX_INVALID_ID);
  return 1;
}
template class Assembly<int>;
template class Assembly<int64_t>;
