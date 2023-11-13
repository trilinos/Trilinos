/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

namespace Excn {
  // Note that these are used as indices into truth table for ELBK, NSET, SSET, EDBLK, FABLK...
  enum class ObjectType {
    EBLK  = 0,
    NSET  = 1,
    SSET  = 2,
    EDBLK = 3,
    FABLK = 4,
    NODE,
    ELEM,
    EDGE,
    FACE,
    GLOBAL,
    ASSM,
    UNSET
  };
} // namespace Excn
