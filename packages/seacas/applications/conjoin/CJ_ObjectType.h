// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef ObjectType_H
#define ObjectType_H

namespace Excn {
  // Note that these are used as indices into truth table for ELBK, NSET, SSET...
  enum ObjectType { EBLK = 0, NSET = 1, SSET = 2, NODE, ELEM, GLOBAL, TIME, DIM };
} // namespace Excn

#endif
