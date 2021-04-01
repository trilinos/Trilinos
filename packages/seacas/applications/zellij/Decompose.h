// Copyright(C) 2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef ZE_Decompose_H
#define ZE_Decompose_H
#include <string>

class Grid;
void decompose_grid(Grid &grid, int ranks, const std::string &method);
#endif
