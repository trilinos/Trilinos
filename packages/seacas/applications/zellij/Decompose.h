// Copyright(C) 2021, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <string>

class Grid;
void decompose_grid(Grid &grid, int ranks, const std::string &method);
