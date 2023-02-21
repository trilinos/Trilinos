// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "EJ_match_xyz.h"

template <typename INT>
void match_node_xyz(RegionVector &part_mesh, double tolerance, std::vector<INT> &global_node_map,
                    std::vector<INT> &local_node_map);
