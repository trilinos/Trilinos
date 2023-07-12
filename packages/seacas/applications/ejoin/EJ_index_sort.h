// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once
#include <vector>

template <typename INT>
void index_coord_sort(const std::vector<double> &xyz, std::vector<INT> &index, int axis);

template <typename INT> void index_sort(const std::vector<INT> &ids, std::vector<INT> &index);
