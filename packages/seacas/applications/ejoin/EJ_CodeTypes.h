// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <Ioss_Region.h>
#include <set>
#include <string>
#include <utility>
#include <vector>

using IntVector = std::vector<int>;
using IdMap     = std::vector<int>;

using RealVector     = std::vector<double>;
using StringIdVector = std::vector<std::pair<std::string, size_t>>;
using StringVector   = std::vector<std::string>;
using RegionVector   = std::vector<Ioss::Region *>;
using Omissions      = std::vector<StringVector>;
