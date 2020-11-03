// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef SEACAS_CodeTypes_H
#define SEACAS_CodeTypes_H

#include <Ioss_Region.h>
#include <set>
#include <string>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
#ifdef _WIN64
#define ssize_t __int64
#else
#define ssize_t long
#endif
#endif

using IntVector = std::vector<int>;
using IdMap     = std::vector<int>;

using RealVector = std::vector<double>;
typedef std::vector<std::pair<std::string, size_t>> StringIdVector;
using StringVector = std::vector<std::string>;
using RegionVector = std::vector<Ioss::Region *>;
using Omissions    = std::vector<StringVector>;
#endif
