// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_ChainGenerator.h>
#include <Ioss_Region.h>
#include <SL_SystemInterface.h>
#include <vector>

#pragma once
template <typename INT>
std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                    const std::vector<float> &weights, IOSS_MAYBE_UNUSED INT dummy);
