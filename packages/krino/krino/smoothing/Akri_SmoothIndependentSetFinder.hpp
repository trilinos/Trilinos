// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_SMOOTHINDEPENDENTSETFINDER_H
#define AKRI_SMOOTHINDEPENDENTSETFINDER_H

#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothInfo.hpp>
#include <stk_emend/independent_set/IndependentSetFinder.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace krino
{

using SmoothInfoIndependentSetFinder = independent_set::IndependentSetFinder<SmoothInfo>;

inline
std::vector<SmoothInfo> find_smooth_info_independent_sets(const std::vector<SmoothInfo> &allSmoothInfos,
    stk::ParallelMachine comm)
{
  SmoothInfo::Comparator comparator;
  return SmoothInfoIndependentSetFinder::find_independent_set(allSmoothInfos, comparator, comm);
}

}

#endif /* AKRI_SNAPINDEPENDENTSETFINDER_H */
