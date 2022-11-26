// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_SNAPINDEPENDENTSETFINDER_H
#define AKRI_SNAPINDEPENDENTSETFINDER_H

#include <Akri_QualityMetric.hpp>
#include <Akri_SnapInfo.hpp>
#include <stk_emend/independent_set/IndependentSetFinder.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace krino
{

typedef independent_set::IndependentSetFinder<SnapInfo> SnapInfoIndependentSetFinder;

inline
std::vector<SnapInfo> find_snap_info_independent_sets(const std::vector<SnapInfo> &allSnapInfos,
    const QualityMetric & qualityMetric,
    stk::ParallelMachine comm)
{
    SnapInfo::Comparator comparator {qualityMetric};
    return SnapInfoIndependentSetFinder::find_independent_set(allSnapInfos, comparator, comm);
}

}

#endif /* AKRI_SNAPINDEPENDENTSETFINDER_H */
