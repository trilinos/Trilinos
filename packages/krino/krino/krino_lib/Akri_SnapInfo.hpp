// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_SNAPINFO_H
#define AKRI_SNAPINFO_H

#include <Akri_QualityMetric.hpp>
#include <stk_math/StkVector.hpp>
#include <vector>

namespace krino
{

class SnapInfo
{
public:
    typedef size_t ExclusionIdentifierType;
    typedef std::array<size_t,2> GlobalId;

    SnapInfo(const size_t nodeGlobalId,
        const size_t intersectionPointIndex,
        const stk::math::Vector3d & nodeLocation,
        const int owner,
        const std::vector<int> & procsThatNeedToKnowAboutThisInfo,
        const std::vector<size_t> & globalIdsOfSnapNodeElems,
        const double& postSnapQuality,
        const stk::math::Vector3d & snapLocation,
        const int snapRank)
    : mUniqueId(make_unique_id(nodeGlobalId, intersectionPointIndex)),
      mOwner{owner},
      mProcsThatNeedToKnowAboutThisInfo{procsThatNeedToKnowAboutThisInfo},
      mGlobalIdsOfSnapNodeElems{globalIdsOfSnapNodeElems},
      mPostWorstQuality(postSnapQuality), // () for intel 17
      mNodeLocation(nodeLocation),
      mSnapLocation(snapLocation),
      mSnapRank{snapRank}
    {
    }

    static GlobalId make_unique_id(const size_t globalNodeId, const size_t intersectionPointIndex) { return GlobalId{globalNodeId, intersectionPointIndex}; }
    const GlobalId &get_unique_id() const { return mUniqueId; }
    size_t get_node_global_id() const { return mUniqueId[0]; }
    size_t get_intersection_point_index() const { return mUniqueId[1]; }
    void set_intersection_point_index(const size_t intPtIndex) { mUniqueId[1] = intPtIndex; }
    int get_owner() const { return mOwner; }
    const std::vector<int> &get_procs_that_need_to_know_about_this_info() const { return mProcsThatNeedToKnowAboutThisInfo; }
    double get_post_worst_quality() const { return mPostWorstQuality; }
    const stk::math::Vector3d & get_snap_location() const { return mSnapLocation; };
    const stk::math::Vector3d & get_node_location() const { return mNodeLocation; };
    int get_snap_rank() const { return mSnapRank; }
    const std::vector<ExclusionIdentifierType> &get_conflicting_ids() const { return mGlobalIdsOfSnapNodeElems; }

    class Comparator
    {
    public:
        Comparator(const QualityMetric &qualityMetric) : mQualityMetric{qualityMetric}{}
        bool is_first_higher_priority_than_second(const SnapInfo& tetSnapInfoA,const SnapInfo& tetSnapInfoB) const;
        bool does_first_win_priority_tie_with_second(const SnapInfo& tetSnapInfoA,const SnapInfo& tetSnapInfoB) const;

    private:
        const QualityMetric &mQualityMetric;
    };

private:
    GlobalId mUniqueId;
    int mOwner{0u};
    std::vector<int> mProcsThatNeedToKnowAboutThisInfo;
    std::vector<size_t> mGlobalIdsOfSnapNodeElems;
    double mPostWorstQuality{0.0};
    stk::math::Vector3d mNodeLocation{};
    stk::math::Vector3d mSnapLocation{};
    int mSnapRank;
};

std::ostream & operator<<(std::ostream & os, const SnapInfo& snapInfo);

}

#endif /* AKRI_SNAPINFO_H */
