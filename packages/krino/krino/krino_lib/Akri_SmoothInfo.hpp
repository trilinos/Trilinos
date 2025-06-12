// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_SMOOTHINFO_H
#define AKRI_SMOOTHINFO_H

#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <vector>

namespace krino
{

class SmoothInfo
{
public:
    typedef stk::mesh::EntityId ExclusionIdentifierType;
    typedef stk::mesh::EntityId GlobalId;

    SmoothInfo(const stk::mesh::EntityId nodeId,
        const int owner,
        const stk::math::Vector3d & preSmoothLocation,
        const stk::math::Vector3d & postSmoothLocation,
        const std::vector<ExclusionIdentifierType> & globalIdsOfCavityElems,
        const std::vector<int> & procsThatNeedToKnowAboutThisInfo)
    : mUniqueId{nodeId},
      mOwner(owner),
      mPreSmoothLocation(preSmoothLocation),
      mPostSmoothLocation(postSmoothLocation),
      mGlobalIdsOfSmoothNodeElems(globalIdsOfCavityElems),
      mProcsThatNeedToKnowAboutThisInfo(procsThatNeedToKnowAboutThisInfo)
    {
    }

    const GlobalId &get_unique_id() const { return mUniqueId; }
    stk::mesh::EntityId get_node_global_id() const { return mUniqueId; }

    int get_owner() const { return mOwner; }
    const std::vector<int> &get_procs_that_need_to_know_about_this_info() const { return mProcsThatNeedToKnowAboutThisInfo; }
    const stk::math::Vector3d & get_pre_smooth_location() const { return mPreSmoothLocation; }
    const stk::math::Vector3d & get_post_smooth_location() const { return mPostSmoothLocation; }
    const std::vector<ExclusionIdentifierType> &get_conflicting_ids() const { return mGlobalIdsOfSmoothNodeElems; }
    const std::vector<stk::mesh::EntityId> &get_ids_of_elements_impacted_by_smoothing() const { return mGlobalIdsOfSmoothNodeElems; }

    class Comparator
    {
    public:
        bool is_first_higher_priority_than_second(const SmoothInfo& a,const SmoothInfo& b) const;
        bool does_first_win_priority_tie_with_second(const SmoothInfo& a,const SmoothInfo& b) const;
    };

private:
  const GlobalId mUniqueId;
  int mOwner;
  stk::math::Vector3d mPreSmoothLocation;
  stk::math::Vector3d mPostSmoothLocation;
  std::vector<ExclusionIdentifierType> mGlobalIdsOfSmoothNodeElems;
  std::vector<int> mProcsThatNeedToKnowAboutThisInfo;
};

std::ostream & operator<<(std::ostream & os, const SmoothInfo& smoothInfo);

}

#endif /* AKRI_SMOOTHINFO_H */
