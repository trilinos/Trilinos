// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SnapInfo.hpp>
#include <Akri_EntityOrder.hpp>

#include <Akri_QualityMetric.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <vector>
#include <iostream>

namespace krino
{

template <typename Container>
std::string
to_string(const Container & container)
{
  std::ostringstream os;
  os << " {";
  for(auto data : container)
      os << " " << data;
  os << " }";
  return os.str();
}

bool SnapInfo::Comparator::is_first_higher_priority_than_second(const SnapInfo& snapInfoA,const SnapInfo& snapInfoB) const
{
    if( snapInfoA.get_snap_rank() > snapInfoB.get_snap_rank() )
        return true;
    else if ( snapInfoA.get_snap_rank() < snapInfoB.get_snap_rank() )
        return false;

    if (mIsPredictable)
    {
      // Closest snap first, then the one whose distance grows slowest with the level set, so that
      // the choice for a small increase of the level set is the one made at the current design.
      const float distA = static_cast<float>(snapInfoA.get_snap_distance());
      const float distB = static_cast<float>(snapInfoB.get_snap_distance());
      if (distA != distB) return distA < distB;
      return static_cast<float>(snapInfoA.get_snap_distance_rate()) <
          static_cast<float>(snapInfoB.get_snap_distance_rate());
    }

    if(mQualityMetric.is_first_quality_metric_better_than_second(snapInfoA.get_post_worst_quality(),snapInfoB.get_post_worst_quality()))
        return true;

    return false;
}

bool SnapInfo::Comparator::does_first_win_priority_tie_with_second(const SnapInfo& snapInfoA,const SnapInfo& snapInfoB) const
{
    if (snapInfoA.get_node_global_id() != snapInfoB.get_node_global_id())
    {
      if ( less_by_coordinates_because_it_impacts_structure(snapInfoA.get_node_location(), snapInfoB.get_node_location()) )
          return true;
      else if ( less_by_coordinates_because_it_impacts_structure(snapInfoB.get_node_location(), snapInfoA.get_node_location()) )
          return false;
    }

    if ( less_by_coordinates_because_it_impacts_structure(snapInfoA.get_snap_location(), snapInfoB.get_snap_location()) )
        return true;
    if (less_by_coordinates_because_it_impacts_structure(
            snapInfoB.get_snap_location(), snapInfoA.get_snap_location()))
      return false;

    // Snaps of the same node to the same location (a node on the interface): the one that will be
    // first in x then y then z after a small increase of the level set.
    if (less_by_coordinates_because_it_impacts_structure(
            snapInfoA.get_snap_direction(), snapInfoB.get_snap_direction()))
      return true;

    return false;
}

std::ostream & operator<<(std::ostream & os, const SnapInfo& snapInfo)
{
    os << "Snap Node: " << snapInfo.get_node_global_id() << std::endl;
    os << "IntptId: " << snapInfo.get_intersection_point_index() << std::endl;
    os << "Owner: " << snapInfo.get_owner() << std::endl;
    os << "ConflictingIds: " << to_string(snapInfo.get_conflicting_ids()) << std::endl;
    os << "SnapLocation: " << snapInfo.get_snap_location().to_string(16) << std::endl;
    os << "PostSnapWorstQuality: " << snapInfo.get_post_worst_quality() << std::endl;
    os << "GetProcsThatNeedToKnow: " << to_string(snapInfo.get_procs_that_need_to_know_about_this_info()) << std::endl;
    os << "Snap Rank: " << snapInfo.get_snap_rank() << std::endl;
    os << "Snap Distance Rate: " << snapInfo.get_snap_distance_rate() << std::endl;
    os << "Snap Direction: " << snapInfo.get_snap_direction() << std::endl;
    return os;
}

}
