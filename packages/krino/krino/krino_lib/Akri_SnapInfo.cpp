// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SnapInfo.hpp>

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

    if(mQualityMetric.is_first_quality_metric_better_than_second(snapInfoA.get_post_worst_quality(),snapInfoB.get_post_worst_quality()))
        return true;

    return false;
}

bool SnapInfo::Comparator::does_first_win_priority_tie_with_second(const SnapInfo& snapInfoA,const SnapInfo& snapInfoB) const
{
    if (snapInfoA.get_node_global_id() != snapInfoB.get_node_global_id())
    {
      if ( is_less_than_in_x_then_y_then_z(snapInfoA.get_node_location(), snapInfoB.get_node_location()) )
          return true;
      else if ( is_less_than_in_x_then_y_then_z(snapInfoB.get_node_location(), snapInfoA.get_node_location()) )
          return false;
    }

    if ( is_less_than_in_x_then_y_then_z(snapInfoA.get_snap_location(), snapInfoB.get_snap_location()) )
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
    return os;
}

}
