// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothInfo.hpp>

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

bool SmoothInfo::Comparator::is_first_higher_priority_than_second(const SmoothInfo& a, const SmoothInfo& b) const
{
  const double infoASqrLen = (a.get_post_smooth_location()-a.get_pre_smooth_location()).length_squared();
  const double infoBSqrLen = (b.get_post_smooth_location()-b.get_pre_smooth_location()).length_squared();
  if (infoASqrLen > infoBSqrLen)
    return true;
  else if (infoBSqrLen > infoASqrLen)
    return false;

  return false;
}

bool SmoothInfo::Comparator::does_first_win_priority_tie_with_second(const SmoothInfo& a, const SmoothInfo& b) const
{
  if (a.get_node_global_id() != b.get_node_global_id())
    return is_less_than_in_x_then_y_then_z(a.get_pre_smooth_location(), b.get_pre_smooth_location());

  throw std::runtime_error("Comparator failed to distinguish between given SmoothInfos for nodes " + std::to_string(a.get_node_global_id()) + " and " + std::to_string(b.get_node_global_id()));
}

std::ostream & operator<<(std::ostream & os, const SmoothInfo& smoothInfo)
{
    os << "  Owner: " << smoothInfo.get_owner() << "\n";
    os << "Smooth Node: " << smoothInfo.get_node_global_id() << "\n";
    os << "Pre-smooth location: " << smoothInfo.get_pre_smooth_location().to_string(16) << "\n";
    os << "Post-smooth location: " << smoothInfo.get_post_smooth_location().to_string(16) << "\n";
    os << "GetProcsThatNeedToKnow: " << to_string(smoothInfo.get_procs_that_need_to_know_about_this_info()) << "\n";
    return os;
}

}
