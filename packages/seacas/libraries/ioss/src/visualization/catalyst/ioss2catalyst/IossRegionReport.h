// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef Ioss_Region_Report_h
#define Ioss_Region_Report_h

#include "Ioss_Region.h"
#include <iostream>
#include <string>
#include <vector>

namespace ioss_region_report {

  using Message = std::string;
  using Key     = std::string;

  struct Messages
  {
    std::string          begin{""};
    std::vector<Message> messages;

    Messages &operator+=(const Message &rhs)
    {
      messages.push_back(begin + rhs);
      return *this;
    }

    Messages &operator+=(const Messages &rhs)
    {
      for (auto msg : rhs.messages)
        messages.push_back(begin + msg);
      return *this;
    }
  };

  std::ostream &operator<<(std::ostream &os, const Messages &messages);
  Messages      region_report(const Ioss::Region &region);

} // namespace ioss_region_report

#endif
