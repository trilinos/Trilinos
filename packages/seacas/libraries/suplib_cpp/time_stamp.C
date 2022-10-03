// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "fmt/chrono.h"
#include <time_stamp.h>

std::string time_stamp(const std::string &format)
{
  if (format == "") {
    return std::string("");
  }

  time_t      calendar_time = std::time(nullptr);
  struct tm  *local_time    = std::localtime(&calendar_time);
  std::string time_string   = fmt::format(format, *local_time);
  return time_string;
}
