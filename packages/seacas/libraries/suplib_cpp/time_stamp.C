// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <fmt/chrono.h>
#include <time_stamp.h>

std::string time_stamp(const std::string &format)
{
  if (format.empty()) {
    return {""};
  }

  std::time_t t           = std::time(nullptr);
  std::string time_string = fmt::format(fmt::runtime(format), fmt::localtime(t));
  return time_string;
}
