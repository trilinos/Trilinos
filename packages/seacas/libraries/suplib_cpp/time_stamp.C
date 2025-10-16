// Copyright(C) 1999-2021, 2023, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <ctime>
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <time_stamp.h>

std::string time_stamp(const std::string &format)
{
  if (format.empty()) {
    return {""};
  }

  auto        now         = std::chrono::system_clock::now();
  std::string time_string = fmt::format(fmt::runtime(format), now);
  return time_string;
}
