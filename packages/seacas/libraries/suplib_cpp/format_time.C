/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <fmt/format.h>
#include <format_time.h>

std::string format_time(double seconds)
{
  std::string suffix("u");
  if (seconds > 0.0 && seconds < 1.0) {
    seconds *= 1000.;
    suffix = "ms";
  }
  else if (seconds > 86400) {
    suffix = "d";
    seconds /= 86400.;
  }
  else if (seconds > 3600) {
    suffix = "h";
    seconds /= 3600.;
  }
  else if (seconds > 60) {
    suffix = "m";
    seconds /= 60.;
  }
  else {
    suffix = "s";
  }
  return fmt::format("{:.3}{}", seconds, suffix);
}
