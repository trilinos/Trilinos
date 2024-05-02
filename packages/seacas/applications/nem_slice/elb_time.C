/*
 * Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <chrono> // for duration, etc
double get_time()
{
  static auto                   start = std::chrono::steady_clock::now();
  auto                          now   = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff  = now - start;
  return diff.count();
}
