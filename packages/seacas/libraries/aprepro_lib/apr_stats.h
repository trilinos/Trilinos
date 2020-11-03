/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "apr_util.h" // for check_valid_var, new_string, etc

namespace SEAMS {
  class Stats
  {
  public:
    Stats() = default;

    void   newsample(int n);
    double mean() const;
    double deviation() const;
    double variance() const;

  private:
    size_t Numnums{0};
    double Mean{0.0};
    double StdDev{0.0};
  };
} // namespace SEAMS
