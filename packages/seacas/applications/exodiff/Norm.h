// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef ED_NORM_H
#define ED_NORM_H

#include <cmath>
class Norm
{
public:
  Norm() = default;

  double diff(int order) const
  {
    if (order == 1) {
      return l1_norm_d;
    }
    if (order == 2) {
      return std::sqrt(l2_norm_d);
    }
    else {
      return 0.0;
    }
  }

  double left(int order) const
  {
    if (order == 1) {
      return l1_norm_1;
    }
    if (order == 2) {
      return std::sqrt(l2_norm_1);
    }
    else {
      return 0.0;
    }
  }

  double right(int order) const
  {
    if (order == 1) {
      return l1_norm_2;
    }
    if (order == 2) {
      return std::sqrt(l2_norm_2);
    }
    else {
      return 0.0;
    }
  }

  double relative(int order) const
  {
    double l      = left(order);
    double r      = right(order);
    double lr_max = l > r ? l : r;
    return diff(order) / lr_max;
  }

  void add_value(double val1, double val2)
  {
    l1_norm_d += std::fabs(val1 - val2);
    l1_norm_1 += std::fabs(val1);
    l1_norm_2 += std::fabs(val2);

    l2_norm_d += (val1 - val2) * (val1 - val2);
    l2_norm_1 += val1 * val1;
    l2_norm_2 += val2 * val2;
  }

  double l1_norm_1{0.0};
  double l1_norm_2{0.0};
  double l1_norm_d{0.0};

  double l2_norm_1{0.0};
  double l2_norm_2{0.0};
  double l2_norm_d{0.0};
};

#endif
