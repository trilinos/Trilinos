// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include <cmath>
class DiffData
{
public:
  enum Type {
    mm_unknown = 0,
    mm_time    = 1, // Only val and step valid.
    mm_global  = 2, // Only val and step valid.
    mm_nodal   = 3, // Only val, step, and id valid.
    mm_element = 4, // All fields valid for the rest.
    mm_sideset = 5,
    mm_nodeset = 6,
    mm_elematt = 7 // step not valid
  };

  DiffData() {}

  void set_max(double d, double val_1, double val_2, size_t id_ = 0, size_t blk_ = 0)
  {
    if (diff < d) {
      diff = d;
      val1 = val_1;
      val2 = val_2;
      id   = id_;
      blk  = blk_;
    }
  }

  double diff{0.0};
  double val1{0.0};
  double val2{0.0};
  size_t id{0};
  size_t blk{0};

  Type type{mm_unknown};
};

class MinMaxData
{
public:
  enum Type {
    mm_unknown = 0,
    mm_time    = 1, // Only val and step valid.
    mm_global  = 2, // Only val and step valid.
    mm_nodal   = 3, // Only val, step, and id valid.
    mm_element = 4, // All fields valid for the rest.
    mm_sideset = 5,
    mm_nodeset = 6,
    mm_elematt = 7 // step not valid
  };
  MinMaxData() : min_val(DBL_MAX) {}

  void spec_min_max(double val, int step, size_t id = 0, size_t blk = 0)
  {
    if (std::fabs(val) < min_val) {
      min_val  = std::fabs(val);
      min_step = step;
      min_id   = id;
      min_blk  = blk;
    }

    if (std::fabs(val) > max_val) {
      max_val  = std::fabs(val);
      max_step = step;
      max_id   = id;
      max_blk  = blk;
    }
  }

  double min_val{};
  int    min_step{0};
  size_t min_id{0};
  size_t min_blk{0};

  double max_val{-1.0};
  int    max_step{0};
  size_t max_id{0};
  size_t max_blk{0};

  Type type{mm_unknown};
};
