// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_CDFEM_SNAPPER_H_
#define KRINO_INCLUDE_AKRI_CDFEM_SNAPPER_H_
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

class CDFEM_Snapper {
public:
  CDFEM_Snapper() : my_cdfem_edge_tol(default_edge_tol) {}

  void set_edge_tolerance(double edge_tol) { my_cdfem_edge_tol = edge_tol; }
  double get_edge_tolerance() const { STK_ThrowRequireMsg(my_cdfem_edge_tol < 0.5, "Should not be requesting tolerance when always snapping."); return my_cdfem_edge_tol; }
  bool always_snap() const { return my_cdfem_edge_tol > 0.5; }
  bool snap_lo(double crossingVal) const { return always_snap() ? (crossingVal <= 0.5) : (crossingVal < get_edge_tolerance()); }
  bool snap_hi(double crossingVal) const { return always_snap() ? (crossingVal > 0.5) : (crossingVal > 1.0-get_edge_tolerance()); }
private:
  static constexpr double default_edge_tol = 1.e-3;
  double my_cdfem_edge_tol;
};

}

#endif /* KRINO_INCLUDE_AKRI_CDFEM_SNAPPER_H_ */
