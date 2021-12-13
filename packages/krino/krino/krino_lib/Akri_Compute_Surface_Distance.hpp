// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Compute_Surface_Distance_h
#define Akri_Compute_Surface_Distance_h

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/diag/Timer.hpp>

namespace krino {

class Compute_Surface_Distance {
public:
  static void calculate(
    const stk::mesh::BulkData & mesh,
    const stk::diag::Timer &parent_timer,
    const stk::mesh::Field<double>& coordinates,
    const stk::mesh::Field<double>& distance,
    const stk::mesh::Selector & surface_selector,
    const double narrowBandSize = 0.0,
    const double farFieldValue = 0.0);

  static void calculate(
    const stk::mesh::BulkData & mesh,
    const stk::diag::Timer &parent_timer,
    const stk::mesh::Field<double>& coordinates,
    const stk::mesh::Field<double>& distance,
    const stk::mesh::Selector & volume_selector,
    const stk::mesh::Selector & surface_selector,
    const double narrowBandSize = 0.0,
    const double farFieldValue = 0.0);
  };

} // namespace krino

#endif // Akri_Compute_Surface_Distance_h
