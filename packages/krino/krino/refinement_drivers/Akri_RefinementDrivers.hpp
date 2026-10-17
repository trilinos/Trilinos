// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTDRIVERS_HPP_
#define KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTDRIVERS_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <functional>

namespace krino {

class RefinementManager;

enum class Refinement_Marker
{
  COARSEN = -1,
  NOTHING = 0,
  REFINE = 1
};

void
perform_multilevel_adaptivity(RefinementManager & refinement,
    stk::mesh::BulkData & mesh,
    const std::function<void(int)> & marker_function,
    const stk::mesh::Selector & do_not_refine_selector,
    const int rebalanceInterval = -1);

}

#endif /* KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTDRIVERS_HPP_ */
