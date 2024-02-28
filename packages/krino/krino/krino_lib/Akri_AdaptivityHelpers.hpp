// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_ADAPTIVITYHELPERS_H_
#define KRINO_INCLUDE_AKRI_ADAPTIVITYHELPERS_H_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <functional>
#include <string>

namespace krino {

class CDFEM_Support;
class RefinementInterface;

stk::mesh::Selector cdfem_do_not_refine_or_unrefine_selector(const CDFEM_Support & cdfem_support);

enum class Refinement_Marker
{
  COARSEN = -1,
  NOTHING = 0,
  REFINE = 1
};

void
perform_multilevel_adaptivity(RefinementInterface & refinement,
    stk::mesh::BulkData & mesh,
    const std::function<void(int)> & marker_function,
    const stk::mesh::Selector & do_not_refine_selector = stk::mesh::Selector());

}

#endif /* KRINO_INCLUDE_AKRI_ADAPTIVITYHELPERS_H_ */
