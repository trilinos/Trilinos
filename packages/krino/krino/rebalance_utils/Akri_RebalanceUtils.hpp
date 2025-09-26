// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_REBALANCE_UTILS_INCLUDE_AKRI_REBALANCEUTILS_H_
#define KRINO_REBALANCE_UTILS_INCLUDE_AKRI_REBALANCEUTILS_H_

#include <stk_mesh/base/Selector.hpp>
#include <string>
#include <vector>

namespace krino { class CDMesh; }
namespace krino { class RefinementManager; }

namespace krino {
namespace rebalance_utils {

bool have_parmetis();

// This function will call STK rebalance to rebalance the specified selections of the
// mesh based on the ELEMENT_RANK field element_weights_field. If cdmesh != nullptr the
// rebalance operation will ensure that all CDFEM child elements are moved to the same
// processor as their parent elements. The weights field will be adjusted
// such that the parent elements are weighted with the sum of their child
// weights, and the child weights are 0 for both CDFEM and adaptivity parents/children.
bool rebalance_mesh(stk::mesh::BulkData & bulk_data,
    const RefinementManager * refinement,
    CDMesh * cdmesh,
    const std::string & element_weights_field_name,
    const std::string & coordinates_field_name,
    const std::vector<stk::mesh::Selector> & selections_to_rebalance_separately,
    const unsigned max_num_nodal_rebal_iters,
    const std::string & decomp_method = "parmetis",
    const double imbalance_threshold = 1.05);

// This version handles multiple criteria rebalancing with different weights for each
// criterion.
bool rebalance_mesh(stk::mesh::BulkData & bulk_data,
    const RefinementManager * refinement,
    CDMesh * cdmesh,
    const std::vector<std::string> & element_weights_field_names,
    const std::string & coordinates_field_name,
    const unsigned max_num_nodal_rebal_iters,
    const std::string & decomp_method = "parmetis",
    const double imbalance_threshold = 1.05);
}
}

#endif /* KRINO_REBALANCE_UTILS_INCLUDE_AKRI_REBALANCEUTILS_H_ */
