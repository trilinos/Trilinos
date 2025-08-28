// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_REBALANCE_UTILS_INCLUDE_AKRI_REBALANCEUTILS_IMPL_H_
#define KRINO_REBALANCE_UTILS_INCLUDE_AKRI_REBALANCEUTILS_IMPL_H_

#include <stk_mesh/base/Field.hpp>

namespace stk { namespace balance { class DecompositionChangeList; } }
namespace stk { namespace mesh { class BulkData; } }
namespace krino { class CDMesh; }
namespace krino { class RefinementManager; }

namespace krino {
namespace rebalance_utils {
namespace impl {

void
update_rebalance_for_adaptivity(stk::balance::DecompositionChangeList & decomp,
    const RefinementManager& refinement,
    const stk::mesh::BulkData & bulk_data);

void
update_rebalance_for_cdfem(stk::balance::DecompositionChangeList & decomp,
    const stk::mesh::BulkData & bulk_data,
    const CDMesh & cdmesh);

void
accumulate_cdfem_child_weights_to_parents(const stk::mesh::BulkData & bulk_data,
    stk::mesh::Field<double> & element_weights_field,
    const CDMesh & cdmesh);

void accumulate_adaptivity_child_weights_to_parents(
    const stk::mesh::BulkData & bulk_data, const RefinementManager& refinement, stk::mesh::Field<double> & element_weights_field);

bool check_family_tree_element_and_side_ownership(const stk::mesh::BulkData & bulk_data);
}
}
}

#endif /* KRINO_REBALANCE_UTILS_INCLUDE_AKRI_REBALANCEUTILS_IMPL_H_ */
