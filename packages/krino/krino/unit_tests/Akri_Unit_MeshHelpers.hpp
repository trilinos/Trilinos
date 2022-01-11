// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_UNIT_TESTS_INCLUDE_AKRI_UNIT_MESHHELPERS_H_
#define KRINO_UNIT_TESTS_INCLUDE_AKRI_UNIT_MESHHELPERS_H_

#include <stk_mesh/base/Types.hpp>
#include <vector>

namespace krino
{

void build_mesh(stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::EntityIdVector> & elem_nodes,
    const std::vector<int> & elem_procs,
    std::vector<stk::mesh::PartVector> & elem_parts);
}

#endif /* KRINO_UNIT_TESTS_INCLUDE_AKRI_UNIT_MESHHELPERS_H_ */
