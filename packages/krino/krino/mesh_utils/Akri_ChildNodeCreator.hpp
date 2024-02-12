// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CHILDNODECREATOR_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CHILDNODECREATOR_HPP_
#include <stk_mesh/base/BulkData.hpp>

namespace krino {

struct ChildNodeRequest
{
    std::vector<stk::mesh::Entity*> parents;
    stk::mesh::Entity* child;

    ChildNodeRequest(const std::vector<stk::mesh::Entity*> & in_parents, stk::mesh::Entity *in_child)
    : parents(in_parents), child(in_child) {}
};

typedef std::function<void(stk::topology::rank_t, size_t, std::vector<stk::mesh::EntityId>&)> GenerateNewEntityIdsFunction;

void batch_create_child_nodes(stk::mesh::BulkData & mesh, const std::vector< ChildNodeRequest > & child_node_requests, const stk::mesh::PartVector & node_parts, const GenerateNewEntityIdsFunction & generate_new_ids);

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CHILDNODECREATOR_HPP_ */
