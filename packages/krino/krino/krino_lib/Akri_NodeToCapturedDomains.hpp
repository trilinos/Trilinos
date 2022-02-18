// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_NODETOCAPTUREDDOMAINS_H_
#define KRINO_INCLUDE_AKRI_NODETOCAPTUREDDOMAINS_H_
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <map>
#include <vector>

namespace krino {

typedef std::map<stk::mesh::Entity, std::vector<int>> NodeToCapturedDomainsMap;

void communicate_node_captured_domains_for_given_nodes(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & nodes,
    NodeToCapturedDomainsMap & nodesToCapturedDomains);

void communicate_node_captured_domains_for_all_nodes(const stk::mesh::BulkData & mesh,
    NodeToCapturedDomainsMap & nodesToCapturedDomains);

typedef std::map<stk::mesh::Entity, int> ElementToDomainMap;

}

#endif /* KRINO_INCLUDE_AKRI_NODETOCAPTUREDDOMAINS_H_ */
