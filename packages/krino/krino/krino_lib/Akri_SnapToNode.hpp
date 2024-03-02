// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_SNAPTONODE_H_
#define KRINO_INCLUDE_AKRI_SNAPTONODE_H_

#include <Akri_Intersection_Points.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace krino {

class CDFEM_Support;
class InterfaceGeometry;
class Phase_Support;

void snap_to_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & interfaceGeometry,
    const CDFEM_Snapper & snapper,
    NodeToCapturedDomainsMap & nodesToCapturedDomains);

}
#endif /* KRINO_INCLUDE_AKRI_SNAPTONODE_H_ */
