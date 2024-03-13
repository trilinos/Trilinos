// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Snapper.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>

namespace krino {

void determine_node_snapping_from_intersection_points(const stk::mesh::BulkData & mesh,
    const std::vector<IntersectionPoint> & intersectionPoints,
    const CDFEM_Snapper & snapper,
    NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  const double snapTol = snapper.get_edge_tolerance();
  for (auto && intersectionPoint : intersectionPoints)
  {
    const auto & nodes = intersectionPoint.get_nodes();
    const auto & weights = intersectionPoint.get_weights();
    for (size_t nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex)
    {
      if (weights[nodeIndex] > 1.-snapTol)
      {
        if (krinolog.shouldPrint(LOG_DEBUG))
        {
          krinolog << "Snapping node " << debug_output(mesh, intersectionPoint) << " to " << mesh.identifier(nodes[nodeIndex]) << stk::diag::dendl;
        }

        const auto & intersectionPointDomains = intersectionPoint.get_sorted_domains();
        auto & nodeCapturedDomains = nodesToCapturedDomains[nodes[nodeIndex]];
        nodeCapturedDomains.insert(nodeCapturedDomains.end(), intersectionPointDomains.begin(), intersectionPointDomains.end());
      }
    }
  }

  for (auto && nodeToCapturedDomains : nodesToCapturedDomains)
    stk::util::sort_and_unique(nodeToCapturedDomains.second);
}

void snap_to_node(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & interfaceGeometry,
    const CDFEM_Snapper & snapper,
    NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
    const std::vector<IntersectionPoint> intersectionPoints = build_uncaptured_intersection_points(mesh, elementSelector, interfaceGeometry, nodesToCapturedDomains);

    determine_node_snapping_from_intersection_points(mesh, intersectionPoints, snapper, nodesToCapturedDomains);
    communicate_node_captured_domains_for_all_nodes(mesh, nodesToCapturedDomains);
}

}


