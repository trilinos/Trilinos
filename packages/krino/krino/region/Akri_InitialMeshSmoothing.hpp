// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_REGION_AKRI_INITIALMESHSMOOTHING_HPP_
#define KRINO_KRINO_REGION_AKRI_INITIALMESHSMOOTHING_HPP_

#include <optional>

namespace krino {

// How the initial mesh is smoothed, after any initial refinement and conforming decomposition (the region
// option smooth_initial_mesh).
struct InitialMeshSmoothing
{
  bool smooth{false};
  // Whether the nodes on the boundaries, sidesets, nodesets and interfaces slide within their surfaces
  // and along their curves (smooth_mesh_respecting_mesh_based_geometry), or are fixed
  // (smooth_mesh_keeping_boundary_nodes_fixed).
  bool boundaryNodesSlide{true};
  // With sliding boundary nodes, the curves and corners also where a surface is sharper than this
  // interior angle (180 is flat), as well as those of the topology.
  std::optional<double> sharpFeatureAngleInDegrees;
};

}

#endif /* KRINO_KRINO_REGION_AKRI_INITIALMESHSMOOTHING_HPP_ */
