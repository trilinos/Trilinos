// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Region_Parser_h
#define Akri_Region_Parser_h

#include <Akri_InitialMeshSmoothing.hpp>

namespace krino { class Simulation; }
namespace krino { namespace Parser { class Node; } }

namespace krino {
namespace Region_Parser {
 void parse(const Parser::Node & node, Simulation & simulation);

 // The region's smooth_initial_mesh: absent or false, no smoothing; true, or a block whose options are
 // all optional (boundary_nodes: slide (the default) or fixed; sharp_feature_angle, in degrees, only
 // with slide), smoothing.
 InitialMeshSmoothing parse_initial_mesh_smoothing(const Parser::Node & regionNode);
}
}

#endif // Akri_Region_Parser_h
