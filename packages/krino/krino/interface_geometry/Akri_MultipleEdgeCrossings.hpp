// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_MULTIPLEEDGECROSSINGS_HPP_
#define KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_MULTIPLEEDGECROSSINGS_HPP_

#include <vector>

#include <stk_math/StkVector.hpp>

namespace krino {

// The multiple crossings of a surface to cut on the edge from x0 to x1, from the positions of its
// crossings along it (from x0): two neighbors within the close-crossing tolerance of each other are
// close crossings and are removed, pairing them from the same end in every element with the edge,
// and a crossing within the end tolerance of an end is that end.  Returns the positions to cut
// (from x0, sorted), 0 or 1 for an end.
std::vector<double> multiple_crossing_positions_to_cut(const std::vector<double> & positionsFromNode0,
    const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const double closeCrossingTol,
    const double atEndTol);

// How far from an end a multiple crossing of a surface must be to be cut, as a fraction of the
// extent of the edge (of its mesh ancestors): the edge tolerance for a surface with sharp features,
// to capture the features beside a node, and 0.33 for one without.  Most multiple crossings of a
// surface without sharp features are beside a node already on the surface, where the surface
// crosses back just past it: cutting them adds little volume and makes slivers.
double multiple_crossing_end_tolerance(const bool surfaceHasSharpFeatures, const double edgeTol);

// Whether an edge whose ends have these domains is a chord of an interface with these domains (both
// ends have them): part of its linear approximation.  Where a curved surface crosses the chord is
// its curvature, as on any edge, and cutting there would only make chords that cross it again,
// without end, so a chord has no multiple crossings to cut.
bool edge_is_chord_of_interface(const std::vector<int> & end0SortedDomains,
    const std::vector<int> & end1SortedDomains,
    const std::vector<int> & interfaceSortedDomains);

// The positions to cut (multiple_crossing_positions_to_cut) without an end (0 or 1) whose node
// already has the interface's domains: the crossing there is that node.
std::vector<double> multiple_crossing_positions_to_insert(const std::vector<double> & positionsFromNode0,
    const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const std::vector<int> & end0SortedDomains,
    const std::vector<int> & end1SortedDomains,
    const std::vector<int> & interfaceSortedDomains,
    const double closeCrossingTol,
    const double atEndTol);

}

#endif /* KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_MULTIPLEEDGECROSSINGS_HPP_ */
