// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_DECOMP_AKRI_SHARPEDGERECOVERYBYSWAPS_HPP_
#define KRINO_KRINO_DECOMP_AKRI_SHARPEDGERECOVERYBYSWAPS_HPP_

#include <array>
#include <functional>
#include <optional>
#include <vector>

#include <stk_math/StkVector.hpp>

namespace krino {

// Swaps of the tets of an element's leaves that make the segments of sharp edges inside it edges of
// the tets, where they can, without making nodes.  The swaps are found one at a time (find_swap),
// until no internal face of the tets has a sharp-edge crossing inside it that the face does not
// capture, and no internal edge one within the edge tolerance of it:
// - a 2-3 swap of the two tets on a crossed face, if a crossing is of the segment joining their
//   opposite nodes: both nodes have its domains, and it is within the edge tolerance of the line
//   through them, in owner coordinates;
// - or else an edge swap of an edge of the face, to the configuration of best quality;
// - or, for a crossing within the edge tolerance of an edge of a face, an edge swap of that edge.
//   The sharp edge crosses the leaf edge there, and which faces around it find the crossing depends
//   on roundoff, so it is a crossing of the edge, whichever face it is found on.
// A swap is made only if none of its new internal faces has a crossing inside it, and none of its
// new internal edges (on no face of only one of its tets) one near it, so each swap removes at least
// one crossed face or edge.  A crossing within the edge tolerance of a node of a face is left to the
// cutting: no swap moves the node.  The crossed faces and edges are tried in the order of their node
// coordinates (less_by_coordinates_because_it_impacts_structure), so the swaps do not depend on the ids or the order
// of the tets.
//
// NODE is a node handle: a const SubElementNode * in the decomposition, an int in the unit tests
// (explicitly instantiated for both).
template <typename NODE>
class SharpEdgeRecoveryBySwaps
{
public:
  using Tet = std::array<NODE, 4>;
  using Face = std::array<NODE, 3>;

  // A crossing of a face by a sharp edge, at these parametric coordinates of the face (with its
  // nodes in canonical order), with the sharp edge's domains.
  struct FaceCrossing
  {
    stk::math::Vector3d faceCoords;
    std::vector<int> sortedDomains;
  };

  // A swap: the tets replaced, by index, and the tets that replace them, positively oriented.
  struct Swap
  {
    std::vector<size_t> oldTets;
    std::vector<Tet> newTets;
  };

  // coordinates: a node's location, for the tets' orientation and quality and for the order in which
  // crossed faces and edges are tried.  ownerCoordinates: its parametric coordinates in the element.
  // nodeHasDomains: whether the node's domains capture a crossing's (node_domains_capture_point_domains).
  // canonicalFace: a face's nodes in an order that does not depend on the tet it is in, so that its
  // crossings do not depend on it, even to roundoff.  uncapturedFaceCrossings: the sharp-edge
  // crossings of a face (with its nodes in canonical order) that it does not capture
  // (tet_face_point_is_captured).
  SharpEdgeRecoveryBySwaps(const double edgeTol,
      std::function<stk::math::Vector3d(NODE)> coordinates,
      std::function<stk::math::Vector3d(NODE)> ownerCoordinates,
      std::function<bool(NODE, const std::vector<int> &)> nodeHasDomains,
      std::function<Face(NODE, NODE, NODE)> canonicalFace,
      std::function<std::vector<FaceCrossing>(const Face &)> uncapturedFaceCrossings);

  // The first swap of the tets, in the order above, or none.
  std::optional<Swap> find_swap(const std::vector<Tet> & tets) const;

  // Make the swaps, one at a time, until there is none: each replaces its old tets, removed, by its
  // new ones, added at the end.
  void make_swaps(std::vector<Tet> & tets) const;

  // Whether an internal face of the tets (in two of them) has a crossing inside it, or an internal
  // edge (on no face in only one of them) has one within the edge tolerance of it.
  bool have_uncaptured_crossing_inside(const std::vector<Tet> & tets) const;

private:
  using Edge = std::array<NODE, 2>;
  // The crossings of the face inside it, and the edges of the face (nodes in coordinate order) that
  // a crossing is within the edge tolerance of.
  std::vector<FaceCrossing> fill_crossings_inside_face_and_crossed_edges(
      const Face & canonicalFaceNodes, std::vector<Edge> & crossedFaceEdges) const;
  std::optional<Swap> find_face_swap(const Face & canonicalFaceNodes,
      const std::vector<FaceCrossing> & crossings,
      std::vector<std::pair<size_t, NODE>> tetsAndOppositeNodes) const;
  std::optional<Swap> find_edge_swap(const std::vector<Tet> & tets, const NODE edgeNode0, const NODE edgeNode1) const;
  bool is_less_in_coordinate_order(const NODE a, const NODE b) const;
  Edge edge_in_coordinate_order(const NODE a, const NODE b) const;
  std::vector<stk::math::Vector3d> coordinates_in_coordinate_order(const std::vector<NODE> & nodes) const;

  double myEdgeTol;
  std::function<stk::math::Vector3d(NODE)> myCoordinates;
  std::function<stk::math::Vector3d(NODE)> myOwnerCoordinates;
  std::function<bool(NODE, const std::vector<int> &)> myNodeHasDomains;
  std::function<Face(NODE, NODE, NODE)> myCanonicalFace;
  std::function<std::vector<FaceCrossing>(const Face &)> myUncapturedFaceCrossings;
};

} // namespace krino

#endif /* KRINO_KRINO_DECOMP_AKRI_SHARPEDGERECOVERYBYSWAPS_HPP_ */
