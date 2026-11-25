// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SharpEdgeRecoveryBySwaps.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <utility>

#include <Akri_EntityOrder.hpp>
#include <Akri_SharpFeatureLeafUtils.hpp>
#include <Akri_TetSwaps.hpp>

namespace krino {

class SubElementNode;

template <typename NODE>
SharpEdgeRecoveryBySwaps<NODE>::SharpEdgeRecoveryBySwaps(const double edgeTol,
    std::function<stk::math::Vector3d(NODE)> coordinates,
    std::function<stk::math::Vector3d(NODE)> ownerCoordinates,
    std::function<bool(NODE, const std::vector<int> &)> nodeHasDomains,
    std::function<Face(NODE, NODE, NODE)> canonicalFace,
    std::function<std::vector<FaceCrossing>(const Face &)> uncapturedFaceCrossings)
    : myEdgeTol(edgeTol),
      myCoordinates(std::move(coordinates)),
      myOwnerCoordinates(std::move(ownerCoordinates)),
      myNodeHasDomains(std::move(nodeHasDomains)),
      myCanonicalFace(std::move(canonicalFace)),
      myUncapturedFaceCrossings(std::move(uncapturedFaceCrossings))
{
}

template <typename NODE>
bool SharpEdgeRecoveryBySwaps<NODE>::is_less_in_coordinate_order(const NODE a, const NODE b) const
{
  return less_by_coordinates_because_it_impacts_structure(myCoordinates(a), myCoordinates(b));
}

// The same for the edge whichever face it is found on: by coordinates, or by handle where they tie.
template <typename NODE>
typename SharpEdgeRecoveryBySwaps<NODE>::Edge
SharpEdgeRecoveryBySwaps<NODE>::edge_in_coordinate_order(const NODE a, const NODE b) const
{
  if (is_less_in_coordinate_order(b, a) || (!is_less_in_coordinate_order(a, b) && std::less<NODE>()(b, a)))
    return {{b, a}};
  return {{a, b}};
}

template <typename NODE>
std::vector<stk::math::Vector3d> SharpEdgeRecoveryBySwaps<NODE>::coordinates_in_coordinate_order(
    const std::vector<NODE> & nodes) const
{
  std::vector<stk::math::Vector3d> coords;
  for (const NODE node : nodes)
    coords.push_back(myCoordinates(node));
  std::sort(coords.begin(), coords.end(), LessByCoordinatesBecauseItImpactsStructure());
  return coords;
}

template <typename NODE>
std::vector<typename SharpEdgeRecoveryBySwaps<NODE>::FaceCrossing>
SharpEdgeRecoveryBySwaps<NODE>::fill_crossings_inside_face_and_crossed_edges(
    const Face & canonicalFaceNodes, std::vector<Edge> & crossedFaceEdges) const
{
  std::vector<FaceCrossing> insideCrossings;
  for (auto & crossing : myUncapturedFaceCrossings(canonicalFaceNodes))
  {
    const LocationOnTetFace location = locate_crossing_on_tet_face(crossing.faceCoords, myEdgeTol);
    if (location.location == TetFaceCrossingLocation::INSIDE)
      insideCrossings.push_back(crossing);
    else if (location.location == TetFaceCrossingLocation::NEAR_EDGE)
      crossedFaceEdges.push_back(edge_in_coordinate_order(
          canonicalFaceNodes[location.edgeNodeOrdinals[0]], canonicalFaceNodes[location.edgeNodeOrdinals[1]]));
  }
  return insideCrossings;
}

template <typename NODE>
bool SharpEdgeRecoveryBySwaps<NODE>::have_uncaptured_crossing_inside(const std::vector<Tet> & tets) const
{
  std::map<Face, int> numTetsWithFace;
  for (auto & tet : tets)
    for (auto & side : tetSideNodeOrdinals)
      ++numTetsWithFace[myCanonicalFace(tet[side[0]], tet[side[1]], tet[side[2]])];

  // A crossing near an edge on the outside of the tets is not changed by replacing them.  Every face
  // with an internal edge is internal, so each internal edge's crossing is found, on some face.
  std::set<Edge> outsideEdges;
  for (auto & [face, numTets] : numTetsWithFace)
    if (numTets == 1)
      for (int i = 0; i < 3; ++i)
        outsideEdges.insert(edge_in_coordinate_order(face[i], face[(i + 1) % 3]));

  for (auto & [face, numTets] : numTetsWithFace)
  {
    if (numTets != 2) continue;
    std::vector<Edge> crossedFaceEdges;
    if (!fill_crossings_inside_face_and_crossed_edges(face, crossedFaceEdges).empty()) return true;
    for (auto & edge : crossedFaceEdges)
      if (!outsideEdges.contains(edge)) return true;
  }
  return false;
}

template <typename NODE>
std::optional<typename SharpEdgeRecoveryBySwaps<NODE>::Swap>
SharpEdgeRecoveryBySwaps<NODE>::find_face_swap(const Face & canonicalFaceNodes,
    const std::vector<FaceCrossing> & crossings,
    std::vector<std::pair<size_t, NODE>> tetsAndOppositeNodes) const
{
  std::sort(tetsAndOppositeNodes.begin(),
      tetsAndOppositeNodes.end(),
      [&](const auto & a, const auto & b) { return is_less_in_coordinate_order(a.second, b.second); });
  const NODE p = tetsAndOppositeNodes[0].second;
  const NODE q = tetsAndOppositeNodes[1].second;
  const stk::math::Vector3d pOwnerCoords = myOwnerCoordinates(p);
  const stk::math::Vector3d pqOwnerCoords = myOwnerCoordinates(q) - pOwnerCoords;
  const auto is_crossing_of_segment = [&](const FaceCrossing & crossing)
  {
    if (!myNodeHasDomains(p, crossing.sortedDomains) || !myNodeHasDomains(q, crossing.sortedDomains))
      return false;
    const stk::math::Vector3d & faceCoords = crossing.faceCoords;
    const stk::math::Vector3d ownerCoords =
        (1. - faceCoords[0] - faceCoords[1]) * myOwnerCoordinates(canonicalFaceNodes[0]) +
        faceCoords[0] * myOwnerCoordinates(canonicalFaceNodes[1]) +
        faceCoords[1] * myOwnerCoordinates(canonicalFaceNodes[2]);
    return Cross(ownerCoords - pOwnerCoords, pqOwnerCoords).length() <= myEdgeTol * pqOwnerCoords.length();
  };
  if (std::none_of(crossings.begin(), crossings.end(), is_crossing_of_segment)) return std::nullopt;

  std::vector<Tet> newTets = tets_of_tet_face_swap(canonicalFaceNodes, p, q);
  if (!orient_tets_positively(newTets, myCoordinates) || have_uncaptured_crossing_inside(newTets))
    return std::nullopt;
  return Swap{{tetsAndOppositeNodes[0].first, tetsAndOppositeNodes[1].first}, std::move(newTets)};
}

template <typename NODE>
std::optional<typename SharpEdgeRecoveryBySwaps<NODE>::Swap> SharpEdgeRecoveryBySwaps<NODE>::find_edge_swap(
    const std::vector<Tet> & tets, const NODE edgeNode0, const NODE edgeNode1) const
{
  std::vector<size_t> tetsAroundEdgeIndices;
  std::vector<Tet> tetsAroundEdge;
  for (size_t iTet = 0; iTet < tets.size(); ++iTet)
  {
    const Tet & tet = tets[iTet];
    if (std::find(tet.begin(), tet.end(), edgeNode0) == tet.end() ||
        std::find(tet.begin(), tet.end(), edgeNode1) == tet.end())
      continue;
    tetsAroundEdgeIndices.push_back(iTet);
    tetsAroundEdge.push_back(tet);
  }

  std::vector<NODE> ring;
  if (!fill_ring_of_nodes_around_tet_edge(tetsAroundEdge, edgeNode0, edgeNode1, ring)) return std::nullopt;
  const int iConfig = best_tet_edge_swap_configuration(ring,
      edgeNode0,
      edgeNode1,
      myCoordinates,
      [&](const std::vector<Tet> & configTets, size_t) { return !have_uncaptured_crossing_inside(configTets); });
  if (iConfig < 0) return std::nullopt;

  std::vector<Tet> newTets = tets_of_tet_edge_swap_configuration(ring, edgeNode0, edgeNode1, iConfig);
  orient_tets_positively(newTets, myCoordinates);
  return Swap{std::move(tetsAroundEdgeIndices), std::move(newTets)};
}

template <typename NODE>
std::optional<typename SharpEdgeRecoveryBySwaps<NODE>::Swap>
SharpEdgeRecoveryBySwaps<NODE>::find_swap(const std::vector<Tet> & tets) const
{
  // An edge of a face on the outside of the tets is on the outside too, so the crossings near
  // internal edges are all found on internal faces.  An edge on the outside has no closed ring of
  // tets, and is not swapped.
  struct CrossedFaceOrEdge
  {
    std::vector<stk::math::Vector3d> sortedCoords;
    std::optional<size_t> iFace;
    Face canonicalFaceNodes;
    std::vector<FaceCrossing> faceCrossings;
    Edge edgeNodes;
  };
  const auto faces = unique_faces_of_tets(tets);
  std::vector<CrossedFaceOrEdge> crossedFacesAndEdges;
  std::set<Edge> crossedEdges;
  for (size_t iFace = 0; iFace < faces.size(); ++iFace)
  {
    if (faces[iFace].tetsAndOppositeNodes.size() != 2) continue;
    const Face & faceNodes = faces[iFace].nodes;
    const Face canonicalFaceNodes = myCanonicalFace(faceNodes[0], faceNodes[1], faceNodes[2]);
    std::vector<Edge> crossedFaceEdges;
    std::vector<FaceCrossing> crossings =
        fill_crossings_inside_face_and_crossed_edges(canonicalFaceNodes, crossedFaceEdges);
    if (!crossings.empty())
      crossedFacesAndEdges.push_back({coordinates_in_coordinate_order(
                                          std::vector<NODE>(canonicalFaceNodes.begin(), canonicalFaceNodes.end())),
          iFace,
          canonicalFaceNodes,
          std::move(crossings),
          {}});
    for (auto & edge : crossedFaceEdges)
      if (crossedEdges.insert(edge).second)
        crossedFacesAndEdges.push_back(
            {coordinates_in_coordinate_order({edge[0], edge[1]}), std::nullopt, {}, {}, edge});
  }
  std::sort(crossedFacesAndEdges.begin(),
      crossedFacesAndEdges.end(),
      [](const CrossedFaceOrEdge & a, const CrossedFaceOrEdge & b)
      {
        return std::lexicographical_compare(a.sortedCoords.begin(),
            a.sortedCoords.end(),
            b.sortedCoords.begin(),
            b.sortedCoords.end(),
            LessByCoordinatesBecauseItImpactsStructure());
      });

  for (const CrossedFaceOrEdge & crossed : crossedFacesAndEdges)
  {
    if (!crossed.iFace)
    {
      if (auto swap = find_edge_swap(tets, crossed.edgeNodes[0], crossed.edgeNodes[1])) return swap;
      continue;
    }
    if (auto swap = find_face_swap(crossed.canonicalFaceNodes,
            crossed.faceCrossings,
            faces[*crossed.iFace].tetsAndOppositeNodes))
      return swap;
    Face edgeOrderNodes = crossed.canonicalFaceNodes;
    std::sort(edgeOrderNodes.begin(),
        edgeOrderNodes.end(),
        [&](const NODE a, const NODE b) { return is_less_in_coordinate_order(a, b); });
    for (auto [i, j] : {std::pair{0, 1}, std::pair{0, 2}, std::pair{1, 2}})
      if (auto swap = find_edge_swap(tets, edgeOrderNodes[i], edgeOrderNodes[j])) return swap;
  }
  return std::nullopt;
}

template <typename NODE>
void SharpEdgeRecoveryBySwaps<NODE>::make_swaps(std::vector<Tet> & tets) const
{
  while (auto swap = find_swap(tets))
  {
    std::vector<size_t> oldTets = swap->oldTets;
    std::sort(oldTets.begin(), oldTets.end(), std::greater<size_t>());
    for (const size_t iTet : oldTets)
      tets.erase(tets.begin() + iTet);
    tets.insert(tets.end(), swap->newTets.begin(), swap->newTets.end());
  }
}

template class SharpEdgeRecoveryBySwaps<const SubElementNode *>;
template class SharpEdgeRecoveryBySwaps<int>;

} // namespace krino
