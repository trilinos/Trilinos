// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_DECOMP_AKRI_SHARPFEATURELEAFUTILS_HPP_
#define KRINO_KRINO_DECOMP_AKRI_SHARPFEATURELEAFUTILS_HPP_

#include <algorithm>
#include <array>
#include <map>
#include <utility>
#include <vector>

#include <stk_math/StkVector.hpp>

namespace krino {

// Helpers for cutting the sharp features of an element's leaves (Mesh_Element), free of the element
// so that they can be tested alone.

// The nodes of the sides of a tet, as in its topology.
inline constexpr std::array<std::array<int, 3>, 4> tetSideNodeOrdinals{
    {{{0, 1, 3}}, {{1, 2, 3}}, {{0, 3, 2}}, {{0, 2, 1}}}};

// The weights of the nodes of a simplex with numNodes nodes (a tet, a tri or an edge) at these
// parametric coordinates: {1 - p0 - p1 - ..., p0, p1, ...}.
std::vector<double> simplex_weights_at_parametric_coordinates(
    const stk::math::Vector3d & parametricCoords, const unsigned numNodes);

// The ordinals of the weights above tol, and those weights scaled to sum to 1: the entity of the
// simplex (a node, an edge, a face or the simplex itself) that a point with these weights is on.
void fill_entity_of_weights_above_tolerance(const std::vector<double> & weights,
    const double tol,
    std::vector<int> & entityOrdinals,
    std::vector<double> & entityWeights);

// Where on a tet face a crossing at these face parametric coordinates is: inside, within the edge
// tolerance in the face's weights of an edge of the face (edgeNodeOrdinals are that edge's nodes in
// the face), or of a node of it.
enum class TetFaceCrossingLocation { INSIDE, NEAR_EDGE, NEAR_NODE };
struct LocationOnTetFace
{
  TetFaceCrossingLocation location;
  std::array<int, 2> edgeNodeOrdinals;
};
LocationOnTetFace locate_crossing_on_tet_face(const stk::math::Vector3d & faceCoords, const double edgeTol);

// Whether one of the nodes, with these domains, captures the point: one of higher rank, whose
// domains contain the point's and are more, counting repeats.  So a corner {0,0,0} captures the
// crossings {0,0} of sharp edges, but a node on a sharp edge captures no crossing, though all sharp
// edges of one surface have the same domains.
bool some_node_captures_point_domains(const std::vector<const std::vector<int> *> & nodesSortedDomains,
    const std::vector<int> & sortedPointDomains);

// Whether a point on a tet face, at these face parametric coordinates, is captured by the face, its
// nodes with these domains: by a node of higher rank (some_node_captures_point_domains), or, within
// the edge tolerance in the face's weights, by a node that already has the point's domains, as where
// a sharp edge crosses the face at a node on the edge, or by an edge whose two nodes have them, a
// chord of the sharp edge.  Near such a node or chord the sharp edge runs nearly in the face, so
// where it crosses the face moves far for a roundoff change of the nodes.
bool tet_face_point_is_captured(const std::array<const std::vector<int> *, 3> & faceNodesSortedDomains,
    const stk::math::Vector3d & faceCoords,
    const std::vector<int> & sortedPointDomains,
    const double edgeTol);

// A face of some tets, with its nodes as in the first tet found with it (tetSideNodeOrdinals), and
// each tet with it, by index, with the node of that tet opposite it.
template <typename NODE>
struct FaceOfTets
{
  std::array<NODE, 3> nodes;
  std::vector<std::pair<size_t, NODE>> tetsAndOppositeNodes;
};

// The faces of the tets (each a container of 4 nodes), each once, in the order found: by tet, then
// by side (tetSideNodeOrdinals).
template <typename TET>
std::vector<FaceOfTets<typename TET::value_type>> unique_faces_of_tets(const std::vector<TET> & tets)
{
  using NODE = typename TET::value_type;
  constexpr std::array<int, 4> oppositeNodeOrdinalOfSide{{2, 0, 1, 3}};
  std::map<std::array<NODE, 3>, size_t> faceIndices;
  std::vector<FaceOfTets<NODE>> faces;
  for (size_t iTet = 0; iTet < tets.size(); ++iTet)
  {
    const TET & tet = tets[iTet];
    for (int iSide = 0; iSide < 4; ++iSide)
    {
      const auto & side = tetSideNodeOrdinals[iSide];
      const std::array<NODE, 3> faceNodes{{tet[side[0]], tet[side[1]], tet[side[2]]}};
      std::array<NODE, 3> sortedFaceNodes = faceNodes;
      std::sort(sortedFaceNodes.begin(), sortedFaceNodes.end());
      const auto [iter, inserted] = faceIndices.emplace(sortedFaceNodes, faces.size());
      if (inserted) faces.push_back({faceNodes, {}});
      faces[iter->second].tetsAndOppositeNodes.emplace_back(iTet, tet[oppositeNodeOrdinalOfSide[iSide]]);
    }
  }
  return faces;
}

} // namespace krino

#endif /* KRINO_KRINO_DECOMP_AKRI_SHARPFEATURELEAFUTILS_HPP_ */
