// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <array>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_SharpEdgeRecoveryBySwaps.hpp>
#include <Akri_SharpFeatureIntersectionPoints.hpp>
#include <Akri_SharpFeatureLeafUtils.hpp>
#include <Akri_TetSwaps.hpp>

namespace krino
{

namespace
{

using Recovery = SharpEdgeRecoveryBySwaps<int>;
using Tet = Recovery::Tet;

const std::vector<int> sharpEdgeDomains{0, 0};

// Tets of nodes with these coordinates (their own owner coordinates), and a straight sharp edge
// from node segmentBegin to node segmentEnd, the nodes with its domains.  A face's crossing is where
// the segment crosses its plane, inside it, unless the face captures it (tet_face_point_is_captured).
class SegmentInTets
{
public:
  SegmentInTets(std::vector<stk::math::Vector3d> coords, const int segmentBegin, const int segmentEnd)
      : myCoords(std::move(coords)), myDomains(myCoords.size()), myBegin(segmentBegin), myEnd(segmentEnd)
  {
    myDomains[segmentBegin] = sharpEdgeDomains;
    myDomains[segmentEnd] = sharpEdgeDomains;
  }

  Recovery recovery() const
  {
    return Recovery(edgeTol,
        [this](const int node) { return myCoords[node]; },
        [this](const int node) { return myCoords[node]; },
        [this](const int node, const std::vector<int> & domains)
        { return node_domains_capture_point_domains(myDomains[node], domains); },
        [](const int a, const int b, const int c)
        {
          Recovery::Face face{{a, b, c}};
          std::sort(face.begin(), face.end());
          return face;
        },
        [this](const Recovery::Face & face) { return crossings(face); });
  }

  // The tets, each as its sorted vertex coordinates, so that two sets of tets can be compared
  // whatever their node numbers and order.
  std::set<std::array<std::array<double, 3>, 4>> geometry(const std::vector<Tet> & tets) const
  {
    std::set<std::array<std::array<double, 3>, 4>> tetGeometry;
    for (const Tet & tet : tets)
    {
      std::array<std::array<double, 3>, 4> vertices;
      for (int i = 0; i < 4; ++i)
        vertices[i] = {{myCoords[tet[i]][0], myCoords[tet[i]][1], myCoords[tet[i]][2]}};
      std::sort(vertices.begin(), vertices.end());
      tetGeometry.insert(vertices);
    }
    return tetGeometry;
  }

  double volume(const std::vector<Tet> & tets) const
  {
    double sum = 0.;
    for (const Tet & tet : tets)
      sum += tet_volume({{myCoords[tet[0]], myCoords[tet[1]], myCoords[tet[2]], myCoords[tet[3]]}});
    return sum;
  }

  static constexpr double edgeTol = 0.01;

private:
  std::vector<Recovery::FaceCrossing> crossings(const Recovery::Face & face) const
  {
    const stk::math::Vector3d & x0 = myCoords[face[0]];
    const stk::math::Vector3d e1 = myCoords[face[1]] - x0;
    const stk::math::Vector3d e2 = myCoords[face[2]] - x0;
    const stk::math::Vector3d a = myCoords[myBegin];
    const stk::math::Vector3d ab = myCoords[myEnd] - a;
    const stk::math::Vector3d normal = Cross(e1, e2);
    const double denom = Dot(normal, ab);
    if (std::abs(denom) <= 1.e-12 * normal.length() * ab.length()) return {}; // in the face's plane
    const double t = Dot(normal, x0 - a) / denom;
    if (t < -1.e-12 || t > 1. + 1.e-12) return {};
    const stk::math::Vector3d d = a + t * ab - x0;
    const double e11 = Dot(e1, e1), e12 = Dot(e1, e2), e22 = Dot(e2, e2);
    const double det = e11 * e22 - e12 * e12;
    const double w1 = (e22 * Dot(d, e1) - e12 * Dot(d, e2)) / det;
    const double w2 = (e11 * Dot(d, e2) - e12 * Dot(d, e1)) / det;
    if (w1 < -1.e-12 || w2 < -1.e-12 || w1 + w2 > 1. + 1.e-12) return {};
    const stk::math::Vector3d faceCoords(w1, w2, 0.);
    if (tet_face_point_is_captured({{&myDomains[face[0]], &myDomains[face[1]], &myDomains[face[2]]}},
            faceCoords,
            sharpEdgeDomains,
            edgeTol))
      return {};
    return {{faceCoords, sharpEdgeDomains}};
  }

  std::vector<stk::math::Vector3d> myCoords;
  std::vector<std::vector<int>> myDomains;
  int myBegin;
  int myEnd;
};

bool have_edge(const std::vector<Tet> & tets, const int node0, const int node1)
{
  return std::any_of(tets.begin(),
      tets.end(),
      [&](const Tet & tet)
      {
        return std::count(tet.begin(), tet.end(), node0) == 1 && std::count(tet.begin(), tet.end(), node1) == 1;
      });
}

// Two tets on the face of nodes 0, 1 and 2, with the segment between their opposite nodes 3 and 4
// crossing it inside.
SegmentInTets segment_across_shared_face()
{
  return SegmentInTets({stk::math::Vector3d(1., 0., 0.),
                           stk::math::Vector3d(-0.5, 0.8, 0.),
                           stk::math::Vector3d(-0.5, -0.8, 0.),
                           stk::math::Vector3d(0., 0., 1.),
                           stk::math::Vector3d(0., 0., -1.)},
      3,
      4);
}

// The four tets around the edge from node 4, at (0,0,-1), to node 5, at (0,0,1), with the ring of
// nodes 0 to 3 at (1,0,0), (0,1,0), (-1,0,0) and (0,-1,0), and the segment from node 1 to node 3
// crossing the edge: it is in the planes of the faces of the edge with nodes 1 and 3, and crosses
// those with nodes 0 and 2 on the edge, within roundoff of it.
SegmentInTets segment_across_internal_edge()
{
  return SegmentInTets({stk::math::Vector3d(1., 0., 0.),
                           stk::math::Vector3d(0., 1., 0.),
                           stk::math::Vector3d(-1., 0., 0.),
                           stk::math::Vector3d(0., -1., 0.),
                           stk::math::Vector3d(0., 0., -1.),
                           stk::math::Vector3d(0., 0., 1.)},
      1,
      3);
}

const std::vector<Tet> tetsAroundInternalEdge{{{0, 1, 4, 5}}, {{1, 2, 4, 5}}, {{2, 3, 4, 5}}, {{3, 0, 4, 5}}};

} // namespace

TEST(SharpEdgeRecoveryBySwaps, segmentAcrossSharedFace_twoThreeSwapMakesItAnEdge)
{
  const SegmentInTets segment = segment_across_shared_face();
  const Recovery recovery = segment.recovery();
  std::vector<Tet> tets{{{0, 1, 2, 3}}, {{0, 2, 1, 4}}};
  const double volume = segment.volume(tets);

  const std::optional<Recovery::Swap> swap = recovery.find_swap(tets);
  ASSERT_TRUE(swap);
  EXPECT_EQ((std::set<size_t>{0, 1}), std::set<size_t>(swap->oldTets.begin(), swap->oldTets.end()));
  ASSERT_EQ(3u, swap->newTets.size());
  for (const Tet & tet : swap->newTets)
  {
    EXPECT_EQ(1, std::count(tet.begin(), tet.end(), 3));
    EXPECT_EQ(1, std::count(tet.begin(), tet.end(), 4));
  }

  recovery.make_swaps(tets);
  EXPECT_TRUE(have_edge(tets, 3, 4));
  EXPECT_FALSE(recovery.find_swap(tets));
  EXPECT_NEAR(volume, segment.volume(tets), 1.e-14);
}

TEST(SharpEdgeRecoveryBySwaps, segmentOutsideTheTets_noSwap)
{
  const SegmentInTets segment({stk::math::Vector3d(1., 0., 0.),
                                  stk::math::Vector3d(-0.5, 0.8, 0.),
                                  stk::math::Vector3d(-0.5, -0.8, 0.),
                                  stk::math::Vector3d(0., 0., 1.),
                                  stk::math::Vector3d(0., 0., -1.),
                                  stk::math::Vector3d(3., 0., 1.),
                                  stk::math::Vector3d(3., 0., -1.)},
      5,
      6);
  EXPECT_FALSE(segment.recovery().find_swap({{{0, 1, 2, 3}}, {{0, 2, 1, 4}}}));
}

TEST(SharpEdgeRecoveryBySwaps, segmentAcrossInternalEdge_edgeSwapMakesItAnEdge)
{
  const SegmentInTets segment = segment_across_internal_edge();
  const Recovery recovery = segment.recovery();
  std::vector<Tet> tets = tetsAroundInternalEdge;
  const double volume = std::abs(segment.volume(tets));
  EXPECT_TRUE(recovery.have_uncaptured_crossing_inside(tets));

  const std::optional<Recovery::Swap> swap = recovery.find_swap(tets);
  ASSERT_TRUE(swap);
  EXPECT_EQ(4u, swap->oldTets.size());

  recovery.make_swaps(tets);
  EXPECT_TRUE(have_edge(tets, 1, 3));
  EXPECT_FALSE(have_edge(tets, 4, 5));
  EXPECT_FALSE(recovery.have_uncaptured_crossing_inside(tets));
  EXPECT_FALSE(recovery.find_swap(tets));
  EXPECT_NEAR(volume, segment.volume(tets), 1.e-14);
}

// The order of the tets and of their nodes does not change the swaps.
TEST(SharpEdgeRecoveryBySwaps, segmentAcrossInternalEdge_sameTetsWhateverTheirOrder)
{
  const SegmentInTets segment = segment_across_internal_edge();
  const Recovery recovery = segment.recovery();
  std::vector<Tet> tets = tetsAroundInternalEdge;
  recovery.make_swaps(tets);

  std::vector<Tet> reorderedTets{{{5, 4, 0, 3}}, {{5, 4, 2, 1}}, {{3, 5, 4, 2}}, {{1, 0, 5, 4}}};
  recovery.make_swaps(reorderedTets);
  EXPECT_EQ(segment.geometry(tets), segment.geometry(reorderedTets));
}

} // namespace krino
