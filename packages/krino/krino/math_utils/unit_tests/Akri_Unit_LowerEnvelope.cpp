// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_LowerEnvelope.hpp>
#include <gtest/gtest.h>


namespace krino
{

TEST(Lower_Envelope, Find_Crossing_Segment)
{
  double s = 0.0;

  ASSERT_TRUE(find_lower_envelope_crossing_point({{{0., 1.}, {0., 0.}}}, s));
  EXPECT_DOUBLE_EQ(1.0, s);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{{0., 1.}, {1., 0.}}}, s));
  EXPECT_DOUBLE_EQ(0.5, s);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{{0., 1.}, {0.5, 0.}}}, s));
  EXPECT_DOUBLE_EQ(2./3., s);
}

TEST(Lower_Envelope, Find_Crossing_Triangle)
{
  std::array<double,2> pt = {{0., 0.}};

  const std::vector<double> phi000 = {0., 0., 0.};
  const std::vector<double> phi011 = {0., 1., 1.};
  const std::vector<double> phi101 = {1., 0., 1.};
  const std::vector<double> phi110 = {1., 1., 0.};

  EXPECT_FALSE(find_lower_envelope_crossing_point({{phi110, phi110, phi101}}, pt));

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi110, phi011, phi101}}, pt));
  std::array<double,2> goldPt = {{1./3., 1./3.}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{{1., 1., 0.75}, {0., 1., 0.75}, {1., 0., 0.75}}}, pt));
  goldPt = {{0.25, 0.25}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi000, phi011, phi101}}, pt));
  goldPt = {{0., 0.}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi011, phi000, phi101}}, pt));
  goldPt = {{1., 0.}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi011, phi101, phi000}}, pt));
  goldPt = {{0., 1.}};
  EXPECT_EQ(goldPt, pt);
}

TEST(Lower_Envelope, Find_Crossing_Tetrahedron)
{
  std::array<double,3> pt = {{0., 0., 0.}};

  const std::vector<double> phi0000 = {0., 0., 0., 0.};
  const std::vector<double> phi1110 = {1., 1., 1., 0.};
  const std::vector<double> phi1101 = {1., 1., 0., 1.};
  const std::vector<double> phi1011 = {1., 0., 1., 1.};
  const std::vector<double> phi0111 = {0., 1., 1., 1.};

  EXPECT_FALSE(find_lower_envelope_crossing_point({{phi1110, phi1110, phi1011, phi0111}}, pt));

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi1110, phi1101, phi1011, phi0111}}, pt));
  std::array<double,3> goldPt = {{0.25, 0.25, 0.25}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi0000, phi0111, phi1011, phi1101}}, pt));
  goldPt = {{0., 0., 0.}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi0111, phi0000, phi1011, phi1101}}, pt));
  goldPt = {{1., 0., 0.}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi0111, phi1011, phi0000, phi1101}}, pt));
  goldPt = {{0., 1., 0.}};
  EXPECT_EQ(goldPt, pt);

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi0111, phi1011, phi1101, phi0000}}, pt));
  goldPt = {{0., 0., 1.}};
  EXPECT_EQ(goldPt, pt);
}


TEST(Lower_Envelope, Two_LS_Bug)
{
  Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope({0., 0.}, {1., 0.});
  ASSERT_EQ(1u, envelope.size());
  EXPECT_EQ(1, envelope[0].ls_index());
  EXPECT_DOUBLE_EQ(0.0, envelope[0].left_endpoint());
  EXPECT_DOUBLE_EQ(1., envelope[0].right_endpoint());
}

TEST(Lower_Envelope, Two_LS_Lower_Envelope)
{
  {
    // 2 level sets, infinitesimal crossing on left by convention that highest level set index wins
    Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope({0., 0.}, {0.25, 1.});
    ASSERT_EQ(2u, envelope.size());
    const LS_Segment & env1 = envelope[0];
    const LS_Segment & env2 = envelope[1];
    EXPECT_EQ(1, env1.ls_index());
    EXPECT_EQ(0, env2.ls_index());
    EXPECT_DOUBLE_EQ(0.0, env1.left_endpoint());
    EXPECT_DOUBLE_EQ(0.0, env1.right_endpoint());
    EXPECT_DOUBLE_EQ(0.0, env2.left_endpoint());
    EXPECT_DOUBLE_EQ(1.0, env2.right_endpoint());
  }

  {
    // 2 level sets, infinitesimal crossing on right by convention that highest level set index wins
    Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope({0.3, 0.6}, {0., 0.});
    ASSERT_EQ(2u, envelope.size());
    const LS_Segment & env1 = envelope[0];
    const LS_Segment & env2 = envelope[1];
    EXPECT_EQ(0, env1.ls_index());
    EXPECT_EQ(1, env2.ls_index());
    EXPECT_DOUBLE_EQ(0.0, env1.left_endpoint());
    EXPECT_DOUBLE_EQ(1.0, env1.right_endpoint());
    EXPECT_DOUBLE_EQ(1.0, env2.left_endpoint());
    EXPECT_DOUBLE_EQ(1.0, env2.right_endpoint());
  }
}

void expect_good_edge(const Segment_Vector & envelope)
{
  for (size_t i=0; i<envelope.size(); ++i)
  {
    if (i==0 || i==envelope.size()-1)
    {
      EXPECT_TRUE(envelope[i].length()==0 || envelope[i].length() >= SegmentLowerEnvelope::MinSize()) << "Unexpected infinitesimal, non-zero segment at end of edge " << envelope;
    }
    else
    {
      EXPECT_GE(envelope[i].length(), SegmentLowerEnvelope::MinSize()) << "Unexpected internal infinitesimal internal segment on edge " << envelope;
    }

    if (i<envelope.size()-1)
    {
      EXPECT_EQ(envelope[i].right_endpoint(), envelope[i+1].left_endpoint()) << "Missing portion of edge between segments " << i << " and " << i+1 << " on edge " << envelope;
      EXPECT_NE(envelope[i].ls_index(), envelope[i+1].ls_index()) << "Identical phase on consecutive segments " << i << " and " << i+1 << " on edge " << envelope;
    }
  }
}

void expect_end_segments_to_match_end_phases(const Segment_Vector & envelope, const std::vector<double> & phi0, const std::vector<double> & phi1)
{
  ASSERT_FALSE(envelope.empty());
  EXPECT_EQ(get_min(phi0), envelope.front().ls_index()) << "Front segment of edge does not match phase at node on edge " << envelope;
  EXPECT_EQ(get_min(phi1), envelope.back().ls_index()) << "Back segment of edge does not match phase at node on edge " << envelope;
}

TEST(Lower_Envelope, Three_LS_Lower_Envelope)
{
  {
    // 3 level sets, edge is B from 0 to 0.5, C from 0.5 to 1.0 and
    // 1st level set intersects both others at 0.5 but is never the lowest
    Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope({0.5, 0., 1.}, {0.5, 1., 0.});
    ASSERT_EQ(2u, envelope.size());
    const LS_Segment & env1 = envelope[0];
    const LS_Segment & env2 = envelope[1];
    EXPECT_EQ(1, env1.ls_index());
    EXPECT_EQ(2, env2.ls_index());
    EXPECT_DOUBLE_EQ(0.0, env1.left_endpoint());
    EXPECT_DOUBLE_EQ(0.5, env1.right_endpoint());
    EXPECT_DOUBLE_EQ(0.5, env2.left_endpoint());
    EXPECT_DOUBLE_EQ(1., env2.right_endpoint());

    expect_good_edge(envelope);
  }

  {
    // 3 level sets, left is A, middle C, right B
    Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope({0., 1., 0.25}, {1., 0., 0.25});
    ASSERT_EQ(3u, envelope.size());
    const LS_Segment & env1 = envelope[0];
    const LS_Segment & env2 = envelope[1];
    const LS_Segment & env3 = envelope[2];
    EXPECT_EQ(0, env1.ls_index());
    EXPECT_EQ(2, env2.ls_index());
    EXPECT_EQ(1, env3.ls_index());
    EXPECT_DOUBLE_EQ(0.0, env1.left_endpoint());
    EXPECT_DOUBLE_EQ(0.25, env1.right_endpoint());
    EXPECT_DOUBLE_EQ(0.25, env2.left_endpoint());
    EXPECT_DOUBLE_EQ(0.75, env2.right_endpoint());
    EXPECT_DOUBLE_EQ(0.75, env3.left_endpoint());
    EXPECT_DOUBLE_EQ(1.0, env3.right_endpoint());

    expect_good_edge(envelope);
  }

  {
    // 3 level sets with infinitesmal transitions near ends
    const std::vector<double> phi0 = {1.56125e-17, -4.77049e-18, 0.};
    const std::vector<double> phi1 = {-0.0417425, 0.0477226, 0.};

    {
      Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(phi0, phi1);

      expect_good_edge(envelope);
      expect_end_segments_to_match_end_phases(envelope, phi0, phi1);
    }
    {
      Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(phi1, phi0);

      expect_good_edge(envelope);
      expect_end_segments_to_match_end_phases(envelope, phi1, phi0);
    }
  }

  {
    // 3 level sets with infinitesmal transitions near middle
    const std::vector<double> phi0 = {-1., (1.+1.e-12), 0.};
    const std::vector<double> phi1 = {(1.+1.e-12), -1., 0.};

    {
      Segment_Vector envelope01 = SegmentLowerEnvelope::find_lower_envelope(phi0, phi1);
      expect_good_edge(envelope01);
      expect_end_segments_to_match_end_phases(envelope01, phi0, phi1);
      EXPECT_EQ(2u, envelope01.size());

      Segment_Vector envelope10 = SegmentLowerEnvelope::find_lower_envelope(phi1, phi0);
      expect_good_edge(envelope10);
      expect_end_segments_to_match_end_phases(envelope10, phi1, phi0);
      EXPECT_EQ(2u, envelope10.size());

      EXPECT_DOUBLE_EQ(envelope01[0].right_endpoint(), envelope10[0].right_endpoint());
    }
  }
}

void expect_segments_lengths(const Segment_Vector & segments, const std::vector<double> goldSegmentLengthsByLS)
{
  for (auto && segment : segments)
  {
    ASSERT_TRUE(segment.ls_index() < (int)goldSegmentLengthsByLS.size());
    EXPECT_DOUBLE_EQ(goldSegmentLengthsByLS[segment.ls_index()], segment.length());
  }
}

TEST(Lower_Envelope,Sensitive_LS)
{
  std::array<std::vector<double>,2> phi = {{ {-1.e-17,-2.e-17,1.0,0}, {2.e-17,1.0,-2.e-17,0} }};
  const double goldSegmentLengthForPhi0 = phi[0][0]/(phi[0][0]-phi[1][0]);
  const std::vector<double> goldSegmentLengthsByLS = { goldSegmentLengthForPhi0, 0., 0., 1.-goldSegmentLengthForPhi0 };

  Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(phi[0], phi[1]);
  expect_segments_lengths(envelope, goldSegmentLengthsByLS);
}

TEST(Lower_Envelope, Find_Crossing_Position_And_Sign)
{
  // 2 nodes at positions {0,1}, 2 level sets. The (lsA, lsB) pair selects which
  // level sets to difference: key_isovar = phi[min(lsA,lsB)] - phi[max(lsA,lsB)].
  const std::vector<double> pos = {0., 1.};

  {
    // isovar goes -1 -> +1: crossing at midpoint, sign taken from the right node (+1).
    const std::vector<std::vector<double>> phi = {{0., 1.}, {1., 0.}};
    double crossingPos = 0.;
    int crossingSign = 0;
    EXPECT_TRUE(find_crossing_position_and_sign(0, 1, pos, phi, crossingPos, crossingSign));
    EXPECT_DOUBLE_EQ(0.5, crossingPos);
    EXPECT_EQ(1, crossingSign);
  }

  {
    // Reversed: isovar goes +1 -> -1, crossing still at midpoint but sign is now -1.
    const std::vector<std::vector<double>> phi = {{1., 0.}, {0., 1.}};
    double crossingPos = 0.;
    int crossingSign = 0;
    EXPECT_TRUE(find_crossing_position_and_sign(0, 1, pos, phi, crossingPos, crossingSign));
    EXPECT_DOUBLE_EQ(0.5, crossingPos);
    EXPECT_EQ(-1, crossingSign);
  }

  {
    // Asymmetric: isovar -1 -> +0.5 crosses at 2/3 of the way.
    const std::vector<std::vector<double>> phi = {{0., 1.}, {0.5, 0.}};
    double crossingPos = 0.;
    int crossingSign = 0;
    EXPECT_TRUE(find_crossing_position_and_sign(0, 1, pos, phi, crossingPos, crossingSign));
    EXPECT_DOUBLE_EQ(2./3., crossingPos);
    EXPECT_EQ(1, crossingSign);
  }

  {
    // No crossing (isovar stays negative): returns false, leaves crossingPos
    // untouched, and reports the common end-node sign.
    const std::vector<std::vector<double>> phi = {{0., 1.}, {0., 2.}};
    double crossingPos = -42.;
    int crossingSign = 0;
    EXPECT_FALSE(find_crossing_position_and_sign(0, 1, pos, phi, crossingPos, crossingSign));
    EXPECT_DOUBLE_EQ(-42., crossingPos) << "crossingPos must be left unchanged when there is no crossing";
    EXPECT_EQ(-1, crossingSign);
  }

  {
    // Non-zero key indices exercise the min/max column indexing:
    // (1,2) differences columns 1 and 2 and ignores column 0.
    const std::vector<std::vector<double>> phi = {{5., 0., 1.}, {5., 1., 0.}};
    double crossingPos = 0.;
    int crossingSign = 0;
    EXPECT_TRUE(find_crossing_position_and_sign(1, 2, pos, phi, crossingPos, crossingSign));
    EXPECT_DOUBLE_EQ(0.5, crossingPos);
    EXPECT_EQ(1, crossingSign);
  }
}

TEST(Lower_Envelope, Find_Crossing_Position_And_Sign_Piecewise)
{
  // 3 nodes: the interior node bends the isovar so the crossing lands at 2/3
  // rather than the linear-endpoint midpoint. key_isovar = {-1, -0.5, +1}, so
  // the crossing is found in the second interval.
  const std::vector<double> pos = {0., 0.5, 1.};
  const std::vector<std::vector<double>> phi = {{0., 1.}, {0., 0.5}, {1., 0.}};

  double crossingPos = 0.;
  int crossingSign = 0;
  EXPECT_TRUE(find_crossing_position_and_sign(0, 1, pos, phi, crossingPos, crossingSign));
  EXPECT_DOUBLE_EQ(2./3., crossingPos);
  EXPECT_EQ(1, crossingSign);
}

TEST(Lower_Envelope, Piecewise_Lower_Envelope_RelocatesCrossing)
{
  // The linear envelope from the two endpoints alone would cross at 0.5, but the
  // interior sample pulls the true crossing to 2/3. This drives
  // adjust_piecewise_linear_positions, which the 2-point overloads never reach.
  const std::vector<double> pos = {0., 0.5, 1.};
  const std::vector<std::vector<double>> phi = {{0., 1.}, {0., 0.5}, {1., 0.}};

  const Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(pos, phi);

  ASSERT_EQ(2u, envelope.size());
  EXPECT_EQ(0, envelope[0].ls_index());
  EXPECT_EQ(1, envelope[1].ls_index());
  EXPECT_DOUBLE_EQ(0.0, envelope[0].left_endpoint());
  EXPECT_DOUBLE_EQ(2./3., envelope[0].right_endpoint());
  EXPECT_DOUBLE_EQ(2./3., envelope[1].left_endpoint());
  EXPECT_DOUBLE_EQ(1.0, envelope[1].right_endpoint());

  expect_good_edge(envelope);
  expect_end_segments_to_match_end_phases(envelope, phi.front(), phi.back());
}

TEST(Lower_Envelope, Piecewise_Lower_Envelope_ThreePhase)
{
  // Three-phase piecewise edge (A, then C dipping lowest in the middle, then B).
  // Exercises the size > 2 branch of adjust_piecewise_linear_positions. The
  // interior node is the exact linear midpoint, so the two crossings stay at
  // 0.25 and 0.75 and remain correctly ordered (no inverted segment).
  const std::vector<double> pos = {0., 0.5, 1.};
  const std::vector<std::vector<double>> phi = {{0., 1., 0.25}, {0.5, 0.5, 0.25}, {1., 0., 0.25}};

  const Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(pos, phi);

  ASSERT_EQ(3u, envelope.size());
  EXPECT_EQ(0, envelope[0].ls_index());
  EXPECT_EQ(2, envelope[1].ls_index());
  EXPECT_EQ(1, envelope[2].ls_index());
  EXPECT_DOUBLE_EQ(0.25, envelope[0].right_endpoint());
  EXPECT_DOUBLE_EQ(0.75, envelope[2].left_endpoint());

  expect_good_edge(envelope);
  expect_end_segments_to_match_end_phases(envelope, phi.front(), phi.back());
}

TEST(Lower_Envelope, Explicit_Endpoints_Lower_Envelope)
{
  // The explicit-endpoint overload places the crossing in absolute coordinates
  // between x0 and x1 rather than on [0,1].
  const Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(2., {0., 1.}, 5., {1., 0.});

  ASSERT_EQ(2u, envelope.size());
  EXPECT_EQ(0, envelope[0].ls_index());
  EXPECT_EQ(1, envelope[1].ls_index());
  EXPECT_DOUBLE_EQ(2.0, envelope[0].left_endpoint());
  EXPECT_DOUBLE_EQ(3.5, envelope[0].right_endpoint());
  EXPECT_DOUBLE_EQ(3.5, envelope[1].left_endpoint());
  EXPECT_DOUBLE_EQ(5.0, envelope[1].right_endpoint());
}

TEST(Lower_Envelope, Explicit_Endpoints_InfinitesimalCuts_UseEdgeEndpoints)
{
  // Regression test for the infinitesimal-cut handling assuming a [0,1] domain.
  // These are the Two_LS_Lower_Envelope cases shifted onto [2,5]: the zero-length
  // end segments must sit at the actual edge endpoints (2 and 5), not at 0 and 1.
  {
    // Infinitesimal crossing on the left: the leading zero-length segment must be
    // stamped at x0 == 2, not at 0.
    const Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(2., {0., 0.}, 5., {0.25, 1.});
    ASSERT_EQ(2u, envelope.size());
    EXPECT_EQ(1, envelope[0].ls_index());
    EXPECT_EQ(0, envelope[1].ls_index());
    EXPECT_DOUBLE_EQ(2.0, envelope[0].left_endpoint());
    EXPECT_DOUBLE_EQ(2.0, envelope[0].right_endpoint());
    EXPECT_DOUBLE_EQ(2.0, envelope[1].left_endpoint());
    EXPECT_DOUBLE_EQ(5.0, envelope[1].right_endpoint());
  }

  {
    // Infinitesimal crossing on the right: the trailing zero-length segment must be
    // stamped at x1 == 5. The old x1==1. guard silently dropped it for xEnd != 1.
    const Segment_Vector envelope = SegmentLowerEnvelope::find_lower_envelope(2., {0.3, 0.6}, 5., {0., 0.});
    ASSERT_EQ(2u, envelope.size());
    EXPECT_EQ(0, envelope[0].ls_index());
    EXPECT_EQ(1, envelope[1].ls_index());
    EXPECT_DOUBLE_EQ(2.0, envelope[0].left_endpoint());
    EXPECT_DOUBLE_EQ(5.0, envelope[0].right_endpoint());
    EXPECT_DOUBLE_EQ(5.0, envelope[1].left_endpoint());
    EXPECT_DOUBLE_EQ(5.0, envelope[1].right_endpoint());
  }
}

TEST(Lower_Envelope, Find_Crossing_Tetrahedron_OutsideSimplex_ReturnsFalse)
{
  // Four distinct per-node minima (so the quadruple-point solve is attempted), but
  // level set 0 is low enough everywhere that the all-phases-equal point lands
  // outside the reference tet (barycentric sum > 1). Must report no valid crossing.
  std::array<double,3> pt = {{-1., -1., -1.}};

  const std::vector<double> phi0 = {0.0, 1.0, 1.0, 1.0}; // min at LS 0
  const std::vector<double> phi1 = {0.4, 0.0, 1.0, 1.0}; // min at LS 1
  const std::vector<double> phi2 = {0.4, 1.0, 0.0, 1.0}; // min at LS 2
  const std::vector<double> phi3 = {0.4, 1.0, 1.0, 0.0}; // min at LS 3

  EXPECT_FALSE(find_lower_envelope_crossing_point({{phi0, phi1, phi2, phi3}}, pt));
}

} // namespace krino
