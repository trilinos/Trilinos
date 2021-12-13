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

  ASSERT_TRUE(find_lower_envelope_crossing_point({{phi011, phi000, phi101}}, pt));
  goldPt = {{1., 0.}};

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

} // namespace krino
