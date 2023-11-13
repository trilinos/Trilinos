// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Phase_Support.hpp>

#include <Akri_Unit_Single_Element_Fixtures.hpp>

namespace krino
{

TEST(CDFEM_Parent_Edge_One_LS, Two_Nodes_No_Snapping)
{
  const InterfaceID iface(0,0);
  std::vector<std::vector<double> > nodes_isovar(2);
  nodes_isovar[0].resize(1);
  nodes_isovar[1].resize(1);
  nodes_isovar[0][0] = 1;
  nodes_isovar[1][0] = 1;

  const bool oneLSPerPhase = false;
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
  EXPECT_FALSE(edge.have_any_crossings());
  EXPECT_EQ(1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = -1;
  nodes_isovar[1][0] = -1;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_FALSE(edge.have_any_crossings());
  EXPECT_EQ(-1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = -1;
  nodes_isovar[1][0] = 2;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_TRUE(edge.have_any_crossings());
  EXPECT_TRUE(edge.have_crossing(iface));
  EXPECT_DOUBLE_EQ(1./3., edge.get_crossing_position(iface));
  EXPECT_EQ(1, edge.get_crossing_sign(iface));
}

TEST(CDFEM_Parent_Edge_One_LS, Three_Nodes_No_Snapping)
{
  std::vector<std::vector<double> > nodes_isovar(3);
  nodes_isovar[0].resize(1);
  nodes_isovar[1].resize(1);
  nodes_isovar[2].resize(1);
  nodes_isovar[0][0] = 1;
  nodes_isovar[1][0] = 1;
  nodes_isovar[2][0] = 1;
  const InterfaceID iface(0,0);

  const bool oneLSPerPhase = false;
  CDFEM_Parent_Edge edge(oneLSPerPhase, {0.0,0.5,1.0}, nodes_isovar);
  EXPECT_FALSE(edge.have_any_crossings());
  EXPECT_EQ(1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = -1;
  nodes_isovar[1][0] = -1;
  nodes_isovar[2][0] = -1;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_FALSE(edge.have_any_crossings());
  EXPECT_EQ(-1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = -1;
  nodes_isovar[1][0] = 1;
  nodes_isovar[2][0] = -1;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_FALSE(edge.have_any_crossings());
  EXPECT_EQ(-1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = 1;
  nodes_isovar[1][0] = -1;
  nodes_isovar[2][0] = 1;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_FALSE(edge.have_any_crossings());
  EXPECT_EQ(1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = -1;
  nodes_isovar[1][0] = 2;
  nodes_isovar[2][0] = 4;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_TRUE(edge.have_any_crossings());
  const InterfaceID interface_id(0,0);
  EXPECT_TRUE(edge.have_crossing(interface_id));
  EXPECT_DOUBLE_EQ(1./6., edge.get_crossing_position(interface_id));
  EXPECT_EQ(1, edge.get_crossing_sign(iface));

  nodes_isovar[0][0] = -1;
  nodes_isovar[1][0] = -1;
  nodes_isovar[2][0] = 1;
  edge.find_crossings(oneLSPerPhase, nodes_isovar);
  EXPECT_TRUE(edge.have_any_crossings());
  EXPECT_TRUE(edge.have_crossing(interface_id));
  EXPECT_DOUBLE_EQ(3./4., edge.get_crossing_position(interface_id));
  EXPECT_EQ(1, edge.get_crossing_sign(iface));
}

TEST(CDFEM_Parent_Edge_Two_LS, Two_Nodes_No_Snapping)
{
  const bool oneLSPerPhase = true;

  const InterfaceID iface(0,1);

  std::vector<std::vector<double> > nodes_isovar(2);
  nodes_isovar[0].resize(2);
  nodes_isovar[1].resize(2);

  // crossing_sign convention, iface.first phase corresponds to (-) for one LS=one interface case
  // iface.first lower everywhere
  {
    nodes_isovar[0][0] = 0;
    nodes_isovar[0][1] = 1;
    nodes_isovar[1][0] = 0;
    nodes_isovar[1][1] = 1;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_EQ(0, edge.get_uncrossed_phase());
  }

  // iface.second lower everywhere
  {
    nodes_isovar[0][0] = 2;
    nodes_isovar[0][1] = 1;
    nodes_isovar[1][0] = 2;
    nodes_isovar[1][1] = 1;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_EQ(1, edge.get_uncrossed_phase());
  }

  // first lower at node 0, second lower at node 1
  {
    nodes_isovar[0][0] = 0.5;
    nodes_isovar[0][1] = 1;
    nodes_isovar[1][0] = 1.5;
    nodes_isovar[1][1] = 1;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface));
    EXPECT_DOUBLE_EQ(0.5, edge.get_crossing_position(iface));
    EXPECT_EQ(1, edge.get_crossing_sign(iface));
  }

  // second lower at node 0, first lower at node 1
  {
    nodes_isovar[0][0] = 1.5;
    nodes_isovar[0][1] = 1;
    nodes_isovar[1][0] = 0.5;
    nodes_isovar[1][1] = 1;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface));
    EXPECT_DOUBLE_EQ(0.5, edge.get_crossing_position(iface));
    EXPECT_EQ(-1, edge.get_crossing_sign(iface));
  }

  // Equal at node 0, first lower at node 1
  {
    nodes_isovar[0][0] = 0.;
    nodes_isovar[0][1] = 0.;
    nodes_isovar[1][0] = 0;
    nodes_isovar[1][1] = 1.;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface));
    EXPECT_DOUBLE_EQ(0., edge.get_crossing_position(iface));
    EXPECT_EQ(-1, edge.get_crossing_sign(iface));
  }

  // Equal at node 0, second lower at node 1
  {
    nodes_isovar[0][0] = 0.;
    nodes_isovar[0][1] = 0.;
    nodes_isovar[1][0] = 0;
    nodes_isovar[1][1] = -1.;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface));
  }
    // Equal at node 0, first lower at node 1
  {
    nodes_isovar[0][0] = 0.;
    nodes_isovar[0][1] = 1.;
    nodes_isovar[1][0] = 0;
    nodes_isovar[1][1] = 0.;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface));
    EXPECT_DOUBLE_EQ(1., edge.get_crossing_position(iface));
    EXPECT_EQ(1, edge.get_crossing_sign(iface));
  }

  // Equal at node 1, second lower at node 0
  {
    nodes_isovar[0][0] = 0.;
    nodes_isovar[0][1] = -1.;
    nodes_isovar[1][0] = 0.;
    nodes_isovar[1][1] = 0.;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface));
  }

  // (0, +delta) at node 0, (0, -) at node 1
  {
    nodes_isovar[0][0] = 0.;
    nodes_isovar[0][1] = 1.e-16;
    nodes_isovar[1][0] = 0.;
    nodes_isovar[1][1] = -1.;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface));
    EXPECT_NEAR(0., edge.get_crossing_position(iface), 1.e-15);
    EXPECT_EQ(1, edge.get_crossing_sign(iface));
  }

  // (0, -) at node 0, (0, +delta) at node 1
  {
    nodes_isovar[0][0] = 0.;
    nodes_isovar[0][1] = -1.;
    nodes_isovar[1][0] = 0.;
    nodes_isovar[1][1] = 1.e-16;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface));
    EXPECT_NEAR(1., edge.get_crossing_position(iface), 1.e-15);
    EXPECT_EQ(-1, edge.get_crossing_sign(iface));
  }
}

TEST(CDFEM_Parent_Edge_Three_LS, Two_Nodes_No_Snapping)
{
  const bool oneLSPerPhase = true;

  const InterfaceID iface01(0,1);
  const InterfaceID iface02(0,2);
  const InterfaceID iface12(1,2);

  std::vector<std::vector<double> > nodes_isovar(2);
  nodes_isovar[0].resize(3);
  nodes_isovar[1].resize(3);

  {
    nodes_isovar[0][0] = -0.01288;
    nodes_isovar[0][1] = 0.021996;
    nodes_isovar[0][2] = 0.01;
    nodes_isovar[1][0] = 0.038335;
    nodes_isovar[1][1] = 0.037250;
    nodes_isovar[1][2] = 0.01;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface01));
    EXPECT_FALSE(edge.have_crossing(iface12));
    EXPECT_TRUE(edge.have_crossing(iface02));
    const double crossing_pos = (nodes_isovar[0][0] - nodes_isovar[0][2]) /
                              (nodes_isovar[0][0] - nodes_isovar[1][0] + nodes_isovar[1][2] - nodes_isovar[0][2]);
    EXPECT_DOUBLE_EQ(crossing_pos, edge.get_crossing_position(iface02));
    EXPECT_EQ(1, edge.get_crossing_sign(iface02));
  }

  {
    nodes_isovar[0][0] = 0.038335;
    nodes_isovar[0][1] = 0.037250;
    nodes_isovar[0][2] = 0.01;
    nodes_isovar[1][0] = -0.01288;
    nodes_isovar[1][1] = 0.021996;
    nodes_isovar[1][2] = 0.01;
    CDFEM_Parent_Edge edge(oneLSPerPhase, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface01));
    EXPECT_FALSE(edge.have_crossing(iface12));
    EXPECT_TRUE(edge.have_crossing(iface02));
    const double crossing_pos = (nodes_isovar[0][0] - nodes_isovar[0][2]) /
                              (nodes_isovar[0][0] - nodes_isovar[1][0] + nodes_isovar[1][2] - nodes_isovar[0][2]);
    EXPECT_DOUBLE_EQ(crossing_pos, edge.get_crossing_position(iface02));
    EXPECT_EQ(-1, edge.get_crossing_sign(iface02));
  }
}

TEST(CDFEM_Parent_Edge_Three_LS, Three_Nodes_No_Snapping)
{
  const bool oneLSPerPhase = true;

  const InterfaceID iface01(0,1);
  const InterfaceID iface02(0,2);
  const InterfaceID iface12(1,2);

  std::vector<std::vector<double> > nodes_isovar(3);
  nodes_isovar[0].resize(3);
  nodes_isovar[1].resize(3);
  nodes_isovar[2].resize(3);

  {
    nodes_isovar[0][0] = 0.2;
    nodes_isovar[0][1] = 0.0025391738062501383;
    nodes_isovar[0][2] = 0.01;
    nodes_isovar[1][0] = 0.19;
    nodes_isovar[1][1] = -0.00048874386459929649;
    nodes_isovar[1][2] = 0.01;
    nodes_isovar[2][0] = 0.17;
    nodes_isovar[2][1] = -0.0052539980592431739;
    nodes_isovar[2][2] = 0.01;
    CDFEM_Parent_Edge edge(oneLSPerPhase, {0.0,0.5,1.0}, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface01));
    EXPECT_FALSE(edge.have_crossing(iface12));
    EXPECT_FALSE(edge.have_crossing(iface02));
  }

  // Check that uncrossed phase is correct for 3 node problem with 2 crossings on the edge.
  {
    nodes_isovar[0][0] = 0.015;
    nodes_isovar[0][1] = 1.;
    nodes_isovar[0][2] = 0.01;
    nodes_isovar[1][0] = 0.0095;
    nodes_isovar[1][1] = 1.;
    nodes_isovar[1][2] = 0.01;
    nodes_isovar[2][0] = 0.015;
    nodes_isovar[2][1] = 1.;
    nodes_isovar[2][2] = 0.01;
    CDFEM_Parent_Edge edge(oneLSPerPhase, {0.0,0.5,1.0}, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface01));
    EXPECT_FALSE(edge.have_crossing(iface12));
    EXPECT_FALSE(edge.have_crossing(iface02));
    EXPECT_EQ(2, edge.get_uncrossed_phase());
  }

  // Test edge that goes from phase 1-0-2 in the piecewise approximation, but there is no 0-2 in the linear version,
  // so result is just 1-2.
  {
    nodes_isovar[0][0] = 0.72535;
    nodes_isovar[0][1] = -0.844886;
    nodes_isovar[0][2] = 0.10576;
    nodes_isovar[1][0] = -0.58386;
    nodes_isovar[1][1] = -0.931365;
    nodes_isovar[1][2] = 0.7754522;
    nodes_isovar[2][0] = -0.28731;
    nodes_isovar[2][1] = 0.711750;
    nodes_isovar[2][2] = -0.5794;
    CDFEM_Parent_Edge edge(oneLSPerPhase, {0.0,0.5,1.0}, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface01));
    EXPECT_TRUE(edge.have_crossing(iface12));
    EXPECT_FALSE(edge.have_crossing(iface02));
  }

  // Test edge that goes from phase 1-0-2 in the piecewise approximation, with different locations than those
  // given by a simple linear approximation.
  {
    nodes_isovar[0][0] = 0.25;
    nodes_isovar[0][1] = 1.0;
    nodes_isovar[0][2] = 0.0;
    nodes_isovar[1][0] = 0.25;
    nodes_isovar[1][1] = 1.0;
    nodes_isovar[1][2] = 0.0;
    nodes_isovar[2][0] = 0.25;
    nodes_isovar[2][1] = 0.0;
    nodes_isovar[2][2] = 1.0;
    CDFEM_Parent_Edge edge(oneLSPerPhase, {0.0,0.5,1.0}, nodes_isovar);
    EXPECT_TRUE(edge.have_any_crossings());
    EXPECT_TRUE(edge.have_crossing(iface01));
    EXPECT_FALSE(edge.have_crossing(iface12));
    EXPECT_TRUE(edge.have_crossing(iface02));
    EXPECT_DOUBLE_EQ(0.625, edge.get_crossing_position(iface02));
    EXPECT_DOUBLE_EQ(0.875, edge.get_crossing_position(iface01));
  }

}

TEST(CDFEM_Parent_Edge_Two_LS, Two_Crossings_Same_Edge)
{
  const bool oneLSPerPhase = true;

  const InterfaceID iface01(0,1);

  std::vector<std::vector<double> > nodes_isovar(3);
  nodes_isovar[0].resize(2);
  nodes_isovar[1].resize(2);
  nodes_isovar[2].resize(2);

  // Interface (0,1) has 2 crossings, one between each parent and the mid node.
  // We will treat this as an uncrossed edge regardless of snapping.
  nodes_isovar[0][0] = 0.02;
  nodes_isovar[0][1] = 0.01;
  nodes_isovar[1][0] = 0.;
  nodes_isovar[1][1] = 0.01;
  nodes_isovar[2][0] = 0.02;
  nodes_isovar[2][1] = 0.01;

  // No snapping
  {
    CDFEM_Parent_Edge edge(oneLSPerPhase, {0.0,0.5,1.0}, nodes_isovar);
    EXPECT_FALSE(edge.have_any_crossings());
    EXPECT_FALSE(edge.have_crossing(iface01));
  }
}

void expect_all_edge_segments_have_finite_or_zero_length(const CDFEM_Parent_Edge & edge, const double snapTol)
{
  std::vector<double> crossingLocations;
  for (auto && crossing : edge.get_crossings())
  {
    crossingLocations.push_back(crossing.second);
  }
  std::sort(crossingLocations.begin(), crossingLocations.end());
  double previousCrossingLocation = 0.0;
  for (auto && crossingLocation : crossingLocations)
  {
    const double intervalSize = crossingLocation- previousCrossingLocation;
    EXPECT_TRUE(intervalSize == 0.0 || intervalSize >= snapTol) << "Found infinitesmal interval " << intervalSize << " on edge " << edge;
    previousCrossingLocation = crossingLocation;
  }
  const double lastIntervalSize = 1.0 - previousCrossingLocation;
  EXPECT_TRUE(lastIntervalSize == 0.0 || lastIntervalSize >= snapTol) << "Found infinitesmal interval " << lastIntervalSize << " on edge " << edge;
}

TEST(CDFEM_Parent_Edge_Two_LS, LSPerPhaseInfinitesimalSnapTo0)
{
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{-CDFEM_Parent_Edge::MinSize(),0.0},{1.0,0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);
  EXPECT_TRUE(edge.have_any_crossings());

  expect_all_edge_segments_have_finite_or_zero_length(edge, CDFEM_Parent_Edge::MinSize());
}

TEST(CDFEM_Parent_Edge_Two_LS, LSPerPhaseInfinitesimalSnapTo1)
{
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{1.0,0.0},{-CDFEM_Parent_Edge::MinSize(),0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);
  EXPECT_TRUE(edge.have_any_crossings());

  expect_all_edge_segments_have_finite_or_zero_length(edge, CDFEM_Parent_Edge::MinSize());
}

void expect_all_fake_crossings_are_really_fake(const CDFEM_Parent_Edge & edge)
{
  EXPECT_TRUE(edge.all_fake_crossings_are_really_fake());
}

TEST(CDFEM_Parent_Edge_Three_LS, LSPerPhaseInfinitesimalSnapInMiddle)
{
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{-1.0,1.+CDFEM_Parent_Edge::MinSize(),0.0},{1.+CDFEM_Parent_Edge::MinSize(),-1.0,0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);
  EXPECT_EQ(1u,edge.get_crossings().size());

  expect_all_edge_segments_have_finite_or_zero_length(edge, CDFEM_Parent_Edge::MinSize());
  expect_all_fake_crossings_are_really_fake(edge);
}

TEST(CDFEM_Parent_Edge_Three_LS, LSPerPhaseTieInMiddle)
{
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{0.,-1.,1.},{0.,1.,-1.}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);
  EXPECT_EQ(1u,edge.get_crossings().size());

  expect_all_edge_segments_have_finite_or_zero_length(edge, CDFEM_Parent_Edge::MinSize());
  expect_all_fake_crossings_are_really_fake(edge);
}

TEST(CDFEM_Parent_Edge_Three_LS, UnderflowAtEndStillResultsInCorrectFakeCrossing)
{
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{ -0.6, -0.6, 0.6 },{ -6e-16, -5.9e-16, -7e-16 }};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  expect_all_edge_segments_have_finite_or_zero_length(edge, CDFEM_Parent_Edge::MinSize());
  expect_all_fake_crossings_are_really_fake(edge);
}

TEST(CDFEM_Parent_Edge_Three_LS, WithScaling_UnderflowAtEndStillResultsInCorrectFakeCrossing)
{
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{ -0.6e8, -0.6e8, 0.6e8 },{ -6e-8, -5.9e-8, -7e-8 }};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  expect_all_edge_segments_have_finite_or_zero_length(edge, CDFEM_Parent_Edge::MinSize());
  expect_all_fake_crossings_are_really_fake(edge);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

TEST(CDFEM_Parent_Edge_Two_LS, SnapTo0)
{
  const double snapTol = 0.1;
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{-0.01,0.0},{1.0,0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);
  EXPECT_TRUE(edge.have_any_crossings());

  edge.collapse_small_segments_while_preserving_topology(snapTol);
  expect_all_edge_segments_have_finite_or_zero_length(edge, snapTol);
}

TEST(CDFEM_Parent_Edge_Two_LS, SnapTo1)
{
  const double snapTol = 0.1;
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{1.0,0.0},{-0.01,0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.collapse_small_segments_while_preserving_topology(snapTol);
  expect_all_edge_segments_have_finite_or_zero_length(edge, snapTol);
}

TEST(CDFEM_Parent_Edge_Three_LS, SnapInMiddle)
{
  const double snapTol = 0.1;
  const bool oneLSPerPhase = true;

  const std::vector<std::vector<double> > nodesIsovar = {{-1.0,1.1,0.0},{1.1,-1.0,0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.collapse_small_segments_while_preserving_topology(snapTol);
  expect_all_edge_segments_have_finite_or_zero_length(edge, snapTol);
}

TEST(CDFEM_Parent_Edge_Three_LS, OneSnapInMiddleMakesOtherSnapUnnecessary)
{
  const double snapTol = 0.1;
  const bool oneLSPerPhase = true;
  std::vector<double> locs = {0.4,0.44,0.53};
  const double phi0At0 = -1.0;
  const double phi0At1 = (0.-(1.-locs[0])*phi0At0)/locs[0];
  const double phi1At1 = -0.1;
  const double phi1At0 = (0.-locs[1]*phi1At1)/(1.-locs[1]);
  const double phi1AtLocs2 = (1.-locs[2])*phi1At0 + locs[2]*phi1At1;
  const double phi2At1 = -1.0;
  const double phi2At0 = (phi1AtLocs2-locs[2]*phi2At1)/(1.-locs[2]);

  const std::vector<std::vector<double> > nodesIsovar = {{phi0At0, phi1At0, phi2At0, 0.0}, {phi0At1, phi1At1, phi2At1, 0.0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.collapse_small_segments_while_preserving_topology(snapTol);
  expect_all_edge_segments_have_finite_or_zero_length(edge, snapTol);

  const std::vector<std::vector<double> > oppositeNodesIsovar = {{phi0At1, phi1At1, phi2At1, 0.0}, {phi0At0, phi1At0, phi2At0, 0.0}};
  CDFEM_Parent_Edge oppositeEdge(oneLSPerPhase, oppositeNodesIsovar);

  oppositeEdge.collapse_small_segments_while_preserving_topology(snapTol);
  expect_all_edge_segments_have_finite_or_zero_length(oppositeEdge, snapTol);

  for (auto && crossing : edge.get_crossings())
    EXPECT_DOUBLE_EQ(crossing.second, 1.-oppositeEdge.get_crossing_position(crossing.first)) << "Snapping opposite edges give different results " << edge << " compared to " << oppositeEdge;
}

void expect_crossing_at_position_with_sign(const CDFEM_Parent_Edge & edge, const InterfaceID & iface, const double pos, const int sign)
{
  EXPECT_TRUE(edge.have_crossing(iface));
  EXPECT_EQ(pos, edge.get_crossing_position(iface));
  EXPECT_EQ(sign, edge.get_crossing_sign(iface));
}

void expect_crossing_at_position_with_sign(const CDFEM_Parent_Edge & edge, const int ls, const double pos, const int sign)
{
  const InterfaceID iface(ls,ls);
  EXPECT_TRUE(edge.have_crossing(iface));
  EXPECT_EQ(pos, edge.get_crossing_position(iface));
  EXPECT_EQ(sign, edge.get_crossing_sign(iface));
}

void expect_fake_crossing_at_position_with_sign(const bool oneLSPerPhase, const CDFEM_Parent_Edge & edge, const InterfaceID & iface, const double pos, const int sign)
{
  const auto result = edge.get_crossing_position_and_sign(oneLSPerPhase, iface);
  EXPECT_LE(0, std::get<0>(result)) << "Did not find expected fake crossing " << iface << " on edge " << edge;
  EXPECT_EQ(pos, std::get<0>(result));
  EXPECT_EQ(sign, std::get<1>(result));
  EXPECT_TRUE(std::get<2>(result)) << "Found real crossing when fake one was expected for " << iface << " on edge " << edge;
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAt0MovesNegCrossingTo0)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{1}, {-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0}, {});

  expect_crossing_at_position_with_sign(edge, 0, 0., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAt0MovesPosCrossingTo0)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{-1}, {1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0}, {});

  expect_crossing_at_position_with_sign(edge, 0, 0., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAt1MovesNegCrossingTo1)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{1}, {-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {}, {0});

  expect_crossing_at_position_with_sign(edge, 0, 1., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAt1MovesPosCrossingTo1)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{-1}, {1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {}, {0});

  expect_crossing_at_position_with_sign(edge, 0, 1., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAtBothEndsMovesNegCrossingToClosestAt0)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{0.9}, {-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0}, {0});

  expect_crossing_at_position_with_sign(edge, 0, 0., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAtBothEndsMovesNegCrossingToClosestAt1)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{1}, {-0.9}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0}, {0});

  expect_crossing_at_position_with_sign(edge, 0, 1., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_One_LS, nodeDomainsAtBothEndsMovesPosCrossingToClosestAt0)
{
  const bool oneLSPerPhase = false;
  const std::vector<std::vector<double> > nodesIsovar = {{-0.9}, {1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0}, {0});

  expect_crossing_at_position_with_sign(edge, 0, 0., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAt0MovesNegCrossingTo0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{1,-1}, {-1,1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1}, {});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 0., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAt0MovesPosCrossingTo0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{-1,1}, {1,-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1}, {});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 0., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAt1MovesNegCrossingTo1)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{1,-1}, {-1,1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {}, {0,1});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 1., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAt1MovesPosCrossingTo1)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{-1,1}, {1,-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {}, {0,1});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 1., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAtBothEndsMovesNegCrossingToClosestAt0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{0.9,-0.9}, {-1,1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1}, {0,1});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 0., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAtBothEndsMovesNegCrossingToClosestAt1)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{1,-1}, {-0.9,0.9}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1}, {0,1});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 1., -1);
}

TEST(CDFEM_Parent_Edge_Snapping_Two_LS_Per_Phase, nodeDomainsAtBothEndsMovesPosCrossingToClosestAt0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{-0.9,0.9}, {1,-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1}, {0,1});

  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 0., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_Three_LS_Per_Phase, nodeDomainsMoveAllCrossingsTo0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{-1,0,2}, {2,0,-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1,2}, {});

  expect_fake_crossing_at_position_with_sign(oneLSPerPhase, edge, InterfaceID(0,1), 0., 1);
  expect_fake_crossing_at_position_with_sign(oneLSPerPhase, edge, InterfaceID(1,2), 0., 1);
  expect_crossing_at_position_with_sign(edge, InterfaceID(0,2), 0., 1);
}

TEST(CDFEM_Parent_Edge_Snapping_Three_LS_Per_Phase, nodeDomainsWithMiddlePhase2MoveAllCrossingsTo0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{-1,2,0}, {2,-1,0}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,1,2}, {});

  expect_fake_crossing_at_position_with_sign(oneLSPerPhase, edge, InterfaceID(0,2), 0., 1);
  expect_fake_crossing_at_position_with_sign(oneLSPerPhase, edge, InterfaceID(1,2), 0., -1);
  expect_crossing_at_position_with_sign(edge, InterfaceID(0,1), 0., 1);
}


TEST(CDFEM_Parent_Edge_Snapping_Three_LS_Per_Phase, nodeDomainsCenterPhaseMoveAllCrossingsTo0)
{
  const bool oneLSPerPhase = true;
  const std::vector<std::vector<double> > nodesIsovar = {{-1,0,2}, {2,0,-1}};
  CDFEM_Parent_Edge edge(oneLSPerPhase, nodesIsovar);

  edge.adjust_crossing_locations_based_on_node_captured_domains(oneLSPerPhase, {0,2}, {});

  expect_fake_crossing_at_position_with_sign(oneLSPerPhase, edge, InterfaceID(0,1), 0., 1);
  expect_crossing_at_position_with_sign(edge, InterfaceID(0,2), 0., 1);
}

} // namespace krino
