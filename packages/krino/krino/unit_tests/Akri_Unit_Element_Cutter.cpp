// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Snapper.hpp>
#include <Akri_Cutting_Surface.hpp>
#include <Akri_Element_Cutter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Plane_Intersections.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_topology/topology.hpp>
#include <Akri_MasterElementDeterminer.hpp>

namespace krino
{

static std::function<bool(const std::array<unsigned,4> &)>
build_always_false_diagonal_picker()
{
  auto diagonalPicker =
  [](const std::array<unsigned,4> & /*faceNodes*/)
  {
    return false;
  };
  return diagonalPicker;
}

static void build_simple_parent_edges(const bool oneLSPerPhase,
    const stk::topology topology,
    const std::vector<stk::mesh::EntityId> & nodeIds,
    const std::vector<std::vector<double>> & nodalIsovars,
    ParentEdgeMap & parentEdges,
    std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    std::vector<bool> & areParentEdgesAreOrientedSameAsElementEdges)
{
  const unsigned numEdges = topology.num_edges();

  elementParentEdges.clear();
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * edgeLNN = get_edge_node_ordinals(topology, i);
    const unsigned i0 = edgeLNN[0];
    const unsigned i1 = edgeLNN[1];
    const ParentEdgeKey edge_key(nodeIds[i0], nodeIds[i1]);
    CDFEM_Parent_Edge & parentEdge = parentEdges[edge_key];

    if(!parentEdge.valid())
      parentEdge = CDFEM_Parent_Edge(oneLSPerPhase, {nodalIsovars[i0], nodalIsovars[i1]});

    elementParentEdges.push_back(&parentEdge);
  }

  areParentEdgesAreOrientedSameAsElementEdges.clear();
  areParentEdgesAreOrientedSameAsElementEdges.resize(numEdges, true);
}

struct ElementWithCutter : public ::testing::Test
{
  ElementWithCutter() {}

  void build_parent_edges_and_cutter(const stk::topology topology,
      const std::vector<stk::mesh::EntityId> & nodeIds,
      const std::vector<std::vector<double> > & nodalIsovars)
  {
    const auto diagonalPicker = build_always_false_diagonal_picker();

    const MasterElement & masterElem = MasterElementDeterminer::getMasterElement(topology);
    const bool oneLSPerPhase = true;

    std::vector<const CDFEM_Parent_Edge *> elementParentEdges;
    std::vector<bool> areParentEdgesAreOrientedSameAsElementEdges;

    build_simple_parent_edges(oneLSPerPhase, topology, nodeIds, nodalIsovars, parentEdges, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges);

    cutter.reset( new One_LS_Per_Phase_Cutter(masterElem, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges, diagonalPicker) );
  }

  ParentEdgeMap parentEdges;
  CDFEM_Snapper snapper;
  std::unique_ptr<Element_Cutter> cutter;
};

struct TriangleWithTriplePoint : public ElementWithCutter
{
  const stk::topology topology{stk::topology::TRIANGLE_3_2D};
  const std::vector<std::vector<double> > nodalIsovars{ {-1., 0., 0.}, {0., -1., 0.}, {1.,1.,0.} };
  const std::vector<stk::mesh::EntityId> nodeIds{1,2,3};

  const InterfaceID iface01{0,1};
  const InterfaceID iface02{0,2};
  const InterfaceID iface12{1,2};
};

struct TriangleWithFakeTriplePoint : public ElementWithCutter
{
  const stk::topology topology{stk::topology::TRIANGLE_3_2D};
  const std::vector<std::vector<double> > nodalIsovars{ {2., 2.,-1., 0.}, {-1.,0.5, 2., 0.}, {0.5, -1., 2., 0.} };
  const std::vector<stk::mesh::EntityId> nodeIds{1,2,3};
};

TEST_F(TriangleWithTriplePoint, givenCutter_haveExpectedInterfaces)
{
  build_parent_edges_and_cutter(topology, nodeIds, nodalIsovars);

  std::vector<InterfaceID> interfacesWithCuttingSurface;
  cutter->fill_interfaces_with_cutting_surface(interfacesWithCuttingSurface);
  EXPECT_EQ(3u, interfacesWithCuttingSurface.size());

  EXPECT_TRUE(cutter->have_cutting_surface(iface01));
  EXPECT_TRUE(cutter->have_cutting_surface(iface02));
  EXPECT_TRUE(cutter->have_cutting_surface(iface12));
}

static bool is_nearly_eq(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1, const double relativeTol=1.e-6)
{
  const double absoluteTol = relativeTol * (v0.length() + v1.length());
  for (int i=0; i<3; ++i)
    if (std::abs(v0[i]-v1[i]) > absoluteTol) return false;
  return true;
}

void expect_to_find_all_first_in_second(const std::vector<ElementIntersection> & vec0, const std::vector<ElementIntersection> & vec1, const std::string & errorMsg)
{
  for (auto && val0 : vec0)
  {
    bool found = false;
    for (auto && val1 : vec1)
    {
      if (is_nearly_eq(val0.parametricCoords, val1.parametricCoords) && val0.sortedDomains == val1.sortedDomains)
      {
        found = true;
        break;
      }
    }
    EXPECT_TRUE(found) << errorMsg << val0;
  }
}

void expect_num_interfaces_with_cutting_surface(size_t goldNumInterfacesWithCuttingSurfaces, const Element_Cutter & cutter)
{
  std::vector<InterfaceID> interfacesWithCuttingSurface;
  cutter.fill_interfaces_with_cutting_surface(interfacesWithCuttingSurface);
  EXPECT_EQ(goldNumInterfacesWithCuttingSurfaces, interfacesWithCuttingSurface.size());
}

void expect_intersections(const std::vector<ElementIntersection> & goldIntersections, const std::vector<ElementIntersection> & actualIntersections)
{
  EXPECT_EQ(goldIntersections.empty(), actualIntersections.empty());
  expect_to_find_all_first_in_second(actualIntersections, goldIntersections, "Actual intersection not found in gold intersections: ");
  expect_to_find_all_first_in_second(goldIntersections, actualIntersections, "Gold intersection not found in actual intersections: ");
}

TEST_F(TriangleWithTriplePoint, whenFindingIntersectionPoints_findPointAtCentroid)
{
  build_parent_edges_and_cutter(topology, nodeIds, nodalIsovars);

  std::vector<ElementIntersection> triangleIntersections;
  cutter->fill_interior_intersections(triangleIntersections);

  const std::vector<ElementIntersection> goldIntersections{ {stk::math::Vector3d{1./3.,1./3.,0.}, std::vector<int>{0,1,2}} };
  expect_intersections(goldIntersections, triangleIntersections);
}

TEST_F(TriangleWithFakeTriplePoint, whenFindingIntersectionPoints_findCorrectPoints)
{
  build_parent_edges_and_cutter(topology, nodeIds, nodalIsovars);

  std::vector<ElementIntersection> triangleIntersections;
  cutter->fill_interior_intersections(triangleIntersections);

  expect_num_interfaces_with_cutting_surface(4u, *cutter);

  const std::vector<ElementIntersection> goldIntersections{ {stk::math::Vector3d{4./9., 4./9., 0.}, std::vector<int>{0,1,3}} };
  expect_intersections(goldIntersections, triangleIntersections);
}


}
