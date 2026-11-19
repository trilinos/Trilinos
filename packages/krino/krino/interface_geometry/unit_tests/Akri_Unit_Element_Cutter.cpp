// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>

#include <gtest/gtest.h>

#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Snapper.hpp>
#include <Akri_Cutting_Surface.hpp>
#include <Akri_Element_Cutter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_ParametricElement.hpp>
#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Plane_Intersections.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_topology/topology.hpp>

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

    const ParametricElement & parametricElem = parametricElemLib.get_parametric_element(topology);
    const bool oneLSPerPhase = true;

    std::vector<const CDFEM_Parent_Edge *> elementParentEdges;
    std::vector<bool> areParentEdgesAreOrientedSameAsElementEdges;

    build_simple_parent_edges(oneLSPerPhase, topology, nodeIds, nodalIsovars, parentEdges, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges);

    cutter.reset( new One_LS_Per_Phase_Cutter(parametricElem, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges, diagonalPicker) );
  }

  ParametricElementLibrary parametricElemLib;
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

// The tet of the cdfem_hidden_interface rtest, with its node ids, and its level sets 2(1-y),
// 2(y-0.2), x and -x.  The sharp edge of phases 1, 2 and 3 is the line x=0, y=0.2.
struct TetWithSharpEdgeInPlaneOfInterface : public ElementWithCutter
{
  const stk::topology topology{stk::topology::TETRAHEDRON_4};
  const std::vector<std::vector<double>> nodalIsovars{
      {2., -0.4, -1., 1.}, {-1., 2.6, 0., 0.}, {2., -0.4, 0., 0.}, {2., -0.4, 1., -1.}};
  const std::vector<stk::mesh::EntityId> nodeIds{3, 4, 1, 2};

  // The sharp-edge crossings of a face, given by the parametric coordinates of its nodes in the
  // tet, from every pair of interfaces.
  std::vector<ElementIntersection>
  sharp_edge_crossings_of_face(const std::array<stk::math::Vector3d, 3> & faceNodes)
  {
    build_parent_edges_and_cutter(topology, nodeIds, nodalIsovars);
    std::vector<InterfaceID> interfaces;
    cutter->fill_interfaces_with_cutting_surface(interfaces);
    const ElementIntersectionPointFilter keepAll = [](const std::vector<int> &) { return true; };
    std::vector<ElementIntersection> crossings;
    for (size_t i1 = 0; i1 < interfaces.size(); ++i1)
      for (size_t i2 = i1 + 1; i2 < interfaces.size(); ++i2)
        cutter->append_tetrahedron_face_interior_intersections(
            faceNodes, interfaces[i1], interfaces[i2], keepAll, crossings);
    return crossings;
  }
};

static bool has_crossing(const std::vector<ElementIntersection> & crossings,
    const stk::math::Vector3d & faceCoords,
    const std::vector<int> & sortedDomains)
{
  for (auto && crossing : crossings)
    if (crossing.sortedDomains == sortedDomains &&
        (crossing.parametricCoords - faceCoords).length() < 1.e-12)
      return true;
  return false;
}

// A face in the plane x=0 of the interface between phases 2 and 3, (0,0,1), (0,1,0) and (0,0.2,0):
// the sharp edge lies in its plane, and crosses it where it crosses its edge from (0,0,1) to
// (0,1,0), at (0,0.2,0.8).  The planes of the interfaces and the face do not meet at a point.
TEST_F(TetWithSharpEdgeInPlaneOfInterface, faceInPlaneOfInterface_crossingOnEdgeOfFace)
{
  const std::array<stk::math::Vector3d, 3> faceNodes{{stk::math::Vector3d{0., 1., 0.},
      stk::math::Vector3d{2. / 3., 0., 1. / 6.},
      stk::math::Vector3d{2. / 15., 0., 13. / 30.}}};
  const std::vector<ElementIntersection> crossings = sharp_edge_crossings_of_face(faceNodes);
  EXPECT_TRUE(has_crossing(crossings, stk::math::Vector3d{0.2, 0., 0.}, {1, 2, 3}));
}

// The face z=0 of the tet, (-1,0,0), (0,1.5,0) and (1,0,0), is not in the plane of an interface:
// the sharp edge crosses it at the one point (0,0.2,0), inside it.
TEST_F(TetWithSharpEdgeInPlaneOfInterface, faceNotInPlaneOfInterface_crossingInsideFace)
{
  const std::array<stk::math::Vector3d, 3> faceNodes{{stk::math::Vector3d{0., 0., 0.},
      stk::math::Vector3d{1., 0., 0.},
      stk::math::Vector3d{0., 0., 1.}}};
  const std::vector<ElementIntersection> crossings = sharp_edge_crossings_of_face(faceNodes);
  EXPECT_TRUE(has_crossing(crossings, stk::math::Vector3d{2. / 15., 13. / 30., 0.}, {1, 2, 3}));
}

// A face with two nodes on the sharp edge, (0,0.2,0) and (0,0.2,13/15), and its third node (1,0,0)
// off the planes of the interfaces: the sharp edge lies in the plane of the face, which lies in no
// plane of an interface.  It crosses the face only at those two nodes, with no point off the line,
// as the near-singular intersection of the three planes gave.
TEST_F(TetWithSharpEdgeInPlaneOfInterface, faceWithTwoNodesOnSharpEdge_crossingsAtThoseNodesOnly)
{
  const std::array<stk::math::Vector3d, 3> faceNodes{{stk::math::Vector3d{2. / 15., 0., 13. / 30.},
      stk::math::Vector3d{2. / 15., 13. / 15., 0.},
      stk::math::Vector3d{0., 0., 1.}}};
  const std::vector<ElementIntersection> crossings = sharp_edge_crossings_of_face(faceNodes);
  EXPECT_TRUE(has_crossing(crossings, stk::math::Vector3d{0., 0., 0.}, {1, 2, 3}));
  EXPECT_TRUE(has_crossing(crossings, stk::math::Vector3d{1., 0., 0.}, {1, 2, 3}));
  for (auto && crossing : crossings)
  {
    if (crossing.sortedDomains == std::vector<int>{1, 2, 3})
    {
      const double distanceToNode0 =
          (crossing.parametricCoords - stk::math::Vector3d{0., 0., 0.}).length();
      const double distanceToNode1 =
          (crossing.parametricCoords - stk::math::Vector3d{1., 0., 0.}).length();
      EXPECT_LT(std::min(distanceToNode0, distanceToNode1), 1.e-12) << crossing.parametricCoords;
    }
  }
}

// Faces in the plane x=0 of the interface between phases 2 and 3, each with two nodes on the sharp
// edge and its third node off it, at positions whose coordinates are not exact.  The distances of
// the nodes from that plane are only roundoff, so that each face lies in the plane of the line,
// which crosses it only at the two nodes.  Taken as a generic face, the point solved from the
// roundoff could be on the face's edge between them (in the cdfem_octagon tests, its midpoint,
// cut again in each round).
TEST_F(TetWithSharpEdgeInPlaneOfInterface,
    facesInPlaneOfInterfaceWithTwoNodesOnSharpEdge_crossingsAtThoseNodesOnly)
{
  // In the parametric coordinates of the tet: two points of the sharp edge, (0,0.2,0) and
  // (0,0.2,13/15), and a point of the plane x=0 off it, (0,0,1).
  const stk::math::Vector3d edge0{2. / 15., 0., 13. / 30.};
  const stk::math::Vector3d edge1{2. / 15., 13. / 15., 0.};
  const stk::math::Vector3d offEdge{0., 1., 0.};
  for (int i = 0; i < 40; ++i)
  {
    const double s0 = 0.1 + 0.37 * std::fmod(0.7071067811865476 * (i + 1), 1.);
    const double s1 = 0.55 + 0.4 * std::fmod(0.4142135623730951 * (i + 3), 1.);
    const double t = 0.2 + 0.7 * std::fmod(0.5772156649015329 * (i + 5), 1.);
    const std::array<stk::math::Vector3d, 3> faceNodes{{(1. - s0) * edge0 + s0 * edge1,
        (1. - s1) * edge0 + s1 * edge1,
        (1. - t) * offEdge + t * edge0}};
    const std::vector<ElementIntersection> crossings = sharp_edge_crossings_of_face(faceNodes);
    for (auto && crossing : crossings)
    {
      if (crossing.sortedDomains == std::vector<int>{1, 2, 3})
      {
        const double distanceToNode0 =
            (crossing.parametricCoords - stk::math::Vector3d{0., 0., 0.}).length();
        const double distanceToNode1 =
            (crossing.parametricCoords - stk::math::Vector3d{1., 0., 0.}).length();
        EXPECT_LT(std::min(distanceToNode0, distanceToNode1), 1.e-12)
            << "face " << i << ": " << crossing.parametricCoords;
      }
    }
  }
}

// A face in the plane x=0.5, parallel to the sharp edge x=0, y=0.2, off it: no crossing of it.
TEST_F(TetWithSharpEdgeInPlaneOfInterface, faceParallelToSharpEdge_noCrossing)
{
  // (0.5,0,0), (0.5,0.75,0) and (0.5,0,0.5), in the parametric coordinates of the tet.
  const std::array<stk::math::Vector3d, 3> faceNodes{{stk::math::Vector3d{0., 0., 0.75},
      stk::math::Vector3d{0.5, 0., 0.5},
      stk::math::Vector3d{0., 0.5, 0.5}}};
  const std::vector<ElementIntersection> crossings = sharp_edge_crossings_of_face(faceNodes);
  for (auto && crossing : crossings)
  {
    EXPECT_NE((std::vector<int>{1, 2, 3}), crossing.sortedDomains) << crossing.parametricCoords;
  }
}

// Three phases with level sets a, -a and 0 (with a = x): the planes of all their interfaces are the
// plane x=0, so there is no sharp edge, and no face, even one crossed by the plane, has a crossing
// of one.
struct TetWithCoincidentInterfacePlanes : public ElementWithCutter
{
  const stk::topology topology{stk::topology::TETRAHEDRON_4};
  const std::vector<std::vector<double>> nodalIsovars{
      {-1., 1., 0.}, {0., 0., 0.}, {0., 0., 0.}, {1., -1., 0.}};
  const std::vector<stk::mesh::EntityId> nodeIds{3, 4, 1, 2};
};

TEST_F(TetWithCoincidentInterfacePlanes, faceCrossedByThePlane_noSharpEdgeCrossing)
{
  build_parent_edges_and_cutter(topology, nodeIds, nodalIsovars);
  std::vector<InterfaceID> interfaces;
  cutter->fill_interfaces_with_cutting_surface(interfaces);
  const ElementIntersectionPointFilter keepAll = [](const std::vector<int> &) { return true; };
  // The face z=0 of the tet, (-1,0,0), (0,1.5,0) and (1,0,0), crossed by the plane x=0.
  const std::array<stk::math::Vector3d, 3> faceNodes{{stk::math::Vector3d{0., 0., 0.},
      stk::math::Vector3d{1., 0., 0.},
      stk::math::Vector3d{0., 0., 1.}}};
  std::vector<ElementIntersection> crossings;
  for (size_t i1 = 0; i1 < interfaces.size(); ++i1)
    for (size_t i2 = i1 + 1; i2 < interfaces.size(); ++i2)
      cutter->append_tetrahedron_face_interior_intersections(
          faceNodes, interfaces[i1], interfaces[i2], keepAll, crossings);
  EXPECT_TRUE(crossings.empty());
}
}
