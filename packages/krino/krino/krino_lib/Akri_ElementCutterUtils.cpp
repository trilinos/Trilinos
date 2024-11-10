// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Element_Intersections.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <vector>

namespace krino {

static void fill_nodal_shape_functions(const MasterElement & masterElement, const stk::math::Vector3d & elementParametricCoords, std::vector<double> & nodalShapeFunctions)
{
  const int npe = masterElement.num_nodes();
  nodalShapeFunctions.resize(npe);
  masterElement.shape_fcn(1, elementParametricCoords.data(), nodalShapeFunctions.data());
}

static ElementIntersectionPointFilter build_element_intersection_filter(const std::vector<stk::mesh::Entity> & nodes,
    const IntersectionPointFilter & intersectionPointFilter)
{
  auto filter =
  [&nodes, &intersectionPointFilter](const std::vector<int> & intersectionPointSortedDomains)
  {
    return intersectionPointFilter(nodes, intersectionPointSortedDomains);
  };
  return filter;
}

static void append_intersection_points_from_interior(const MasterElement & masterElement,
    const std::vector<stk::mesh::Entity> & intersectionPointNodes,
    const std::vector<ElementIntersection> & interiorIntersections,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  if (interiorIntersections.empty())
    return;

  const bool intersectionPointIsOwned = true;

  std::vector<double> intersectionPointWeights;
  for (auto & interiorIntersection : interiorIntersections)
  {
    fill_nodal_shape_functions(masterElement, interiorIntersection.parametricCoords, intersectionPointWeights);
    intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, intersectionPointWeights, interiorIntersection.sortedDomains);
  }
}

static void append_intersection_points_from_element_face(const MasterElement & elementMasterElement,
    const std::vector<stk::mesh::Entity> & elementNodes,
    const int iFace,
    const ElementCutter & elementCutter,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  std::vector<ElementIntersection> interiorIntersections;

  const std::vector<InterfaceID> interfaces = elementCutter.get_sorted_cutting_interfaces();

  std::array<int,3> faceNodeOrdinals;
  elementMasterElement.get_topology().face_node_ordinals(iFace, faceNodeOrdinals.data());
  const std::vector<stk::mesh::Entity> faceNodes{elementNodes[faceNodeOrdinals[0]], elementNodes[faceNodeOrdinals[1]], elementNodes[faceNodeOrdinals[2]]};
  const double * elemNodeParamCoords = elementMasterElement.nodal_parametric_coordinates();
  const std::array<stk::math::Vector3d,3> faceNodeCoordinates = {{stk::math::Vector3d(elemNodeParamCoords+3*faceNodeOrdinals[0]), stk::math::Vector3d(elemNodeParamCoords+3*faceNodeOrdinals[1]), stk::math::Vector3d(elemNodeParamCoords+3*faceNodeOrdinals[2])}};

  const MasterElement & faceMasterElement = MasterElementDeterminer::getMasterElement(elementMasterElement.get_topology().face_topology(iFace));
  const ElementIntersectionPointFilter faceIntersectionPointFilter = build_element_intersection_filter(faceNodes, intersectionPointFilter);

  std::vector<int> sortedDomains;
  for (size_t i1=0; i1<interfaces.size(); ++i1)
  {
    for (size_t i2=i1+1; i2<interfaces.size(); ++i2)
    {
      elementCutter.fill_tetrahedron_face_interior_intersections(faceNodeCoordinates, interfaces[i1], interfaces[i2], faceIntersectionPointFilter, interiorIntersections);
      append_intersection_points_from_interior(faceMasterElement, faceNodes, interiorIntersections, intersectionPoints);
    }
  }
}

static bool element_owns_face(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const stk::topology & elementTopology,
    const stk::mesh::Entity element,
    const std::vector<stk::mesh::Entity> & elementNodes,
    const int iFace)
{
  stk::mesh::EntityId elementId = mesh.identifier(element);
  std::vector<stk::mesh::Entity> faceNodes(elementTopology.face_topology(iFace).num_nodes());
  std::vector<stk::mesh::Entity> faceElements;
  elementTopology.face_nodes(elementNodes, iFace, faceNodes.data());
  stk::mesh::get_entities_through_relations(mesh, faceNodes, stk::topology::ELEMENT_RANK, faceElements);

  for (auto faceElement : faceElements)
    if (mesh.identifier(faceElement) < elementId && elementSelector(mesh.bucket(faceElement)))
      return false;
  return true;
}

static void append_intersection_points_from_within_element_and_owned_faces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    const stk::mesh::Entity element,
    const InterfaceGeometry & geometry,
    const std::function<bool(const std::array<unsigned,4> &)> & diagonalPicker,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  stk::topology elementTopology = mesh.bucket(element).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  const std::vector<stk::mesh::Entity> elementNodes(mesh.begin_nodes(element), mesh.end_nodes(element));

  std::unique_ptr<ElementCutter> elementCutter = geometry.build_element_cutter(mesh, element, diagonalPicker);
  if (elementCutter->might_have_interior_or_face_intersections())
  {
    append_intersection_points_from_element_interior(masterElement, elementNodes, *elementCutter, intersectionPointFilter, intersectionPoints);

    const int numFaces = masterElement.get_topology().num_faces();
    if (numFaces > 0)
    {
      for (int iFace=0; iFace<numFaces; ++iFace)
      {
        if (element_owns_face(mesh, parentElementSelector, elementTopology, element, elementNodes, iFace))
        {
          append_intersection_points_from_element_face(masterElement, elementNodes, iFace, *elementCutter, intersectionPointFilter, intersectionPoints);
        }
      }
    }
  }
}

static std::function<bool(const std::array<unsigned,4> &)>
temporary_build_always_true_diagonal_picker()
{
  auto diagonalPicker =
  [](const std::array<unsigned,4> & faceNodes)
  {
    return true;
  };
  return diagonalPicker;
}

void append_intersection_points_from_within_elements_and_owned_faces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    const std::vector<stk::mesh::Entity> & elements,
    const InterfaceGeometry & geometry,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  const auto diagonalPicker = temporary_build_always_true_diagonal_picker();

  for (auto element : elements)
  {
    if (parentElementSelector(mesh.bucket(element)))
    {
      append_intersection_points_from_within_element_and_owned_faces(mesh,
          parentElementSelector,
          element,
          geometry,
          diagonalPicker,
          intersectionPointFilter,
          intersectionPoints);
    }
  }
}

void append_intersection_points_from_element_interior(const MasterElement & masterElement,
    const std::vector<stk::mesh::Entity> & nodes,
    const ElementCutter & elementCutter,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  std::vector<ElementIntersection> interiorIntersections;
  const ElementIntersectionPointFilter elementIntersectionPointFilter = build_element_intersection_filter(nodes, intersectionPointFilter);
  elementCutter.fill_interior_intersections(elementIntersectionPointFilter, interiorIntersections);
  append_intersection_points_from_interior(masterElement, nodes, interiorIntersections, intersectionPoints);
}

}
