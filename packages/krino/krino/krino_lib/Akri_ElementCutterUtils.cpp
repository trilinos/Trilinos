// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <krino/geometry/Akri_Triangle.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/MetaData.hpp>
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
  std::array<int,3> faceNodeOrdinals;
  elementMasterElement.get_topology().face_node_ordinals(iFace, faceNodeOrdinals.data());

  const std::vector<stk::mesh::Entity> faceNodes{elementNodes[faceNodeOrdinals[0]], elementNodes[faceNodeOrdinals[1]], elementNodes[faceNodeOrdinals[2]]};
  const ElementIntersectionPointFilter faceIntersectionPointFilter = build_element_intersection_filter(faceNodes, intersectionPointFilter);

  std::vector<ElementIntersection> faceIntersections;
  elementCutter.fill_tetrahedron_face_interior_intersections(faceNodeOrdinals, faceIntersectionPointFilter, faceIntersections);

  const MasterElement & faceMasterElement = MasterElementDeterminer::getMasterElement(elementMasterElement.get_topology().face_topology(iFace));
  append_intersection_points_from_interior(faceMasterElement, faceNodes, faceIntersections, intersectionPoints);
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

void append_intersection_points_from_within_element_and_owned_faces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    const stk::mesh::Entity element,
    const ElementCutter & elementCutter,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  stk::topology elementTopology = mesh.bucket(element).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  const std::vector<stk::mesh::Entity> elementNodes(mesh.begin_nodes(element), mesh.end_nodes(element));

  if (elementCutter.might_have_interior_or_face_intersections())
  {
    append_intersection_points_from_element_interior(masterElement, elementNodes, elementCutter, intersectionPointFilter, intersectionPoints);

    const int numFaces = masterElement.get_topology().num_faces();
    if (numFaces > 0)
    {
      for (int iFace=0; iFace<numFaces; ++iFace)
      {
        if (element_owns_face(mesh, parentElementSelector, elementTopology, element, elementNodes, iFace))
        {
          append_intersection_points_from_element_face(masterElement, elementNodes, iFace, elementCutter, intersectionPointFilter, intersectionPoints);
        }
      }
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

std::map<stk::mesh::Entity, std::vector<unsigned>> map_intersections_to_elements(
  const std::vector<IntersectionPoint> & intersectionPoints,
  const stk::mesh::BulkData & mesh,
  const stk::mesh::Selector & parentElementSelector)
{
  std::map<stk::mesh::Entity, std::vector<unsigned>> elementToIntersectionMap;

  for(unsigned i=0; i<intersectionPoints.size(); i++)
  {
    const auto & intersectionPt = intersectionPoints[i];
    if(intersectionPt.get_nodes().size() != 2) continue;
    stk::mesh::Entity node0 = intersectionPt.get_nodes()[0];
    stk::mesh::Entity node1 = intersectionPt.get_nodes()[1];

    for (auto elem : StkMeshEntities{mesh.begin_elements(node0), mesh.end_elements(node0)})
    {
      if(!parentElementSelector(mesh.bucket(elem))) continue;
      for (auto node : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
      {
        if(node == node1) 
        {
          auto [iter, didInsert] = elementToIntersectionMap.insert(std::make_pair(elem, std::vector<unsigned>()));
          iter->second.push_back(i);
          break;
        }
      }
    }
  }
  return elementToIntersectionMap;
}

std::vector<stk::math::Vector3d> get_element_intersection_points(const std::vector<IntersectionPoint> & intersectionPoints,
  const stk::mesh::BulkData & /*mesh*/, FieldRef coordsField, stk::mesh::Entity /*elem*/, const std::vector<unsigned> & intersectionIds, 
  int dim)
{
  std::vector<stk::math::Vector3d> elemIntersections;
  for(auto && interId : intersectionIds)
  {
    const auto & intersection = intersectionPoints[interId];
    elemIntersections.push_back(
      stk::math::Vector3d(field_data<double>(coordsField, intersection.get_nodes()[0]), dim)*intersection.get_weights()[0] +
      stk::math::Vector3d(field_data<double>(coordsField, intersection.get_nodes()[1]), dim)*intersection.get_weights()[1]);
  }
  return elemIntersections;
}

double get_element_length_scale_squared(const std::vector<IntersectionPoint> & intersectionPoints,
  const stk::mesh::BulkData & /*mesh*/, FieldRef coordsField, stk::mesh::Entity /*elem*/, const std::vector<unsigned> & intersectionIds, 
  int dim)
{
  double minLenSquared = std::numeric_limits<double>::max();
  for(auto && interId : intersectionIds)
  {
    const auto & intersection = intersectionPoints[interId];
    const double edgeLenSq = (stk::math::Vector3d(field_data<double>(coordsField, intersection.get_nodes()[0]), dim) -
      stk::math::Vector3d(field_data<double>(coordsField, intersection.get_nodes()[1]), dim)).length_squared();
    minLenSquared = std::min(edgeLenSq, minLenSquared);
  }
  return minLenSquared;
}

void append_closest_point_intersections(std::vector<IntersectionPoint> & intersectionPoints,
  const stk::mesh::BulkData & mesh,
  const stk::mesh::Selector & parentElementSelector,
  const FieldRef coordsField,
  const IntersectionPointFilter & intersectionPointFilter)
{
  const bool intersectionPointIsOwned = true;
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  constexpr double tol = 1e-14;

  auto elemToIntersectionMap = map_intersections_to_elements(intersectionPoints, mesh, parentElementSelector);

  for(auto && elemWithIntersection : elemToIntersectionMap)
  {
    if(elemWithIntersection.second.size() < 3) continue;
    const auto elem = elemWithIntersection.first;
    const auto intersectIds = elemWithIntersection.second;

    auto elemIntersectionPoints = get_element_intersection_points(intersectionPoints,
      mesh, coordsField, elem, elemWithIntersection.second, dim);
    const double lenScaleSq = get_element_length_scale_squared(intersectionPoints,
      mesh, coordsField, elem, elemWithIntersection.second, dim);
    const double areaSqZero = 0.25 * tol * lenScaleSq * lenScaleSq;
    STK_ThrowRequire(elemIntersectionPoints.size() == 3 || elemIntersectionPoints.size() == 4);

    const auto intersectionPointSortedDomains = intersectionPoints[intersectIds[0]].get_sorted_domains();

    for(unsigned i=1; i<intersectIds.size(); i++)
    {
      STK_ThrowRequire(intersectionPointSortedDomains == intersectionPoints[intersectIds[i]].get_sorted_domains());
    }

    std::vector<stk::mesh::Entity> elemNodes;
    std::vector<stk::math::Vector3d> nodeCoords;
    for (auto node : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
    {
      elemNodes.push_back(node);
      nodeCoords.emplace_back(field_data<double>(coordsField, node), dim);
    }

    const std::vector<std::array<stk::math::Vector3d, 3>> tris = make_tris_from_intersections(elemIntersectionPoints);

    for(auto && tri : tris)
    {
      const double areaSq = krino::CalcTriangle3<double>::area_vector(tri[0], tri[1], tri[2]).length_squared();
      if(areaSq < areaSqZero) continue;
      for(unsigned n=0; n<elemNodes.size(); n++)
      {
        const auto closestPoint = krino::CalcTriangle3<double>::closest_point(tri[0], 
          tri[1], tri[2], nodeCoords[n]);
        auto parCoords = get_parametric_coordinates_of_point(nodeCoords, closestPoint);
        std::ostringstream oss;
        oss << parCoords;
        STK_ThrowRequireMsg(parCoords[0] > -tol && parCoords[1] > -tol && parCoords[2] > -tol && 
          parCoords[0] + parCoords[1] + parCoords[2] < 1. + tol, oss.str());
        unsigned nzero = 0;
        double pCoordSum = 0.;
        for(unsigned p=0; p<parCoords.size(); p++)
        {
          if(std::fabs(parCoords[p]) < tol)
          {
            parCoords[p] = 0.;
            nzero++;
          }
          pCoordSum += parCoords[p];
        }
        if(1. - pCoordSum < tol) nzero++;
        if(nzero >= 2) continue;
        if (intersectionPointFilter(elemNodes, intersectionPointSortedDomains))
        {
          intersectionPoints.emplace_back(intersectionPointIsOwned, elemNodes, 
            std::vector<double>{1. - parCoords[0] - parCoords[1] - parCoords[2], 
              parCoords[0], parCoords[1], parCoords[2]}, intersectionPointSortedDomains);
        }
      }
    }
  }
}

}
