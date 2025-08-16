// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Edge.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Surface.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include "Akri_DiagWriter.hpp"
#include "Akri_Phase_Support.hpp"
#include "Akri_PhaseTag.hpp"
#include "Akri_SnapInfo.hpp"
#include <Akri_Sign.hpp>

namespace krino {

static int surface_sign_at_position(const Surface & surface, const stk::math::Vector3d & pt)
{
  const double phi = surface.point_signed_distance(pt);
  return sign(phi);
}

static int compute_uncrossed_element_sign(const Surface & surface, const std::vector<stk::math::Vector3d> & elemNodesCoords)
{
  // Can't just test centroid due to curved penetrating surfaces
  // Do we need to add centroid check to handle case where all nodes are on surface but centroid is not?
  double signedDistWithMaxMag = 0.;
  for (auto & nodeCoords : elemNodesCoords)
  {
    const double nodeSignedDist = surface.point_signed_distance(nodeCoords);
    if (std::abs(nodeSignedDist) > std::abs(signedDistWithMaxMag))
      signedDistWithMaxMag = nodeSignedDist;
  }
  return signedDistWithMaxMag < 0. ? -1 : 1;
}

static stk::math::Vector3d get_centroid(const std::vector<stk::math::Vector3d> & elemNodesCoords)
{
  stk::math::Vector3d centroid = stk::math::Vector3d::ZERO;
  for(auto && nodeCoords : elemNodesCoords)
  {
    centroid += nodeCoords;
  }
  centroid *= 1./elemNodesCoords.size();
  return centroid;
}

static bool does_any_surface_have_sharp_features(const std::vector<const Surface *> & surfaces)
{
  for (auto * surf : surfaces)
    if (surf->has_sharp_features())
      return true;
  return false;
}

SurfaceElementCutter::SurfaceElementCutter(const stk::mesh::BulkData & mesh,
  const FieldRef coordsField,
  stk::mesh::Entity element,
  const std::vector<const Surface *> & surfaces,
  const std::vector<int8_t> & elementSigns,
  const bool mightHaveInteriorOrFaceCrossings,
  const double edgeTol)
: myMasterElem(MasterElementDeterminer::getMasterElement(mesh.bucket(element).topology())),
  mySurfaces(surfaces),
  myElementSigns(elementSigns),
  myMightHaveInteriorOrFaceCrossings(mightHaveInteriorOrFaceCrossings),
  myEdgeCrossingTol(edgeTol)
{
  fill_element_node_coordinates(mesh, element, coordsField, myElementNodeCoords);
}

std::vector<InterfaceID> SurfaceElementCutter::get_sorted_cutting_interfaces() const
{
  std::vector<InterfaceID> interfaces;
  for (size_t i=0; i<myElementSigns.size(); ++i)
    if (0 == myElementSigns[i])
      interfaces.push_back(InterfaceID(i,i));
  return interfaces;
}

std::vector<int> SurfaceElementCutter::get_interface_signs_based_on_crossings(const std::vector<stk::math::Vector3d> & /*elemNodesCoords*/,
    const std::vector<const std::vector<int> *> & /*elemNodesSnappedDomains*/) const
{
  unsigned numInterfaces = 0;
  for (size_t i=0; i<myElementSigns.size(); ++i)
    if (0 == myElementSigns[i])
      ++numInterfaces;

  std::vector<int> interfaceSigns(numInterfaces, 0);
  return interfaceSigns;
}


const Surface & SurfaceElementCutter::get_surface(const InterfaceID interface) const
{
  STK_ThrowAssert(interface.is_single_ls());
  const int lsIndex = interface.first_ls();
  STK_ThrowAssert(lsIndex < (int)mySurfaces.size());
  return *mySurfaces[lsIndex];
}

int SurfaceElementCutter::interface_sign_for_uncrossed_element(const InterfaceID interface, const std::vector<stk::math::Vector3d> & elemNodesParamCoords) const
{
  const std::vector<stk::math::Vector3d> elemNodesCoords = parametric_to_global_coordinates(elemNodesParamCoords);
  return compute_uncrossed_element_sign(get_surface(interface), elemNodesCoords);
}

static void append_intersections(const int intPtRank, const int iSurf, const std::vector<stk::math::Vector3d> & intParamCoords, std::vector<ElementIntersection> & intersections)
{
  const std::vector<int> sortedDomains(intPtRank, iSurf);
  for (auto & intPt : intParamCoords)
    intersections.emplace_back(intPt, sortedDomains);
}

static void fill_tetrahedron_intersections(const std::vector<stk::math::Vector3d> & elementNodeCoords, const std::vector<const Surface*> & surfaces, std::vector<ElementIntersection> & intersections)
{
  STK_ThrowAssert(elementNodeCoords.size() >= 4);
  const std::array<stk::math::Vector3d,4> tetNodeCoordinates = {{elementNodeCoords[0], elementNodeCoords[1], elementNodeCoords[2], elementNodeCoords[3]}};
  std::vector<stk::math::Vector3d> intParamCoords;
  for (size_t i=0; i<surfaces.size(); ++i)
  {
    surfaces[i]->fill_tetrahedon_intersection_parametric_coordinates(tetNodeCoordinates, intParamCoords);
    append_intersections(3, i, intParamCoords, intersections);
  }
}

static void fill_triangle_intersections(const std::vector<stk::math::Vector3d> & elementNodeCoords, const std::vector<const Surface*> & surfaces, std::vector<ElementIntersection> & intersections)
{
  STK_ThrowAssert(elementNodeCoords.size() >= 3);
  const std::array<stk::math::Vector3d,3> triNodeCoordinates = {{elementNodeCoords[0], elementNodeCoords[1], elementNodeCoords[2]}};
  std::vector<stk::math::Vector3d> intParamCoords;
  for (size_t i=0; i<surfaces.size(); ++i)
  {
    surfaces[i]->fill_triangle_intersection_parametric_coordinates(triNodeCoordinates, intParamCoords);
    append_intersections(2, i, intParamCoords, intersections);
  }
}

void SurfaceElementCutter::fill_interior_intersections(const ElementIntersectionPointFilter & /*intersectionPointFilter*/, std::vector<ElementIntersection> & intersections) const
{
  if (myMightHaveInteriorOrFaceCrossings)
  {
    const stk::topology baseTopology = myMasterElem.get_topology().base();
    if (baseTopology == stk::topology::TETRAHEDRON_4)
      fill_tetrahedron_intersections(myElementNodeCoords, mySurfaces, intersections);
    else if (baseTopology == stk::topology::TRIANGLE_3_2D || baseTopology == stk::topology::TRIANGLE_3)
      fill_triangle_intersections(myElementNodeCoords, mySurfaces, intersections);
    else
      STK_ThrowErrorMsg("Unexepected topology " << baseTopology.name());
  }
}

std::vector<stk::math::Vector3d> SurfaceElementCutter::parametric_to_global_coordinates(const std::vector<stk::math::Vector3d> & nodesParamCoords) const
{
  std::vector<stk::math::Vector3d> nodesCoords;
  nodesCoords.reserve(nodesParamCoords.size());
  for (auto && nodeParamCoords : nodesParamCoords)
    nodesCoords.push_back(parametric_to_global_coordinates(nodeParamCoords));
  return nodesCoords;
}

void SurfaceElementCutter::append_triangle_intersections(const std::array<stk::math::Vector3d,3> & triCoords,
    const ElementIntersectionPointFilter & /*intersectionPointFilter*/,
    std::vector<ElementIntersection> & intersections) const
{
  std::vector<stk::math::Vector3d> intParamCoords;
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    mySurfaces[i]->fill_triangle_intersection_parametric_coordinates(triCoords, intParamCoords);
    append_intersections(2, i, intParamCoords, intersections);
  }
}

void SurfaceElementCutter::fill_tetrahedron_face_interior_intersections(const std::array<stk::math::Vector3d,3> & faceNodeParamCoords,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & faceIntersections) const
{
  if (myMightHaveInteriorOrFaceCrossings)
  {
    const std::array<stk::math::Vector3d,3> faceNodesCoords{{ parametric_to_global_coordinates(faceNodeParamCoords[0]), parametric_to_global_coordinates(faceNodeParamCoords[1]), parametric_to_global_coordinates(faceNodeParamCoords[2]) }};
    append_triangle_intersections(faceNodesCoords, intersectionPointFilter, faceIntersections);
  }
}

void SurfaceElementCutter::fill_tetrahedron_face_interior_intersections(const std::array<int,3> & faceNodeOrdinals,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & faceIntersections) const
{
  if (myMightHaveInteriorOrFaceCrossings)
  {
    const std::array<stk::math::Vector3d,3> faceNodesCoords{{ myElementNodeCoords[faceNodeOrdinals[0]], myElementNodeCoords[faceNodeOrdinals[1]], myElementNodeCoords[faceNodeOrdinals[2]] }};
    append_triangle_intersections(faceNodesCoords, intersectionPointFilter, faceIntersections);
  }
}

std::pair<int, double> SurfaceElementCutter::interface_edge_crossing_sign_and_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const Surface & surface = get_surface(interface);
  return surface.compute_intersection_with_segment(parametric_to_global_coordinates(edgeNodeCoords[0]), parametric_to_global_coordinates(edgeNodeCoords[1]), myEdgeCrossingTol);
}

stk::math::Vector3d SurfaceElementCutter::parametric_to_global_coordinates(const stk::math::Vector3d & pCoords) const
{
  std::vector<double> nodalShapeFunctions(myMasterElem.num_nodes());
  myMasterElem.shape_fcn(1, pCoords.data(), nodalShapeFunctions.data());
  stk::math::Vector3d pt(stk::math::Vector3d::ZERO);
  for (unsigned n=0; n<myMasterElem.num_nodes(); ++n)
    pt += nodalShapeFunctions[n] * myElementNodeCoords[n];
  return pt;
}

static void append_surface_edge_intersection_points(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<Edge> & edgesToIntersect,
    const InterfaceID interface,
    const Surface & surface,
    const double edgeCrossingTol,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  const bool intersectionPointIsOwned = true;
  std::vector<int> intersectionPointSortedDomains;
  interface.fill_sorted_domains(intersectionPointSortedDomains);
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::mesh::Entity> intersectionPointNodes;
  for (auto & edge : edgesToIntersect)
  {
    fill_edge_nodes(edge, intersectionPointNodes);
    if (intersectionPointFilter(intersectionPointNodes, intersectionPointSortedDomains))
    {
      const auto [crossingSign, interfaceLoc] = surface.compute_intersection_with_segment(get_vector_field(mesh, coordsField, intersectionPointNodes[0], dim), get_vector_field(mesh, coordsField, intersectionPointNodes[1], dim), edgeCrossingTol);
      if (crossingSign != 0)
        intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, std::vector<double>{1.-interfaceLoc, interfaceLoc}, intersectionPointSortedDomains);
    }
  }
}

FieldRef AnalyticSurfaceInterfaceGeometry::get_coordinates_field(const stk::mesh::BulkData & mesh) const
{
  FieldRef coordsField = myCdfemSupport.get_coords_field();
  if (!coordsField.valid())
  {
    coordsField = mesh.mesh_meta_data().coordinate_field();
    STK_ThrowRequireMsg(coordsField.valid(), "No valid coordinates field.");
  }
  return coordsField;
}

stk::mesh::Selector AnalyticSurfaceInterfaceGeometry::get_mesh_parent_element_selector() const
{
  const stk::mesh::Selector parentElementSelector = (myPhaseSupport.is_cdfem_use_case()) ?
      get_decomposed_cdfem_parent_element_selector(myActivePart, myCdfemSupport, myPhaseSupport) :
      stk::mesh::Selector(myActivePart);
  return parentElementSelector;
}

std::vector<stk::mesh::Entity> AnalyticSurfaceInterfaceGeometry::get_mesh_parent_elements(const stk::mesh::BulkData & mesh) const
{
  return get_owned_parent_elements(mesh, get_mesh_parent_element_selector());
}

void AnalyticSurfaceInterfaceGeometry::set_elements_to_intersect_and_prepare_to_compute_with_surfaces(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect) const
{
  myElementsToIntersect = elementsToIntersect;

  BoundingBox nodeBbox = compute_nodal_bbox(mesh, myActivePart, get_coordinates_field(mesh));
  for (auto && surface : mySurfaces)
  {
    Surface * nonConstSurface = const_cast<Surface*>(surface);
    nonConstSurface->prepare_to_compute(0.0, nodeBbox, 0.); // Setup including communication of facets that are within this processors narrow band
  }
}

void AnalyticSurfaceInterfaceGeometry::set_node_signs(NodeToSignsMap & nodesToSigns) const
{
  myNodesToSigns.swap(nodesToSigns);
}

void AnalyticSurfaceInterfaceGeometry::set_node_signs_from_surfaces(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Selector> & perSurfaceElementSelector,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  myNodesToSigns = determine_node_signs(mesh, get_coordinates_field(mesh), get_mesh_parent_element_selector(), perSurfaceElementSelector, mySurfaces, nodesToCapturedDomains);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, get_mesh_parent_elements(mesh));
  set_node_signs_from_surfaces(mesh, mySurfaceElementSelectors, nodesToCapturedDomains);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh) const
{
  const NodeToCapturedDomainsMap emptyNodesToCapturedDomains;
  prepare_to_intersect_elements(mesh, emptyNodesToCapturedDomains);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & /*nodesToCapturedDomains*/) const
{
  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, get_mesh_parent_elements(mesh));
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
  const std::vector<stk::mesh::Entity> & elementsToIntersect,
  const NodeToCapturedDomainsMap & /*nodesToCapturedDomains*/) const
{
  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, elementsToIntersect);
}

static bool edge_is_possibly_cut(const std::array<stk::math::Vector3d,2> & edgeNodeCoords, const std::array<double,2> & edgeNodeDist)
{
  const double edgeLength = (edgeNodeCoords[1]-edgeNodeCoords[0]).length();

  return std::abs(edgeNodeDist[0]) < edgeLength && std::abs(edgeNodeDist[1]) < edgeLength;
}

static bool element_has_possibly_cut_edge(stk::topology elemTopology, const std::vector<stk::math::Vector3d> & elemNodeCoords, const std::vector<double> & elemNodeDist)
{
  const unsigned numEdges = elemTopology.num_edges();
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, i);
    if (edge_is_possibly_cut({{elemNodeCoords[edgeNodeOrdinals[0]], elemNodeCoords[edgeNodeOrdinals[1]]}}, {{elemNodeDist[edgeNodeOrdinals[0]], elemNodeDist[edgeNodeOrdinals[1]]}}))
      return true;
  }
  return false;
}

static void fill_point_distances(const Surface & surface, const std::vector<stk::math::Vector3d> & points, std::vector<double> & pointDist)
{
  pointDist.clear();
  for (auto && point : points)
    pointDist.push_back(surface.point_signed_distance(point));
}

bool AnalyticSurfaceInterfaceGeometry::does_element_possibly_have_cut_edge(const stk::mesh::BulkData & mesh, stk::topology elemTopology, const stk::mesh::Entity elem, const std::vector<stk::math::Vector3d> & elemNodeCoords, std::vector<double> & elemNodeDistWorkspace) const
{
  for (unsigned s=0; s<mySurfaces.size(); ++s)
  {
    if (is_entity_selected(mesh, mySurfaceElementSelectors[s], elem))
    {
      fill_point_distances(*mySurfaces[s], elemNodeCoords, elemNodeDistWorkspace);
      if (element_has_possibly_cut_edge(elemTopology, elemNodeCoords, elemNodeDistWorkspace))
        return true;
    }
  }
  return false;
}

std::vector<stk::mesh::Entity> AnalyticSurfaceInterfaceGeometry::get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const
{
  prepare_to_intersect_elements(mesh);

  std::vector<stk::mesh::Entity> possibleCutElements;
  std::vector<stk::math::Vector3d> elementNodeCoords;
  std::vector<double> elementNodeDist;
  const FieldRef coordsField = get_coordinates_field(mesh);

  const stk::mesh::Selector activeLocallyOwned = myActivePart & mesh.mesh_meta_data().locally_owned_part();

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
  {
    for(const auto & elem : *bucketPtr)
    {
      fill_element_node_coordinates(mesh, elem, coordsField, elementNodeCoords);
      if (does_element_possibly_have_cut_edge(mesh, bucketPtr->topology(), elem, elementNodeCoords, elementNodeDist))
        possibleCutElements.push_back(elem);
    }
  }

  return possibleCutElements;
}

static bool element_intersects_distance_interval(const Surface & surface, const std::vector<stk::math::Vector3d> & elemNodeCoords, const std::array<double,2> & loAndHi, std::vector<double> & elemNodeDistWorkspace)
{
  fill_point_distances(surface, elemNodeCoords, elemNodeDistWorkspace);
  return InterfaceGeometry::element_with_nodal_distance_intersects_distance_interval(elemNodeDistWorkspace, loAndHi);
}

bool AnalyticSurfaceInterfaceGeometry::have_enough_surfaces_to_have_interior_intersections_or_multiple_crossings() const
{
  const unsigned minNumLSForInteriorIntersectionsOrMultipleElementCrossings = myPhaseSupport.has_one_levelset_per_phase() ? 3 : 2;
  return mySurfaces.size() >= minNumLSForInteriorIntersectionsOrMultipleElementCrossings;
}

bool AnalyticSurfaceInterfaceGeometry::snapped_elements_may_have_new_intersections() const
{
  return have_enough_surfaces_to_have_interior_intersections_or_multiple_crossings();
}

unsigned AnalyticSurfaceInterfaceGeometry::get_index_of_surface_with_identifer(const Surface_Identifier surfaceIdentifier) const
{
  auto iter = std::find(mySurfaceIdentifiers.begin(), mySurfaceIdentifiers.end(), surfaceIdentifier);
  STK_ThrowRequire(iter != mySurfaceIdentifiers.end());
  return std::distance(mySurfaceIdentifiers.begin(), iter);
}

void AnalyticSurfaceInterfaceGeometry::fill_elements_that_intersect_distance_interval(const stk::mesh::BulkData & mesh, const Surface_Identifier surfaceIdentifier, const std::array<double,2> loAndHi, std::vector<stk::mesh::Entity> & elementsThaIntersectInterval) const
{
  std::vector<stk::math::Vector3d> elementNodeCoords;
  std::vector<double> elementNodeDist;
  const FieldRef coordsField = get_coordinates_field(mesh);

  const stk::mesh::Selector activeLocallyOwned = myActivePart & (mesh.mesh_meta_data().locally_owned_part());
  const unsigned surfIndex = get_index_of_surface_with_identifer(surfaceIdentifier);

  elementsThaIntersectInterval.clear();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
  {
    if (are_entities_selected(mySurfaceElementSelectors[surfIndex], *bucketPtr))
    {
      for(const auto & elem : *bucketPtr)
      {
	fill_element_node_coordinates(mesh, elem, coordsField, elementNodeCoords);
	if (element_intersects_distance_interval(*mySurfaces[surfIndex], elementNodeCoords, loAndHi, elementNodeDist))
	  elementsThaIntersectInterval.push_back(elem);
      }
    }
  }
}

static bool all_nodes_of_element_will_be_snapped(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  for (auto node : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
    if (node != snapNode && nodesToCapturedDomains.find(node) == nodesToCapturedDomains.end())
      return false;
  return true;
}

static void set_domains_for_element_if_it_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const Surface & surface,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    ElementToDomainMap & elementsToDomain )
{
  auto iter = elementsToDomain.lower_bound(element);
  if (iter == elementsToDomain.end() || iter->first != element)
  {
    if (all_nodes_of_element_will_be_snapped(mesh, element, snapNode, nodesToCapturedDomains))
    {
      std::vector<stk::math::Vector3d> elemNodesCoords;
      fill_element_node_coordinates(mesh, element, coordsField, elemNodesCoords);
      const int elementSign = surface_sign_at_position(surface, get_centroid(elemNodesCoords));

      elementsToDomain.emplace_hint(iter, element, elementSign);
    }
  }
}

AnalyticSurfaceInterfaceGeometry::AnalyticSurfaceInterfaceGeometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
  : myMightHaveInteriorOrFaceCrossings(false),
    myActivePart(activePart),
    myCdfemSupport(cdfemSupport),
    myPhaseSupport(phaseSupport)
{
  myEdgeCrossingTol = std::min(1.e-6, 0.1*cdfemSupport.get_snapper().get_edge_tolerance());
  STK_ThrowRequireMsg(myEdgeCrossingTol > 0., "Invalid minimum edge crossing tolerance " << myEdgeCrossingTol);
}

AnalyticSurfaceInterfaceGeometry::AnalyticSurfaceInterfaceGeometry(const std::vector<Surface_Identifier> & surfaceIdentifiers,
    const std::vector<const Surface*> & surfaces,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
  : AnalyticSurfaceInterfaceGeometry(activePart, cdfemSupport, phaseSupport)
{
  mySurfaceIdentifiers = surfaceIdentifiers;
  mySurfaces = surfaces;
  for (auto & surfId : surfaceIdentifiers)
    mySurfaceElementSelectors.push_back(myPhaseSupport.get_levelset_decomposed_blocks_selector(surfId));
  myMightHaveInteriorOrFaceCrossings = does_any_surface_have_sharp_features(mySurfaces);
}

void AnalyticSurfaceInterfaceGeometry::add_surface(const Surface_Identifier surfaceIdentifier, const Surface & surface, const stk::mesh::Selector & surfaceElementSelector)
{
  mySurfaceIdentifiers.push_back(surfaceIdentifier);
  mySurfaces.push_back(&surface);
  mySurfaceElementSelectors.push_back(surfaceElementSelector);
  myMightHaveInteriorOrFaceCrossings = does_any_surface_have_sharp_features(mySurfaces);
}

void AnalyticSurfaceInterfaceGeometry::store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & /*intersectionPoints*/,
      const std::vector<SnapInfo> & snapInfos,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  const bool oneLSPerPhase = mySurfaceIdentifiers.size() > 1 && myPhaseSupport.has_one_levelset_per_phase();
  if (!oneLSPerPhase && mySurfaceIdentifiers.size() > 1)
    return; //FIXME: Fix for more than one ls per interface
  STK_ThrowAssert(mySurfaces.size() == 1);
  const Surface & surface = *mySurfaces[0];
  const FieldRef coordsField = get_coordinates_field(mesh);

  for (auto && snapInfo : snapInfos)
  {
    stk::mesh::Entity snapNode = mesh.get_entity(stk::topology::NODE_RANK, snapInfo.get_node_global_id());
    for (auto elem : StkMeshEntities{mesh.begin_elements(snapNode), mesh.end_elements(snapNode)})
      if (mesh.bucket(elem).owned() && mesh.bucket(elem).member(myActivePart))
        set_domains_for_element_if_it_will_be_uncut_after_snapping(mesh, coordsField, surface, elem, snapNode, nodesToCapturedDomains, myUncutElementPhases);
  }
}

std::vector<IntersectionPoint> AnalyticSurfaceInterfaceGeometry::get_edge_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & /*nodesToCapturedDomains*/) const
{
  NodeToCapturedDomainsMap nodesToSnappedDomains;
  prepare_to_intersect_elements(mesh, nodesToSnappedDomains);
  const FieldRef coordsField = get_coordinates_field(mesh);

  const IntersectionPointFilter intersectionPointFilter = keep_all_intersection_points_filter();
  std::vector<IntersectionPoint> intersectionPoints;
  std::vector<Edge> edgesToIntersect;
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    InterfaceID interface(i,i);
    if (i==0 || mySurfaceElementSelectors[i] != mySurfaceElementSelectors[i-1])
      edgesToIntersect = get_edges_of_selected_elements(mesh, mySurfaceElementSelectors[i], myElementsToIntersect);
    append_surface_edge_intersection_points(mesh, coordsField, edgesToIntersect, interface, *mySurfaces[i], myEdgeCrossingTol, intersectionPointFilter, intersectionPoints);
  }
  return intersectionPoints;
}

static void append_intersection_points_from_within_elements_and_owned_faces(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const std::vector<const Surface*> & surfaces,
    const std::vector<stk::mesh::Selector> & surfaceElementSelectors,
    const bool mightHaveInteriorOrFaceCrossings,
    const double edgeCrossingTol,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  const std::vector<int8_t> emptyElementSigns; // Is it ok to omit signs for intersection points?
  for (auto element : elementsToIntersect)
  {
    const SurfaceElementCutter elementCutter(mesh, coordsField, element, surfaces, emptyElementSigns, mightHaveInteriorOrFaceCrossings, edgeCrossingTol);
    for (size_t i=0; i<surfaces.size(); ++i)
    {
      if (surfaceElementSelectors[i](mesh.bucket(element)))
      {
        append_intersection_points_from_within_element_and_owned_faces(mesh,
            surfaceElementSelectors[i],
            element,
            elementCutter,
            intersectionPointFilter,
            intersectionPoints);
      }
    }
  }
}

void AnalyticSurfaceInterfaceGeometry::append_element_intersection_points(const stk::mesh::BulkData & mesh,
  const NodeToCapturedDomainsMap & nodesToCapturedDomains,
  const std::vector<stk::mesh::Entity> & elementsToIntersect,
  const IntersectionPointFilter & intersectionPointFilter,
  std::vector<IntersectionPoint> & intersectionPoints) const
{
  prepare_to_intersect_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
  const FieldRef coordsField = get_coordinates_field(mesh);
  std::vector<Edge> edgesToIntersect;
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    InterfaceID interface(i,i);
    if (i==0 || mySurfaceElementSelectors[i] != mySurfaceElementSelectors[i-1])
      edgesToIntersect = get_edges_of_selected_elements(mesh, mySurfaceElementSelectors[i], myElementsToIntersect);
    append_surface_edge_intersection_points(mesh, coordsField, edgesToIntersect, interface, *mySurfaces[i], myEdgeCrossingTol, intersectionPointFilter, intersectionPoints);
  }

  append_intersection_points_from_within_elements_and_owned_faces(mesh,
      coordsField,
      myElementsToIntersect,
      mySurfaces,
      mySurfaceElementSelectors,
      myMightHaveInteriorOrFaceCrossings,
      myEdgeCrossingTol,
      intersectionPointFilter,
      intersectionPoints);
}

static int determine_element_sign_from_node_signs(const StkMeshEntities & elemNodes,
    const NodeToSignsMap & nodesToSigns,
    const int surfIndex)
{
  bool hasNeg = false;
  bool hasPos = false;

  for (auto & node : elemNodes)
  {
    const int8_t nodeSign = nodesToSigns.at(node)[surfIndex];
    if (nodeSign < 0)
      hasNeg = true;
    else if (nodeSign > 0)
      hasPos = true;
  }

  if (hasNeg && !hasPos)
    return -1;
  if (!hasNeg && hasPos)
    return 1;
  return 0;
}

std::vector<int8_t> compute_element_signs_from_node_signs(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::vector<stk::mesh::Selector> & surfaceElementSelectors,
    const NodeToSignsMap & nodesToSigns)
{
  const StkMeshEntities elemNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  std::vector<int8_t> elementSigns(surfaceElementSelectors.size(), -2);
  for (size_t i=0; i<surfaceElementSelectors.size(); ++i)
    if (is_entity_selected(mesh, surfaceElementSelectors[i], element))
      elementSigns[i] = determine_element_sign_from_node_signs(elemNodes, nodesToSigns, i);

  return elementSigns;
}

std::unique_ptr<ElementCutter> AnalyticSurfaceInterfaceGeometry::build_element_cutter(const stk::mesh::BulkData & mesh,
  stk::mesh::Entity element,
  const std::function<bool(const std::array<unsigned,4> &)> & /*intersectingPlanesDiagonalPicker*/) const
{
  std::unique_ptr<ElementCutter> cutter;
  const FieldRef coordsField = get_coordinates_field(mesh);
  const std::vector<int8_t> & elementSigns = compute_element_signs_from_node_signs(mesh, element, mySurfaceElementSelectors, myNodesToSigns);
  cutter.reset( new SurfaceElementCutter(mesh, coordsField, element, mySurfaces, elementSigns, myMightHaveInteriorOrFaceCrossings, myEdgeCrossingTol) );
  return cutter;
}

PhaseTag AnalyticSurfaceInterfaceGeometry::get_starting_phase(const ElementCutter * cutter) const
{
  const SurfaceElementCutter * surfaceCutter = dynamic_cast<const SurfaceElementCutter *>(cutter);
  STK_ThrowRequire(surfaceCutter);

  const auto & elementSigns = surfaceCutter->get_element_signs();

  PhaseTag phase;
  for (size_t i=0; i<mySurfaceIdentifiers.size(); ++i)
    phase.add(mySurfaceIdentifiers[i], elementSigns[i]);
  return phase;
}

} // namespace krino
