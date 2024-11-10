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

SurfaceElementCutter::SurfaceElementCutter(const stk::mesh::BulkData & mesh,
  stk::mesh::Entity element,
  const std::vector<const Surface *> & surfaces,
  const std::vector<int8_t> & elementSigns,
  const double edgeTol)
: myMasterElem(MasterElementDeterminer::getMasterElement(mesh.bucket(element).topology())),
  mySurfaces(surfaces),
  myElementSigns(elementSigns),
  myEdgeCrossingTol(edgeTol)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
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

std::vector<int> SurfaceElementCutter::get_interface_signs_based_on_crossings(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains) const
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
  std::vector<stk::math::Vector3d> elemNodesCoords;
  elemNodesCoords.reserve(elemNodesParamCoords.size());
  for (auto && nodeParamCoords : elemNodesParamCoords)
    elemNodesCoords.push_back(parametric_to_global_coordinates(nodeParamCoords));
  return compute_uncrossed_element_sign(get_surface(interface), elemNodesCoords);
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
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const InterfaceID interface,
    const Surface & surface,
    const double edgeCrossingTol,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  const bool intersectionPointIsOwned = true;
  std::vector<int> intersectionPointSortedDomains;
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  std::set<std::array<stk::mesh::EntityId,2>> edgesAlreadyChecked;
  for (auto && elem : elementsToIntersect)
  {
    const stk::topology topology = mesh.bucket(elem).topology();
    const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(elem);
    const unsigned numEdges = topology.num_edges();

    for (unsigned iedge = 0; iedge < numEdges; ++iedge)
    {
      const unsigned * edge_node_ordinals = get_edge_node_ordinals(topology, iedge);
      const stk::mesh::Entity node0 = elem_nodes[edge_node_ordinals[0]];
      const stk::mesh::Entity node1 = elem_nodes[edge_node_ordinals[1]];
      const stk::mesh::EntityId node0Id = mesh.identifier(node0);
      const stk::mesh::EntityId node1Id = mesh.identifier(node1);
      const std::array<stk::mesh::EntityId,2> edgeNodeIds = (node0Id < node1Id) ? std::array<stk::mesh::EntityId,2>{node0Id, node1Id} : std::array<stk::mesh::EntityId,2>{node1Id, node0Id};
      auto iter = edgesAlreadyChecked.lower_bound(edgeNodeIds);
      if (iter == edgesAlreadyChecked.end() || edgeNodeIds != *iter)
      {
        edgesAlreadyChecked.insert(iter, edgeNodeIds);
        const auto [crossingSign, interfaceLoc] = surface.compute_intersection_with_segment(get_vector_field(mesh, coordsField, node0, dim), get_vector_field(mesh, coordsField, node1, dim), edgeCrossingTol);
        if (crossingSign != 0)
        {
          interface.fill_sorted_domains(intersectionPointSortedDomains);
          const std::vector<stk::mesh::Entity> intersectionPointNodes{node0,node1};
          if (intersectionPointFilter(intersectionPointNodes, intersectionPointSortedDomains))
            intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, std::vector<double>{1.-interfaceLoc, interfaceLoc}, intersectionPointSortedDomains);
        }
      }
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
  const stk::mesh::Selector parentElementSelector = (is_cdfem_use_case(myPhaseSupport)) ?
      get_cdfem_parent_element_selector(myActivePart, myCdfemSupport, myPhaseSupport) :
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

void AnalyticSurfaceInterfaceGeometry::set_element_signs(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Selector> & perSurfaceElementSelector) const
{
  myElementsToSigns = determine_element_signs(mesh, get_coordinates_field(mesh), get_mesh_parent_element_selector(), perSurfaceElementSelector, mySurfaces);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, get_mesh_parent_elements(mesh));

  const std::vector<stk::mesh::Selector> selectAllPerSurfaceElementSelector(mySurfaces.size(), mesh.mesh_meta_data().universal_part());
  set_element_signs(mesh, selectAllPerSurfaceElementSelector);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh) const
{
  const NodeToCapturedDomainsMap emptyNodesToCapturedDomains;
  prepare_to_intersect_elements(mesh, emptyNodesToCapturedDomains);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, get_mesh_parent_elements(mesh));
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
  const std::vector<stk::mesh::Entity> & elementsToIntersect,
  const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
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

static bool element_has_possibly_cut_edge(stk::topology elemTopology, const std::vector<const Surface*> surfaces, const std::vector<stk::math::Vector3d> & elemNodeCoords, std::vector<double> & elemNodeDistWorkspace)
{
  for (auto && surface : surfaces)
  {
    fill_point_distances(*surface, elemNodeCoords, elemNodeDistWorkspace);
    if (element_has_possibly_cut_edge(elemTopology, elemNodeCoords, elemNodeDistWorkspace))
      return true;
  }
  return false;
}

std::vector<stk::mesh::Entity> AnalyticSurfaceInterfaceGeometry::get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const
{
  prepare_to_intersect_elements(mesh);

  std::vector<stk::mesh::Entity> possibleCutElements;
  std::vector<stk::math::Vector3d> elementNodeCoords;
  std::vector<double> elementNodeDist;
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  const stk::mesh::Selector activeLocallyOwned = myActivePart & mesh.mesh_meta_data().locally_owned_part();

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
  {
    for(const auto & elem : *bucketPtr)
    {
      fill_element_node_coordinates(mesh, elem, coordsField, elementNodeCoords);
      if (element_has_possibly_cut_edge(bucketPtr->topology(), mySurfaces, elementNodeCoords, elementNodeDist))
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

const Surface & AnalyticSurfaceInterfaceGeometry::get_surface_with_identifer(const Surface_Identifier surfaceIdentifier) const
{
  auto iter = std::find(mySurfaceIdentifiers.begin(), mySurfaceIdentifiers.end(), surfaceIdentifier);
  STK_ThrowRequire(iter != mySurfaceIdentifiers.end());
  return *mySurfaces[std::distance(mySurfaceIdentifiers.begin(), iter)];
}

void AnalyticSurfaceInterfaceGeometry::fill_elements_that_intersect_distance_interval(const stk::mesh::BulkData & mesh, const Surface_Identifier surfaceIdentifier, const std::array<double,2> loAndHi, std::vector<stk::mesh::Entity> & elementsThaIntersectInterval) const
{
  std::vector<stk::math::Vector3d> elementNodeCoords;
  std::vector<double> elementNodeDist;
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  const stk::mesh::Selector activeLocallyOwned = myActivePart & (mesh.mesh_meta_data().locally_owned_part());

  elementsThaIntersectInterval.clear();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
  {
    for(const auto & elem : *bucketPtr)
    {
      fill_element_node_coordinates(mesh, elem, coordsField, elementNodeCoords);
      if (element_intersects_distance_interval(get_surface_with_identifer(surfaceIdentifier), elementNodeCoords, loAndHi, elementNodeDist))
        elementsThaIntersectInterval.push_back(elem);
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
      const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
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
  : myActivePart(activePart),
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
}

void AnalyticSurfaceInterfaceGeometry::add_surface(const Surface_Identifier surfaceIdentifier, const Surface & surface)
{
  mySurfaceIdentifiers.push_back(surfaceIdentifier);
  mySurfaces.push_back(&surface);
}

void AnalyticSurfaceInterfaceGeometry::store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<SnapInfo> & snapInfos,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  const bool oneLSPerPhase = mySurfaceIdentifiers.size() > 1 && myPhaseSupport.has_one_levelset_per_phase();
  if (!oneLSPerPhase && mySurfaceIdentifiers.size() > 1)
    return; //FIXME: Fix for more than one ls per interface
  STK_ThrowAssert(mySurfaces.size() == 1);
  const Surface & surface = *mySurfaces[0];

  for (auto && snapInfo : snapInfos)
  {
    stk::mesh::Entity snapNode = mesh.get_entity(stk::topology::NODE_RANK, snapInfo.get_node_global_id());
    for (auto elem : StkMeshEntities{mesh.begin_elements(snapNode), mesh.end_elements(snapNode)})
      if (mesh.bucket(elem).owned() && mesh.bucket(elem).member(myActivePart))
        set_domains_for_element_if_it_will_be_uncut_after_snapping(mesh, surface, elem, snapNode, nodesToCapturedDomains, myUncutElementPhases);
  }
}

std::vector<IntersectionPoint> AnalyticSurfaceInterfaceGeometry::get_edge_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  NodeToCapturedDomainsMap nodesToSnappedDomains;
  prepare_to_intersect_elements(mesh, nodesToSnappedDomains);

  const IntersectionPointFilter intersectionPointFilter = keep_all_intersection_points_filter();
  std::vector<IntersectionPoint> intersectionPoints;
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    InterfaceID interface(i,i);
    append_surface_edge_intersection_points(mesh, myElementsToIntersect, interface, *mySurfaces[i], myEdgeCrossingTol, intersectionPointFilter, intersectionPoints);
  }
  return intersectionPoints;
}

void AnalyticSurfaceInterfaceGeometry::append_element_intersection_points(const stk::mesh::BulkData & mesh,
  const NodeToCapturedDomainsMap & nodesToCapturedDomains,
  const std::vector<stk::mesh::Entity> & elementsToIntersect,
  const IntersectionPointFilter & intersectionPointFilter,
  std::vector<IntersectionPoint> & intersectionPoints) const
{
  prepare_to_intersect_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    InterfaceID interface(i,i);
    append_surface_edge_intersection_points(mesh, myElementsToIntersect, interface, *mySurfaces[i], myEdgeCrossingTol, intersectionPointFilter, intersectionPoints);
  }
}

std::unique_ptr<ElementCutter> AnalyticSurfaceInterfaceGeometry::build_element_cutter(const stk::mesh::BulkData & mesh,
  stk::mesh::Entity element,
  const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const
{
  std::unique_ptr<ElementCutter> cutter;
  cutter.reset( new SurfaceElementCutter(mesh, element, mySurfaces, myElementsToSigns.at(element), myEdgeCrossingTol) );
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
