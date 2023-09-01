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
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Surface.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include "Akri_DiagWriter.hpp"
#include "Akri_Phase_Support.hpp"
#include "Akri_PhaseTag.hpp"
#include "Akri_SnapInfo.hpp"

namespace krino {

static int surface_sign_at_position(const Surface & surface, const stk::math::Vector3d & pt)
{
  const double phi = surface.point_signed_distance(pt);
  return ( (phi < 0.) ? -1 : 1 ); // GOMA sign convention
}

static std::function<double(const double)> build_edge_distance_function(const Surface & surface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords)
{
  std::function<double(const double)> distanceFunction =
    [&surface, &edgeNodeCoords](const double x)
    {
      return surface.point_signed_distance((1.-x)*edgeNodeCoords[0] + x*edgeNodeCoords[1]);
    };
  return distanceFunction;
}

static double find_crossing_position(const Surface & surface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords, const double edgeTol)
{
  const double phi0 = surface.point_signed_distance(edgeNodeCoords[0]);
  const double phi1 = surface.point_signed_distance(edgeNodeCoords[1]);
  const int maxIters = 100;
  const auto result = find_root(build_edge_distance_function(surface, edgeNodeCoords), 0., 1., phi0, phi1, maxIters, edgeTol);
  STK_ThrowRequire(result.first);
  return result.second;
}

static int compute_element_sign(const Surface & surface, const std::vector<stk::math::Vector3d> & elemNodesCoords)
{
  int crossingState = 0;
  for(auto && nodeCoords : elemNodesCoords)
  {
    crossingState = crossingState | ((surface.point_signed_distance(nodeCoords) < 0.) ? 1 : 2);
    if (crossingState == 3) return 0;
  }
  STK_ThrowAssert(crossingState == 1 || crossingState == 2);
  return (crossingState == 1) ? -1 : 1;
}

static void compute_element_signs(const std::vector<const Surface *> & surfaces, const std::vector<stk::math::Vector3d> & elemNodesCoords, std::vector<int> & elementSigns)
{
  elementSigns.clear();
  elementSigns.reserve(surfaces.size());

  for(auto && surface : surfaces)
    elementSigns.push_back(compute_element_sign(*surface, elemNodesCoords));
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
  const double edgeTol)
: myMasterElem(MasterElementDeterminer::getMasterElement(mesh.bucket(element).topology())),
  mySurfaces(surfaces),
  myEdgeCrossingTol(edgeTol)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  fill_element_node_coordinates(mesh, element, coordsField, myElementNodeCoords);
  compute_element_signs(surfaces, myElementNodeCoords, myElementSigns);
}

std::vector<InterfaceID> SurfaceElementCutter::get_sorted_cutting_interfaces() const
{
  std::vector<InterfaceID> interfaces;
  for (size_t i=0; i<myElementSigns.size(); ++i)
    if (0 == myElementSigns[i])
      interfaces.push_back(InterfaceID(i,i));
  return interfaces;
}

const Surface & SurfaceElementCutter::get_surface(const InterfaceID interface) const
{
  STK_ThrowAssert(interface.is_single_ls());
  const int lsIndex = interface.first_ls();
  STK_ThrowAssert(lsIndex < (int)mySurfaces.size());
  return *mySurfaces[lsIndex];
}

bool SurfaceElementCutter::have_crossing(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const Surface & surface = get_surface(interface);
  return surface_sign_at_position(surface, parametric_to_global_coordinates(edgeNodeCoords[0])) !=
         surface_sign_at_position(surface, parametric_to_global_coordinates(edgeNodeCoords[1]));
}

double SurfaceElementCutter::interface_crossing_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const Surface & surface = get_surface(interface);
  const std::array<stk::math::Vector3d,2> globalEdge{parametric_to_global_coordinates(edgeNodeCoords[0]), parametric_to_global_coordinates(edgeNodeCoords[1])};
  return find_crossing_position(surface, globalEdge, myEdgeCrossingTol);
}

int SurfaceElementCutter::sign_at_position(const InterfaceID interface, const stk::math::Vector3d & paramCoords) const
{
  const Surface & surface = get_surface(interface);
  return surface_sign_at_position(surface, parametric_to_global_coordinates(paramCoords));
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
        const stk::math::Vector3d node0Coords(field_data<double>(coordsField, node0), dim);
        const stk::math::Vector3d node1Coords(field_data<double>(coordsField, node1), dim);
        const double phi0 = surface.point_signed_distance(node0Coords);
        const double phi1 = surface.point_signed_distance(node1Coords);
        const bool haveCrossing = (phi0 < 0.) ? (phi1 >= 0.) : (phi1 < 0.);
        if (haveCrossing)
        {
          const std::array<stk::math::Vector3d,2> edgeNodeCoords{node0Coords, node1Coords};
          const double location = find_crossing_position(surface, edgeNodeCoords, edgeCrossingTol);
          interface.fill_sorted_domains(intersectionPointSortedDomains);
          const std::vector<stk::mesh::Entity> intersectionPointNodes{node0,node1};
          if (intersectionPointFilter(intersectionPointNodes, intersectionPointSortedDomains))
            intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, std::vector<double>{1.-location, location}, intersectionPointSortedDomains);
        }
      }
    }
  }
}

static BoundingBox compute_nodal_bounding_box(const stk::mesh::BulkData & mesh)
{
  const int nDim = mesh.mesh_meta_data().spatial_dimension();
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  BoundingBox nodeBbox;
  for ( auto && bucket : mesh.buckets(stk::topology::NODE_RANK) )
  {
    double *coord = field_data<double>(coordsField, *bucket);
    for (size_t n = 0; n < bucket->size(); ++n)
      nodeBbox.accommodate( stk::math::Vector3d(coord+n*nDim, nDim) );
  }

  return nodeBbox;
}

static void prepare_to_compute_with_surface(const stk::mesh::BulkData & mesh, const std::vector<const Surface*> & surfaces)
{
  const BoundingBox nodeBbox = compute_nodal_bounding_box(mesh);
  for (auto && surface : surfaces)
  {
    Surface * nonConstSurface = const_cast<Surface*>(surface);
    nonConstSurface->prepare_to_compute(0.0, nodeBbox, 0.); // Setup including communication of facets that are within this processors narrow band
  }
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  const stk::mesh::Selector parentElementSelector = (is_cdfem_use_case(myPhaseSupport)) ?
      get_cdfem_parent_element_selector(myActivePart, myCdfemSupport, myPhaseSupport) :
      stk::mesh::Selector(myActivePart);

  myElementsToIntersect = get_owned_parent_elements(mesh, parentElementSelector);
  prepare_to_compute_with_surface(mesh, mySurfaces);
}

void AnalyticSurfaceInterfaceGeometry::prepare_to_process_elements(const stk::mesh::BulkData & mesh,
  const std::vector<stk::mesh::Entity> & elementsToIntersect,
  const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  myElementsToIntersect = elementsToIntersect;
  prepare_to_compute_with_surface(mesh, mySurfaces);
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
  NodeToCapturedDomainsMap nodesToSnappedDomains;
  prepare_to_process_elements(mesh, nodesToSnappedDomains);

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

static bool element_intersects_interval(const std::vector<const Surface*> surfaces, const std::vector<stk::math::Vector3d> & elemNodeCoords, const std::array<double,2> & loAndHi, std::vector<double> & elemNodeDistWorkspace)
{
  for (auto && surface : surfaces)
  {
    fill_point_distances(*surface, elemNodeCoords, elemNodeDistWorkspace);
    if (InterfaceGeometry::element_with_nodal_distance_intersects_interval(elemNodeDistWorkspace, loAndHi))
      return true;
  }
  return false;
}

std::vector<stk::mesh::Entity> AnalyticSurfaceInterfaceGeometry::get_elements_that_intersect_interval(const stk::mesh::BulkData & mesh, const std::array<double,2> loAndHi) const
{
  NodeToCapturedDomainsMap nodesToSnappedDomains;
  prepare_to_process_elements(mesh, nodesToSnappedDomains);

  std::vector<stk::mesh::Entity> elementsThaIntersectInterval;
  std::vector<stk::math::Vector3d> elementNodeCoords;
  std::vector<double> elementNodeDist;
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  const stk::mesh::Selector activeLocallyOwned = myActivePart & (mesh.mesh_meta_data().locally_owned_part());

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
  {
    for(const auto & elem : *bucketPtr)
    {
      fill_element_node_coordinates(mesh, elem, coordsField, elementNodeCoords);
      if (element_intersects_interval(mySurfaces, elementNodeCoords, loAndHi, elementNodeDist))
        elementsThaIntersectInterval.push_back(elem);
    }
  }

  return elementsThaIntersectInterval;
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
  prepare_to_process_elements(mesh, nodesToSnappedDomains);

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
  prepare_to_process_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
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
  cutter.reset( new SurfaceElementCutter(mesh, element, mySurfaces, myEdgeCrossingTol) );
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
