// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <map>

#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_MultipleEdgeCrossings.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_SharpFeatureIntersectionPoints.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Edge.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Sign.hpp>
#include <Akri_SimplexCurvature.hpp>
#include <Akri_SimplexGradient.hpp>
#include <Akri_SnapInfo.hpp>
#include <Akri_Surface.hpp>

namespace krino {

static int surface_sign_at_position(const Surface & surface, const stk::math::Vector3d & pt)
{
  const double phi = surface.point_signed_distance(pt);
  return sign(phi);
}

static int compute_uncrossed_element_sign(const Surface & surface, const std::vector<stk::math::Vector3d> & elemNodesCoords)
{
  // Can't just test centroid due to curved penetrating surfaces
  double signedDistWithMaxMag = 0.;
  for (auto & nodeCoords : elemNodesCoords)
  {
    const double nodeSignedDist = surface.point_signed_distance(nodeCoords);
    if (std::abs(nodeSignedDist) > std::abs(signedDistWithMaxMag))
      signedDistWithMaxMag = nodeSignedDist;
  }
  const double centroidSignedDist = surface.point_signed_distance(compute_element_centroid(elemNodesCoords));
  if (std::abs(centroidSignedDist) > std::abs(signedDistWithMaxMag))
    signedDistWithMaxMag = centroidSignedDist;
  return signedDistWithMaxMag < 0. ? -1 : 1;
}

// The vertices of the element: all its nodes if it is linear.  A quadratic element is cut as the
// linear element of its vertices.
static StkMeshEntities element_vertex_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element)
{
  const stk::mesh::Entity * elemNodes = mesh.begin_nodes(element);
  return StkMeshEntities{elemNodes, elemNodes + mesh.bucket(element).topology().base().num_nodes()};
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
: myParametricElem(AuxMetaData::get(mesh.mesh_meta_data()).parametric_element_library().get_parametric_element(mesh.bucket(element).topology())),
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

// The domains of the points of rank intPtRank (2 on a face, 3 inside a tet) of surface iSurf: the
// surface's index, once for each of its sharp features that meet at the point.
static std::vector<int> intersection_domains(const int intPtRank, const int iSurf)
{
  return std::vector<int>(intPtRank, iSurf);
}

static void append_intersections(const std::vector<int> & sortedDomains,
    const std::vector<stk::math::Vector3d> & intParamCoords,
    std::vector<ElementIntersection> & intersections)
{
  for (auto & intPt : intParamCoords)
    intersections.emplace_back(intPt, sortedDomains);
}

// The domains of a point of one surface (that surface once per codimension), with each other surface
// that the point at location is on, within distanceTol, added once: so its codimension
// (sharp_feature_codimension) stays the same for two surfaces.  The distances are truncated at
// distanceTol, so a faceted surface only searches its facets that near.
static std::vector<int> domains_with_every_surface_point_is_on(const std::vector<int> & sortedDomains,
    const stk::math::Vector3d & location,
    const double distanceTol,
    const std::vector<const Surface *> & surfaces)
{
  std::vector<int> domains = sortedDomains;
  const double farFieldDistance = 2. * distanceTol;
  for (int iSurf = 0; iSurf < int(surfaces.size()); ++iSurf)
    if (iSurf != sortedDomains.front() &&
        std::abs(surfaces[iSurf]->possibly_truncated_point_signed_distance(location, distanceTol, farFieldDistance)) <= distanceTol)
      domains.push_back(iSurf);
  std::sort(domains.begin(), domains.end());
  return domains;
}

// Whether the points with these domains are on the same surfaces, with the same codimension.
static bool have_same_surfaces_and_codimension(const std::vector<int> & sortedDomains0, const std::vector<int> & sortedDomains1)
{
  if (sharp_feature_codimension(sortedDomains0) != sharp_feature_codimension(sortedDomains1)) return false;
  std::vector<int> surfaces0 = sortedDomains0;
  std::vector<int> surfaces1 = sortedDomains1;
  surfaces0.erase(std::unique(surfaces0.begin(), surfaces0.end()), surfaces0.end());
  surfaces1.erase(std::unique(surfaces1.begin(), surfaces1.end()), surfaces1.end());
  return surfaces0 == surfaces1;
}

// As put_intersection_points_on_every_surface_they_are_on does for the points of the mesh, for the
// points of one surface each in the simplex with these vertices (from firstIntersection on, with
// parametric coordinates in the simplex): each also gets every other surface it is on within
// onSurfaceTol times the size of the simplex, and of those that are then the same, only the first
// is kept.
template <size_t NVERTICES>
static void put_simplex_intersections_on_every_surface_they_are_on(const std::array<stk::math::Vector3d, NVERTICES> & simplexCoords,
    const std::vector<const Surface *> & surfaces,
    const double onSurfaceTol,
    const size_t firstIntersection,
    std::vector<ElementIntersection> & intersections)
{
  if (surfaces.size() < 2 || intersections.size() <= firstIntersection) return;
  double simplexSize = 0.;
  for (auto & x0 : simplexCoords)
    for (auto & x1 : simplexCoords)
      simplexSize = std::max(simplexSize, (x1 - x0).length());

  for (size_t i = firstIntersection; i < intersections.size(); ++i)
  {
    std::vector<int> & domains = intersections[i].sortedDomains;
    if (domains.front() != domains.back()) continue;
    const stk::math::Vector3d & pCoords = intersections[i].parametricCoords;
    stk::math::Vector3d location = simplexCoords[0];
    for (size_t n = 1; n < NVERTICES; ++n)
      location += pCoords[n - 1] * (simplexCoords[n] - simplexCoords[0]);
    domains = domains_with_every_surface_point_is_on(domains, location, onSurfaceTol * simplexSize, surfaces);
  }

  size_t numKept = firstIntersection;
  for (size_t i = firstIntersection; i < intersections.size(); ++i)
  {
    bool isRepeat = false;
    for (size_t k = firstIntersection; !isRepeat && k < numKept; ++k)
    {
      isRepeat = have_same_surfaces_and_codimension(intersections[i].sortedDomains, intersections[k].sortedDomains);
      for (size_t n = 0; isRepeat && n + 1 < NVERTICES; ++n)
        isRepeat = std::abs(intersections[i].parametricCoords[n] - intersections[k].parametricCoords[n]) <= onSurfaceTol;
    }
    if (!isRepeat)
      intersections[numKept++] = intersections[i];
  }
  intersections.erase(intersections.begin() + numKept, intersections.end());
}

static void
fill_tetrahedron_intersections(const std::vector<stk::math::Vector3d> & elementNodeCoords,
    const std::vector<const Surface *> & surfaces,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    const double onSurfaceTol,
    std::vector<ElementIntersection> & intersections)
{
  STK_ThrowAssert(elementNodeCoords.size() >= 4);
  const std::array<stk::math::Vector3d,4> tetNodeCoordinates = {{elementNodeCoords[0], elementNodeCoords[1], elementNodeCoords[2], elementNodeCoords[3]}};
  std::vector<stk::math::Vector3d> intParamCoords;
  const size_t firstIntersection = intersections.size();
  for (size_t i=0; i<surfaces.size(); ++i)
  {
    const std::vector<int> sortedDomains = intersection_domains(3, i);
    if (!intersectionPointFilter(sortedDomains)) continue;
    surfaces[i]->fill_tetrahedon_intersection_parametric_coordinates(tetNodeCoordinates, intParamCoords);
    append_intersections(sortedDomains, intParamCoords, intersections);
  }
  put_simplex_intersections_on_every_surface_they_are_on(tetNodeCoordinates, surfaces, onSurfaceTol, firstIntersection, intersections);
}

static void fill_triangle_intersections(const std::vector<stk::math::Vector3d> & elementNodeCoords,
    const std::vector<const Surface *> & surfaces,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    const double onSurfaceTol,
    std::vector<ElementIntersection> & intersections)
{
  STK_ThrowAssert(elementNodeCoords.size() >= 3);
  const std::array<stk::math::Vector3d,3> triNodeCoordinates = {{elementNodeCoords[0], elementNodeCoords[1], elementNodeCoords[2]}};
  std::vector<stk::math::Vector3d> intParamCoords;
  const size_t firstIntersection = intersections.size();
  for (size_t i=0; i<surfaces.size(); ++i)
  {
    const std::vector<int> sortedDomains = intersection_domains(2, i);
    if (!intersectionPointFilter(sortedDomains)) continue;
    surfaces[i]->fill_triangle_intersection_parametric_coordinates(triNodeCoordinates, intParamCoords);
    append_intersections(sortedDomains, intParamCoords, intersections);
  }
  put_simplex_intersections_on_every_surface_they_are_on(triNodeCoordinates, surfaces, onSurfaceTol, firstIntersection, intersections);
}

void SurfaceElementCutter::fill_interior_intersections(
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const
{
  if (myMightHaveInteriorOrFaceCrossings)
  {
    const stk::topology baseTopology = myParametricElem.get_topology().base();
    if (baseTopology == stk::topology::TETRAHEDRON_4)
      fill_tetrahedron_intersections(
          myElementNodeCoords, mySurfaces, intersectionPointFilter, myEdgeCrossingTol, intersections);
    else if (baseTopology == stk::topology::TRIANGLE_3_2D || baseTopology == stk::topology::TRIANGLE_3)
      fill_triangle_intersections(
          myElementNodeCoords, mySurfaces, intersectionPointFilter, myEdgeCrossingTol, intersections);
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

void SurfaceElementCutter::append_triangle_intersections(
    const std::array<stk::math::Vector3d, 3> & triCoords,
    const ElementIntersectionPointFilter & intersectionPointFilter,
    std::vector<ElementIntersection> & intersections) const
{
  std::vector<stk::math::Vector3d> intParamCoords;
  const size_t firstIntersection = intersections.size();
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    const std::vector<int> sortedDomains = intersection_domains(2, i);
    if (!intersectionPointFilter(sortedDomains)) continue;
    mySurfaces[i]->fill_triangle_intersection_parametric_coordinates(triCoords, intParamCoords);
    append_intersections(sortedDomains, intParamCoords, intersections);
  }
  put_simplex_intersections_on_every_surface_they_are_on(triCoords, mySurfaces, myEdgeCrossingTol, firstIntersection, intersections);
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

bool SurfaceElementCutter::interface_has_sharp_features(const InterfaceID interface) const
{
  return get_surface(interface).has_sharp_features();
}

void SurfaceElementCutter::fill_multiple_edge_crossing_positions(const InterfaceID interface,
    const std::array<stk::math::Vector3d, 2> & edgeNodeGlobalCoords,
    std::vector<double> & positions) const
{
  positions.clear();

  fill_multiple_segment_crossings(get_surface(interface),
      edgeNodeGlobalCoords[0],
      edgeNodeGlobalCoords[1],
      myEdgeCrossingTol,
      positions);
}

stk::math::Vector3d SurfaceElementCutter::parametric_to_global_coordinates(const stk::math::Vector3d & pCoords) const
{
  std::vector<double> nodalShapeFunctions(myParametricElem.num_nodes());
  myParametricElem.shape_fcn(1, pCoords.data(), nodalShapeFunctions.data());
  stk::math::Vector3d pt(stk::math::Vector3d::ZERO);
  for (unsigned n=0; n<myParametricElem.num_nodes(); ++n)
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

// Where surfaces mate, a mesh edge crosses them at the same place, and their sharp features on the
// boundary of the patch they share are at the same places in the faces and elements.  Each point of
// one surface (from firstPoint on) also gets every other surface that it is within onSurfaceTol times
// the size of its stencil of (domains_with_every_surface_point_is_on), so a node snapped there is on
// each of them.  Of the points that are then the same (the same stencil nodes, surfaces and
// codimension, with weights within onSurfaceTol), only the first is kept.  This depends only on where
// the surfaces are, not on whether their facets match.
static void put_intersection_points_on_every_surface_they_are_on(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<const Surface *> & surfaces,
    const size_t firstPoint,
    const double onSurfaceTol,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  if (surfaces.size() < 2) return;
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  for (size_t i = firstPoint; i < intersectionPoints.size(); ++i)
  {
    const IntersectionPoint & point = intersectionPoints[i];
    const std::vector<int> & domains = point.get_sorted_domains();
    if (domains.front() != domains.back()) continue;

    stk::math::Vector3d location = stk::math::Vector3d::ZERO;
    for (size_t n = 0; n < point.get_nodes().size(); ++n)
      location += point.get_weights()[n] * get_vector_field(mesh, coordsField, point.get_nodes()[n], dim);
    double stencilSize = 0.;
    for (auto node0 : point.get_nodes())
      for (auto node1 : point.get_nodes())
        stencilSize = std::max(stencilSize,
            (get_vector_field(mesh, coordsField, node0, dim) - get_vector_field(mesh, coordsField, node1, dim)).length());

    const std::vector<int> newDomains = domains_with_every_surface_point_is_on(domains, location, onSurfaceTol * stencilSize, surfaces);
    if (newDomains.size() > domains.size())
      intersectionPoints[i] = IntersectionPoint(point.is_owned(), point.get_nodes(), point.get_weights(), newDomains);
  }

  // The stencil nodes, sorted, with their weights.
  auto sorted_stencil = [](const IntersectionPoint & point)
  {
    std::vector<std::pair<stk::mesh::Entity, double>> nodesAndWeights;
    for (size_t n = 0; n < point.get_nodes().size(); ++n)
      nodesAndWeights.emplace_back(point.get_nodes()[n], point.get_weights()[n]);
    std::sort(nodesAndWeights.begin(), nodesAndWeights.end(),
        [](const auto & a, const auto & b) { return a.first < b.first; });
    return nodesAndWeights;
  };

  std::map<std::vector<stk::mesh::Entity>, std::vector<size_t>> pointsByStencilNodes;
  for (size_t i = firstPoint; i < intersectionPoints.size(); ++i)
  {
    std::vector<stk::mesh::Entity> nodes = intersectionPoints[i].get_nodes();
    std::sort(nodes.begin(), nodes.end());
    pointsByStencilNodes[nodes].push_back(i);
  }

  std::vector<bool> isRepeat(intersectionPoints.size(), false);
  bool anyRepeat = false;
  for (auto & [stencilNodes, points] : pointsByStencilNodes)
  {
    for (size_t a = 0; a < points.size(); ++a)
    {
      if (isRepeat[points[a]]) continue;
      const auto stencilA = sorted_stencil(intersectionPoints[points[a]]);
      const std::vector<int> & domainsA = intersectionPoints[points[a]].get_sorted_domains();
      for (size_t b = a + 1; b < points.size(); ++b)
      {
        if (isRepeat[points[b]] ||
            !have_same_surfaces_and_codimension(intersectionPoints[points[b]].get_sorted_domains(), domainsA))
          continue;
        const auto stencilB = sorted_stencil(intersectionPoints[points[b]]);
        bool isSamePlace = true;
        for (size_t n = 0; isSamePlace && n < stencilA.size(); ++n)
          isSamePlace = std::abs(stencilA[n].second - stencilB[n].second) <= onSurfaceTol;
        if (isSamePlace)
          isRepeat[points[b]] = anyRepeat = true;
      }
    }
  }

  if (anyRepeat)
  {
    size_t numKept = firstPoint;
    for (size_t i = firstPoint; i < intersectionPoints.size(); ++i)
      if (!isRepeat[i])
        intersectionPoints[numKept++] = intersectionPoints[i];
    intersectionPoints.erase(intersectionPoints.begin() + numKept, intersectionPoints.end());
  }
}

// The multiple crossings of the surface with the mesh edges, as cutting leaves them
// (multiple_crossing_positions_to_cut, with the edge as its own extent): without the close crossings,
// and with those near an end at that end, with the end tolerance of a surface with sharp features
// only if they are cut.
static void append_surface_multiple_edge_intersection_points(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::vector<Edge> & edgesToIntersect,
    const InterfaceID interface,
    const Surface & surface,
    const double edgeCrossingTol,
    const double edgeTol,
    const bool sharpFeaturesAreCut,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  const double endTol =
      multiple_crossing_end_tolerance(sharpFeaturesAreCut && surface.has_sharp_features(), edgeTol);
  const bool intersectionPointIsOwned = true;
  std::vector<int> intersectionPointSortedDomains;
  interface.fill_sorted_domains(intersectionPointSortedDomains);
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::mesh::Entity> intersectionPointNodes;
  std::vector<double> crossingPositions;
  for (auto & edge : edgesToIntersect)
  {
    fill_edge_nodes(edge, intersectionPointNodes);
    const stk::math::Vector3d x0 = get_vector_field(mesh, coordsField, intersectionPointNodes[0], dim);
    const stk::math::Vector3d x1 = get_vector_field(mesh, coordsField, intersectionPointNodes[1], dim);
    fill_multiple_segment_crossings(surface, x0, x1, edgeCrossingTol, crossingPositions);
    if (crossingPositions.empty()) continue;
    const double length = (x1 - x0).length();
    for (const double position :
        multiple_crossing_positions_to_cut(crossingPositions, x0, x1, edgeTol * length, endTol * length))
      intersectionPoints.emplace_back(intersectionPointIsOwned,
          intersectionPointNodes,
          std::vector<double>{1. - position, position},
          intersectionPointSortedDomains);
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
    nonConstSurface->prepare_to_compute(myTime, nodeBbox, 0.); // Setup including communication of facets that are within this processors narrow band
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
  for (auto node : element_vertex_nodes(mesh, element))
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
      for (auto node : element_vertex_nodes(mesh, element))
        elemNodesCoords.push_back(get_vector_field(mesh, coordsField, node, mesh.mesh_meta_data().spatial_dimension()));
      const int elementSign = surface_sign_at_position(surface, compute_element_centroid(elemNodesCoords));

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

std::vector<IntersectionPoint> AnalyticSurfaceInterfaceGeometry::get_all_edge_intersection_points(
    const stk::mesh::BulkData & mesh) const
{
  std::vector<IntersectionPoint> intersectionPoints =
      get_edge_intersection_points(mesh, NodeToCapturedDomainsMap{});

  const FieldRef coordsField = get_coordinates_field(mesh);
  std::vector<Edge> edgesToIntersect;
  for (size_t i = 0; i < mySurfaces.size(); ++i)
  {
    InterfaceID interface(i, i);
    if (i == 0 || mySurfaceElementSelectors[i] != mySurfaceElementSelectors[i - 1])
      edgesToIntersect =
          get_edges_of_selected_elements(mesh, mySurfaceElementSelectors[i], myElementsToIntersect);
    append_surface_multiple_edge_intersection_points(mesh,
        coordsField,
        edgesToIntersect,
        interface,
        *mySurfaces[i],
        myEdgeCrossingTol,
        myCdfemSupport.get_snapper().get_edge_tolerance(),
        myCdfemSupport.cuts_sharp_features(),
        intersectionPoints);
  }
  return intersectionPoints;
}

bool AnalyticSurfaceInterfaceGeometry::node_is_on_sharp_featured_interface(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity node,
    const InterfaceID interface,
    const double tol) const
{
  if (!myMightHaveInteriorOrFaceCrossings || interface.first_ls() != interface.second_ls())
    return false;
  const Surface & surface = *mySurfaces[interface.first_ls()];
  const stk::math::Vector3d x = get_vector_field(
      mesh, get_coordinates_field(mesh), node, mesh.mesh_meta_data().spatial_dimension());
  return std::abs(surface.point_signed_distance(x)) <= tol;
}

static void append_intersection_points_from_within_elements_and_owned_faces(ParametricElementLibrary & parametricElemLibrary,
    const stk::mesh::BulkData & mesh,
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
        append_intersection_points_from_within_element_and_owned_faces(parametricElemLibrary,
            mesh,
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
  ParametricElementLibrary & parametricElemLibrary = AuxMetaData::get(mesh.mesh_meta_data()).parametric_element_library();
  std::vector<Edge> edgesToIntersect;
  const size_t firstPoint = intersectionPoints.size();
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    InterfaceID interface(i,i);
    if (i==0 || mySurfaceElementSelectors[i] != mySurfaceElementSelectors[i-1])
      edgesToIntersect = get_edges_of_selected_elements(mesh, mySurfaceElementSelectors[i], myElementsToIntersect);
    append_surface_edge_intersection_points(mesh, coordsField, edgesToIntersect, interface, *mySurfaces[i], myEdgeCrossingTol, intersectionPointFilter, intersectionPoints);
  }

  append_intersection_points_from_within_elements_and_owned_faces(parametricElemLibrary,
      mesh,
      coordsField,
      myElementsToIntersect,
      mySurfaces,
      mySurfaceElementSelectors,
      myMightHaveInteriorOrFaceCrossings,
      myEdgeCrossingTol,
      intersectionPointFilter,
      intersectionPoints);
  put_intersection_points_on_every_surface_they_are_on(mesh, coordsField, mySurfaces, firstPoint, myEdgeCrossingTol, intersectionPoints);
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
  const StkMeshEntities elemNodes = element_vertex_nodes(mesh, element);
  std::vector<int8_t> elementSigns(surfaceElementSelectors.size(), -2);
  for (size_t i=0; i<surfaceElementSelectors.size(); ++i)
    if (is_entity_selected(mesh, surfaceElementSelectors[i], element))
      elementSigns[i] = determine_element_sign_from_node_signs(elemNodes, nodesToSigns, i);

  return elementSigns;
}

void fill_multiple_segment_crossings(const Surface & surface,
    const stk::math::Vector3d & x0,
    const stk::math::Vector3d & x1,
    const double edgeCrossingTol,
    std::vector<double> & crossingPositions)
{
  // Every element sharing the edge must get the same crossings, whatever the order of its ends.
  const bool reversed =
      std::lexicographical_compare(x1.data(), x1.data() + 3, x0.data(), x0.data() + 3);
  const stk::math::Vector3d & xBegin = reversed ? x1 : x0;
  const stk::math::Vector3d & xEnd = reversed ? x0 : x1;

  std::vector<double> candidatePositions;
  surface.fill_segment_crossings(xBegin, xEnd, edgeCrossingTol, candidatePositions);
  crossingPositions.clear();
  if (candidatePositions.empty()) return;

  // A point within the tolerance of the surface is on it.
  const double onSurfaceTol = edgeCrossingTol * (xEnd - xBegin).length();
  auto sign_at = [&](const double position)
  {
    const double dist = surface.point_signed_distance((1. - position) * xBegin + position * xEnd);
    return (std::abs(dist) <= onSurfaceTol) ? 0 : sign(dist);
  };

  // Keep only the crossings where the edge goes from one side to the other, not those of an edge
  // that runs along the surface, such as one in the plane of a face.
  candidatePositions.insert(candidatePositions.begin(), 0.);
  candidatePositions.push_back(1.);
  std::vector<int> intervalSigns;
  for (size_t i = 0; i + 1 < candidatePositions.size(); ++i)
    intervalSigns.push_back(sign_at(0.5 * (candidatePositions[i] + candidatePositions[i + 1])));
  for (size_t i = 1; i + 1 < candidatePositions.size(); ++i)
    if (intervalSigns[i - 1] * intervalSigns[i] < 0)
      crossingPositions.push_back(candidatePositions[i]);

  const size_t numCrossingsShownBySigns = (sign_at(0.) * sign_at(1.) < 0) ? 1 : 0;
  if (crossingPositions.size() <= numCrossingsShownBySigns)
  {
    crossingPositions.clear();
    return;
  }

  if (reversed)
  {
    std::reverse(crossingPositions.begin(), crossingPositions.end());
    for (auto & position : crossingPositions)
      position = 1. - position;
  }
}

static bool element_has_multiple_edge_crossings(const stk::topology topology,
    const std::vector<stk::math::Vector3d> & elemNodeCoords,
    const Surface & surface,
    const double edgeCrossingTol)
{
  std::array<unsigned, 2> edgeNodeOrdinals;
  std::vector<double> crossingPositions;
  for (unsigned iEdge = 0; iEdge < topology.num_edges(); ++iEdge)
  {
    topology.edge_node_ordinals(iEdge, edgeNodeOrdinals.data());
    fill_multiple_segment_crossings(surface,
        elemNodeCoords[edgeNodeOrdinals[0]],
        elemNodeCoords[edgeNodeOrdinals[1]],
        edgeCrossingTol,
        crossingPositions);
    if (!crossingPositions.empty()) return true;
  }
  return false;
}

static bool any_point_is_away_from_the_vertices(
    const std::vector<stk::math::Vector3d> & paramCoords, const int dim, const double tol)
{
  for (auto & pt : paramCoords)
  {
    double maxWeight = 1. - pt[0] - pt[1] - ((dim == 3) ? pt[2] : 0.);
    for (int d = 0; d < dim; ++d)
      maxWeight = std::max(maxWeight, pt[d]);
    if (maxWeight < 1. - tol) return true;
  }
  return false;
}

static bool element_has_sharp_feature_points(const stk::topology topology,
    const std::vector<stk::math::Vector3d> & elemNodeCoords,
    const Surface & surface,
    const double tol)
{
  std::vector<stk::math::Vector3d> paramCoords;
  if (topology.base() == stk::topology::TETRAHEDRON_4)
  {
    surface.fill_tetrahedon_intersection_parametric_coordinates(
        {{elemNodeCoords[0], elemNodeCoords[1], elemNodeCoords[2], elemNodeCoords[3]}},
        paramCoords);
    if (any_point_is_away_from_the_vertices(paramCoords, 3, tol)) return true;
    std::array<unsigned, 3> faceNodeOrdinals;
    for (unsigned iFace = 0; iFace < topology.num_faces(); ++iFace)
    {
      topology.face_node_ordinals(iFace, faceNodeOrdinals.data());
      surface.fill_triangle_intersection_parametric_coordinates(
          {{elemNodeCoords[faceNodeOrdinals[0]],
              elemNodeCoords[faceNodeOrdinals[1]],
              elemNodeCoords[faceNodeOrdinals[2]]}},
          paramCoords);
      if (any_point_is_away_from_the_vertices(paramCoords, 2, tol)) return true;
    }
    return false;
  }
  surface.fill_triangle_intersection_parametric_coordinates(
      {{elemNodeCoords[0], elemNodeCoords[1], elemNodeCoords[2]}}, paramCoords);
  return any_point_is_away_from_the_vertices(paramCoords, 2, tol);
}

// A surface can cross an element where the signs of its nodes do not show it, wherever an edge
// crosses it more than the signs at its ends show, with or without sharp features, or, with sharp
// features, where a sharp edge or corner passes through a face.
static bool element_has_crossings_not_shown_by_node_signs(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    stk::mesh::Entity element,
    const Surface & surface,
    const double edgeCrossingTol)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  const stk::topology topology = mesh.bucket(element).topology().base();
  std::vector<stk::math::Vector3d> elemNodeCoords;
  for (auto node : element_vertex_nodes(mesh, element))
    elemNodeCoords.push_back(get_vector_field(mesh, coordsField, node, dim));
  return element_has_multiple_edge_crossings(topology, elemNodeCoords, surface, edgeCrossingTol) ||
      element_has_sharp_feature_points(topology, elemNodeCoords, surface, edgeCrossingTol);
}

// Whether the element reaches the other side of the surface from elementSign, the side its nodes off the
// surface are on, between two of its nodes on the surface: the nodes with sign 0, which snapping captured
// on it (their node domains, not their distances, put them there).  The midpoint of the edge joining them is
// strictly on the other side, as where the edge cuts across a convex edge of the surface, though no sign
// changes at its ends.
static bool does_element_reach_other_side_between_its_nodes_on_surface(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Entity element,
    const NodeToSignsMap & nodesToSigns,
    const unsigned surfIndex,
    const Surface & surface,
    const int elementSign)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> nodesOnSurface;
  for (auto node : element_vertex_nodes(mesh, element))
    if (0 == nodesToSigns.at(node)[surfIndex])
      nodesOnSurface.push_back(get_vector_field(mesh, coordsField, node, dim));
  for (size_t i = 0; i < nodesOnSurface.size(); ++i)
    for (size_t j = i + 1; j < nodesOnSurface.size(); ++j)
      if (elementSign * surface.point_signed_distance(0.5 * (nodesOnSurface[i] + nodesOnSurface[j])) < 0.)
        return true;
  return false;
}

std::unique_ptr<ElementCutter> AnalyticSurfaceInterfaceGeometry::build_element_cutter(const stk::mesh::BulkData & mesh,
  stk::mesh::Entity element,
  const std::function<bool(const std::array<unsigned,4> &)> & /*intersectingPlanesDiagonalPicker*/) const
{
  std::unique_ptr<ElementCutter> cutter;
  const FieldRef coordsField = get_coordinates_field(mesh);
  std::vector<int8_t> elementSigns = compute_element_signs_from_node_signs(
      mesh, element, mySurfaceElementSelectors, myNodesToSigns);
  for (size_t i = 0; i < elementSigns.size(); ++i)
    if ((elementSigns[i] == -1 || elementSigns[i] == 1) &&
        (element_has_crossings_not_shown_by_node_signs(mesh, coordsField, element, *mySurfaces[i], myEdgeCrossingTol) ||
            does_element_reach_other_side_between_its_nodes_on_surface(
                mesh, coordsField, element, myNodesToSigns, i, *mySurfaces[i], elementSigns[i])))
      elementSigns[i] = 0;
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

static double tri_closest_point_normal_curvature_magnitude(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<stk::math::Vector3d> & elemNodesClosestPtNormals)
{
  const auto & [vol, normalDivergence] = tri_area_and_divergence_of_node_normals(
      elemNodesCoords[0], elemNodesCoords[1], elemNodesCoords[2],
      elemNodesClosestPtNormals[0], elemNodesClosestPtNormals[1], elemNodesClosestPtNormals[2]);
  return std::abs(normalDivergence);
}

static double tet_closest_point_normal_mean_curvature_magnitude(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<stk::math::Vector3d> & elemNodesClosestPtNormals)
{
  const auto & [vol, normalDivergence] = tet_volume_and_divergence_of_node_normals(
      elemNodesCoords[0], elemNodesCoords[1], elemNodesCoords[2], elemNodesCoords[3],
      elemNodesClosestPtNormals[0], elemNodesClosestPtNormals[1], elemNodesClosestPtNormals[2], elemNodesClosestPtNormals[3]);
  return 0.5*std::abs(normalDivergence);
}

static double finite_difference_curvature_magnitude_2d(const Surface & surf, const stk::math::Vector3d & x0, const double delta)
{
  auto f = [&surf](const stk::math::Vector3d & x) { return surf.point_signed_distance(x); };
  const stk::math::Vector3d dx(delta,0,0);
  const stk::math::Vector3d dy(0,delta,0);

  const double normalDivergence = finite_difference_divergence_of_normal_2d(delta, f(x0),
      f(x0-dx), f(x0+dx), f(x0-dy), f(x0+dy),
      f(x0-dx-dy), f(x0-dx+dy), f(x0+dx-dy), f(x0+dx+dy));
  return std::abs(normalDivergence);
}

static double finite_difference_mean_curvature_magnitude_3d(const Surface & surf, const stk::math::Vector3d & x0, const double delta)
{
  auto f = [&surf](const stk::math::Vector3d & x) { return surf.point_signed_distance(x); };
  const stk::math::Vector3d dx(delta,0,0);
  const stk::math::Vector3d dy(0,delta,0);
  const stk::math::Vector3d dz(0,0,delta);

  const double normalDivergence = finite_difference_divergence_of_normal_3d(delta, f(x0),
      f(x0-dx), f(x0+dx), f(x0-dy), f(x0+dy), f(x0-dz), f(x0+dz),
      f(x0-dx-dy), f(x0-dx+dy), f(x0+dx-dy), f(x0+dx+dy),
      f(x0-dx-dz), f(x0-dx+dz), f(x0+dx-dz), f(x0+dx+dz),
      f(x0-dy-dz), f(x0-dy+dz), f(x0+dy-dz), f(x0+dy+dz));
  return 0.5*std::abs(normalDivergence);
}

std::vector<stk::math::Vector3d> compute_closest_point_normal_for_points(const Surface & surface, const std::vector<stk::math::Vector3d> & coords)
{
  std::vector<stk::math::Vector3d> pointClosestPtNormals;
  pointClosestPtNormals.reserve(coords.size());
  for (auto & ptCoords : coords)
    pointClosestPtNormals.push_back(surface.closest_point_normal(ptCoords));
  return pointClosestPtNormals;
}

static double calculate_surface_curvature_magnitude(const unsigned dim, const Surface & surface, const std::vector<stk::math::Vector3d> & elemNodesCoords, const double elemSize)
{
  if (surface.can_approximate_closest_point_normal())
  {
    const auto elemNodesClosestPtNormals = compute_closest_point_normal_for_points(surface, elemNodesCoords);
    if (3 == dim)
      return tet_closest_point_normal_mean_curvature_magnitude(elemNodesCoords, elemNodesClosestPtNormals);
    return tri_closest_point_normal_curvature_magnitude(elemNodesCoords, elemNodesClosestPtNormals);
  }
  const stk::math::Vector3d centroid = compute_element_centroid(elemNodesCoords);
  if (3 == dim)
    return finite_difference_mean_curvature_magnitude_3d(surface, centroid, 0.5*elemSize);
  return finite_difference_curvature_magnitude_2d(surface, centroid, 0.5*elemSize);
}

double AnalyticSurfaceInterfaceGeometry::estimate_curvature_times_element_size(const stk::mesh::BulkData & mesh, stk::mesh::Entity element) const
{
  std::vector<stk::math::Vector3d> elemNodesCoords;
  fill_element_node_coordinates(mesh, element, get_coordinates_field(mesh), elemNodesCoords);
  const double elemSize = compute_simplex_RMS_edge_lengths(elemNodesCoords);

  double err = 0;
  for (auto * surf : mySurfaces)
  {
    err += elemSize * calculate_surface_curvature_magnitude(mesh.mesh_meta_data().spatial_dimension(), *surf, elemNodesCoords, elemSize);
  }
  return err;
}

static stk::math::Vector3d calculate_surface_normal(const unsigned dim, const Surface & surface, const std::vector<stk::math::Vector3d> & elemNodesCoords)
{
  if (surface.can_approximate_closest_point_normal())
  {
    const stk::math::Vector3d centroid = compute_element_centroid(elemNodesCoords);
    return surface.closest_point_normal(centroid);
  }
  std::vector<double> elemNodesDist;
  fill_point_distances(surface, elemNodesCoords, elemNodesDist);
  stk::math::Vector3d grad = calculate_simplex_gradient(elemNodesCoords, elemNodesDist);
  grad.unitize();
  return grad;
}

stk::math::Vector3d AnalyticSurfaceInterfaceGeometry::compute_interface_normal(const stk::mesh::BulkData & mesh, const Surface_Identifier surfaceIdentifier, const stk::mesh::Entity element) const
{
  std::vector<stk::math::Vector3d> elemNodesCoords;
  fill_element_node_coordinates(mesh, element, get_coordinates_field(mesh), elemNodesCoords);
  const unsigned surfIndex = get_index_of_surface_with_identifer(surfaceIdentifier);
  return calculate_surface_normal(mesh.mesh_meta_data().spatial_dimension(), *mySurfaces[surfIndex], elemNodesCoords);
}

void AnalyticSurfaceInterfaceGeometry::add_surfaces_to_journal(CubitJournal & journal,
    const BoundingBox & region,
    const std::function<std::string(const Surface_Identifier)> & surfaceGroupName) const
{
  static const std::array<std::string,6> colors{{"white", "yellow", "cyan", "pink", "khaki", "lightblue"}};
  // Each surface draws only what is in the region (not every surface implements does_intersect).
  for (size_t i=0; i<mySurfaces.size(); ++i)
  {
    journal.set_group(surfaceGroupName(mySurfaceIdentifiers[i]), colors[i % colors.size()]);
    mySurfaces[i]->add_to_journal(journal, region);
  }
}

} // namespace krino
