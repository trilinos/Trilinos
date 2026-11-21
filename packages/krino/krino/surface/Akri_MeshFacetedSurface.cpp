// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetMeshChecks.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_IntersectionUtils.hpp>
#include <Akri_Sign.hpp>
#include <Akri_SurfaceIntersectionFromSignedDistance.hpp>
#include <Akri_Triangle.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <string>

namespace krino
{

static std::vector<std::pair<FacetVertexId, stk::math::Vector3d>>
vertex_ids_and_coords_from_indices(const std::vector<stk::math::Vector3d> & vertices)
{
  std::vector<std::pair<FacetVertexId, stk::math::Vector3d>> vertexIdsAndCoords;
  vertexIdsAndCoords.reserve(vertices.size());
  for (size_t i = 0; i < vertices.size(); ++i)
    vertexIdsAndCoords.emplace_back(i + 1, vertices[i]);
  return vertexIdsAndCoords;
}

static std::vector<std::array<FacetVertexId, 3>>
facet_vertex_ids_from_indices(const std::vector<std::array<unsigned, 3>> & facetVertexIndices)
{
  std::vector<std::array<FacetVertexId, 3>> facetVertexIds;
  facetVertexIds.reserve(facetVertexIndices.size());
  for (auto & facet : facetVertexIndices)
    facetVertexIds.push_back({{facet[0] + 1, facet[1] + 1, facet[2] + 1}});
  return facetVertexIds;
}

MeshFacetedSurface::MeshFacetedSurface(const stk::ParallelMachine comm,
    const unsigned dim,
    const std::vector<stk::math::Vector3d> & vertices,
    const std::vector<std::array<unsigned, 3>> & facetVertexIndices,
    const bool haveConnectedSurface)
    : MeshFacetedSurface(comm,
          dim,
          vertex_ids_and_coords_from_indices(vertices),
          facet_vertex_ids_from_indices(facetVertexIndices),
          haveConnectedSurface)
{
}

MeshFacetedSurface::MeshFacetedSurface(const stk::ParallelMachine comm,
    const unsigned dim,
    const std::vector<std::pair<FacetVertexId, stk::math::Vector3d>> & vertexIdsAndCoords,
    const std::vector<std::array<FacetVertexId, 3>> & facetVertexIds,
    const bool haveConnectedSurface)
    : MeshFacetedSurface(comm, dim, haveConnectedSurface)
{
  myGivenVertexIdsAndCoords = vertexIdsAndCoords;
  myGivenFacetVertexIds = facetVertexIds;
}

MeshFacetedSurface::MeshFacetedSurface(const stk::ParallelMachine comm,
    const unsigned dim,
    const std::vector<std::pair<FacetVertexId, stk::math::Vector3d>> & vertexIdsAndCoords,
    const std::vector<std::array<FacetVertexId, 3>> & facetVertexIds,
    const std::vector<stk::math::Vector3d> & vertexVelocities,
    const bool haveConnectedSurface)
    : MeshFacetedSurface(
          comm, dim, haveConnectedSurface, FacetMeshOptions{.withVertexVelocities = true})
{
  myGivenVertexIdsAndCoords = vertexIdsAndCoords;
  myGivenFacetVertexIds = facetVertexIds;
  myGivenVertexVelocities = vertexVelocities;
}

MeshFacetedSurface::MeshFacetedSurface(const stk::ParallelMachine comm,
    const unsigned dim,
    const std::vector<stk::math::Vector3d> & facetVertexCoords)
    : MeshFacetedSurface(comm, dim, false, FacetMeshOptions{.sharedVertices = false})
{
  myMesh.set_owned_facets_from_coordinates(facetVertexCoords);
  given_facets_were_set();
}

MeshFacetedSurface::MeshFacetedSurface(const stk::ParallelMachine comm,
    const unsigned dim,
    const bool haveConnectedSurface,
    const FacetMeshOptions & meshOptions,
    const bool rebuildOwnedFacetsEachPrepare)
    : myMesh(comm, dim, meshOptions),
      myHaveConnectedSurface(haveConnectedSurface),
      myRebuildOwnedFacetsEachPrepare(rebuildOwnedFacetsEachPrepare)
{
}

void MeshFacetedSurface::set_facets(const FacetCoordinates & facets)
{
  const unsigned dim = myMesh.spatial_dimension();
  STK_ThrowRequireMsg(facets.spatial_dimension() == dim,
      "MeshFacetedSurface of dimension " << dim << " given facets of dimension "
                                         << facets.spatial_dimension() << ".");
  if (myMesh.has_shared_vertices() || myMesh.has_vertex_velocities())
    myMesh = FacetMesh(myMesh.comm(), dim, FacetMeshOptions{.sharedVertices = false});
  myHaveConnectedSurface = false;
  myMesh.set_owned_facets_from_coordinates(facets.coordinates());
  given_facets_were_set();
}

void MeshFacetedSurface::set_facets(const FacetMesh::VertexIdsAndCoords & vertexIdsAndCoords,
    const FacetMesh::FacetVertexIds & facetVertexIds,
    const bool haveConnectedSurface)
{
  if (!myMesh.has_shared_vertices() || myMesh.has_vertex_velocities())
    myMesh = FacetMesh(myMesh.comm(), myMesh.spatial_dimension());
  myHaveConnectedSurface = haveConnectedSurface;
  myMesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds);
  given_facets_were_set();
}

void MeshFacetedSurface::given_facets_were_set()
{
  myGivenVertexIdsAndCoords.clear();
  myGivenFacetVertexIds.clear();
  myGivenVertexVelocities.clear();
  myHaveBuiltOwnedFacets = true;
  myHaveCheckedOwnedFacets = false;
  myPreparedQuery.reset();
}

bool MeshFacetedSurface::is_last_preparation_valid_for(const BoundingBox & pointBbox, const double truncationLength) const
{
  const bool isValidHere = myPreparedQuery && myPreparedQuery->truncationLength == truncationLength &&
      (!pointBbox.valid() || (myPreparedQuery->pointBbox.valid() && myPreparedQuery->pointBbox.contains(pointBbox)));
  return stk::is_true_on_all_procs(parallel_comm(), isValidHere);
}

FacetCoordinates MeshFacetedSurface::owned_facet_coordinates() const
{
  const unsigned dim = myMesh.spatial_dimension();
  FacetCoordinates facets(dim);
  for (unsigned f = 0; f < myMesh.num_owned_facets(); ++f)
  {
    const auto v = myMesh.facet_vertices(f);
    if (dim == 3)
      facets.emplace_back_3d(myMesh.vertex_coordinates(v[0]),
          myMesh.vertex_coordinates(v[1]),
          myMesh.vertex_coordinates(v[2]));
    else
      facets.emplace_back_2d(myMesh.vertex_coordinates(v[0]), myMesh.vertex_coordinates(v[1]));
  }
  return facets;
}

template <class FACET>
static void accumulate_owned_facet_areas(
    const FacetMesh & mesh, double & minArea, double & maxArea, double & totalArea)
{
  for (unsigned f = 0; f < mesh.num_owned_facets(); ++f)
  {
    const double area = mesh.make_facet<FACET>(f).facet_area();
    totalArea += area;
    minArea = std::min(area, minArea);
    maxArea = std::max(area, maxArea);
  }
}

std::string MeshFacetedSurface::print_sizes() const
{
  const stk::ParallelMachine comm = myMesh.comm();
  const unsigned localNumFacets = myMesh.num_owned_facets();
  unsigned numFacets = 0;
  stk::all_reduce_sum(comm, &localNumFacets, &numFacets, 1);

  double localTotalArea = 0.0;
  double localMaxArea = -1.0;
  double localMinArea = std::numeric_limits<double>::max();
  if (myMesh.spatial_dimension() == 3)
    accumulate_owned_facet_areas<Facet3d>(myMesh, localMinArea, localMaxArea, localTotalArea);
  else
    accumulate_owned_facet_areas<Facet2d>(myMesh, localMinArea, localMaxArea, localTotalArea);
  double minArea = 0., maxArea = 0., totalArea = 0.;
  stk::all_reduce_min(comm, &localMinArea, &minArea, 1);
  stk::all_reduce_max(comm, &localMaxArea, &maxArea, 1);
  stk::all_reduce_sum(comm, &localTotalArea, &totalArea, 1);

  std::ostringstream out;
  out << "\t " << "Global sizes: { " << numFacets << " facets on }" << std::endl;
  out << "\t " << "Local sizes: { " << localNumFacets << " facets on }" << std::endl;
  out << "\t Global areas: { Min = " << minArea << ", Max = " << maxArea
      << ", Total = " << totalArea << " }" << std::endl;
  return out.str();
}

void MeshFacetedSurface::build_local_facets(const std::vector<BoundingBox> & /*procQueryBboxes*/,
    VertexIdsAndCoords & vertexIdsAndCoords,
    FacetVertexIds & facetVertexIds)
{
  vertexIdsAndCoords.swap(myGivenVertexIdsAndCoords);
  facetVertexIds.swap(myGivenFacetVertexIds);
}

MeshFacetedSurface::~MeshFacetedSurface() = default;

void MeshFacetedSurface::build_search_tree(const BoundingBox & localPaddedQueryBbox)
{
  std::vector<unsigned> facetsInBand;
  for (unsigned facet = 0; facet < myMesh.num_facets(); ++facet)
    if (myMesh.facet_bounding_box(facet).intersects(localPaddedQueryBbox))
      facetsInBand.push_back(facet);
  myNumSearchedFacets = facetsInBand.size();

  myTree = std::make_unique<SearchTree<unsigned>>(
      facetsInBand,
      [&](const unsigned facet) { return myMesh.facet_centroid(facet); },
      [&](const unsigned facet, BoundingBox & bbox)
      {
        for (unsigned n = 0; n < myMesh.spatial_dimension(); ++n)
          bbox.accommodate(myMesh.vertex_coordinates(myMesh.facet_vertices(facet)[n]));
      });
}

void MeshFacetedSurface::prepare_to_compute(
    const double /*time*/, const BoundingBox & pointBbox, const double truncationLength)
{
  if (!myHaveBuiltOwnedFacets || myRebuildOwnedFacetsEachPrepare)
  {
    std::vector<BoundingBox> procQueryBboxes;
    BoundingBox::gather_bboxes(pointBbox, procQueryBboxes);
    VertexIdsAndCoords vertexIdsAndCoords;
    FacetVertexIds facetVertexIds;
    build_local_facets(procQueryBboxes, vertexIdsAndCoords, facetVertexIds);
    myMesh.set_owned_facets(vertexIdsAndCoords, facetVertexIds, myGivenVertexVelocities);
    myGivenVertexVelocities.clear();
    myHaveBuiltOwnedFacets = true;
    myHaveCheckedOwnedFacets = false;
  }

  if (!myHaveCheckedOwnedFacets)
  {
    myPreparedQuery.reset();
    // With connected facets (merged vertices) facets and edges are matched by their vertex ids,
    // otherwise by location.
    if (myRemoveCoincidentFacets)
    {
      const size_t numRemoved = myHaveConnectedSurface ?
          remove_coincident_facets_by_vertex_ids(myMesh) :
          remove_coincident_facets_by_location(myMesh);
      if (numRemoved > 0)
        krinolog << "Removed " << numRemoved << " redundant coincident facets." << stk::diag::dendl;
    }
    if (myCheckWatertight)
      report_non_watertight_edges(myMesh.comm(),
          myHaveConnectedSurface ?
              find_non_watertight_edges_by_vertex_ids(myMesh, myWatertightBbox) :
              find_non_watertight_edges_by_location(myMesh, myWatertightBbox));
    if (myCosSharpFeatureAngle)
      mySharpFeaturesFoundHere = find_sharp_features_by_vertex_ids(
          myMesh, *myCosSharpFeatureAngle, myAreOpenBoundaryEdgesFeatures);
    myHaveCheckedOwnedFacets = true;
  }

  if (is_last_preparation_valid_for(pointBbox, truncationLength)) return;

  const std::vector<BoundingBox> procPaddedQueryBboxes = fill_processor_bounding_boxes(
      myMesh.owned_facets_bounding_box(), pointBbox, truncationLength);
  myMesh.update_ghost_facets(procPaddedQueryBboxes);
  build_search_tree(procPaddedQueryBboxes[stk::parallel_machine_rank(myMesh.comm())]);
  if (myCosSharpFeatureAngle) distribute_sharp_features(procPaddedQueryBboxes);
  myPreparedQuery = PreparedQuery{pointBbox, truncationLength};
}

void MeshFacetedSurface::set_sharp_feature_angle(
    const double featureAngleInDegrees, const bool areOpenBoundaryEdgesFeatures)
{
  STK_ThrowRequireMsg(myHaveConnectedSurface && myMesh.has_shared_vertices(),
      "Sharp features are only supported for connected faceted surfaces (with merged vertices).");
  myCosSharpFeatureAngle = std::cos(featureAngleInDegrees * std::numbers::pi / 180.);
  myAreOpenBoundaryEdgesFeatures = areOpenBoundaryEdgesFeatures;
  myHaveCheckedOwnedFacets = false; // find the features at the next prepare_to_compute
}

static BoundingBox bounding_box_of_points(const std::initializer_list<stk::math::Vector3d> & points)
{
  BoundingBox bbox;
  for (auto & point : points)
    bbox.accommodate(point);
  return bbox;
}

void MeshFacetedSurface::distribute_sharp_features(const std::vector<BoundingBox> & procBboxes)
{
  const int numProcs = procBboxes.size();
  stk::CommSparse commSparse(myMesh.comm());
  stk::pack_and_communicate(commSparse,
      [&]()
      {
        for (int proc = 0; proc < numProcs; ++proc)
        {
          stk::CommBuffer & buffer = commSparse.send_buffer(proc);
          for (auto & edge : mySharpFeaturesFoundHere.edges)
            if (procBboxes[proc].intersects(
                    bounding_box_of_points({edge.coords[0], edge.coords[1]})))
            {
              buffer.pack(true);
              buffer.pack(edge);
            }
          for (auto & corner : mySharpFeaturesFoundHere.corners)
            if (procBboxes[proc].contains(corner.coords))
            {
              buffer.pack(false);
              buffer.pack(corner);
            }
        }
      });

  myLocalSharpFeatures.edges.clear();
  myLocalSharpFeatures.corners.clear();
  stk::unpack_communications(commSparse,
      [&](int proc)
      {
        stk::CommBuffer & buffer = commSparse.recv_buffer(proc);
        while (buffer.remaining())
        {
          bool isEdge = false;
          buffer.unpack(isEdge);
          if (isEdge)
            buffer.unpack(myLocalSharpFeatures.edges.emplace_back());
          else
            buffer.unpack(myLocalSharpFeatures.corners.emplace_back());
        }
      });

  const auto & edges = myLocalSharpFeatures.edges;
  std::vector<unsigned> edgeIndices(edges.size());
  std::iota(edgeIndices.begin(), edgeIndices.end(), 0);
  myFeatureEdgeTree = std::make_unique<SearchTree<unsigned>>(
      edgeIndices,
      [&](const unsigned e) { return 0.5 * (edges[e].coords[0] + edges[e].coords[1]); },
      [&](const unsigned e, BoundingBox & bbox)
      {
        bbox.accommodate(edges[e].coords[0]);
        bbox.accommodate(edges[e].coords[1]);
      });
  const auto & corners = myLocalSharpFeatures.corners;
  std::vector<unsigned> cornerIndices(corners.size());
  std::iota(cornerIndices.begin(), cornerIndices.end(), 0);
  myCornerTree = std::make_unique<SearchTree<unsigned>>(
      cornerIndices,
      [&](const unsigned c) { return corners[c].coords; },
      [&](const unsigned c, BoundingBox & bbox) { bbox.accommodate(corners[c].coords); });
}

void MeshFacetedSurface::fill_triangle_intersection_parametric_coordinates(
    const std::array<stk::math::Vector3d, 3> & triNodes,
    std::vector<stk::math::Vector3d> & intParamCoords) const
{
  intParamCoords.clear();
  if (!myCosSharpFeatureAngle) return;
  STK_ThrowRequireMsg(myCornerTree,
      "MeshFacetedSurface: call prepare_to_compute before finding sharp feature intersections.");
  BoundingBox triBbox = bounding_box_of_points({triNodes[0], triNodes[1], triNodes[2]});
  triBbox.pad_epsilon();
  std::vector<unsigned> candidates;
  if (spatial_dimension() == 3)
  {
    // The feature edges crossing the triangle (a face of a tetrahedron).
    myFeatureEdgeTree->get_intersecting_entities(triBbox, candidates);
    const auto [triArea, triNormal] = CalcTriangle3<double>::area_and_normal(triNodes);
    for (const unsigned e : candidates)
      append_triangle_edge_intersection_parametric_coordinates(triNodes,
          triArea,
          triNormal,
          myLocalSharpFeatures.edges[e].coords[0],
          myLocalSharpFeatures.edges[e].coords[1],
          intParamCoords);
  }
  else
  {
    // The corners in the triangle (a 2D element).
    myCornerTree->get_intersecting_entities(triBbox, candidates);
    for (const unsigned c : candidates)
      append_triangle_intersection_parametric_coordinates(
          triNodes, myLocalSharpFeatures.corners[c].coords, intParamCoords);
  }
}

void MeshFacetedSurface::fill_tetrahedon_intersection_parametric_coordinates(
    const std::array<stk::math::Vector3d, 4> & tetNodes,
    std::vector<stk::math::Vector3d> & intParamCoords) const
{
  intParamCoords.clear();
  if (!myCosSharpFeatureAngle) return;
  STK_ThrowRequireMsg(myCornerTree,
      "MeshFacetedSurface: call prepare_to_compute before finding sharp feature intersections.");
  BoundingBox tetBbox =
      bounding_box_of_points({tetNodes[0], tetNodes[1], tetNodes[2], tetNodes[3]});
  tetBbox.pad_epsilon();
  std::vector<unsigned> candidates;
  myCornerTree->get_intersecting_entities(tetBbox, candidates);
  for (const unsigned c : candidates)
    append_tetrahedron_intersection_parametric_coordinates(
        tetNodes, myLocalSharpFeatures.corners[c].coords, intParamCoords);
}

std::vector<unsigned> MeshFacetedSurface::find_nearest_facets(
    const stk::math::Vector3d & x, const double maxSearchRadius) const
{
  STK_ThrowRequireMsg(
      myTree, "MeshFacetedSurface: call prepare_to_compute before computing distances.");
  std::vector<unsigned> nearestFacets;
  myTree->find_closest_entities(x, nearestFacets, maxSearchRadius);
  return nearestFacets;
}

template <class FACET>
std::vector<FACET> MeshFacetedSurface::build_facets(const std::vector<unsigned> & facets) const
{
  std::vector<FACET> result;
  result.reserve(facets.size());
  for (auto facet : facets)
    result.push_back(myMesh.make_facet<FACET>(facet));
  return result;
}

void MeshFacetedSurface::add_to_journal(CubitJournal & journal, const BoundingBox & region) const
{
  if (!myTree) return;
  std::vector<unsigned> facetsInRegion;
  myTree->get_intersecting_entities(region, facetsInRegion);
  journal.add_comment(std::to_string(facetsInRegion.size()) + " facets of a faceted surface");
  if (spatial_dimension() == 3)
  {
    for (const Facet3d & facet : build_facets<Facet3d>(facetsInRegion))
      journal.add_tri({{facet.facet_vertex(0), facet.facet_vertex(1), facet.facet_vertex(2)}});
  }
  else
  {
    for (const Facet2d & facet : build_facets<Facet2d>(facetsInRegion))
      journal.add_segment(facet.facet_vertex(0), facet.facet_vertex(1));
  }
}

template <class FACET>
static std::vector<const FACET *> pointers_to(const std::vector<FACET> & facets)
{
  std::vector<const FACET *> facetPtrs;
  facetPtrs.reserve(facets.size());
  for (auto & facet : facets)
    facetPtrs.push_back(&facet);
  return facetPtrs;
}

template <class FACET>
stk::math::Vector3d MeshFacetedSurface::compute_connected_pseudo_normal(
    const unsigned closestFacet, const stk::math::Vector3d & closestPtWeights) const
{
  // The closest point is on the vertex or edge of closestFacet given by its vertices with
  // non-negligible weight (the tolerance of is_closest_point_on_facet_edge_or_vertex).
  constexpr double edgeTol = 1.e-6;
  const auto & closestFacetVertices = myMesh.facet_vertices(closestFacet);
  std::array<unsigned, 2> featureVertices{{0, 0}};
  unsigned numFeatureVertices = 0;
  for (int n = 0; n < FACET::DIM; ++n)
    if (closestPtWeights[n] >= edgeTol)
      featureVertices[numFeatureVertices++] = closestFacetVertices[n];
  STK_ThrowRequire(numFeatureVertices == 1 || numFeatureVertices == 2);
  const bool isVertex = numFeatureVertices == 1;

  // Every facet sharing the vertex or edge is on this processor: it contains the closest point,
  // which is inside the padded query box, so it is owned or ghosted here.
  stk::math::Vector3d pseudoNormal = stk::math::Vector3d::ZERO;
  const unsigned vertex = featureVertices[0];
  for (const unsigned facetIndex : myMesh.vertex_facets(vertex))
  {
    const auto facetVertices = myMesh.facet_vertices(facetIndex);
    if (!isVertex &&
        std::find(facetVertices.begin(), facetVertices.end(), featureVertices[1]) ==
            facetVertices.end())
      continue;

    const FACET facet = myMesh.make_facet<FACET>(facetIndex);
    if (facet.degenerate()) continue;

    if constexpr (3 == FACET::DIM)
    {
      if (isVertex)
      {
        // Weight by the facet's angle at the vertex.
        const int n0 =
            std::find(facetVertices.begin(), facetVertices.end(), vertex) - facetVertices.begin();
        const stk::math::Vector3d edge0 = facet.facet_vertex((n0 + 1) % 3) - facet.facet_vertex(n0);
        const stk::math::Vector3d edge1 = facet.facet_vertex((n0 + 2) % 3) - facet.facet_vertex(n0);
        const double cosAngle =
            std::clamp(Dot(edge0, edge1) / (edge0.length() * edge1.length()), -1., 1.);
        pseudoNormal += facet.facet_normal() * std::acos(cosAngle);
        continue;
      }
    }
    pseudoNormal += facet.facet_normal();
  }
  return pseudoNormal;
}

template <class FACET>
stk::math::Vector3d MeshFacetedSurface::compute_pseudo_normal_at_closest_point(
    const std::vector<FacetDistanceQuery<FACET>> & queries,
    const FacetDistanceQuery<FACET> & closestQuery,
    const unsigned closestFacet) const
{
  if (myHaveConnectedSurface)
    return compute_connected_pseudo_normal<FACET>(
        closestFacet, closestQuery.closest_point_weights());
  return compute_pseudo_normal(queries, closestQuery);
}

template <class FACET>
double MeshFacetedSurface::point_distance(const stk::math::Vector3d & x,
    const double truncationLength,
    const double farFieldValue,
    const bool computeSignedDistance) const
{
  const std::vector<unsigned> nearestFacetIndices = find_nearest_facets(x, truncationLength);
  const std::vector<FACET> nearestFacets = build_facets<FACET>(nearestFacetIndices);
  if (!computeSignedDistance || nearestFacets.empty())
    return point_distance_given_nearest_facets(
        x, pointers_to(nearestFacets), truncationLength, farFieldValue, computeSignedDistance);

  const std::vector<FacetDistanceQuery<FACET>> queries =
      build_distance_queries(x, pointers_to(nearestFacets));
  STK_ThrowRequireMsg(
      !queries.empty(), "All facets are degenerate in MeshFacetedSurface::point_distance.");
  const FacetDistanceQuery<FACET> & closestQuery = find_closest_facet(queries);
  const unsigned closestFacet = nearestFacetIndices[&closestQuery.facet() - nearestFacets.data()];

  const double dist = signed_distance_given_closest_facet(x,
      closestQuery,
      [&]()
      {
        return compute_pseudo_normal_at_closest_point<FACET>(queries, closestQuery, closestFacet);
      });
  return (0.0 != truncationLength && std::abs(dist) > truncationLength) ? farFieldValue : dist;
}

double MeshFacetedSurface::possibly_truncated_point_signed_distance(
    const stk::math::Vector3d & x, const double truncationLength, const double farFieldValue) const
{
  return (myMesh.spatial_dimension() == 3) ?
      point_distance<Facet3d>(x, truncationLength, farFieldValue, true) :
      point_distance<Facet2d>(x, truncationLength, farFieldValue, true);
}

double MeshFacetedSurface::point_unsigned_distance(
    const stk::math::Vector3d & x, const double truncationLength, const double farFieldValue) const
{
  return (myMesh.spatial_dimension() == 3) ?
      point_distance<Facet3d>(x, truncationLength, farFieldValue, false) :
      point_distance<Facet2d>(x, truncationLength, farFieldValue, false);
}

template <class FACET>
stk::math::Vector3d MeshFacetedSurface::closest_point_or_normal(
    const stk::math::Vector3d & x, const bool computeNormal) const
{
  const std::vector<unsigned> nearestFacetIndices = find_nearest_facets(x, 0.);
  const std::vector<FACET> nearestFacets = build_facets<FACET>(nearestFacetIndices);
  if (nearestFacets.empty()) return stk::math::Vector3d::ZERO;
  if (!computeNormal) return compute_closest_point(x, pointers_to(nearestFacets));

  const std::vector<FacetDistanceQuery<FACET>> queries =
      build_distance_queries(x, pointers_to(nearestFacets));
  STK_ThrowRequireMsg(!queries.empty(),
      "All facets are degenerate in MeshFacetedSurface::closest_point_or_normal.");
  const FacetDistanceQuery<FACET> & closestQuery = find_closest_facet(queries);
  const unsigned closestFacet = nearestFacetIndices[&closestQuery.facet() - nearestFacets.data()];

  return closest_point_normal_given_closest_facet(x,
      closestQuery,
      [&]()
      {
        return compute_pseudo_normal_at_closest_point<FACET>(queries, closestQuery, closestFacet);
      });
}

stk::math::Vector3d MeshFacetedSurface::closest_point(const stk::math::Vector3d & x) const
{
  return (myMesh.spatial_dimension() == 3) ? closest_point_or_normal<Facet3d>(x, false) :
                                             closest_point_or_normal<Facet2d>(x, false);
}

template <class FACET>
stk::math::Vector3d
MeshFacetedSurface::velocity_at_closest_point_for_facets(const stk::math::Vector3d & x) const
{
  const std::vector<unsigned> nearestFacetIndices = find_nearest_facets(x, 0.);
  const std::vector<FACET> nearestFacets = build_facets<FACET>(nearestFacetIndices);
  if (nearestFacets.empty()) return stk::math::Vector3d::ZERO;
  const std::vector<FacetDistanceQuery<FACET>> queries =
      build_distance_queries(x, pointers_to(nearestFacets));
  STK_ThrowRequireMsg(!queries.empty(),
      "All facets are degenerate in MeshFacetedSurface::velocity_at_closest_point.");
  return find_closest_facet(queries).facet().velocity_at_closest_point(x);
}

stk::math::Vector3d
MeshFacetedSurface::velocity_at_closest_point(const stk::math::Vector3d & x) const
{
  STK_ThrowRequireMsg(myMesh.has_vertex_velocities(),
      "MeshFacetedSurface::velocity_at_closest_point requires vertex velocities.");
  return (myMesh.spatial_dimension() == 3) ?
      velocity_at_closest_point_for_facets<FacetWithVelocity3d>(x) :
      velocity_at_closest_point_for_facets<FacetWithVelocity2d>(x);
}

stk::math::Vector3d MeshFacetedSurface::closest_point_normal(const stk::math::Vector3d & x) const
{
  return (myMesh.spatial_dimension() == 3) ? closest_point_or_normal<Facet3d>(x, true) :
                                             closest_point_or_normal<Facet2d>(x, true);
}

template <class FACET>
std::pair<int, double> MeshFacetedSurface::compute_intersection_with_segment_for_facets(
    const stk::math::Vector3d & pt0, const stk::math::Vector3d & pt1) const
{
  const double narrowBand = (pt0 - pt1).length();

  const double dist0 = possibly_truncated_point_signed_distance(pt0, narrowBand, narrowBand);
  if (dist0 == narrowBand) return {0, -1.};
  const double dist1 = possibly_truncated_point_signed_distance(pt1, narrowBand, narrowBand);
  if (dist1 == narrowBand) return {0, -1.};

  if (sign_change(dist0, dist1))
  {
    constexpr double tightEnoughTolToCalcDistanceAtSignedDistCrossing = 0.001;
    constexpr double looseEnoughTolToHandleSmallGaps = 0.05;

    BoundingBox segmentBbox;
    segmentBbox.accommodate(pt0);
    segmentBbox.accommodate(pt1);
    segmentBbox.pad_epsilon();
    std::vector<unsigned> candidateFacets;
    myTree->get_intersecting_entities(segmentBbox, candidateFacets);
    const std::vector<FACET> candidates = build_facets<FACET>(candidateFacets);

    const double intersectionLoc =
        compute_intersection_between_surface_facets_and_edge(pointers_to(candidates), pt0, pt1);
    if (intersectionLoc >= 0.) return {sign(dist1), intersectionLoc};
    auto [crossingSign, linearCrossingLoc] =
        compute_surface_intersection_with_crossed_segment_from_signed_distance(
            *this, pt0, pt1, dist0, dist1, tightEnoughTolToCalcDistanceAtSignedDistCrossing);
    const stk::math::Vector3d linearCrossingCoords =
        (1. - linearCrossingLoc) * pt0 + linearCrossingLoc * pt1;
    const double crossingDist = point_signed_distance(linearCrossingCoords);
    if (std::abs(crossingDist) < looseEnoughTolToHandleSmallGaps * (pt1 - pt0).length())
      return {sign(dist1), linearCrossingLoc};
  }
  return {0, -1.};
}

std::pair<int, double>
MeshFacetedSurface::compute_intersection_with_segment(const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const double /*edgeCrossingTol*/) const
{
  return (myMesh.spatial_dimension() == 3) ?
      compute_intersection_with_segment_for_facets<Facet3d>(pt0, pt1) :
      compute_intersection_with_segment_for_facets<Facet2d>(pt0, pt1);
}

template <class FACET>
void MeshFacetedSurface::fill_segment_crossings_for_facets(const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const double edgeCrossingTol,
    std::vector<double> & crossingPositions) const
{
  crossingPositions.clear();
  if (!myTree) return;

  BoundingBox segmentBbox;
  segmentBbox.accommodate(pt0);
  segmentBbox.accommodate(pt1);
  segmentBbox.pad_epsilon();
  std::vector<unsigned> candidateFacets;
  myTree->get_intersecting_entities(segmentBbox, candidateFacets);

  std::vector<std::pair<double, int>> positionsAndSigns;
  for (const FACET & facet : build_facets<FACET>(candidateFacets))
  {
    const auto [crossingSign, position] = compute_facet_edge_intersection(facet, pt0, pt1);
    if (crossingSign != 0 && position > edgeCrossingTol && position < 1. - edgeCrossingTol)
      positionsAndSigns.emplace_back(position, crossingSign);
  }
  std::sort(positionsAndSigns.begin(), positionsAndSigns.end());

  // Along the segment the crossings of a closed surface alternate between entering and leaving.
  // The facets are enlarged, so those sharing an edge where the segment crosses both report the
  // crossing: keep the first of each run in the same direction.
  int previousSign = 0;
  for (const auto & [position, crossingSign] : positionsAndSigns)
  {
    if (crossingSign != previousSign) crossingPositions.push_back(position);
    previousSign = crossingSign;
  }
}

void MeshFacetedSurface::fill_segment_crossings(const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const double edgeCrossingTol,
    std::vector<double> & crossingPositions) const
{
  if (myMesh.spatial_dimension() == 3)
    fill_segment_crossings_for_facets<Facet3d>(pt0, pt1, edgeCrossingTol, crossingPositions);
  else
    fill_segment_crossings_for_facets<Facet2d>(pt0, pt1, edgeCrossingTol, crossingPositions);
}

size_t MeshFacetedSurface::storage_size() const
{
  return sizeof(MeshFacetedSurface) + myMesh.storage_size() + (myTree ? myTree->storage_size() : 0);
}

void default_watertight_bounding_box_from_domain(
    const std::vector<std::shared_ptr<Surface>> & surfaces, const BoundingBox & bbox)
{
  for (const std::shared_ptr<Surface> & surf : surfaces)
  {
    MeshFacetedSurface * facetedSurface = dynamic_cast<MeshFacetedSurface *>(surf.get());
    if (facetedSurface && facetedSurface->get_check_watertight() &&
        !facetedSurface->get_watertight_bounding_box().has_value())
      facetedSurface->set_watertight_bounding_box(bbox);
  }
}

} // namespace krino
