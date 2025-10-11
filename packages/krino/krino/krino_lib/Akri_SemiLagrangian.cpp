#include <Akri_AdaptiveContourTet.hpp>
#include <Akri_AdaptiveContourTri.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_Composite_Surface.hpp>
#include <Akri_ContourElement.hpp>
#include <Akri_SemiLagrangian.hpp>

#include <Akri_Faceted_Surface.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_Sign.hpp>
#include <Akri_String_Function_Expression.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_string_expressions(const std::vector<String_Function_Expression> & interfaceVelocity)
{
  auto extv = [&](const double time, const stk::math::Vector3d & closestPt)
  {
    return evaluate_vector_at_location(time, closestPt, interfaceVelocity);
  };
  return extv;
}

template <typename FACET>
ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_facets_with_velocity(const FacetedSurfaceBase & facets)
{
  const auto & facetsWithVelocity = facets.as_derived_type<FACET>();
  auto extv = [&](const double /*time*/, const stk::math::Vector3d & closestPt)
  {
    const auto * nearest = facetsWithVelocity.get_closest_facet(closestPt);
    return nearest->velocity_at_closest_point(closestPt);
  };
  return extv;
}

ExtensionVelocityFunction build_extension_velocity_at_closest_point_using_facets_with_velocity(const int dim, const FacetedSurfaceBase & facets)
{
  if (dim == 2)
      return build_extension_velocity_at_closest_point_using_facets_with_velocity<FacetWithVelocity2d>(facets);
  return build_extension_velocity_at_closest_point_using_facets_with_velocity<FacetWithVelocity3d>(facets);
}

ExtensionVelocityFunction build_extension_velocity_using_velocity_at_closest_point(const FacetedSurfaceBase & facets, const ExtensionVelocityFunction & velAtClosestPt)
{
  auto extv = [&](const double time, const stk::math::Vector3d & pt)
  {
    const stk::math::Vector3d closestPt = facets.closest_point(pt);
    return velAtClosestPt(time, closestPt);
  };
  return extv;
}

template <typename FACET>
double compute_max_facet_velocity_magnitude(const double time,
    const std::vector<FACET> & facets,
    const std::vector<String_Function_Expression> & interfaceVelocity)
{
  constexpr int NUMFACETNODES = FACET::DIM;
  double maxSqrMag = 0.;
  for (auto & facet : facets)
  {
    for (int n=0; n<NUMFACETNODES; ++n)
    {
      const double velSqrMag = (evaluate_vector_at_location(time, facet.facet_vertex(n), interfaceVelocity)).length_squared();
      if (velSqrMag > maxSqrMag)
        maxSqrMag = velSqrMag;
    }
  }
  return std::sqrt(maxSqrMag);
}

template <typename FACET>
double compute_max_facet_velocity_magnitude_for_facets_with_velocity(const FacetedSurfaceBase & facets)
{
  const auto & facetsWithVelocity = facets.as_derived_type<FACET>();
  double maxSqrMag = 0.;
  for (auto & facet : facetsWithVelocity.get_facets())
  {
    for (auto & vel : facet.get_velocity())
    {
      const double velSqrMag = vel.length_squared();
      if (velSqrMag > maxSqrMag)
        maxSqrMag = velSqrMag;
    }
  }
  return std::sqrt(maxSqrMag);
}

static void pad_bounding_box_using_time_step_and_velocity_magnitude(BoundingBox & bbox, const double dt, const double velocityMagnitude)
{
  const double paddingFactorOfSafety = 1.5;
  bbox.pad(paddingFactorOfSafety*velocityMagnitude*dt);  // Need something better?
}

BoundingBox compute_padded_node_bounding_box_for_semilagrangian_using_string_velocity_expressions(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const std::vector<String_Function_Expression> & interfaceVelocity,
    const FacetedSurfaceBase & facets)
{
  BoundingBox nodeBBox = krino::compute_nodal_bbox(mesh, activeFieldSelector, coordsField);

  const double timeMid = 0.5*(timeN+timeNp1);
  const double velMag = (2==mesh.mesh_meta_data().spatial_dimension()) ?
      compute_max_facet_velocity_magnitude(timeMid, facets.get_facets_2d(), interfaceVelocity) :
      compute_max_facet_velocity_magnitude(timeMid, facets.get_facets_3d(), interfaceVelocity);
  pad_bounding_box_using_time_step_and_velocity_magnitude(nodeBBox, timeNp1-timeN, velMag);
  return nodeBBox;
}

BoundingBox compute_padded_node_bounding_box_for_semilagrangian_using_facets_with_velocity(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double dt,
    const FieldRef coordsField,
    const FacetedSurfaceBase & facets)
{
  BoundingBox nodeBBox = krino::compute_nodal_bbox(mesh, activeFieldSelector, coordsField);

  const double velMag = (2==mesh.mesh_meta_data().spatial_dimension()) ?
      compute_max_facet_velocity_magnitude_for_facets_with_velocity<FacetWithVelocity2d>(facets) :
      compute_max_facet_velocity_magnitude_for_facets_with_velocity<FacetWithVelocity3d>(facets);
  pad_bounding_box_using_time_step_and_velocity_magnitude(nodeBBox, dt, velMag);
  return nodeBBox;
}

static stk::math::Vector3d compute_semilagrangian_departure_point(const int /*dim*/,
    const double timeN,
    const double timeNp1,
    const stk::math::Vector3d & pt,
    const ExtensionVelocityFunction & extension_velocity)
{
  const double dt = timeNp1 - timeN;
  const double tMid = 0.5*(timeN+timeNp1);

#if 1
  // Midpoint
  const auto velN = extension_velocity(timeN, pt);

  const stk::math::Vector3d coordsHalf = pt - 0.5*dt*velN;
  const auto velHalf = extension_velocity(tMid, coordsHalf);

  const stk::math::Vector3d coordsNp1 = pt - dt*velHalf;
#endif
#if 0
  // Trapezoidal
  const auto velN = extension_velocity(timeN, pt);

  const stk::math::Vector3d coordsPred = pt - dt*velN;
  const auto velPred = extension_velocity(timeNp1, coordsPred);

  const stk::math::Vector3d coordsNp1 = pt - 0.5*dt*(velN+velPred);
#endif
#if 0
  const auto velN = extension_velocity(timeN, pt);

  const stk::math::Vector3d coords2 = pt - 0.5*dt*velN;
  const auto vel2 = extension_velocity(tMid, coords2);

  const stk::math::Vector3d coords3 = pt - 0.5*dt*vel2;
  const auto vel3 = extension_velocity(tMid, coords3);

  const stk::math::Vector3d coords4 = pt - dt*vel3;
  const auto vel4 = extension_velocity(timeNp1, coords4);

  const stk::math::Vector3d coordsNp1 = pt - dt/6.*(velN + 2.*vel2 + 2.*vel3 + vel4);
#endif
#if 0
// use local velocity instead of extension velocity (20240916: This would now be accomplished using a different extension velocity)
  const auto velN = evaluate_vector_at_location(dim, timeN, pt, interfaceVelocityExpr);
  const stk::math::Vector3d coordsHalf = pt - 0.5*dt*velN;
  const auto velHalf = evaluate_vector_at_location(dim, tMid, coordsHalf, interfaceVelocityExpr);

  const stk::math::Vector3d coordsNp1 = pt - dt*velHalf;
#endif
#if 0
  const auto velN = extension_velocity(timeN, pt);

  const stk::math::Vector3d coordsHalf = pt - 0.5*dt*velN;
  const auto velHalf = extension_velocity(tMid, coordsHalf);

  const stk::math::Vector3d coordsPred = pt - dt*velN;
  const auto velPred = extension_velocity(tMid, coordsPred);

  const stk::math::Vector3d coordsNp1 = pt - 0.25*dt*(velN+2*velHalf+velPred);
#endif
#if 0
  //BFECC
  const auto velN = extension_velocity(timeN, pt);
  const stk::math::Vector3d coordsBack = pt - dt*velN;

  const auto velBack = extension_velocity(timeNp1, coordsBack);
  const stk::math::Vector3d coordsForth = coordsBack + dt*velBack;

  const stk::math::Vector3d corrected = pt - 0.5*(coordsForth-pt);

  const auto velCorr = extension_velocity(timeN, corrected);
  const stk::math::Vector3d coordsNp1 = corrected - dt*velCorr;
#endif
  return coordsNp1;
}

static stk::math::Vector3d compute_predicted_departure_point(const int /*dim*/,
    const double timeN,
    const double timeNp1,
    const stk::math::Vector3d & pt,
    const ExtensionVelocityFunction & extension_velocity)
{
  const double dt = timeNp1 - timeN;
  const auto velN = extension_velocity(timeN, pt);
  const stk::math::Vector3d coordsTilde = pt - dt*velN;

  return coordsTilde;
}

static stk::math::Vector3d compute_corrected_departure_point(const int /*dim*/,
    const double timeN,
    const double timeNp1,
    const stk::math::Vector3d & pt,
    const ExtensionVelocityFunction & extension_velocity_old,
    const ExtensionVelocityFunction & extension_velocity_pred)
{
  const double dt = timeNp1 - timeN;
  const auto velN = extension_velocity_old(timeN, pt);
  const stk::math::Vector3d coordsTilde = pt - dt*velN;
  const auto vel1 = extension_velocity_old(timeN, coordsTilde);
  const auto vel2 = extension_velocity_pred(timeNp1, pt);
  const auto velCorr = 0.5*(vel1+vel2);
  const stk::math::Vector3d coordsCorr = pt - dt*velCorr;

  return coordsCorr;
}

static std::function<double(const stk::math::Vector3d & pt)> build_initial_distance_at_point(const Composite_Surface & initSurfaces, const double narrowBandSize)
{
  auto fn = [&initSurfaces, narrowBandSize](const stk::math::Vector3d & pt)
    {
      return initSurfaces.point_signed_distance_with_narrow_band(pt, narrowBandSize);
    };
  return fn;
}

static std::function<double(const stk::math::Vector3d & pt)> build_distance_at_departure_point(const FacetedSurfaceBase & facets)
{
  auto fn = [&facets](const stk::math::Vector3d & departurePt)
    {
      constexpr double zeroNarrowBandSize = 0.;
      return facets.truncated_point_signed_distance(departurePt, zeroNarrowBandSize, zeroNarrowBandSize);
    };
  return fn;
}

static std::function<double(const stk::math::Vector3d & pt, const int sign)> build_narrow_band_distance_at_departure_point(const FacetedSurfaceBase & facets, const double narrowBandSize)
{
  auto fn = [&facets, narrowBandSize](const stk::math::Vector3d & departurePt, const int sign)
    {
      return facets.truncated_point_signed_distance(departurePt, narrowBandSize, sign*narrowBandSize);
    };
  return fn;
}

static std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> build_semilagrangian_departure_point_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const ExtensionVelocityFunction & extension_velocity)
{
  auto fn = [dim, timeN, timeNp1, &extension_velocity](const stk::math::Vector3d & pt)
    {
      return compute_semilagrangian_departure_point(dim, timeN, timeNp1, pt, extension_velocity);
    };
  return fn;
}

static std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> build_predicted_departure_point_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const ExtensionVelocityFunction & extension_velocity)
{
  auto fn = [dim, timeN, timeNp1, &extension_velocity](const stk::math::Vector3d & pt)
    {
      return compute_predicted_departure_point(dim, timeN, timeNp1, pt, extension_velocity);
    };
  return fn;
}

static std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> build_corrected_departure_point_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const ExtensionVelocityFunction & extension_velocity_old,
    const ExtensionVelocityFunction & extension_velocity_pred)
{
  auto fn = [dim, timeN, timeNp1, &extension_velocity_old, &extension_velocity_pred](const stk::math::Vector3d & pt)
    {
      return compute_corrected_departure_point(dim, timeN, timeNp1, pt, extension_velocity_old, extension_velocity_pred);
    };
  return fn;
}

static std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> build_departure_point_as_point()
{
  auto fn = [](const stk::math::Vector3d & pt)
    {
      return pt;
    };
  return fn;
}

void build_nonadaptive_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const FieldRef coordsField,
    const FieldRef distField,
    const double lengthScale,
    FacetedSurfaceBase & facets)
{
  facets.clear();
  for ( auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeFieldSelector & mesh.mesh_meta_data().locally_owned_part()) )
  {
    for (auto elem : *bucketPtr)
    {
      ContourElement lsElem( mesh, elem, coordsField, distField );
      lsElem.compute_subelement_decomposition(lengthScale);

      lsElem.build_subelement_facets( facets );
    }
  }
}

static void calc_semilagrangian_nodal_distance(const int dim,
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> & departure_point_at_point,
    const std::function<double(const stk::math::Vector3d & pt, const int sign)> & distance_at_departure_point)
{
  for ( auto && bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, activeFieldSelector) )
  {
    const double * coordsData = field_data<double>(coordsField , *bucketPtr);
    double * distData = field_data<double>(distField, *bucketPtr);

    for (size_t i = 0; i < bucketPtr->size(); ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*dim, dim);
      const int previousSign = sign(distData[i]);
      distData[i] = distance_at_departure_point(departure_point_at_point(nodeCoords), previousSign);
    }
  }
}

void calc_single_step_nonadaptive_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const auto departure_point_at_point = build_semilagrangian_departure_point_at_point(dim, timeN, timeNp1, extension_velocity);
  const auto distance_at_departure_point = build_distance_at_departure_point(facetsN);
  const auto narrow_band_distance_at_departure_point = build_narrow_band_distance_at_departure_point(facetsN, narrowBandSize);

  calc_semilagrangian_nodal_distance(dim, mesh, activeFieldSelector, coordsField, distField, departure_point_at_point, narrow_band_distance_at_departure_point);
  build_nonadaptive_facets(mesh, activeFieldSelector, coordsField, distField, avgEdgeLength, facetsNp1);
}

template <size_t NVERT>
bool has_any_chance_of_cut_edge(const std::array<stk::math::Vector3d,NVERT> & coords,
    const std::array<double,NVERT> & dist)
{
  for (size_t n1=0; n1<NVERT; ++n1)
  {
    for (size_t n2=n1; n2<NVERT; ++n2)
    {
      const double sqrLen = (coords[n1] - coords[n2]).length_squared();
      if (dist[n1]*dist[n1] <= sqrLen || dist[n2]*dist[n2] <= sqrLen)
        return true;
    }
  }
  return false;
}

static void adaptively_append_facets_for_tri_mesh_element_using_semilagrangian_distance(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef isoField,
    const stk::mesh::Entity elem,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> & departure_point_at_point,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_departure_point,
    const double lengthScale,
    const int minDepth,
    const int maxDepth,
    FacetedSurfaceBase & facets)
{
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
  const std::array<stk::math::Vector3d,3> nodeCoords = get_triangle_vector(mesh, coordsField, elemNodes, 2);
  const std::array<double,3> nodeDist = get_triangle_scalar(mesh, isoField, elemNodes);

  if (has_any_chance_of_cut_edge(nodeCoords, nodeDist))
  {
    const std::array<stk::math::Vector3d,3> nodalDepartureCoords = {departure_point_at_point(nodeCoords[0]), departure_point_at_point(nodeCoords[1]), departure_point_at_point(nodeCoords[2])};
    adaptively_append_facets_for_tri_using_semilagrangian_distance(nodeCoords, nodalDepartureCoords, nodeDist, distance_at_departure_point, lengthScale, facets, 0, minDepth, maxDepth);
  }
}

static void adaptively_append_facets_for_tet_mesh_element_using_semilagrangian_distance(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef isoField,
    const stk::mesh::Entity elem,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> & departure_point_at_point,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_departure_point,
    const double lengthScale,
    const int minDepth,
    const int maxDepth,
    FacetedSurfaceBase & facets)
{
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
  const std::array<stk::math::Vector3d,4> nodeCoords = get_tetrahedron_vector(mesh, coordsField, elemNodes);
  const std::array<double,4> nodeDist = get_tetrahedron_scalar(mesh, isoField, elemNodes);

  if (has_any_chance_of_cut_edge(nodeCoords, nodeDist))
  {
    const std::array<stk::math::Vector3d,4> nodalDepartureCoords = {departure_point_at_point(nodeCoords[0]), departure_point_at_point(nodeCoords[1]), departure_point_at_point(nodeCoords[2]), departure_point_at_point(nodeCoords[3])};
    adaptively_append_facets_for_tet_using_semilagrangian_distance(nodeCoords, nodalDepartureCoords, nodeDist, distance_at_departure_point, lengthScale, facets, 0, minDepth, maxDepth);
  }
}

static void adaptively_append_facets_for_mesh_element_using_semilagrangian_distance(const int dim,
    const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef isoField,
    const stk::mesh::Entity elem,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> & departure_point_at_point,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_departure_point,
    const double lengthScale,
    const int minDepth,
    const int maxDepth,
    FacetedSurfaceBase & facets)
{
  if (2 == dim)
  {
    adaptively_append_facets_for_tri_mesh_element_using_semilagrangian_distance(mesh, coordsField, isoField, elem, departure_point_at_point, distance_at_departure_point, lengthScale, minDepth, maxDepth, facets);
  }
  else
  {
    adaptively_append_facets_for_tet_mesh_element_using_semilagrangian_distance(mesh, coordsField, isoField, elem, departure_point_at_point, distance_at_departure_point, lengthScale, minDepth, maxDepth, facets);
  }
}

static void build_adaptive_facets_using_semilagrangian_distance(const int dim,
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d & pt)> & departure_point_at_point,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_departure_point,
    const double avgEdgeLength,
    const int minDepth,
    const int maxDepth,
    FacetedSurfaceBase & facets)
{
  facets.clear();
  for ( auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeFieldSelector) )
  {
    STK_ThrowRequireMsg(bucketPtr->topology() == stk::topology::TRIANGLE_3_2D || bucketPtr->topology() == stk::topology::TETRAHEDRON_4, "Only Tri3 and Tet4 elements currently supported.");
    for (auto elem : *bucketPtr)
      adaptively_append_facets_for_mesh_element_using_semilagrangian_distance(dim, mesh, coordsField, distField, elem, departure_point_at_point, distance_at_departure_point, avgEdgeLength, minDepth, maxDepth, facets);
  }
}

void build_initial_adaptive_facets_after_nodal_distance_is_initialized_from_initial_surfaces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double time,
    const FieldRef coordsField,
    const FieldRef distField,
    const double avgEdgeLength,
    const Composite_Surface & initSurfaces,
    FacetedSurfaceBase & facets)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const int minDepth = (2 == dim) ? 5 : 4;
  const int maxDepth = 5;

  const auto departure_point_at_point = build_departure_point_as_point();
  const auto initial_distance_at_point = build_initial_distance_at_point(initSurfaces, time);

  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, departure_point_at_point, initial_distance_at_point, avgEdgeLength, minDepth, maxDepth, facets);
}

void calc_single_step_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1)
{
  const int minDepth = 2;
  const int maxDepth = 5;

  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const auto distance_at_departure_point = build_distance_at_departure_point(facetsN);
  const auto narrow_band_distance_at_departure_point = build_narrow_band_distance_at_departure_point(facetsN, narrowBandSize);
  const auto departure_point_at_point = build_semilagrangian_departure_point_at_point(dim, timeN, timeNp1, extension_velocity);

  calc_semilagrangian_nodal_distance(dim, mesh, activeFieldSelector, coordsField, distField, departure_point_at_point, narrow_band_distance_at_departure_point);
  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, departure_point_at_point, distance_at_departure_point, avgEdgeLength, minDepth, maxDepth, facetsNp1);
}

void predict_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsPred)
{
  const int minDepth = 1;
  const int maxDepth = 2;

  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const auto distance_at_departure_point = build_distance_at_departure_point(facetsN);
  const auto narrow_band_distance_at_departure_point = build_narrow_band_distance_at_departure_point(facetsN, narrowBandSize);
  const auto predicted_departure_point_at_point = build_predicted_departure_point_at_point(dim, timeN, timeNp1, extension_velocity);

  calc_semilagrangian_nodal_distance(dim, mesh, activeFieldSelector, coordsField, distField, predicted_departure_point_at_point, narrow_band_distance_at_departure_point);
  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, predicted_departure_point_at_point, distance_at_departure_point, avgEdgeLength, minDepth, maxDepth, facetsPred);
}

void correct_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const ExtensionVelocityFunction & extension_velocity_old,
    const ExtensionVelocityFunction & extension_velocity_pred,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1)
{
  const int minDepth = 2;
  const int maxDepth = 5;

  const int dim = mesh.mesh_meta_data().spatial_dimension();

  const auto distance_at_departure_point = build_distance_at_departure_point(facetsN);
  const auto narrow_band_distance_at_departure_point = build_narrow_band_distance_at_departure_point(facetsN, narrowBandSize);
  const auto corrected_departure_point_at_point = build_corrected_departure_point_at_point(dim, timeN, timeNp1, extension_velocity_old, extension_velocity_pred);

  calc_semilagrangian_nodal_distance(dim, mesh, activeFieldSelector, coordsField, distField, corrected_departure_point_at_point, narrow_band_distance_at_departure_point);
  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, corrected_departure_point_at_point, distance_at_departure_point, avgEdgeLength, minDepth, maxDepth, facetsNp1);
}

}


