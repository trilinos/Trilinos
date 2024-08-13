#include <Akri_AdaptiveElementContour.hpp>
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

static stk::math::Vector3d compute_interface_velocity_at_point(const int dim, const double time, const stk::math::Vector3d & coords, const std::vector<String_Function_Expression> & interfaceVelocityExpr)
{
  if (2 == dim)
    return stk::math::Vector3d(interfaceVelocityExpr[0].evaluate(time, coords), interfaceVelocityExpr[1].evaluate(time, coords), 0.0);
  return stk::math::Vector3d(interfaceVelocityExpr[0].evaluate(time, coords), interfaceVelocityExpr[1].evaluate(time, coords), interfaceVelocityExpr[2].evaluate(time, coords));
}

static stk::math::Vector3d compute_semilagrangian_evaluation_point(const int dim,
    const double timeN,
    const double timeNp1,
    const FacetedSurfaceBase & facets,
    const stk::math::Vector3d & pt,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr)
{
  const double dt = timeNp1 - timeN;
  const double tMid = 0.5*(timeN+timeNp1);

#if 0
  const stk::math::Vector3d closestPtN = facets.closest_point(pt);
  const auto velN = compute_interface_velocity_at_point(dim, timeN, closestPtN, interfaceVelocityExpr);

  const stk::math::Vector3d coords2 = pt - 0.5*dt*velN;
  const stk::math::Vector3d closestPt2 = facets.closest_point(coords2);
  const auto vel2 = compute_interface_velocity_at_point(dim, tMid, closestPt2, interfaceVelocityExpr);

  const stk::math::Vector3d coords3 = pt - 0.5*dt*vel2;
  const stk::math::Vector3d closestPt3 = facets.closest_point(coords3);
  const auto vel3 = compute_interface_velocity_at_point(dim, tMid, closestPt3, interfaceVelocityExpr);

  const stk::math::Vector3d coords4 = pt - dt*vel3;
  const stk::math::Vector3d closestPt4 = facets.closest_point(coords4);
  const auto vel4 = compute_interface_velocity_at_point(dim, timeNp1, closestPt4, interfaceVelocityExpr);

  const stk::math::Vector3d coordsNp1 = pt - dt/6.*(velN + 2.*vel2 + 2.*vel3 + vel4);
#endif
#if 1
  const stk::math::Vector3d closestPtN = facets.closest_point(pt);
  const auto velN = compute_interface_velocity_at_point(dim, timeN, closestPtN, interfaceVelocityExpr);

  const stk::math::Vector3d coordsHalf = pt - 0.5*dt*velN;
  const stk::math::Vector3d closestPtHalf = facets.closest_point(coordsHalf);
  const auto velHalf = compute_interface_velocity_at_point(dim, tMid, closestPtHalf, interfaceVelocityExpr);

  const stk::math::Vector3d coordsNp1 = pt - dt*velHalf;
#endif
#if 0
// use local velocity instead of extension velocity
  const auto velN = compute_interface_velocity_at_point(dim, timeN, pt, interfaceVelocityExpr);
  const stk::math::Vector3d coordsHalf = pt - 0.5*dt*velN;
  const auto velHalf = compute_interface_velocity_at_point(dim, tMid, coordsHalf, interfaceVelocityExpr);

  const stk::math::Vector3d coordsNp1 = pt - dt*velHalf;
#endif
#if 0
  const stk::math::Vector3d closestPtN = facets.closest_point(pt);
  const auto velN = compute_interface_velocity_at_point(dim, timeN, closestPtN, interfaceVelocityExpr);

  const stk::math::Vector3d coordsPred = pt - dt*velN;
  const stk::math::Vector3d closestPtPred = facets.closest_point(coordsPred);
  const auto velPred = compute_interface_velocity_at_point(dim, tMid, closestPtPred, interfaceVelocityExpr);

  const stk::math::Vector3d coordsNp1 = pt - 0.5*dt*(velN+velPred);
#endif
#if 0
  const stk::math::Vector3d closestPtN = facets.closest_point(pt);
  const auto velN = compute_interface_velocity_at_point(dim, timeN, closestPtN, interfaceVelocityExpr);

  const stk::math::Vector3d coordsHalf = pt - 0.5*dt*velN;
  const stk::math::Vector3d closestPtHalf = facets.closest_point(coordsHalf);
  const auto velHalf = compute_interface_velocity_at_point(dim, tMid, closestPtHalf, interfaceVelocityExpr);

  const stk::math::Vector3d coordsPred = pt - dt*velN;
  const stk::math::Vector3d closestPtPred = facets.closest_point(coordsPred);
  const auto velPred = compute_interface_velocity_at_point(dim, tMid, closestPtPred, interfaceVelocityExpr);

  const stk::math::Vector3d coordsNp1 = pt - 0.25*dt*(velN+2*velHalf+velPred);
#endif
  return coordsNp1;
}

static double compute_semilagrangian_distance_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const FacetedSurfaceBase & facets,
    const stk::math::Vector3d & pt,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double farFieldValue)
{
  const auto prevCoords = compute_semilagrangian_evaluation_point(dim, timeN, timeNp1, facets, pt, interfaceVelocityExpr);
  return facets.truncated_point_signed_distance(prevCoords, narrowBandSize, farFieldValue);
}

static stk::math::Vector3d compute_semilagrangian_predicted_evaluation_point(const int dim,
    const double timeN,
    const double timeNp1,
    const FacetedSurfaceBase & facetsN,
    const stk::math::Vector3d & pt,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr)
{
  const double dt = timeNp1 - timeN;
  const stk::math::Vector3d closestPtN = facetsN.closest_point(pt);
  const auto velN = compute_interface_velocity_at_point(dim, timeN, closestPtN, interfaceVelocityExpr);
  const stk::math::Vector3d coordsTilde = pt - dt*velN;

  return coordsTilde;
}

static stk::math::Vector3d compute_semilagrangian_corrected_evaluation_point(const int dim,
    const double timeN,
    const double timeNp1,
    const FacetedSurfaceBase & facetsN,
    const FacetedSurfaceBase & facetsPred,
    const stk::math::Vector3d & pt,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr)
{
  const double dt = timeNp1 - timeN;
  const stk::math::Vector3d facetsNClosestPt = facetsN.closest_point(pt);
  const auto velN = compute_interface_velocity_at_point(dim, timeN, facetsNClosestPt, interfaceVelocityExpr);
  const stk::math::Vector3d coordsTilde = pt - dt*velN;
  const stk::math::Vector3d facetsNClosestPtTilde = facetsN.closest_point(coordsTilde);
  const auto vel1 = compute_interface_velocity_at_point(dim, timeN, facetsNClosestPtTilde, interfaceVelocityExpr);
  const stk::math::Vector3d facetsPredClosestPt = facetsPred.closest_point(pt);
  const auto vel2 = compute_interface_velocity_at_point(dim, timeNp1, facetsPredClosestPt, interfaceVelocityExpr);
  const auto velCorr = 0.5*(vel1+vel2);
  const stk::math::Vector3d coordsCorr = pt - dt*velCorr;

  return coordsCorr;
}

static double compute_semilagrangian_distance_prediction_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const FacetedSurfaceBase & facetsN,
    const stk::math::Vector3d & pt,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double farFieldValue)
{
  const stk::math::Vector3d coordsTilde = compute_semilagrangian_predicted_evaluation_point(dim, timeN, timeNp1, facetsN, pt, interfaceVelocityExpr);
  return facetsN.truncated_point_signed_distance(coordsTilde, narrowBandSize, farFieldValue);
}

static double compute_semilagrangian_distance_correction_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const FacetedSurfaceBase & facetsN,
    const FacetedSurfaceBase & facetsPred,
    const stk::math::Vector3d & pt,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double farFieldValue)
{
  const stk::math::Vector3d coordsCorr = compute_semilagrangian_corrected_evaluation_point(dim, timeN, timeNp1, facetsN, facetsPred, pt, interfaceVelocityExpr);
  return facetsN.truncated_point_signed_distance(coordsCorr, narrowBandSize, farFieldValue);
}

static std::function<double(const stk::math::Vector3d & pt)> build_initial_distance_at_point(const Composite_Surface & initSurfaces, const double narrowBandSize)
{
  auto fn = [&initSurfaces, narrowBandSize](const stk::math::Vector3d & pt)
    {
      return initSurfaces.point_signed_distance_with_narrow_band(pt, narrowBandSize);
    };
  return fn;
}

static std::function<double(const stk::math::Vector3d & pt)> build_semilagrangian_distance_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const FacetedSurfaceBase & facets)
{
  auto fn = [dim, timeN, timeNp1, &interfaceVelocityExpr, &facets](const stk::math::Vector3d & pt)
    {
      constexpr double zeroNarrowBandSize = 0.;
      return compute_semilagrangian_distance_at_point(dim, timeN, timeNp1, facets, pt, interfaceVelocityExpr, zeroNarrowBandSize, zeroNarrowBandSize);
    };
  return fn;
}

static std::function<double(const stk::math::Vector3d & pt)> build_semilagrangian_distance_predictor_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const FacetedSurfaceBase & facets)
{
  auto fn = [dim, timeN, timeNp1, &interfaceVelocityExpr, &facets](const stk::math::Vector3d & pt)
    {
      constexpr double zeroNarrowBandSize = 0.;
      return compute_semilagrangian_distance_prediction_at_point(dim, timeN, timeNp1, facets, pt, interfaceVelocityExpr, zeroNarrowBandSize, zeroNarrowBandSize);
    };
  return fn;
}

static std::function<double(const stk::math::Vector3d & pt)> build_semilagrangian_distance_corrector_at_point(const int dim,
    const double timeN,
    const double timeNp1,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const FacetedSurfaceBase & facetsN,
    const FacetedSurfaceBase & facetsPred)
{
  auto fn = [dim, timeN, timeNp1, &interfaceVelocityExpr, &facetsN, &facetsPred](const stk::math::Vector3d & pt)
    {
      constexpr double zeroNarrowBandSize = 0.;
      return compute_semilagrangian_distance_correction_at_point(dim, timeN, timeNp1, facetsN, facetsPred, pt, interfaceVelocityExpr, zeroNarrowBandSize, zeroNarrowBandSize);
    };
  return fn;
}

template <typename FACET>
double compute_max_facet_velocity_magnitude(const double time,
    const std::vector<FACET> & facets,
    const std::vector<String_Function_Expression> & interfaceVelocity)
{
  double maxSqrMag = 0.;
  for (auto & facet : facets)
  {
    for (int n=0; n<FACET::DIM; ++n)
    {
      const double velSqrMag = (compute_interface_velocity_at_point(FACET::DIM, time, facet.facet_vertex(n), interfaceVelocity)).length_squared();
      if (velSqrMag > maxSqrMag)
        maxSqrMag = velSqrMag;
    }
  }
  return std::sqrt(maxSqrMag);
}

BoundingBox compute_padded_node_bounding_box_for_semilagrangian(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const FacetedSurfaceBase & facets)
{
  BoundingBox nodeBBox = krino::compute_nodal_bbox(mesh, activeFieldSelector, coordsField);

  const double timeMid = 0.5*(timeN+timeNp1);
  const double velMag = (2==mesh.mesh_meta_data().spatial_dimension()) ?
      compute_max_facet_velocity_magnitude(timeMid, facets.get_facets_2d(), interfaceVelocityExpr) :
      compute_max_facet_velocity_magnitude(timeMid, facets.get_facets_3d(), interfaceVelocityExpr);
  const double paddingFactorOfSafety = 1.5;
  nodeBBox.pad(paddingFactorOfSafety*velMag*(timeNp1-timeN));  // Need something better?
  return nodeBBox;
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

static void calc_single_step_semilagrangian_nodal_distance(const int dim,
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const FacetedSurfaceBase & facetsN)
{
  for ( auto && bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, activeFieldSelector) )
  {
    const double * coordsData = field_data<double>(coordsField , *bucketPtr);
    double * distData = field_data<double>(distField, *bucketPtr);

    for (size_t i = 0; i < bucketPtr->size(); ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*dim, dim);
      const int previousSign = sign(distData[i]);
      distData[i] = compute_semilagrangian_distance_at_point(dim, timeN, timeNp1, facetsN, nodeCoords, interfaceVelocityExpr, narrowBandSize, previousSign*narrowBandSize);
    }
  }
}

void calc_single_step_nonadaptive_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  calc_single_step_semilagrangian_nodal_distance(dim, mesh, activeFieldSelector, timeN, timeNp1, coordsField, distField, interfaceVelocityExpr, narrowBandSize, facetsN);

  build_nonadaptive_facets(mesh, activeFieldSelector, coordsField, distField, avgEdgeLength, facetsNp1);
}

static void adaptively_append_facets_for_mesh_element_using_semilagrangian_distance(const int dim,
    const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef isoField,
    const stk::mesh::Entity elem,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
    const double lengthScale,
    const int minDepth,
    const int maxDepth,
    FacetedSurfaceBase & facets)
{
  STK_ThrowRequire(dim == 2);
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
  const std::array<stk::math::Vector3d,3> nodeCoords = get_triangle_vector(mesh, coordsField, elemNodes, 2);
  const std::array<double,3> nodeDist = get_triangle_scalar(mesh, isoField, elemNodes);
  adaptively_append_facets_for_tri_using_semilagrangian_distance(nodeCoords, nodeDist, distance_at_point, lengthScale, facets, 0, minDepth, maxDepth);
}

static void build_adaptive_facets_using_semilagrangian_distance(const int dim,
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::function<double(const stk::math::Vector3d & pt)> & distance_at_point,
    const double avgEdgeLength,
    const int minDepth,
    const int maxDepth,
    FacetedSurfaceBase & facets)
{
  facets.clear();
  for ( auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeFieldSelector) )
  {
    STK_ThrowRequireMsg(bucketPtr->topology() == stk::topology::TRIANGLE_3_2D, "Only Tri3d elements currently supported.");
    for (auto elem : *bucketPtr)
      adaptively_append_facets_for_mesh_element_using_semilagrangian_distance(dim, mesh, coordsField, distField, elem, distance_at_point, avgEdgeLength, minDepth, maxDepth, facets);
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
  const int minDepth = 5;
  const int maxDepth = 5;

  const auto initial_distance_at_point = build_initial_distance_at_point(initSurfaces, time);

  build_adaptive_facets_using_semilagrangian_distance(mesh.mesh_meta_data().spatial_dimension(), mesh, activeFieldSelector, coordsField, distField, initial_distance_at_point, avgEdgeLength, minDepth, maxDepth, facets);
}

void calc_single_step_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsNp1)
{
  const int minDepth = 2;
  const int maxDepth = 5;

  const int dim = mesh.mesh_meta_data().spatial_dimension();
  calc_single_step_semilagrangian_nodal_distance(dim, mesh, activeFieldSelector, timeN, timeNp1, coordsField, distField, interfaceVelocityExpr, narrowBandSize, facetsN);

  const auto distance_at_point = build_semilagrangian_distance_at_point(dim, timeN, timeNp1, interfaceVelocityExpr, facetsN);

  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, distance_at_point, avgEdgeLength, minDepth, maxDepth, facetsNp1);
}

void predict_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    FacetedSurfaceBase & facetsPred)
{
  const int minDepth = 1;
  const int maxDepth = 2;

  const int dim = mesh.mesh_meta_data().spatial_dimension();
  for ( auto && bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, activeFieldSelector) )
  {
    const double * coordsData = field_data<double>(coordsField , *bucketPtr);
    double * distData = field_data<double>(distField , *bucketPtr);

    for (size_t i = 0; i < bucketPtr->size(); ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*dim, dim);
      const int previousSign = sign(distData[i]);
      distData[i] = compute_semilagrangian_distance_prediction_at_point(dim, timeN, timeNp1, facetsN, nodeCoords, interfaceVelocityExpr, narrowBandSize, previousSign*narrowBandSize);
    }
  }

  const auto predict_distance_at_point = build_semilagrangian_distance_predictor_at_point(dim, timeN, timeNp1, interfaceVelocityExpr, facetsN);

  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, predict_distance_at_point, avgEdgeLength, minDepth, maxDepth, facetsPred);
}

void correct_semilagrangian_nodal_distance_and_build_facets(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeFieldSelector,
    const double timeN,
    const double timeNp1,
    const FieldRef coordsField,
    const FieldRef distField,
    const std::vector<String_Function_Expression> & interfaceVelocityExpr,
    const double narrowBandSize,
    const double avgEdgeLength,
    const FacetedSurfaceBase & facetsN,
    const FacetedSurfaceBase & facetsPred,
    FacetedSurfaceBase & facetsNp1)
{
  const int minDepth = 2;
  const int maxDepth = 5;

  const int dim = mesh.mesh_meta_data().spatial_dimension();
  for ( auto && bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, activeFieldSelector) )
  {
    const double * coordsData = field_data<double>(coordsField , *bucketPtr);
    double * distData = field_data<double>(distField , *bucketPtr);

    for (size_t i = 0; i < bucketPtr->size(); ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*dim, dim);
      const int previousSign = sign(distData[i]);
      distData[i] = compute_semilagrangian_distance_correction_at_point(dim, timeN, timeNp1, facetsN, facetsPred, nodeCoords, interfaceVelocityExpr, narrowBandSize, previousSign*narrowBandSize);
    }
  }

  const auto correct_distance_at_point = build_semilagrangian_distance_corrector_at_point(dim, timeN, timeNp1, interfaceVelocityExpr, facetsN, facetsPred);

  build_adaptive_facets_using_semilagrangian_distance(dim, mesh, activeFieldSelector, coordsField, distField, correct_distance_at_point, avgEdgeLength, minDepth, maxDepth, facetsNp1);
}

}


