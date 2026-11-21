#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OrientedSideNodes.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SegmentWithSensitivities.hpp>
#include <Akri_TriangleWithSensitivities.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_Unit_DecomposeWithSensitivities_TestUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/EnvData.hpp>

// Compares the shape sensitivities of the area (length in 2D) of the interface and the volume (area
// in 2D) of the negative phase of a circle or sphere with forward finite differences of the
// re-decomposed mesh, for a perturbation of the level set phi -> phi + h*delta over a sweep of h:
// a uniform shift (delta = 1), or a random perturbation of each node (what gradient checkers use).

namespace
{

struct AreaAndVolume
{
  double area{0.};
  double volume{0.};
};

// The volume enclosed by the oriented facets (divergence theorem) and its derivatives with respect
// to the facet node coordinates.
double facet_volume_and_sensitivities(const std::array<stk::math::Vector3d,3> & x, double * dVolume)
{
  const stk::math::Vector3d dV0 = Cross(x[1], x[2]) / 6.;
  const stk::math::Vector3d dV1 = Cross(x[2], x[0]) / 6.;
  const stk::math::Vector3d dV2 = Cross(x[0], x[1]) / 6.;
  for (unsigned d = 0; d < 3; ++d)
  {
    dVolume[d] = dV0[d];
    dVolume[3+d] = dV1[d];
    dVolume[6+d] = dV2[d];
  }
  return Dot(x[0], dV0);
}

double facet_volume_and_sensitivities(const std::array<stk::math::Vector3d,2> & x, double * dVolume)
{
  dVolume[0] = 0.5*x[1][1];
  dVolume[1] = -0.5*x[1][0];
  dVolume[2] = 0.;
  dVolume[3] = -0.5*x[0][1];
  dVolume[4] = 0.5*x[0][0];
  dVolume[5] = 0.;
  return 0.5*(x[0][0]*x[1][1] - x[0][1]*x[1][0]);
}

double facet_area_and_sensitivities(const std::array<stk::math::Vector3d,3> & x, double * dArea)
{
  return krino::TriangleWithSens::area_and_optional_sensitivities(x[0], x[1], x[2], dArea);
}

double facet_area_and_sensitivities(const std::array<stk::math::Vector3d,2> & x, double * dArea)
{
  return krino::SegmentWithSens::length_and_optional_sensitivities(x[0], x[1], dArea);
}

double compute_negative_phase_element_volume(const stk::mesh::BulkData & mesh, const krino::LS_Field & lsField)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const stk::mesh::Selector negativeOwnedActive = krino::Phase_Support::get(meta).get_negative_levelset_block_selector(lsField.identifier) &
      meta.locally_owned_part() & krino::AuxMetaData::get(meta).active_part();
  double volume = krino::compute_volume_or_surface_area(mesh, stk::topology::ELEMENT_RANK, negativeOwnedActive);
  krino::all_reduce_sum(mesh.parallel(), volume);
  return volume;
}

enum class Perturbation
{
  UNIFORM,
  RANDOM
};

// The perturbation direction at a node, a function of its id only, so the same on every processor.
double perturbation(const Perturbation perturbation, const stk::mesh::EntityId nodeId)
{
  if (perturbation == Perturbation::UNIFORM) return 1.;
  std::uint64_t x = nodeId * 0x9E3779B97F4A7C15ull;
  x ^= x >> 31;
  x *= 0xBF58476D1CE4E5B9ull;
  x ^= x >> 29;
  return 2. * (static_cast<double>(x >> 11) / static_cast<double>(1ull << 53)) - 1.;
}

void add_perturbation_to_levelset(const stk::mesh::BulkData & mesh,
    const krino::LS_Field & lsField,
    const Perturbation direction,
    const double h)
{
  auto lsData = lsField.isovar.data<double, stk::mesh::ReadWrite>();
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(lsField.isovar)))
  {
    auto ls = lsData.bucket_values(*bucketPtr);
    for (stk::mesh::EntityIdx n : bucketPtr->entities())
      ls(n) += h * perturbation(direction, mesh.identifier((*bucketPtr)[n]));
  }
}

void set_levelset_to_perturbed_sphere_distance(const stk::mesh::BulkData & mesh,
    const krino::LS_Field & lsField,
    const double radius,
    const Perturbation direction,
    const double h)
{
  const krino::FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  compute_nodal_distance_from_spheres(mesh, coordsField, lsField.isovar, {{stk::math::Vector3d::ZERO, radius}});
  add_perturbation_to_levelset(mesh, lsField, direction, h);
}

using LevelSetFunction = std::function<double(const stk::math::Vector3d &)>;

void set_levelset_to_perturbed_function(const stk::mesh::BulkData & mesh,
    const krino::LS_Field & lsField,
    const LevelSetFunction & levelSet,
    const Perturbation direction,
    const double h)
{
  const krino::FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  auto lsData = lsField.isovar.data<double, stk::mesh::ReadWrite>();
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(lsField.isovar)))
  {
    auto ls = lsData.bucket_values(*bucketPtr);
    for (stk::mesh::EntityIdx n : bucketPtr->entities())
      ls(n) = levelSet(get_vector_field(mesh, coordsField, (*bucketPtr)[n], dim));
  }
  add_perturbation_to_levelset(mesh, lsField, direction, h);
}

const char * method_name(const krino::ShapeSensitivityMethod method)
{
  switch (method)
  {
  case krino::ShapeSensitivityMethod::EDGE_DIRECTION_USING_EDGE_GRADIENT:
    return "edge direction using edge gradient";
  case krino::ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT:
    return "surface normal using edge gradient";
  }
  return "unknown";
}
}

constexpr auto EDGE_DIRECTION = krino::ShapeSensitivityMethod::EDGE_DIRECTION_USING_EDGE_GRADIENT;
constexpr auto SURFACE_NORMAL = krino::ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT;

// A method to check and its tolerances on the area (length in 2D) and volume (area in 2D) at the
// checked step.  A zero tolerance means the exact tolerance, for the edge direction form and the
// surface normal form's volume.
struct MethodCheck
{
  krino::ShapeSensitivityMethod method;
  AreaAndVolume tol{};
};

class DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences : public ::testing::Test
{
protected:
  // The area and volume of the interface facets, and the directional derivatives of each for the
  // perturbation from the method's shape sensitivities.
  template <unsigned DIM>
  void compute_area_and_volume_and_sensitivities(
      const krino::ShapeSensitivityMethod method, AreaAndVolume & value, AreaAndVolume & sens)
  {
    constexpr unsigned numFacetNodes = DIM;
    const stk::mesh::BulkData & mesh = myDecomposedMesh.mesh->bulk_data();

    const krino::InterfaceFacetsWithSensitivities<numFacetNodes> facets =
        myDecomposedMesh.decomposition->facets_and_sensitivities<numFacetNodes>(0, method);

    value = AreaAndVolume{};
    sens = AreaAndVolume{};
    double dArea[3*numFacetNodes];
    double dVolume[3*numFacetNodes];
    for (size_t facet = 0; facet < facets.num_facets(); ++facet)
    {
      const std::array<const krino::LevelSetShapeSensitivity *, numFacetNodes>
          facetNodeSensitivities = facets.node_sensitivities(facet);
      const std::array<stk::mesh::Entity, numFacetNodes> facetNodes = facets.nodes(mesh, facet);
      const std::array<stk::math::Vector3d,numFacetNodes> facetCoords = krino::get_nodes_vector(mesh, mesh.mesh_meta_data().coordinate_field(), facetNodes, DIM);

      value.area += facet_area_and_sensitivities(facetCoords, dArea);
      value.volume += facet_volume_and_sensitivities(facetCoords, dVolume);

      for (unsigned n=0; n<numFacetNodes; ++n)
      {
        // The motion of the node for the perturbation of its parents' level set values.
        stk::math::Vector3d dX_dShift = stk::math::Vector3d::ZERO;
        const krino::LevelSetShapeSensitivity & nodeSens = *facetNodeSensitivities[n];
        for (size_t i = 0; i < nodeSens.parentNodeIds.size(); ++i)
          dX_dShift += perturbation(myPerturbation, nodeSens.parentNodeIds[i]) *
              nodeSens.dCoordsdParentLevelSets[i];
        sens.area += Dot(stk::math::Vector3d(dArea+3*n), dX_dShift);
        sens.volume += Dot(stk::math::Vector3d(dVolume+3*n), dX_dShift);
      }
    }

    krino::all_reduce_sum(mesh.parallel(), value.area, value.volume, sens.area, sens.volume);
  }

  // The area and volume only, for the perturbed decompositions: the same facets, in the same order,
  // as facets_and_sensitivities, without computing the sensitivities.
  template <unsigned DIM> AreaAndVolume compute_area_and_volume()
  {
    const stk::mesh::BulkData & mesh = myDecomposedMesh.mesh->bulk_data();
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    const krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);
    const stk::mesh::Selector interfaceSelector =
        phaseSupport.get_negative_levelset_interface_selector(levelset_field().identifier);
    const stk::mesh::Selector negativeSideElementSelector =
        phaseSupport.get_negative_levelset_block_selector(levelset_field().identifier);

    AreaAndVolume value;
    double dArea[3 * DIM];
    double dVolume[3 * DIM];
    for (auto * bucketPtr : mesh.get_buckets(meta.side_rank(),
             meta.locally_owned_part() & krino::AuxMetaData::get(meta).active_part() &
                 interfaceSelector))
    {
      for (auto side : *bucketPtr)
      {
        std::array<stk::mesh::Entity, DIM> facetNodes;
        if constexpr (DIM == 2)
          facetNodes = krino::get_oriented_line_side_nodes(mesh, negativeSideElementSelector, side);
        else
          facetNodes =
              krino::get_oriented_triangle_side_nodes(mesh, negativeSideElementSelector, side);
        const std::array<stk::math::Vector3d, DIM> facetCoords =
            krino::get_nodes_vector(mesh, meta.coordinate_field(), facetNodes, DIM);
        value.area += facet_area_and_sensitivities(facetCoords, dArea);
        value.volume += facet_volume_and_sensitivities(facetCoords, dVolume);
      }
    }

    krino::all_reduce_sum(mesh.parallel(), value.area, value.volume);
    return value;
  }

  // Removes the domain wall side parts, so that snapping (and the sensitivities) see only the
  // exposed boundary.
  void remove_domain_side_parts()
  {
    stk::mesh::BulkData & mesh = myDecomposedMesh.mesh->bulk_data();
    const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::PartVector sideParts;
    for (const char * partName : {"Xminus", "Xplus", "Yminus", "Yplus", "Zminus", "Zplus"})
      if (stk::mesh::Part * part = meta.get_part(partName)) sideParts.push_back(part);
    std::vector<stk::mesh::Entity> sides;
    for (auto * bucketPtr : mesh.get_buckets(
             meta.side_rank(), stk::mesh::selectUnion(sideParts) & meta.locally_owned_part()))
      sides.insert(sides.end(), bucketPtr->begin(), bucketPtr->end());
    mesh.batch_change_entity_parts(sides,
        std::vector<stk::mesh::PartVector>(sides.size()),
        std::vector<stk::mesh::PartVector>(sides.size(), sideParts));
  }

  template <unsigned DIM>
  AreaAndVolume decompose_with_shift_and_compute_area_and_volume(const double shift)
  {
    stk::mesh::BulkData & mesh = myDecomposedMesh.mesh->bulk_data();
    krino::CDMesh::reset_mesh_to_original_undecomposed_state(mesh);
    if (myDoRemoveDomainSideParts)
    {
      remove_domain_side_parts();
      myDoRemoveDomainSideParts = false;
    }
    if (myLevelSet)
      set_levelset_to_perturbed_function(mesh, levelset_field(), myLevelSet, myPerturbation, shift);
    else
      set_levelset_to_perturbed_sphere_distance(
          mesh, levelset_field(), myRadius, myPerturbation, shift);
    myDecomposedMesh.decomposition->decompose();
    return compute_area_and_volume<DIM>();
  }

  // One sweep of forward differences, compared with the sensitivities of each method on the same
  // decomposition (it does not depend on the method).  The checked step is h = 1e-7 for a uniform
  // shift and 1e-8 for a random direction, whose small directional derivatives (much of the motion
  // cancels) give a larger relative O(h) truncation error.
  template <unsigned DIM>
  void sweep_finite_difference_step_and_compare_with_sensitivities(const double meshSize,
      const bool doSnapping,
      const Perturbation direction,
      const std::vector<MethodCheck> & checks)
  {
    myMethods.clear();
    for (auto & check : checks)
      myMethods.push_back(check.method);
    myPerturbation = direction;
    const bool doWriteMesh = false;
    myDecomposedMesh = build_circle_or_sphere_conforming_bounding_box_mesh(DIM,
        {-1., -1., -1.},
        {1., 1., 1.},
        meshSize,
        {0., 0., 0.},
        myRadius,
        doSnapping,
        doWriteMesh);
    krino::CDFEM_Support::get(myDecomposedMesh.mesh->meta_data())
        .set_snapping_sharp_feature_angle_in_degrees(mySharpFeatureAngleInDegrees);

    const AreaAndVolume value = myLevelSet ?
        decompose_with_shift_and_compute_area_and_volume<DIM>(0.) :
        compute_area_and_volume<DIM>();
    std::vector<AreaAndVolume> sens(checks.size());
    for (size_t m = 0; m < checks.size(); ++m)
    {
      AreaAndVolume methodValue;
      compute_area_and_volume_and_sensitivities<DIM>(checks[m].method, methodValue, sens[m]);
      EXPECT_EQ(value.area, methodValue.area);
      EXPECT_EQ(value.volume, methodValue.volume);
    }

    // The facets are oriented by the negative-phase element, so they enclose the negative phase.
    if (myIsInterfaceClosed)
    {
      const double elementVolume = compute_negative_phase_element_volume(
          myDecomposedMesh.mesh->bulk_data(), levelset_field());
      EXPECT_NEAR(elementVolume, value.volume, 1.e-10 * std::abs(elementVolume));
    }

    // The surface normal keeps each node's motion along the derivative of the volume, so the volume
    // sensitivity is the edge direction one.
    const auto edgeIter = std::find(myMethods.begin(), myMethods.end(), EDGE_DIRECTION);
    const auto surfaceIter = std::find(myMethods.begin(), myMethods.end(), SURFACE_NORMAL);
    if (edgeIter != myMethods.end() && surfaceIter != myMethods.end())
    {
      const double edgeVolumeSens = sens[edgeIter - myMethods.begin()].volume;
      EXPECT_NEAR(edgeVolumeSens,
          sens[surfaceIter - myMethods.begin()].volume,
          1.e-12 * std::abs(edgeVolumeSens));
    }

    const double pi = M_PI;
    const AreaAndVolume analyticSens = (2 == DIM) ?
        AreaAndVolume{-2. * pi, -2. * pi * myRadius} :
        AreaAndVolume{-8. * pi * myRadius, -4. * pi * myRadius * myRadius};

    const bool isRoot = 0 == stk::EnvData::parallel_rank();
    if (isRoot)
    {
      std::printf("FD sweep: dim %u, mesh size %g, %s%s\n",
          DIM,
          meshSize,
          doSnapping ? "snapping" : "no snapping",
          direction == Perturbation::RANDOM ? ", random direction" : "");
      std::printf("  area %.12g, volume %.12g\n", value.area, value.volume);
      for (size_t m = 0; m < checks.size(); ++m)
        std::printf("  sens %-34s dArea %.10g, dVolume %.10g\n",
            method_name(checks[m].method),
            sens[m].area,
            sens[m].volume);
      if (direction == Perturbation::UNIFORM && !myLevelSet)
        std::printf("  sens %-34s dArea %.10g, dVolume %.10g\n",
            "analytic",
            analyticSens.area,
            analyticSens.volume);
      std::printf("  %8s %18s %18s", "h", "FD dArea", "FD dVolume");
      for (size_t m = 0; m < checks.size(); ++m)
        std::printf("   relErrArea[%zu]  relErrVol[%zu]", m, m);
      std::printf("\n");
    }

    for (int p = 1; p <= 12; ++p)
    {
      const double h = std::pow(10., -p);
      const AreaAndVolume plus = decompose_with_shift_and_compute_area_and_volume<DIM>(h);
      const AreaAndVolume fd{(plus.area - value.area) / h, (plus.volume - value.volume) / h};
      EXPECT_TRUE(std::isfinite(fd.area) && std::isfinite(fd.volume));
      if (isRoot) std::printf("  %8.0e %18.10g %18.10g", h, fd.area, fd.volume);

      const int checkedPower = (direction == Perturbation::UNIFORM) ? 7 : 8;
      for (size_t m = 0; m < checks.size(); ++m)
      {
        const double relErrArea = std::abs(fd.area - sens[m].area) / std::abs(sens[m].area);
        const double relErrVolume = std::abs(fd.volume - sens[m].volume) / std::abs(sens[m].volume);
        if (isRoot) std::printf("  %13.3e  %13.3e", relErrArea, relErrVolume);

        if (checkedPower == p)
        {
          // For an exact derivative the error is the forward difference's O(h) truncation and
          // roundoff (observed up to 1.4e-7 at h=1e-7 for a uniform shift; 5.1e-5 at h=1e-8 for
          // the random direction on the sphere, whose directional derivative mostly cancels).
          const krino::ShapeSensitivityMethod method = checks[m].method;
          const bool isExactArea = method == EDGE_DIRECTION;
          const bool isExactVolume = isExactArea || method == SURFACE_NORMAL;
          const double exactTol = (direction == Perturbation::UNIFORM) ? 1.e-6 : 2.e-4;
          const AreaAndVolume & tol = checks[m].tol;
          EXPECT_LT(relErrArea, (isExactArea && tol.area == 0.) ? exactTol : tol.area)
              << method_name(method);
          if (myIsInterfaceClosed)
          {
            EXPECT_LT(relErrVolume, (isExactVolume && tol.volume == 0.) ? exactTol : tol.volume)
                << method_name(method);
          }
        }
      }
      if (isRoot) std::printf("\n");
    }
  }

  const krino::LS_Field & levelset_field() const
  {
    return myDecomposedMesh.decomposition->levelset_fields()[0];
  }

  double myRadius = 0.76;
  LevelSetFunction myLevelSet; // the sphere distance of myRadius if empty
  bool myIsInterfaceClosed = true; // otherwise only the area is checked
  bool myDoRemoveDomainSideParts = false; // before the first decomposition of myLevelSet
  double mySharpFeatureAngleInDegrees = 0.;
  std::vector<krino::ShapeSensitivityMethod> myMethods;
  Perturbation myPerturbation = Perturbation::UNIFORM;
  DecomposedBoundingBoxMesh myDecomposedMesh;
};

// Surface normal area tolerances, about 3 times the errors observed at the checked step.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, circle)
{
  sweep_finite_difference_step_and_compare_with_sensitivities<2>(
      0.1, false, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {5.e-4, 0.}}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, circle_snapping)
{
  sweep_finite_difference_step_and_compare_with_sensitivities<2>(
      0.1, true, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {5.e-3, 0.}}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere)
{
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {1.e-2, 0.}}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_snapping)
{
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., true, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {5.e-3, 0.}}});
}

// The uniform shift above moves every parent's level set by the same amount, so it checks only each
// node's total motion, the sum of its dCoordsdParentLevelSets, and per-node errors can cancel over
// the interface.  A random direction perturbs each parent independently (delta in [-1,1], a hash of
// the node id), so it checks each parent's contribution separately, as a gradient checker comparing
// a directional derivative with a finite difference does.  The perturbed shape is no longer a
// circle or sphere, and much of the motion cancels, so the directional derivatives are small and
// the relative finite difference error larger.
//
// Random directions are checked without snapping only.  With snapping, a snapped node's position
// can jump by most of an edge for a perturbation of 1e-6 (observed 2.4e-2 in 2D, 7.1e-2 in 3D, with
// the set of snapped nodes unchanged), so the snapped objective is discontinuous in the level set
// close to these designs, and a finite difference across a jump agrees with no sensitivity.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, circle_randomDirection)
{
  sweep_finite_difference_step_and_compare_with_sensitivities<2>(
      0.1, false, Perturbation::RANDOM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {5.e-2, 0.}}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_randomDirection)
{
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::RANDOM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {1.e-1, 0.}}});
}

// Designs whose interface passes through background nodes, which have zero level set and no cut
// edge.  Without snapping, the edge direction form is exact there for the volume but not the area
// (see ShapeSensitivityMethod; observed 3.0e-3 for the uniform shift, 1.9e-3 for the random
// direction), so the area is checked with a tolerance.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_throughNodes)
{
  myRadius = std::sqrt(2. / 3.);
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(1. / 3.,
      false,
      Perturbation::UNIFORM,
      {{EDGE_DIRECTION, {1.e-2, 0.}}, {SURFACE_NORMAL, {3.e-2, 0.}}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_throughOtherNodes)
{
  myRadius = std::sqrt(5.) / 3.;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(1. / 3.,
      false,
      Perturbation::UNIFORM,
      {{EDGE_DIRECTION, {1.e-3, 0.}}, {SURFACE_NORMAL, {3.e-3, 0.}}});
}

// Through cell-centre nodes of the mesh, each with a neighbouring cell-centre node (one mesh edge
// away) also on the interface, so the interface lies along those mesh edges and its area is not
// differentiable there (observed forward and backward derivatives -0.2375 and -0.2433 at one node).
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_throughNeighboringNodes)
{
  myRadius = std::sqrt(2.75) / 3.;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(1. / 3.,
      false,
      Perturbation::UNIFORM,
      {{EDGE_DIRECTION, {3.e-2, 0.}}, {SURFACE_NORMAL, {5.e-2, 0.}}});
}

// With snapping, predictable snap selection moves a background node on the interface along one
// edge as the level set increases, so the edge direction form is exact, for the area too.  This is
// the forward derivative: decreasing the level set moves it along another edge.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_throughNodes_snapping)
{
  myRadius = std::sqrt(2. / 3.);
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., true, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {2.e-2, 0.}}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_throughFewerNodes_snapping)
{
  myRadius = 2. / 3.;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., true, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {5.e-3, 0.}}});
}

TEST_F(
    DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, sphere_throughNodes_randomDirection)
{
  myRadius = 2. / 3.;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(1. / 3.,
      false,
      Perturbation::RANDOM,
      {{EDGE_DIRECTION, {1.e-2, 0.}}, {SURFACE_NORMAL, {1.e-1, 0.}}});
}

// A rhombic dodecahedron, max(|x|+|y|, |y|+|z|, |x|+|z|) = 2/3, whose faces lie on faces of the
// mesh tetrahedra, so every interface side has all its nodes on zero level set and no cut edge. The
// level set is scaled by 0.5 inside, so its gradient jumps across the interface and the volume is
// not differentiable: the interface moves at the inverse of the inside gradient for increasing
// level set and of the outside one for decreasing.  The sensitivities give the derivative for
// increasing level set, which the forward differences of the uniform shift measure.  The faces are
// planar, so the area is exact too.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    rhombicDodecahedronOnMeshFaces_gradientJump)
{
  myLevelSet = [](const stk::math::Vector3d & x)
  {
    const double ax = std::abs(x[0]), ay = std::abs(x[1]), az = std::abs(x[2]);
    const double phi = std::max({ax + ay, ay + az, ax + az}) - 2. / 3.;
    return (std::abs(phi) < 1.e-12) ? 0. : ((phi < 0.) ? 0.5 * phi : phi);
  };
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::UNIFORM, {{EDGE_DIRECTION}, {SURFACE_NORMAL, {1.e-6, 0.}}});
}

// A plane through the origin and mesh nodes, several on the domain walls, at an angle to the walls.
// Snapping moves a wall node only along edges on that wall, so the sensitivities must too.  The
// interior nodes contribute nothing to the area derivative, so it all comes from the walls.  Random
// direction: the uniform shift's area derivative is zero.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    planeThroughNodesAtAngleToWalls_snapping_randomDirection)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 2., 4.).unit_vector();
  myLevelSet = [normal](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., true, Perturbation::RANDOM, {{EDGE_DIRECTION}});
}

// Without side parts on the domain walls, snapping constrains nodes on the domain's sharp edges
// (with the parser's default sharp feature angle) to move along those edges, so the sensitivities
// must too.  The plane passes through nodes on three box edges near a corner.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    planeThroughNodesOnBoxEdgesWithoutSideParts_snappingWithSharpFeatures)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 1., 1.).unit_vector();
  const stk::math::Vector3d pointOnPlane(1., 1., 2. / 3.);
  myLevelSet = [normal, pointOnPlane](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - pointOnPlane);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  myDoRemoveDomainSideParts = true;
  mySharpFeatureAngleInDegrees = 135.;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., true, Perturbation::UNIFORM, {{EDGE_DIRECTION}});
}

// Without snapping, a background node with zero level set on a domain wall moves within the wall
// (projected from the interface normal, keeping its normal speed), as snapping would let it move.
// For a planar interface on flat walls, all the node's sides move along the line where the plane
// meets the wall, so this is exact.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    planeThroughNodesAtAngleToWalls_randomDirection)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 2., 4.).unit_vector();
  myLevelSet = [normal](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::RANDOM, {{EDGE_DIRECTION}});
}

// A node on a box edge (two wall side parts) moves along the edge.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, planeThroughNodesOnBoxEdges)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 1., 1.).unit_vector();
  const stk::math::Vector3d pointOnPlane(1., 1., 2. / 3.);
  myLevelSet = [normal, pointOnPlane](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - pointOnPlane);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::UNIFORM, {{EDGE_DIRECTION}});
}

// Without side parts, the sharp edges of the domain constrain the nodes on them instead.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    planeThroughNodesOnBoxEdgesWithoutSidePartsWithSharpFeatures)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 1., 1.).unit_vector();
  const stk::math::Vector3d pointOnPlane(1., 1., 2. / 3.);
  myLevelSet = [normal, pointOnPlane](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - pointOnPlane);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  myDoRemoveDomainSideParts = true;
  mySharpFeatureAngleInDegrees = 135.;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::UNIFORM, {{EDGE_DIRECTION}});
}

// A plane through a corner of the box, meeting it only there: the corner moves along the box edge
// where the two walls the plane crosses meet (the edge whose neighbor's level set sign differs from
// the other two), even when that points out of the domain.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, planeThroughBoxCorner)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 1., -1.).unit_vector();
  const stk::math::Vector3d corner(1., 1., 1.);
  myLevelSet = [normal, corner](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - corner);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::UNIFORM, {{EDGE_DIRECTION}});
}

TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    planeThroughBoxCornerWithoutSidePartsWithSharpFeatures)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(3., -1., -1.).unit_vector();
  const stk::math::Vector3d corner(1., 1., 1.);
  myLevelSet = [normal, corner](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - corner);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  myDoRemoveDomainSideParts = true;
  mySharpFeatureAngleInDegrees = 135.;
  // random direction: the uniform shift's area derivative is zero for this plane
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::RANDOM, {{EDGE_DIRECTION}});
}

// A plane containing a box edge: the area is not differentiable there (the interface leaves the
// edge into one wall or the other), and the nodes on the edge move into the wall on the negative
// side, giving the derivative for increasing level set that the forward differences measure.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, planeContainingBoxEdge)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., -2., 0.).unit_vector();
  const stk::math::Vector3d pointOnEdge(1., 1., 0.);
  myLevelSet = [normal, pointOnEdge](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - pointOnEdge);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<3>(
      1. / 3., false, Perturbation::UNIFORM, {{EDGE_DIRECTION}});
}

// 2D: a line through a corner of the square moves along the wall with the negative neighbor.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences, lineThroughSquareCorner)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., -2., 0.).unit_vector();
  const stk::math::Vector3d corner(1., 1., 0.);
  myLevelSet = [normal, corner](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - corner);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<2>(
      0.25, false, Perturbation::UNIFORM, {{EDGE_DIRECTION}});
}

// 2D: a line through a node on a wall, at an angle to it.
TEST_F(DecomposeMeshAndCompareSensitivitiesWithFiniteDifferences,
    lineThroughWallNodeAtAngle_randomDirection)
{
  const stk::math::Vector3d normal = stk::math::Vector3d(1., 2., 0.).unit_vector();
  const stk::math::Vector3d wallNode(1., -1. / 9., 0.); // the lattice rows are at y = -1 + 2k/9
  myLevelSet = [normal, wallNode](const stk::math::Vector3d & x)
  {
    const double phi = Dot(normal, x - wallNode);
    return (std::abs(phi) < 1.e-12) ? 0. : phi;
  };
  myIsInterfaceClosed = false;
  sweep_finite_difference_step_and_compare_with_sensitivities<2>(
      0.25, false, Perturbation::RANDOM, {{EDGE_DIRECTION}});
}

static std::map<stk::mesh::EntityId, stk::math::Vector3d> get_node_coordinates_by_id(const stk::mesh::BulkData & mesh)
{
  const krino::FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::map<stk::mesh::EntityId, stk::math::Vector3d> coordsById;
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, mesh.mesh_meta_data().locally_owned_part()))
    for (auto node : *bucketPtr)
      coordsById[mesh.identifier(node)] = get_vector_field(mesh, coordsField, node, dim);
  return coordsById;
}

TEST(ResetSnappedDecomposition, resetMeshToOriginalUndecomposedState_restoresCoordinatesOfSnappedNodes)
{
  std::vector<krino::LS_Field> lsFields;
  std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(stk::topology::TRIANGLE_3_2D, MPI_COMM_WORLD);
  lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh->meta_data(), 1);
  const bool doSnapping = true;
  krino::ConformingDecompositionWithSensitivities decomposition(
      bboxMesh->meta_data(), lsFields, doSnapping);
  krino::populate_bounding_box_mesh_and_activate(*bboxMesh, {-1.,-1.,0.}, {1.,1.,0.}, 0.05);
  stk::mesh::BulkData & mesh = bboxMesh->bulk_data();

  const std::map<stk::mesh::EntityId, stk::math::Vector3d> origCoords = get_node_coordinates_by_id(mesh);

  for (int iter = 0; iter < 2; ++iter)
  {
    set_levelset_to_perturbed_sphere_distance(mesh, lsFields[0], 0.76, Perturbation::UNIFORM, 0.);
    decomposition.decompose();
    krino::CDMesh::reset_mesh_to_original_undecomposed_state(mesh);

    const std::map<stk::mesh::EntityId, stk::math::Vector3d> resetCoords = get_node_coordinates_by_id(mesh);
    EXPECT_EQ(origCoords.size(), resetCoords.size());
    double maxDisplacement = 0.;
    for (auto & [id, coords] : origCoords)
    {
      const auto iterResetCoords = resetCoords.find(id);
      ASSERT_TRUE(iterResetCoords != resetCoords.end());
      maxDisplacement = std::max(maxDisplacement, (iterResetCoords->second - coords).length());
    }
    EXPECT_LT(maxDisplacement, 1.e-14) << "on reset " << iter;
  }
}
