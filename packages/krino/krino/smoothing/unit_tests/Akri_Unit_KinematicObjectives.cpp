// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <MiniTensor_Tensor.h>
#include <Akri_MeshHelpers.hpp>
#include <Akri_KinematicElement.hpp>
#include <Akri_KinematicUtils.hpp>
#include <Akri_MetricTensorUtils.hpp>
#include <Akri_ModifiedNeoHookeanModel.hpp>
#include <Akri_Objective_KinematicModel.hpp>
#include <Akri_SethHillConfig.hpp>
#include <Akri_SethHillModel.hpp>
#include <Akri_SizeShapeModel.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_UnitTestUtils.hpp>

namespace krino {

namespace KinematicObjectiveGoldData {
constexpr std::array<double, 12> reference_coords
{{
  -2.500000000000000e-001, -1.443375672974064e-001, 0.000000000000000e+000,
  2.500000000000000e-001, -1.443375672974064e-001, 0.000000000000000e+000,
  0.000000000000000e+000, 2.886751345948129e-001, 0.000000000000000e+000,
  0.000000000000000e+000, 0.000000000000000e+000, 4.082482904638630e-001
}};

constexpr std::array<double, 12> current_coords
{{
  -2.952397829717832e-001, -3.236702717044001e-001, 6.686989148589159e-002,
  5.827602170282168e-001, 1.233297282955999e-001, -2.106301085141084e-001,
  -2.875204340564336e-001, 2.003405434088001e-001, 1.437602170282168e-001,
  7.711810206862372e-001, -1.318641978198277e-001, 1.291289342737199e+000
}};

constexpr std::array<double, 9> deformation_gradient
{{
  1.756000000000000e+000, -9.960000000000000e-001, 1.889000000000000e+000,
  8.940000000000000e-001, 6.939999999999999e-001, -3.230000000000000e-001,
  -5.550000000000001e-001, 4.980000000000000e-001, 3.163000000000000e+000
}};

constexpr double element_size = 5.000000000000000e-001;

SethHillParams SH_params{1.0, 1.0, 1, 1};

constexpr double integrated_SH_total_energy = 2.548065224351188e-001;

constexpr std::array<double, 12> SH_sensitivities
{{
  -2.015813702621814e-001, -5.871391833167285e-001, -9.778777789019070e-002,
  3.820905496798706e-001, 8.001783948079290e-002, -3.575788562005628e-001,
  -3.386234721220264e-001, 5.548127054195153e-001, 1.040208235140742e-001,
  1.581142927043372e-001, -4.769136158357974e-002, 3.513458105766793e-001
}};

constexpr double integrated_eta_total_energy = 0.046446210231283848;

constexpr std::array<double, 12> eta_sensitivities
{{
  -0.02406292951103901, -0.012373725517176337, -0.012464090987250461,
  0.018503582957359491, 0.0076960346961769469, -0.025782156428927314,
  -0.02305938377694287, 0.0095102578970487767, -0.0094301199684372304,
  0.028618730330622393, -0.0048325670760493841, 0.047676367384615011
}};
}

namespace KinematicObjectiveAnisotropicGoldData {
constexpr std::array<double, 9> M
{{
    3.406314606115276e+000, 4.084888615312371e-001, -1.039800361583855e+000,
    4.084888615312371e-001, 2.117208513999345e+000, -3.452311194608259e-001,
    -1.039800361583855e+000, -3.452311194608259e-001, 1.476476879885379e+000
}};

// The metric mapping F_M of M is not unique, so the gold data is given in terms of quantities that do
// not depend on the choice: |det(F_M)| = sqrt(det(M))
constexpr double abs_det_F_M = 2.828427124746190e+000;

constexpr std::array<double, 12> current_coordinates
{{
    4.997793476506873e-001, 3.435655182579872e-001, 8.279808348862010e-001,
    1.081618379051203e+000, -2.459267253890746e-001, 1.128829398816410e+000,
    -1.268402805984682e+000, -1.419395994848755e-001, -5.576408745583139e-001,
    -3.129949207172080e-001, 4.430080661596296e-002, -1.399169359144297e+000
}};

// Fu is the unit-space deformation gradient for one choice of F_M; any other choice rotates it.
constexpr std::array<double, 9> Fu
{{
    1.396568207719426e+000, -2.064229515591467e+000, 1.799473067309108e+000,
    6.636422011094258e-001, 8.290479128384540e-001, -8.888435323847280e-001,
    -5.615682518364306e-001, -1.078692701128013e-002, 3.387383879442120e+000
}};

SethHillParams SH_params{1.0, 1.0, 1, 1};

constexpr double integrated_SH_total_energy = 8.782105076342680e-001;

constexpr std::array<double, 12> SH_sensitivities
{{
  6.918983688689938e-002, 7.590994464568388e-001, 2.778284382738434e-001,
  5.127780869200135e-001, -6.250067689558163e-001, 1.935189040964040e-001,
  -1.161471570043401e+000, -3.822165824401763e-001, 4.364359398792074e-001,
  5.795036462364882e-001, 2.481239049391538e-001, -9.077832822494550e-001
}};
}

namespace KinematicObjectiveGoldData2d {
constexpr std::array<double, 6> reference_coords
{{
   0.0, 0.86602540378444,
  -0.5, 0.0,
   0.5, 0.0
}};

constexpr std::array<double, 6> current_coords
{{
   0.0, 0.0,
   1.0, 0.0,
   0.0, 1.0
}};

constexpr std::array<double, 4> deformation_gradient
{{
   -1, -0.57735026918962484,
    1, -0.57735026918962484
}};

constexpr double integrated_eta_total_energy = 0.072168783648703244;
}

namespace KinematicObjectiveAnisotropicGoldData2d {
constexpr std::array<double, 4> M
{{
    2.658795048015561e+000, 8.070256022619632e-001,
    8.070256022619632e-001, 1.091204951984439e+000
}};

// The metric mapping F_M of M is not unique, so the gold data is given in terms of quantities that do
// not depend on the choice: |det(F_M)| = sqrt(det(M))
constexpr double abs_det_F_M = 1.500000000000000e+000;

constexpr std::array<double, 6> current_coordinates
{{
    4.587101103778379e-001, 4.436626516400468e-001,
    -1.587518717905141e-001, 5.632034085810095e-001,
    5.587518717905142e-001, -7.632034085810094e-001
}};

// Fu is the unit-space deformation gradient for one choice of F_M; any other choice rotates it.
constexpr std::array<double, 4> Fu
{{
    1.300000000000000e+000, -4.000000000000000e-001,
    2.500000000000000e-001, 9.000000000000000e-001
}};

SethHillParams SH_params{1.0, 1.0, 1, 1};

constexpr double integrated_SH_total_energy = 5.298639836327011e-002;

constexpr std::array<double, 6> SH_sensitivities
{{
    -3.660672249864304e-001, -6.449460477699333e-002,
    2.178800503408262e-001, 2.599919607854630e-001,
    1.481871746456042e-001, -1.954973560084697e-001
}};
}

static constexpr double tol = 1.e-6;
using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

static std::array<stk::math::Vector3d,4> tet_node_coords(const std::array<double, 12> & x)
{
  return {{ stk::math::Vector3d(&x[0]), stk::math::Vector3d(&x[3]), stk::math::Vector3d(&x[6]), stk::math::Vector3d(&x[9]) }};
}
static std::array<stk::math::Vector2d,3> tri_node_coords(const std::array<double, 6> & x)
{
  return {{ stk::math::Vector2d(&x[0]), stk::math::Vector2d(&x[2]), stk::math::Vector2d(&x[4]) }};
}
template <size_t NDOF>
static void expect_sensitivities(const std::array<double, NDOF> & gold, const std::array<double, NDOF> & sens)
{
  for (unsigned i=0; i<NDOF; ++i)
  {
    EXPECT_NEAR(gold[i], sens[i], tol);
  }
}
template <unsigned DIM>
static void expect_tensor(const minitensor::Tensor<double, DIM> & gold, const minitensor::Tensor<double, DIM> & A)
{
  for (unsigned i=0; i<DIM; ++i)
  {
    for (unsigned j=0; j<DIM; ++j)
    {
      EXPECT_NEAR(gold(i,j), A(i,j), tol);
    }
  }
}
template <unsigned DIM>
static minitensor::Tensor<double, DIM> tensor_from_array(const std::array<double, DIM*DIM> & a)
{
  minitensor::Tensor<double, DIM> A;
  for (unsigned i=0; i<DIM; ++i)
    for (unsigned j=0; j<DIM; ++j)
      A(i,j) = a[i*DIM+j];
  return A;
}

TEST(KinematicObjective, testEnergyAndSensitivities)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  const auto & refCoords = tet_node_coords(KinematicObjectiveGoldData::reference_coords);

  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refCoords, currentCoords);
  expect_tensor(tensor_from_array<3>(KinematicObjectiveGoldData::deformation_gradient), F);

  const double sethHillEnergy = tet4_element_energy(SethHillModel(KinematicObjectiveGoldData::SH_params), refCoords, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_SH_total_energy, sethHillEnergy, tol);
  const auto sethHillForces = tet4_element_forces(SethHillModel(KinematicObjectiveGoldData::SH_params), refCoords, currentCoords);
  expect_sensitivities(KinematicObjectiveGoldData::SH_sensitivities, sethHillForces);

  const double sizeShapePseudoEnergy = tet4_element_energy(SizeShapeModel(), refCoords, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_eta_total_energy, sizeShapePseudoEnergy, tol);
  const auto sizeShapePseudoForces = tet4_element_forces(SizeShapeModel(), refCoords, currentCoords);
  expect_sensitivities(KinematicObjectiveGoldData::eta_sensitivities, sizeShapePseudoForces);
}

TEST(KinematicObjective, testEnergyAndSensitivitiesForConstantH)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);

  const Tensor F_refSize = kinematicUtils::compute_deformation_gradient_tet4(KinematicObjectiveGoldData::element_size, currentCoords);
  const Tensor F_ref     = kinematicUtils::compute_deformation_gradient_tet4(kinematicUtils::compute_ideal_tet4_reference(KinematicObjectiveGoldData::element_size), currentCoords);
  expect_tensor(F_ref, F_refSize);

  const double sethHillEnergy = tet4_element_energy(SethHillModel(KinematicObjectiveGoldData::SH_params), KinematicObjectiveGoldData::element_size, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_SH_total_energy, sethHillEnergy, tol);
  const auto sethHillForces = tet4_element_forces(SethHillModel(KinematicObjectiveGoldData::SH_params), KinematicObjectiveGoldData::element_size, currentCoords);
  expect_sensitivities(KinematicObjectiveGoldData::SH_sensitivities, sethHillForces);

  const double sizeShapePseudoEnergy = tet4_element_energy(SizeShapeModel(), KinematicObjectiveGoldData::element_size, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_eta_total_energy, sizeShapePseudoEnergy, tol);
  const auto sizeShapePseudoForces = tet4_element_forces(SizeShapeModel(), KinematicObjectiveGoldData::element_size, currentCoords);
  expect_sensitivities(KinematicObjectiveGoldData::eta_sensitivities, sizeShapePseudoForces);
}

TEST(KinematicObjective, testEnergyAndSensitivitiesFor2DConstantH)
{
  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  const auto & refCoords = tri_node_coords(KinematicObjectiveGoldData2d::reference_coords);
  const double refSize = 1.;

  const Tensor2d F_refSize = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, currentCoords);
  const Tensor2d F_ref     = kinematicUtils::compute_deformation_gradient_tri3_2d(kinematicUtils::compute_ideal_tri3_2d_reference(refSize), currentCoords);
  expect_tensor(F_ref, F_refSize);

  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refCoords, currentCoords);
  expect_tensor(tensor_from_array<2>(KinematicObjectiveGoldData2d::deformation_gradient), F);

  const double sizeShapePseudoEnergy = tri3_2d_element_energy(SizeShapeModel(), refSize, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData2d::integrated_eta_total_energy, sizeShapePseudoEnergy, tol);

  const auto sizeShapePseudoForcesFromSize = tri3_2d_element_forces(SizeShapeModel(), refSize, currentCoords);
  const auto sizeShapePseudoForcesFromRefCoords = tri3_2d_element_forces(SizeShapeModel(), refCoords, currentCoords);
  expect_sensitivities(sizeShapePseudoForcesFromSize, sizeShapePseudoForcesFromRefCoords);

  const SethHillParams params{1.0, 10.0, 3, 3};
  const auto sethHillForcesFromSize = tri3_2d_element_forces(SethHillModel(params), refSize, currentCoords);
  const auto sethHillForcesFromRefCoords = tri3_2d_element_forces(SethHillModel(params), refCoords, currentCoords);
  expect_sensitivities(sethHillForcesFromSize, sethHillForcesFromRefCoords);
}

template <unsigned DIM, size_t NPE, class ObjFn, class SensFn>
static void test_sensitivity_matches_finite_difference(const std::array<stk::math::Vec<double,DIM>,NPE> & coords,
    const ObjFn & objFn,
    const SensFn & sensFn)
{
  const double obj = objFn(coords);
  const auto sens = sensFn(coords);

  const double delta = 1.e-4;
  for (unsigned n=0; n<NPE; ++n)
  {
    for (unsigned i=0; i<DIM; ++i)
    {
      auto perturbCoords = coords;
      perturbCoords[n][i] += delta;
      const double perturbObj = objFn(perturbCoords);
      const auto perturbSens = sensFn(perturbCoords);
      const double fdSens = (perturbObj-obj)/delta;
      const double deltaFDAnalytic = fdSens-sens[n*DIM+i];
      const double deltaAnalytic = perturbSens[n*DIM+i]-sens[n*DIM+i];
      const bool isInRange = (deltaFDAnalytic > 0 && deltaAnalytic > deltaFDAnalytic) || (deltaFDAnalytic < 0 && deltaAnalytic < deltaFDAnalytic);
      EXPECT_TRUE(isInRange);
    }
  }
}

TEST(KinematicObjective, FiniteDifferenceCheck_SethHill3D)
{
  const double refSize = 1.;
  const SethHillParams params{1.0, 10.0, 3, 3};
  const auto objFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_energy(SethHillModel(params), refSize, coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_forces(SethHillModel(params), refSize, coords); };

  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_SethHill2D)
{
  const double refSize = 1.;
  const SethHillParams params{1.0, 10.0, 3, 3};
  const auto objFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_energy(SethHillModel(params), refSize, coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_forces(SethHillModel(params), refSize, coords); };

  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_SizeShape3D)
{
  const double refSize = 1.;
  const auto objFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_energy(SizeShapeModel(), refSize, coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_forces(SizeShapeModel(), refSize, coords); };

  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_SizeShape2D)
{
  const double refSize = 1.;
  const auto objFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_energy(SizeShapeModel(), refSize, coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_forces(SizeShapeModel(), refSize, coords); };

  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_ModifiedNeoHookean3D)
{
  const double refSize = 1.;
  const auto objFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_energy(ModifiedNeoHookeanModel(), refSize, coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_forces(ModifiedNeoHookeanModel(), refSize, coords); };

  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_ModifiedNeoHookean2D)
{
  const double refSize = 1.;
  const auto objFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_energy(ModifiedNeoHookeanModel(), refSize, coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_forces(ModifiedNeoHookeanModel(), refSize, coords); };

  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, testSimplifiedSethHill3DConstantH)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  const double refSize = 1.;

  const double sethHillEnergy = tet4_element_energy(SethHillModel(KinematicObjectiveGoldData::SH_params), refSize, currentCoords);
  const double simplifiedSethHillEnergy = tet4_element_energy(SimplifiedSethHillModel(), refSize, currentCoords);
  EXPECT_NEAR(sethHillEnergy, simplifiedSethHillEnergy, tol);
  const auto sethHillForcesFromSize = tet4_element_forces(SethHillModel(KinematicObjectiveGoldData::SH_params), refSize, currentCoords);
  const auto simplifiedSethHillForcesFromSize = tet4_element_forces(SimplifiedSethHillModel(), refSize, currentCoords);
  expect_sensitivities(sethHillForcesFromSize, simplifiedSethHillForcesFromSize);
}

TEST(KinematicObjective, testSimplifiedSethHill2DConstantH)
{
  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  const double refSize = 1.;

  const double sethHillEnergy = tri3_2d_element_energy(SethHillModel(KinematicObjectiveGoldData::SH_params), refSize, currentCoords);
  const double simplifiedSethHillEnergy = tri3_2d_element_energy(SimplifiedSethHillModel(), refSize, currentCoords);
  EXPECT_NEAR(sethHillEnergy, simplifiedSethHillEnergy, tol);
  const auto sethHillForcesFromSize = tri3_2d_element_forces(SethHillModel(KinematicObjectiveGoldData::SH_params), refSize, currentCoords);
  const auto simplifiedSethHillForcesFromSize = tri3_2d_element_forces(SimplifiedSethHillModel(), refSize, currentCoords);
  expect_sensitivities(sethHillForcesFromSize, simplifiedSethHillForcesFromSize);
}

TEST(KinematicObjective, testAnisotropicDeformationGradientForConstantH)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);

  const Tensor F_refSize = kinematicUtils::compute_deformation_gradient_tet4(KinematicObjectiveGoldData::element_size, currentCoords);
  const Tensor metricMapping = metric::compute_metric_mapping(metric::metric_tensor_for_size(KinematicObjectiveGoldData::element_size));
  const Tensor F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(metricMapping, currentCoords);
  expect_tensor(F_anisotropic, F_refSize);
}

TEST(KinematicObjective, testAnisotropicDeformationGradientMatchesIsotropicForNonUnitSize)
{
  // With metric I/h^2, the anisotropic deformation gradient must match the isotropic one for size h.
  // A unit size is not enough to test this, since the metric is then the identity.
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  for (const double refSize : {0.5, 1.0, 2.0})
  {
    const Tensor F_refSize = kinematicUtils::compute_deformation_gradient_tet4(refSize, currentCoords);
    const Tensor metricMapping = metric::compute_metric_mapping(metric::metric_tensor_for_size(refSize));
    const Tensor F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(metricMapping, currentCoords);
    expect_tensor(F_refSize, F_anisotropic);
  }
}

TEST(KinematicObjective, testAnisotropicDeformationGradientForGoldMetric)
{
  // The metric mapping is unique only up to a rotation on the left (the eigensolver's choice of
  // eigenvector order and sign), which rotates F_U the same way.  The energy depends on F_U only
  // through C_U = F_U^T F_U and J_U = det(F_U), which do not depend on that choice, so compare those.
  const auto & currentCoords = tet_node_coords(KinematicObjectiveAnisotropicGoldData::current_coordinates);
  const Tensor metricMapping = metric::compute_metric_mapping(tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M));

  const Tensor F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(metricMapping, currentCoords);
  const Tensor goldF_anisotropic = tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::Fu);

  expect_tensor(minitensor::transpose(goldF_anisotropic) * goldF_anisotropic, minitensor::transpose(F_anisotropic) * F_anisotropic);
  EXPECT_NEAR(minitensor::det(goldF_anisotropic), minitensor::det(F_anisotropic), tol);
}

TEST(KinematicObjective, testAnisotropicSethHillEnergyAndSensitivitiesForGoldMetric)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveAnisotropicGoldData::current_coordinates);
  const Tensor metric = tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M);
  const Tensor metricMapping = metric::compute_metric_mapping(metric);
  const auto metricMappingArray = metric::to_array(metricMapping);
  expect_tensor(metric, minitensor::transpose(metricMapping) * metricMapping);
  EXPECT_NEAR(KinematicObjectiveAnisotropicGoldData::abs_det_F_M, std::abs(minitensor::det(metricMapping)), tol);

  const SethHillModel model(SethHillParams{1.0, 1.0, 1, 1});
  EXPECT_NEAR(KinematicObjectiveAnisotropicGoldData::integrated_SH_total_energy, tet4_element_energy(model, metricMappingArray.data(), currentCoords), tol);
  expect_sensitivities(KinematicObjectiveAnisotropicGoldData::SH_sensitivities, tet4_element_forces(model, metricMappingArray.data(), currentCoords));
}

TEST(KinematicObjective, testMetricMappingFactorsAnisotropicMetric)
{
  const Tensor metric = tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M);
  const Tensor metricMapping = metric::compute_metric_mapping(metric);

  expect_tensor(metric, minitensor::transpose(metricMapping) * metricMapping);
  EXPECT_NEAR(std::sqrt(minitensor::det(metric)), minitensor::det(metricMapping), tol);
}

TEST(KinematicObjective, testAnisotropicDeformationGradientIsMetricMappingTimesUnitDeformationGradient)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  const Tensor metricMapping = metric::compute_metric_mapping(tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M));

  const Tensor F_unitSize = kinematicUtils::compute_deformation_gradient_tet4(1., currentCoords);
  const Tensor F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(metricMapping, currentCoords);
  expect_tensor(metricMapping * F_unitSize, F_anisotropic);
}

TEST(KinematicObjective, testAnisotropicDeformationGradientIsIdentityForElementIdealInMetric)
{
  // The unit ideal tet mapped by the inverse metric mapping is ideal in the metric, so it is undeformed.
  const Tensor metricMapping = metric::compute_metric_mapping(tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M));
  const Tensor inverseMetricMapping = minitensor::inverse(metricMapping);

  std::array<stk::math::Vector3d, 4> idealInMetricCoords = kinematicUtils::compute_ideal_tet4_reference(1.);
  for (auto & nodeCoords : idealInMetricCoords)
  {
    const minitensor::Vector<double, 3> mapped = inverseMetricMapping * minitensor::Vector<double, 3>(nodeCoords[0], nodeCoords[1], nodeCoords[2]);
    nodeCoords = stk::math::Vector3d(mapped(0), mapped(1), mapped(2));
  }

  const Tensor F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(metricMapping, idealInMetricCoords);
  expect_tensor(minitensor::identity<double, 3>(), F_anisotropic);
}

template<class MODEL>
static void test_anisotropic_sensitivity_matches_finite_difference(const MODEL & model)
{
  const Tensor metricMapping = metric::compute_metric_mapping(tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M));
  const auto metricMappingArray = metric::to_array(metricMapping);
  const auto objFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_energy(model, metricMappingArray.data(), coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector3d,4> & coords)
      { return tet4_element_forces(model, metricMappingArray.data(), coords); };

  const auto & currentCoords = tet_node_coords(KinematicObjectiveAnisotropicGoldData::current_coordinates);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicModifiedNeoHookean3D)
{
  test_anisotropic_sensitivity_matches_finite_difference(ModifiedNeoHookeanModel());
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicSizeShape3D)
{
  test_anisotropic_sensitivity_matches_finite_difference(SizeShapeModel());
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicSethHill3D)
{
  test_anisotropic_sensitivity_matches_finite_difference(SethHillModel(SethHillParams{1.0, 10.0, 3, 3}));
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicSimplifiedSethHill3D)
{
  test_anisotropic_sensitivity_matches_finite_difference(SimplifiedSethHillModel());
}

template<class MODEL>
static void expect_anisotropic_matches_isotropic_for_size_metric(const MODEL & model)
{
  // With metric I/h^2, the anisotropic energy and forces must match the isotropic ones for size h.
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  for (const double refSize : {0.5, 1.0, 2.0})
  {
    const Tensor metricMapping = metric::compute_metric_mapping(metric::metric_tensor_for_size(refSize));
    const auto metricMappingArray = metric::to_array(metricMapping);
    EXPECT_NEAR(tet4_element_energy(model, refSize, currentCoords), tet4_element_energy(model, metricMappingArray.data(), currentCoords), tol);
    expect_sensitivities(tet4_element_forces(model, refSize, currentCoords), tet4_element_forces(model, metricMappingArray.data(), currentCoords));
  }
}

TEST(KinematicObjective, testAnisotropicEnergyAndForcesMatchIsotropicForSizeMetric)
{
  expect_anisotropic_matches_isotropic_for_size_metric(ModifiedNeoHookeanModel());
  expect_anisotropic_matches_isotropic_for_size_metric(SizeShapeModel());
  expect_anisotropic_matches_isotropic_for_size_metric(SethHillModel(SethHillParams{1.0, 10.0, 3, 3}));
  expect_anisotropic_matches_isotropic_for_size_metric(SimplifiedSethHillModel());
}

TEST(KinematicObjective, testAnisotropicElementObjectiveInterface)
{
  const KinematicModelObjective<ModifiedNeoHookeanModel> objective;
  const KinematicElementObjective & elemObjective = objective;
  const Tensor metricMapping = metric::compute_metric_mapping(tensor_from_array<3>(KinematicObjectiveAnisotropicGoldData::M));
  const auto metricMappingArray = metric::to_array(metricMapping);

  const auto & currentCoords = tet_node_coords(KinematicObjectiveAnisotropicGoldData::current_coordinates);
  const std::array<const double *, 4> nodeCoords{{currentCoords[0].data(), currentCoords[1].data(), currentCoords[2].data(), currentCoords[3].data()}};

  EXPECT_NEAR(tet4_element_energy(ModifiedNeoHookeanModel(), metricMappingArray.data(), currentCoords),
      elemObjective.compute_element_objective(3, 4, metricMappingArray.data(), nodeCoords.data()), tol);

  std::vector<double> elemGradContrib;
  elemObjective.fill_element_sensitivity(3, 4, metricMappingArray.data(), nodeCoords.data(), elemGradContrib);
  const auto goldForces = tet4_element_forces(ModifiedNeoHookeanModel(), metricMappingArray.data(), currentCoords);
  ASSERT_EQ(goldForces.size(), elemGradContrib.size());
  for (size_t i=0; i<goldForces.size(); ++i)
    EXPECT_NEAR(goldForces[i], elemGradContrib[i], tol);

  EXPECT_ANY_THROW(elemObjective.compute_element_objective(3, 3, metricMappingArray.data(), nodeCoords.data()));
}

TEST(KinematicObjective, testAnisotropicDeformationGradient2DMatchesIsotropicForSizeMetric)
{
  // With metric I/h^2, the anisotropic deformation gradient must match the isotropic one for size h.
  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  for (const double refSize : {0.5, 1.0, 2.0})
  {
    const Tensor2d F_refSize = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, currentCoords);
    const Tensor2d metricMapping = metric::compute_metric_mapping_2d(metric::metric_tensor_for_size_2d(refSize));
    const Tensor2d F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tri3_2d(metricMapping, currentCoords);
    expect_tensor(F_refSize, F_anisotropic);
  }
}

TEST(KinematicObjective, testMetricMapping2DFactorsAnisotropicMetric)
{
  const Tensor2d metric = tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M);
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(metric);

  expect_tensor(metric, minitensor::transpose(metricMapping) * metricMapping);
  EXPECT_NEAR(std::sqrt(minitensor::det(metric)), minitensor::det(metricMapping), tol);
}

TEST(KinematicObjective, testAnisotropicDeformationGradient2DForGoldMetric)
{
  // The metric mapping is unique only up to a rotation on the left, which rotates F_U the same way,
  // so compare C_U = F_U^T F_U and J_U = det(F_U), which do not depend on that choice.
  const auto & currentCoords = tri_node_coords(KinematicObjectiveAnisotropicGoldData2d::current_coordinates);
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M));

  const Tensor2d F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tri3_2d(metricMapping, currentCoords);
  const Tensor2d goldF_anisotropic = tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::Fu);

  expect_tensor(minitensor::transpose(goldF_anisotropic) * goldF_anisotropic, minitensor::transpose(F_anisotropic) * F_anisotropic);
  EXPECT_NEAR(minitensor::det(goldF_anisotropic), minitensor::det(F_anisotropic), tol);
}

TEST(KinematicObjective, testAnisotropicSethHill2DEnergyAndSensitivitiesForGoldMetric)
{
  const auto & currentCoords = tri_node_coords(KinematicObjectiveAnisotropicGoldData2d::current_coordinates);
  const Tensor2d metric = tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M);
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(metric);
  const auto metricMappingArray = metric::to_array(metricMapping);
  expect_tensor(metric, minitensor::transpose(metricMapping) * metricMapping);
  EXPECT_NEAR(KinematicObjectiveAnisotropicGoldData2d::abs_det_F_M, std::abs(minitensor::det(metricMapping)), tol);

  const SethHillModel model(KinematicObjectiveAnisotropicGoldData2d::SH_params);
  EXPECT_NEAR(KinematicObjectiveAnisotropicGoldData2d::integrated_SH_total_energy, tri3_2d_element_energy(model, metricMappingArray.data(), currentCoords), tol);
  expect_sensitivities(KinematicObjectiveAnisotropicGoldData2d::SH_sensitivities, tri3_2d_element_forces(model, metricMappingArray.data(), currentCoords));
}

TEST(KinematicObjective, testAnisotropicDeformationGradient2DIsMetricMappingTimesUnitDeformationGradient)
{
  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M));

  const Tensor2d F_unitSize = kinematicUtils::compute_deformation_gradient_tri3_2d(1., currentCoords);
  const Tensor2d F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tri3_2d(metricMapping, currentCoords);
  expect_tensor(metricMapping * F_unitSize, F_anisotropic);
}

TEST(KinematicObjective, testAnisotropicDeformationGradient2DIsIdentityForElementIdealInMetric)
{
  // The unit ideal triangle mapped by the inverse metric mapping is ideal in the metric, so it is undeformed.
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M));
  const Tensor2d inverseMetricMapping = minitensor::inverse(metricMapping);

  std::array<stk::math::Vector2d, 3> idealInMetricCoords = kinematicUtils::compute_ideal_tri3_2d_reference(1.);
  for (auto & nodeCoords : idealInMetricCoords)
  {
    const minitensor::Vector<double, 2> mapped = inverseMetricMapping * minitensor::Vector<double, 2>(nodeCoords[0], nodeCoords[1]);
    nodeCoords = stk::math::Vector2d(mapped(0), mapped(1));
  }

  const Tensor2d F_anisotropic = kinematicUtils::compute_anisotropic_deformation_gradient_tri3_2d(metricMapping, idealInMetricCoords);
  expect_tensor(minitensor::identity<double, 2>(), F_anisotropic);
}

template<class MODEL>
static void test_anisotropic_2d_sensitivity_matches_finite_difference(const MODEL & model)
{
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M));
  const auto metricMappingArray = metric::to_array(metricMapping);
  const auto objFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_energy(model, metricMappingArray.data(), coords); };
  const auto sensFn = [&](const std::array<stk::math::Vector2d,3> & coords)
      { return tri3_2d_element_forces(model, metricMappingArray.data(), coords); };

  const auto & currentCoords = tri_node_coords(KinematicObjectiveAnisotropicGoldData2d::current_coordinates);
  test_sensitivity_matches_finite_difference(currentCoords, objFn, sensFn);
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicModifiedNeoHookean2D)
{
  test_anisotropic_2d_sensitivity_matches_finite_difference(ModifiedNeoHookeanModel());
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicSizeShape2D)
{
  test_anisotropic_2d_sensitivity_matches_finite_difference(SizeShapeModel());
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicSethHill2D)
{
  test_anisotropic_2d_sensitivity_matches_finite_difference(SethHillModel(SethHillParams{1.0, 10.0, 3, 3}));
}

TEST(KinematicObjective, FiniteDifferenceCheck_AnisotropicSimplifiedSethHill2D)
{
  test_anisotropic_2d_sensitivity_matches_finite_difference(SimplifiedSethHillModel());
}

template<class MODEL>
static void expect_anisotropic_2d_matches_isotropic_for_size_metric(const MODEL & model)
{
  // With metric I/h^2, the anisotropic energy and forces must match the isotropic ones for size h.
  const auto & currentCoords = tri_node_coords(KinematicObjectiveGoldData2d::current_coords);
  for (const double refSize : {0.5, 1.0, 2.0})
  {
    const Tensor2d metricMapping = metric::compute_metric_mapping_2d(metric::metric_tensor_for_size_2d(refSize));
    const auto metricMappingArray = metric::to_array(metricMapping);
    EXPECT_NEAR(tri3_2d_element_energy(model, refSize, currentCoords), tri3_2d_element_energy(model, metricMappingArray.data(), currentCoords), tol);
    expect_sensitivities(tri3_2d_element_forces(model, refSize, currentCoords), tri3_2d_element_forces(model, metricMappingArray.data(), currentCoords));
  }
}

TEST(KinematicObjective, testAnisotropicEnergyAndForces2DMatchIsotropicForSizeMetric)
{
  expect_anisotropic_2d_matches_isotropic_for_size_metric(ModifiedNeoHookeanModel());
  expect_anisotropic_2d_matches_isotropic_for_size_metric(SizeShapeModel());
  expect_anisotropic_2d_matches_isotropic_for_size_metric(SethHillModel(SethHillParams{1.0, 10.0, 3, 3}));
  expect_anisotropic_2d_matches_isotropic_for_size_metric(SimplifiedSethHillModel());
}

TEST(KinematicObjective, testAnisotropic2DElementObjectiveInterface)
{
  const KinematicModelObjective<ModifiedNeoHookeanModel> objective;
  const KinematicElementObjective & elemObjective = objective;
  const Tensor2d metricMapping = metric::compute_metric_mapping_2d(tensor_from_array<2>(KinematicObjectiveAnisotropicGoldData2d::M));
  const auto metricMappingArray = metric::to_array(metricMapping);

  const auto & currentCoords = tri_node_coords(KinematicObjectiveAnisotropicGoldData2d::current_coordinates);
  const std::array<const double *, 3> nodeCoords{{currentCoords[0].data(), currentCoords[1].data(), currentCoords[2].data()}};

  EXPECT_NEAR(tri3_2d_element_energy(ModifiedNeoHookeanModel(), metricMappingArray.data(), currentCoords),
      elemObjective.compute_element_objective(2, 3, metricMappingArray.data(), nodeCoords.data()), tol);

  std::vector<double> elemGradContrib;
  elemObjective.fill_element_sensitivity(2, 3, metricMappingArray.data(), nodeCoords.data(), elemGradContrib);
  const auto goldForces = tri3_2d_element_forces(ModifiedNeoHookeanModel(), metricMappingArray.data(), currentCoords);
  ASSERT_EQ(goldForces.size(), elemGradContrib.size());
  for (size_t i=0; i<goldForces.size(); ++i)
    EXPECT_NEAR(goldForces[i], elemGradContrib[i], tol);
}

}
