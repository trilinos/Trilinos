// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Objective_SethHill.hpp>
#include <MiniTensor_Tensor.h>
#include <Akri_MeshHelpers.hpp>
#include <Akri_KinematicUtils.hpp>
#include <Akri_Objective_SizeShape.hpp>
#include <Akri_SethHillConfig.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_UnitTestUtils.hpp>

namespace krino {

namespace KinematicObjectiveGoldData {
constexpr std::array<double, 12> reference_coords
{{
  -0.5, -0.2886751345948129, 0.0,
   0.5, -0.2886751345948129, 0.0,
   0.0,  0.5773502691896258, 0.0,
   0.0,  0.0,                0.816496580927726
}};

constexpr std::array<double, 12> current_coords
{{
  -0.5,                 -0.2886751345948129,  0.0,
   0.8779644730092272,   0.158538460905145,  -0.2773500981126146,
  -0.24253562503633297,  0.6682593600987168,  0.07692307692307693,
   0.816496580927726,    0.047619047619047616, 1.632993161855452
}};

constexpr std::array<double, 9> deformation_gradient
{{
  1.3779644730092273,  -0.49827390704159435,  0.9447114043355755,
  0.4472135954999579,   0.8467738864691344,  -0.1613664824343352,
 -0.2773500981126146,   0.24895127214195267, 2.0818239887633934
}};

constexpr double integrated_eta_total_energy = 0.104381041186412;

constexpr std::array<double, 12> eta_sensitivities
{{
 -0.06631816521021236, -0.02218050023669047, -0.03592089463382403,
  0.05019083546471202,  0.004821870442533769,-0.053762036979147144,
 -0.041357965034233404, 0.0257742984284479,  -0.030958626328387587,
  0.05748529477973371, -0.0084156686342912,  0.12064155794135877
}};

constexpr double integrated_SH_total_energy = 37.78876999854464;

constexpr std::array<double, 12> SH_sensitivities
{{
  -101.05024973329421, -201.6675808194271,    -16.307961024769618,
   150.32087696806474,  -29.507632313702175, -115.44947720635962,
   -80.85413528146395,  248.53458507085506,    10.111444241766453,
    31.58350804669342,  -17.35937193772578,   121.64599398936278,
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
static void expect_deformation_gradient(const std::array<double, 9> & gold, const Tensor & F)
{
  for (unsigned i=0; i<3; ++i)
  {
    for (unsigned j=0; j<3; ++j)
    {
      EXPECT_NEAR(gold[i*3+j], F(i,j), tol);
    }
  }
}
static void expect_deformation_gradient(const std::array<double, 4> & gold, const Tensor2d & F)
{
  for (unsigned i=0; i<2; ++i)
  {
    for (unsigned j=0; j<2; ++j)
    {
      EXPECT_NEAR(gold[i*2+j], F(i,j), tol);
    }
  }
}
static void expect_sensitivities(const std::array<double, 12> & gold, const std::array<stk::math::Vector3d, 4> & sens)
{
  for (unsigned i=0; i<4; ++i)
  {
    for (unsigned j=0; j<3; ++j)
    {
      EXPECT_NEAR(gold[i*3+j], sens[i][j], tol);
    }
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

TEST(KinematicObjective, testEnergyAndSensitivities)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  const auto & refCoords = tet_node_coords(KinematicObjectiveGoldData::reference_coords);

  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refCoords, currentCoords);
  expect_deformation_gradient(KinematicObjectiveGoldData::deformation_gradient, F);

  const SethHillParams params{1.0, 10.0, 3, 3};
  const double sethHillEnergy = SethHillSmoothingObjective::compute_tet4_element_energy(refCoords, currentCoords, params);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_SH_total_energy, sethHillEnergy, tol);
  const auto sethHillForces = SethHillSmoothingObjective::compute_tet4_element_forces(refCoords, currentCoords, params);
  expect_sensitivities(KinematicObjectiveGoldData::SH_sensitivities, sethHillForces);

  const double sizeShapePseudoEnergy = SizeShapeObjective::compute_tet4_element_pseudo_energy(refCoords, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_eta_total_energy, sizeShapePseudoEnergy, tol);
  const auto sizeShapePseudoForces = SizeShapeObjective::compute_tet4_element_pseudo_forces(refCoords, currentCoords);
  expect_sensitivities(KinematicObjectiveGoldData::eta_sensitivities, sizeShapePseudoForces);
}

TEST(KinematicObjective, testEnergyAndSensitivitiesForConstantH)
{
  const auto & currentCoords = tet_node_coords(KinematicObjectiveGoldData::current_coords);
  const double refSize = 1.;

  const Tensor F_refSize = kinematicUtils::compute_deformation_gradient_tet4(refSize, currentCoords);
  const Tensor F_ref     = kinematicUtils::compute_deformation_gradient_tet4(kinematicUtils::compute_ideal_tet4_reference(refSize), currentCoords);
  expect_tensor(F_ref, F_refSize);

  const SethHillParams params{1.0, 10.0, 3, 3};
  const double sethHillEnergy = SethHillSmoothingObjective::compute_tet4_element_energy(refSize, currentCoords, params);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_SH_total_energy, sethHillEnergy, tol);
  const auto sethHillForces = SethHillSmoothingObjective::compute_tet4_element_forces(refSize, currentCoords, params);
  expect_sensitivities(KinematicObjectiveGoldData::SH_sensitivities, sethHillForces);

  const double sizeShapePseudoEnergy = SizeShapeObjective::compute_tet4_element_pseudo_energy(refSize, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData::integrated_eta_total_energy, sizeShapePseudoEnergy, tol);
  const auto sizeShapePseudoForces = SizeShapeObjective::compute_tet4_element_pseudo_forces(refSize, currentCoords);
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
  expect_deformation_gradient(KinematicObjectiveGoldData2d::deformation_gradient, F);

  const double sizeShapePseudoEnergy = SizeShapeObjective::compute_tri3_2d_element_pseudo_energy(refSize, currentCoords);
  EXPECT_NEAR(KinematicObjectiveGoldData2d::integrated_eta_total_energy, sizeShapePseudoEnergy, tol);

  const auto sizeShapePseudoForcesFromSize = SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(refSize, currentCoords);
  const auto sizeShapePseudoForcesFromRefCoords = SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(refCoords, currentCoords);

  for (int i=0; i<6; ++i)
    EXPECT_NEAR(sizeShapePseudoForcesFromSize[i], sizeShapePseudoForcesFromRefCoords[i], tol);

  const SethHillParams params{1.0, 10.0, 3, 3};
  const auto sethHillForcesFromSize = SethHillSmoothingObjective::compute_tri3_2d_element_forces(refSize, currentCoords, params);
  const auto sethHillForcesFromRefCoords = SethHillSmoothingObjective::compute_tri3_2d_element_forces(refCoords, currentCoords, params);

  for (int i=0; i<6; ++i)
    EXPECT_NEAR(sethHillForcesFromSize[i], sethHillForcesFromRefCoords[i], tol);
}

}


