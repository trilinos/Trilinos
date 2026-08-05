// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_SETHHILL_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_SETHHILL_HPP_
#include <vector>

#include <stk_math/StkVector.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_SethHillConfig.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_Objective_KinematicMesh.hpp>

namespace krino {

class SethHillSmoothingObjective : public KinematicElementObjective
{
public:
  SethHillSmoothingObjective(const SethHillParams & params)
  : mySethHillParams(params) {}

  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords) const override;

  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const override;

  template<typename ELEMCOORDS>
  static double compute_tet4_element_energy(
      const double refSize,
      const ELEMCOORDS & current_coords,
      const SethHillParams & params);
  static double compute_tet4_element_energy(
      const std::array<stk::math::Vector3d, 4> & ref,
      const std::array<stk::math::Vector3d, 4> & current_coords,
      const SethHillParams & params);
  static std::array<double, 12> compute_tet4_element_forces(
      const std::array<stk::math::Vector3d, 4> & ref,
      const std::array<stk::math::Vector3d, 4> & current_coords,
      const SethHillParams & params);
  template<typename ELEMCOORDS>
  static std::array<double, 12> compute_tet4_element_forces(
      const double refSize,
      const ELEMCOORDS & current_coords,
      const SethHillParams & params);

  template<typename ELEMCOORDS>
  static double compute_tri3_2d_element_energy(
      const double refSize,
      const ELEMCOORDS & current_coords,
      const SethHillParams & params);

  template<typename ELEMCOORDS>
  static std::array<double, 6> compute_tri3_2d_element_forces(const double refSize,
      const ELEMCOORDS & current_coords,
      const SethHillParams & params);

  static std::array<double, 6> compute_tri3_2d_element_forces(const std::array<stk::math::Vector2d, 3> & ref,
      const std::array<stk::math::Vector2d, 3> & current_coords,
      const SethHillParams & params);

private:
  SethHillParams mySethHillParams;
};

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_SETHHILL_HPP_ */
