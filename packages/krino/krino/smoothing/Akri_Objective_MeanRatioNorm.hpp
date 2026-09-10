// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MEANRATIONORM_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MEANRATIONORM_HPP_

#include <array>
#include <Akri_Objective_Mesh.hpp>

namespace krino {

class MeanRatioNormObjective : public ElementObjective
{
public:
  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const double * const * elemNodeCoords) const override;

  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const override;

  static double qualityEpsilon;

  static void set_quality_epsilon_for_handling_inversion(const double eps);

  static double global_objective_element_contribution(const double elemQuality);
  static double global_objective_sensitivity_wrt_quality(const double elemQuality);
};

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MEANRATIONORM_HPP_ */
