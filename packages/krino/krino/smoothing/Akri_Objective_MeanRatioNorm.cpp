// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Objective_MeanRatioNorm.hpp>

#include <Akri_Smoothing_Utils.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_QualityMetricWithSensitivities.hpp>

namespace krino {

double MeanRatioNormObjective::qualityEpsilon = 0; // default is not to use limiter

static double limited_global_objective_element_contribution(const double qualityEpsilon, const double elemQuality)
{
  if (elemQuality >= qualityEpsilon)
    return 1./elemQuality;
  return 2./qualityEpsilon - elemQuality/(qualityEpsilon*qualityEpsilon);
}

static double limited_global_objective_sensitivity_wrt_quality(const double qualityEpsilon, const double elemQuality)
{
  if (elemQuality >= qualityEpsilon)
    return -1./(elemQuality*elemQuality);
  return -1./(qualityEpsilon*qualityEpsilon);
}

static double unlimited_global_objective_element_contribution(const double elemQuality)
{
  if (elemQuality > 0.)
    return 1./elemQuality;
  return std::numeric_limits<double>::infinity();
}

static double unlimited_global_objective_sensitivity_wrt_quality(const double elemQuality)
{
  if (elemQuality > 0.)
    return -1./(elemQuality*elemQuality);
  return -std::numeric_limits<double>::infinity();
}

void MeanRatioNormObjective::set_quality_epsilon_for_handling_inversion(const double eps)
{
  STK_ThrowRequire(eps >= 0.);
  qualityEpsilon = eps;
}

double MeanRatioNormObjective::global_objective_element_contribution(const double elemQuality)
{
  return (qualityEpsilon == 0.) ?
      unlimited_global_objective_element_contribution(elemQuality) :
      limited_global_objective_element_contribution(qualityEpsilon, elemQuality);
}

double MeanRatioNormObjective::global_objective_sensitivity_wrt_quality(const double elemQuality)
{
  return (qualityEpsilon == 0.) ?
      unlimited_global_objective_sensitivity_wrt_quality(elemQuality) :
      limited_global_objective_sensitivity_wrt_quality(qualityEpsilon, elemQuality);
}

double MeanRatioNormObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const double * const * elemNodeCoords) const
{
  STK_ThrowRequire((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4));
  MeanRatioQualityMetric qualityMetric;
  const double quality = (spatialDim == 2) ?
      qualityMetric.tri2d_mean_ratio(elemNodeCoords) :
      qualityMetric.tet_mean_ratio(elemNodeCoords);
  return global_objective_element_contribution(quality);
}

void MeanRatioNormObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  STK_ThrowRequire((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4));
  if (spatialDim == 2)
  {
    const auto & [quality, elemSens] = MeanRatioQualityMetricWithSensitivities::tri2d_quality_and_sensitivities(elemNodeCoords[0], elemNodeCoords[1], elemNodeCoords[2]);
    const double dObj_dQual = global_objective_sensitivity_wrt_quality(quality);
    elemGradContrib.resize(npe*spatialDim);
    for (unsigned i = 0; i < npe*spatialDim; ++i)
      elemGradContrib[i] = dObj_dQual * elemSens[i];
  }
  else
  {
    const std::array<stk::math::Vector3d,4> elemNodeCoordVecs
      {{
        stk::math::Vector3d(elemNodeCoords[0]), stk::math::Vector3d(elemNodeCoords[1]), stk::math::Vector3d(elemNodeCoords[2]), stk::math::Vector3d(elemNodeCoords[3])
      }};
    const auto & [quality, elemSens] = MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivities(elemNodeCoordVecs);
    const double dObj_dQual = global_objective_sensitivity_wrt_quality(quality);
    elemGradContrib.resize(npe*spatialDim);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < spatialDim; ++j)
          elemGradContrib[spatialDim*i + j] = dObj_dQual * elemSens[i][j];
  }
}

}


