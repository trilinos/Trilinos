// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_QUALITY_METRIC_WITH_SENSITIVITIES_H
#define AKRI_QUALITY_METRIC_WITH_SENSITIVITIES_H
#include "Sacado.hpp"
#include <stk_math/StkVector.hpp>
#include <vector>

namespace krino
{

std::tuple<double, std::array<double,3>> tet_volume_and_sensitivity_to_nth_vertex(const std::array<stk::math::Vector3d, 4> & verts, const int n);

class ScaledJacobianQualityMetricWithSensitivities
{
public:
  ScaledJacobianQualityMetricWithSensitivities() = default;
  virtual ~ScaledJacobianQualityMetricWithSensitivities() {}

  static std::tuple<double, std::array<double,3>> tet_quality_and_sensitivity_to_nth_vertex(const std::array<stk::math::Vector3d, 4> & verts, const int n);
};

class MeanRatioQualityMetricWithSensitivities
{
public:
  MeanRatioQualityMetricWithSensitivities() = default;
  virtual ~MeanRatioQualityMetricWithSensitivities() {}

  static std::tuple<double, std::array<double,3>> tet_quality_and_sensitivity_to_nth_vertex(const std::array<stk::math::Vector3d, 4> & verts, const int n);
};

}

#endif
