// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_QUALITY_METRIC_H
#define AKRI_QUALITY_METRIC_H
#include <stk_math/StkVector.hpp>
#include <vector>

namespace krino
{

class QualityMetric
{
public:
    virtual ~QualityMetric() {}

    virtual bool is_first_quality_metric_better_than_second(const double firstValue, const double secondValue) const = 0;
    virtual double get_element_quality_metric(const std::vector<stk::math::Vector3d> &nodeLocations) const = 0;
    virtual double get_best_value_for_metric() const = 0;
    virtual double get_acceptable_value_for_metric() const = 0;
};

class MeanRatioQualityMetric : public QualityMetric
{
public:
    virtual ~MeanRatioQualityMetric() {}

    virtual double get_best_value_for_metric() const override { return 1.0; }
    virtual bool is_first_quality_metric_better_than_second(const double firstValue, const double secondValue) const override
    {
      return static_cast<float>(firstValue) > static_cast<float>(secondValue);
    }

    double get_element_quality_metric(const std::vector<stk::math::Vector3d> &nodeLocations) const override
    {
      return tet_mean_ratio(nodeLocations);
    }

    static double tet_mean_ratio(const std::vector<stk::math::Vector3d> &nodeLocations);

    virtual double get_acceptable_value_for_metric() const override { return 0.2; }
};

class ScaledJacobianQualityMetric : public QualityMetric
{
public:
  ScaledJacobianQualityMetric() = default;
  virtual ~ScaledJacobianQualityMetric() {}

  virtual double get_best_value_for_metric() const override { return 1.0; }
  virtual bool is_first_quality_metric_better_than_second(
      const double firstValue, const double secondValue) const override
  {
    return static_cast<float>(firstValue) > static_cast<float>(secondValue);
  }

  double get_element_quality_metric(const std::vector<stk::math::Vector3d> &nodeLocations) const override
  {
    return (nodeLocations.size() == 4 || nodeLocations.size() == 10) ? tet_scaled_jacobian(nodeLocations) : tri2d_scaled_jacobian(nodeLocations);
  }

  static double tet_scaled_jacobian(const std::vector<stk::math::Vector3d> &nodeLocations);
  static double tri2d_scaled_jacobian(const std::vector<stk::math::Vector3d> &nodeLocations);
  static double tri3d_scaled_jacobian(const std::vector<stk::math::Vector3d> &nodeLocations);

  virtual double get_acceptable_value_for_metric() const override { return 0.1; }
};

}

#endif
