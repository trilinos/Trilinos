// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MEANRATIOABOUTNODE_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MEANRATIOABOUTNODE_HPP_

#include <vector>

#include <stk_math/StkVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_Objective_Mesh.hpp>

namespace krino {

class MeanRatioAboutNodeObjective : public NodeCoordinateObjectiveInterface
{
public:
  virtual double compute_value(const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const stk::mesh::Entity node,
      const stk::math::Vector3d & x) const override;

  virtual void fill_gradient(const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const GradientFilter & gradientFilter,
      const stk::mesh::Entity node,
      const stk::math::Vector3d & x,
      stk::math::Vector3d & g) const override;
};

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MEANRATIOABOUTNODE_HPP_ */
