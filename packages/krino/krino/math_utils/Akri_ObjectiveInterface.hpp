// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVEINTERFACE_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVEINTERFACE_HPP_

#include <array>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <Akri_DistributedVector.hpp>

namespace krino {

class ObjectiveInterface
{
public:
  using VectorType = DistributedVector;
  virtual ~ObjectiveInterface() {}
  virtual double compute_value(const DistributedVector & x) const = 0;
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const = 0;
  virtual void project_solution(DistributedVector & x) const {} // empty for unconstrained problems
};

class Objective3DInterface
{
public:
  using VectorType = stk::math::Vector3d;
  virtual ~Objective3DInterface() {}
  virtual double compute_value(const stk::math::Vector3d & x) const = 0;
  virtual void fill_gradient(const stk::math::Vector3d & x, stk::math::Vector3d & g) const = 0;
  virtual void project_solution(stk::math::Vector3d & x) const {} // empty for unconstrained problems
};

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVEINTERFACE_HPP_ */
