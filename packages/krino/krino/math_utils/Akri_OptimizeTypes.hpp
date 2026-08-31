// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZETYPES_HPP_
#define KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZETYPES_HPP_
#include <stk_math/StkVector.hpp>
#include <functional>
#include <Akri_DistributedVector.hpp>

namespace krino {

enum class OptimizerType { NormaSteepestDescent, ROL, SteepestDescent, LBFGS };

class ObjectiveInterface;
using OptimizerFunction = std::function<void(const ObjectiveInterface &objFn, DistributedVector &x)>;

class Objective3DInterface;
using Optimizer3DFunction = std::function<void(const Objective3DInterface &objFn, stk::math::Vector3d &x)>;

struct SolverConfig
{
  // Time stepping
  int numTimeSteps = 0;

  OptimizerType solverType = OptimizerType::NormaSteepestDescent;
  int minIterations = 0;
  int maxIterations = 0;
  double relativeTolerance = 0.0;
  double absoluteTolerance = 0.0;
  double stepLength = 0.0;
  double lineSearchBacktrackFactor = 0.0;
  double lineSearchDecreaseFactor = 0.0;
  int lineSearchMaxIterations = 0;
};

}



#endif /* KRINO_KRINO_MATH_UTILS_AKRI_OPTIMIZETYPES_HPP_ */
