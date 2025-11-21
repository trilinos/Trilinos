// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYROL_ETI
#define PYROL_ETI

#include <PyROL_ETI_helper.hpp>

#include <ROL_BoundConstraint_SimOpt.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_DynamicConstraintCheck.hpp>
#include <ROL_DynamicObjectiveCheck.hpp>
#include <ROL_MonteCarloGenerator.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Objective_SimOpt.hpp>
#include <ROL_OED_Factory.hpp>
#include <ROL_PrimalDualRisk.hpp>
#include <ROL_Problem.hpp>
#include <ROL_ReducedDynamicObjective.hpp>
#include <ROL_Reduced_Objective_SimOpt.hpp>
#include <ROL_SampleGenerator.hpp>
#include <ROL_SerialConstraint.hpp>
#include <ROL_SerialObjective.hpp>
#include <ROL_SimConstraint.hpp>
#include <ROL_Solver.hpp>
#include <ROL_StochasticProblem.hpp>
#include <ROL_UserInputGenerator.hpp>
#include <ROL_ValidateFunction.hpp>
#include <ROL_Vector.hpp>
#include <ROL_Vector_SimOpt.hpp>

#define BINDER_ETI_ABSTRACT(CLASS_NAME) \
  template class CLASS_NAME;

// #define BINDER_ETI_WITH_FOO(CLASS_NAME) \
//   template class CLASS_NAME; \
//   template <> inline void PyROL::foo(CLASS_NAME a){}

#define BINDER_ROL_CORE(SCALAR) \
  BINDER_ETI_ABSTRACT(Constraint<SCALAR>) \
  BINDER_ETI_ABSTRACT(Objective<SCALAR>) \
  BINDER_ETI_ABSTRACT(Problem<SCALAR>) \
  BINDER_ETI_ABSTRACT(Solver<SCALAR>) \
  BINDER_ETI_ABSTRACT(Vector<SCALAR>)

#define BINDER_ROL_SIMOPT(SCALAR) \
  BINDER_ETI_ABSTRACT(BoundConstraint_SimOpt<SCALAR>) \
  BINDER_ETI_ABSTRACT(Reduced_Objective_SimOpt<SCALAR>) \
  BINDER_ETI_ABSTRACT(SimConstraint<SCALAR>) \
  BINDER_ETI_ABSTRACT(Vector_SimOpt<SCALAR>)

#define BINDER_ROL_DYNAMIC(SCALAR) \
  BINDER_ETI_ABSTRACT(DynamicConstraintCheck<SCALAR>) \
  BINDER_ETI_ABSTRACT(DynamicObjectiveCheck<SCALAR>) \
  BINDER_ETI_ABSTRACT(ReducedDynamicObjective<SCALAR>) \
  BINDER_ETI_ABSTRACT(SerialConstraint<SCALAR>) \
  BINDER_ETI_ABSTRACT(SerialObjective<SCALAR>)

#define BINDER_ROL_UTILS(SCALAR) \
  BINDER_ETI_ABSTRACT(ValidateFunction<SCALAR>)

#define BINDER_ROL_STOCHASTIC(SCALAR) \
  BINDER_ETI_ABSTRACT(MonteCarloGenerator<SCALAR>) \
  BINDER_ETI_ABSTRACT(PrimalDualRisk<SCALAR>) \
  BINDER_ETI_ABSTRACT(RiskNeutralObjective<SCALAR>) \
  BINDER_ETI_ABSTRACT(SampleGenerator<SCALAR>) \
  BINDER_ETI_ABSTRACT(StochasticProblem<SCALAR>) \
  BINDER_ETI_ABSTRACT(UserInputGenerator<SCALAR>)

#define BINDER_ROL_OED(SCALAR) \
  BINDER_ETI_ABSTRACT(Factory<SCALAR>)

namespace ROL {

  BINDER_ROL_CORE(double)
  BINDER_ROL_SIMOPT(double)
  BINDER_ROL_DYNAMIC(double)
  BINDER_ROL_STOCHASTIC(double)

namespace details {
  BINDER_ROL_UTILS(double)
}

namespace OED {
  BINDER_ROL_OED(double)
}

}

#endif // PYROL_ETI
