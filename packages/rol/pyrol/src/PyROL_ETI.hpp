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

#include <ROL_Vector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Solver.hpp>
#include <ROL_Problem.hpp>

#include <ROL_Vector_SimOpt.hpp>
#include <ROL_Objective_SimOpt.hpp>
#include <ROL_Reduced_Objective_SimOpt.hpp>
#include <ROL_SimConstraint.hpp>
#include <ROL_BoundConstraint_SimOpt.hpp>
#include <ROL_SerialConstraint.hpp>

#include <ROL_OED_Factory.hpp>
#include <ROL_ReducedDynamicObjective.hpp>

#include <PyROL_ETI_helper.hpp>

#define BINDER_ETI_ABSTRACT(CLASS_NAME) \
  template class CLASS_NAME;

#define BINDER_ETI_WITH_FOO(CLASS_NAME) \
  template class CLASS_NAME; \
  template <> inline void PyROL::foo(CLASS_NAME a){}

#define BINDER_ROL_VECTOR(SCALAR) \
  BINDER_ETI_ABSTRACT(Vector<SCALAR>) \
  BINDER_ETI_ABSTRACT(Vector_SimOpt<SCALAR>)

#define BINDER_ROL_OBJECTIVE(SCALAR) \
  BINDER_ETI_ABSTRACT(Objective<SCALAR>) \
  BINDER_ETI_ABSTRACT(Objective_SimOpt<SCALAR>) \
  BINDER_ETI_ABSTRACT(Reduced_Objective_SimOpt<SCALAR>) \
  BINDER_ETI_ABSTRACT(ReducedDynamicObjective<SCALAR>)

#define BINDER_ROL_CONSTRAINT(SCALAR) \
  BINDER_ETI_ABSTRACT(Constraint<SCALAR>) \
  BINDER_ETI_ABSTRACT(SimConstraint<SCALAR>) \
  BINDER_ETI_ABSTRACT(BoundConstraint_SimOpt<SCALAR>) \
  BINDER_ETI_ABSTRACT(SerialConstraint<SCALAR>)

#define BINDER_ROL_SOLVER(SCALAR) \
  BINDER_ETI_ABSTRACT(Solver<SCALAR>)

#define BINDER_ROL_PROBLEM(SCALAR) \
  BINDER_ETI_ABSTRACT(Problem<SCALAR>)

#define BINDER_ROL_OED(SCALAR) \
  BINDER_ETI_ABSTRACT(Factory<SCALAR>)

namespace ROL {

  BINDER_ROL_VECTOR(double)
  BINDER_ROL_OBJECTIVE(double)
  BINDER_ROL_CONSTRAINT(double)
  BINDER_ROL_SOLVER(double)
  BINDER_ROL_PROBLEM(double)

namespace OED {
  BINDER_ROL_OED(double)
}

}

#endif // PYROL_ETI
