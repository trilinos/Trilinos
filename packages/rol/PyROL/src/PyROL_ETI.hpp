#ifndef PYROL_ETI
#define PYROL_ETI

#include <ROL_Vector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_QuadraticObjective.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Solver.hpp>
#include <ROL_Algorithm.hpp>
#include <ROL_TrustRegionStep.hpp>
#include <ROL_Problem.hpp>

#include <PyROL_ETI_helper.hpp>

#define BINDER_ETI_ABSTRACT(CLASS_NAME) \
  template class CLASS_NAME;

#define BINDER_ETI_WITH_FOO(CLASS_NAME) \
  template class CLASS_NAME; \
  template <> inline void PyROL::foo(CLASS_NAME a){}

#define BINDER_ROL_VECTOR(SCALAR) \
  BINDER_ETI_ABSTRACT(Vector<SCALAR>)

#define BINDER_ROL_OBJECTIVE(SCALAR) \
  BINDER_ETI_ABSTRACT(Objective<SCALAR>) \
  BINDER_ETI_WITH_FOO(QuadraticObjective<SCALAR>)

#define BINDER_ROL_CONSTRAINT(SCALAR) \
  BINDER_ETI_ABSTRACT(Constraint<SCALAR>)

#define BINDER_ROL_SOLVER(SCALAR) \
  BINDER_ETI_ABSTRACT(Solver<SCALAR>) \
  BINDER_ETI_WITH_FOO(Algorithm<SCALAR>) \
  BINDER_ETI_WITH_FOO(TrustRegionStep<SCALAR>)

#define BINDER_ROL_PROBLEM(SCALAR) \
  BINDER_ETI_ABSTRACT(Problem<SCALAR>)

namespace ROL {

  BINDER_ROL_VECTOR(double)
  BINDER_ROL_OBJECTIVE(double)
  BINDER_ROL_CONSTRAINT(double)
  BINDER_ROL_SOLVER(double)
  BINDER_ROL_PROBLEM(double)

}

#endif // PYROL_ETI
