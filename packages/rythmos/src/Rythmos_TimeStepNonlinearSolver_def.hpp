//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_DEF_HPP
#define RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_DEF_HPP

#include "Rythmos_TimeStepNonlinearSolver_decl.hpp"

#include "Thyra_TestingTools.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_as.hpp"

namespace Rythmos {


// ////////////////////////
// Defintions


// Static members


template<class Scalar>
const std::string
TimeStepNonlinearSolver<Scalar>::DefaultTol_name_ = "Default Tol";

template<class Scalar>
const double
TimeStepNonlinearSolver<Scalar>::DefaultTol_default_ = 1e-2;


template<class Scalar>
const std::string
TimeStepNonlinearSolver<Scalar>::DefaultMaxIters_name_ = "Default Max Iters";

template<class Scalar>
const int
TimeStepNonlinearSolver<Scalar>::DefaultMaxIters_default_ = 3;


template<class Scalar>
const std::string
TimeStepNonlinearSolver<Scalar>::NonlinearSafetyFactor_name_
= "Nonlinear Safety Factor";

template<class Scalar>
const double
TimeStepNonlinearSolver<Scalar>::NonlinearSafetyFactor_default_ = 0.1;


template<class Scalar>
const std::string
TimeStepNonlinearSolver<Scalar>::LinearSafetyFactor_name_ = "Linear Safety Factor";

template<class Scalar>
const double
TimeStepNonlinearSolver<Scalar>::LinearSafetyFactor_default_ = 0.05;


template<class Scalar>
const std::string
TimeStepNonlinearSolver<Scalar>::RMinFraction_name_ = "R Min Fraction";

template<class Scalar>
const double
TimeStepNonlinearSolver<Scalar>::RMinFraction_default_ = 0.3;


template<class Scalar>
const std::string
TimeStepNonlinearSolver<Scalar>::ThrownOnLinearSolveFailure_name_
= "Thrown on Linear Solve Failure";

template<class Scalar>
const bool
TimeStepNonlinearSolver<Scalar>::ThrownOnLinearSolveFailure_default_ = false;


// Constructors/Intializers/Misc


template <class Scalar>
TimeStepNonlinearSolver<Scalar>::TimeStepNonlinearSolver()
  :J_is_current_(false),
   defaultTol_(DefaultTol_default_),
   defaultMaxIters_(DefaultMaxIters_default_),
   nonlinearSafetyFactor_(NonlinearSafetyFactor_default_),
   linearSafetyFactor_(LinearSafetyFactor_default_),
   RMinFraction_(RMinFraction_default_),
   throwOnLinearSolveFailure_(ThrownOnLinearSolveFailure_default_)
{}


// Overridden from ParameterListAcceptor


template<class Scalar>
void TimeStepNonlinearSolver<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  using Teuchos::get;
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*getValidParameters(),0);
  paramList_ = paramList;
  defaultTol_ = get<double>(*paramList_,DefaultTol_name_);
  defaultMaxIters_ = get<int>(*paramList_,DefaultMaxIters_name_);
  nonlinearSafetyFactor_ = get<double>(*paramList_,NonlinearSafetyFactor_name_);
  linearSafetyFactor_ = get<double>(*paramList_,LinearSafetyFactor_name_);
  RMinFraction_ = get<double>(*paramList_,RMinFraction_name_);
  throwOnLinearSolveFailure_ = get<bool>(
    *paramList_,ThrownOnLinearSolveFailure_name_);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef RYTHMOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // RYTHMOS_DEBUG
}


template<class Scalar>
RCP<ParameterList>
TimeStepNonlinearSolver<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<ParameterList>
TimeStepNonlinearSolver<Scalar>::unsetParameterList()
{
  RCP<ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const ParameterList>
TimeStepNonlinearSolver<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const ParameterList>
TimeStepNonlinearSolver<Scalar>::getValidParameters() const
{
  using Teuchos::setDoubleParameter; using Teuchos::setIntParameter;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    setDoubleParameter(
      DefaultTol_name_, DefaultTol_default_,
      "The default base tolerance for the nonlinear timestep solve.\n"
      "This tolerance can be overridden ???",
      &*pl );
    setIntParameter(
      DefaultMaxIters_name_, DefaultMaxIters_default_,
      "The default maximum number of Newton iterations to perform.\n"
      "This default can be overridden ???",
      &*pl );
    setDoubleParameter(
      NonlinearSafetyFactor_name_, NonlinearSafetyFactor_default_,
      "The factor (< 1.0) to multiply tol to bound R*||dx|||.\n"
      "The exact nonlinear convergence test is:\n"
      "  R*||dx|| <= \"" + NonlinearSafetyFactor_name_ + "\" * tol.",
      &*pl );
    setDoubleParameter(
      LinearSafetyFactor_name_, LinearSafetyFactor_default_,
      "This factor multiplies the nonlinear safety factor which multiplies\n"
      "tol when determining the linear solve tolerence.\n"
      "The exact linear convergence tolerance is:\n"
      "  ||J*dx+f||/||f|| <= \"" + LinearSafetyFactor_name_ + "\" * "
      "\"" + NonlinearSafetyFactor_name_ + "\" * tol.",
      &*pl );
    setDoubleParameter(
      RMinFraction_name_, RMinFraction_default_,
      "The faction below which the R factor is not allowed to drop\n"
      "below each Newton iteration.  The R factor is related to the\n"
      "ratio of ||dx||/||dx_last|| between nonlinear iterations.",
      &*pl );
    pl->set(
      ThrownOnLinearSolveFailure_name_, ThrownOnLinearSolveFailure_default_,
      "If set to true (\"1\"), then an Thyra::CatastrophicSolveFailure\n"
      "exception will be thrown when a linear solve fails to meet it's tolerance."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from NonlinearSolverBase


template <class Scalar>
void TimeStepNonlinearSolver<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
  current_x_ = Teuchos::null;
  J_is_current_ = false;
}


template <class Scalar>
RCP<const Thyra::ModelEvaluator<Scalar> >
TimeStepNonlinearSolver<Scalar>::getModel() const
{
  return model_;
}

template <class Scalar>
Thyra::SolveStatus<Scalar>
TimeStepNonlinearSolver<Scalar>::solve(
  Thyra::VectorBase<Scalar> *x,
  const Thyra::SolveCriteria<Scalar> *solveCriteria,
  Thyra::VectorBase<Scalar> *delta
  )
{
  
#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos:TimeStepNonlinearSolver::solve");
#endif

  using std::endl;
  using Teuchos::incrVerbLevel;
  using Teuchos::describe;
  using Teuchos::as;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  using Teuchos::getFancyOStream;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;
  typedef Thyra::LinearOpWithSolveBase<Scalar> LOWSB;
  typedef Teuchos::VerboseObjectTempState<LOWSB> VOTSLOWSB;

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(0==x);
  THYRA_ASSERT_VEC_SPACES(
    "TimeStepNonlinearSolver<Scalar>::solve(...)",
    *x->space(),*model_->get_x_space() );
  TEST_FOR_EXCEPT(
    0!=solveCriteria && "ToDo: Support passed in solve criteria!" );
#endif

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool showNewtonDetails =
    (!is_null(out) && (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)));
  const bool dumpAll =
    (!is_null(out) && (as<int>(verbLevel) == as<int>(Teuchos::VERB_EXTREME))); 
  TEUCHOS_OSTAB;
  VOTSME stateModel_outputTempState(model_,out,incrVerbLevel(verbLevel,-1));

  if (showNewtonDetails)
    *out
      << "\nEntering TimeStepNonlinearSolver::solve(...) ...\n"
      << "\nmodel = " << Teuchos::describe(*model_,verbLevel);

  if(dumpAll) {
    *out << "\nInitial guess:\n";
    *out << "\nx = " << *x;
  }

  // Initialize storage for algorithm
  if(!J_.get()) J_ = model_->create_W();
  TEST_FOR_EXCEPTION( Teuchos::is_null(J_), std::logic_error,
      "Error!  model->create_W() returned a null pointer!\n"
      );
  RCP<Thyra::VectorBase<Scalar> > f = createMember(model_->get_f_space());
  RCP<Thyra::VectorBase<Scalar> > dx = createMember(model_->get_x_space());
  RCP<Thyra::VectorBase<Scalar> > dx_last = createMember(model_->get_x_space());
  RCP<Thyra::VectorBase<Scalar> > x_curr = createMember(model_->get_x_space());
  if (delta != NULL)
    Thyra::V_S(ptr(delta),ST::zero()); // delta stores the cumulative update to x over the whole Newton solve.
  Thyra::assign(x_curr.ptr(),*x);
  J_is_current_ = false;
  current_x_ = Teuchos::null;

  // Initialize convergence criteria
  ScalarMag R = SMT::one();
  ScalarMag linearTolSafety = linearSafetyFactor_ * nonlinearSafetyFactor_;
  int maxIters = defaultMaxIters_;
  ScalarMag tol = defaultTol_;
  // ToDo: Get above from solveCriteria!

  // Do the undampened Newton iterations
  bool converged = false;
  bool sawFailedLinearSolve = false;
  Thyra::SolveStatus<Scalar> failedLinearSolveStatus;
  ScalarMag nrm_dx = SMT::nan();
  ScalarMag nrm_dx_last = SMT::nan();
  int iter = 1;
  for( ; iter <= maxIters; ++iter ) {
    if (showNewtonDetails)
      *out << "\n*** newtonIter = " << iter << endl;
    if (showNewtonDetails)
      *out << "\nEvaluating the model f and W ...\n";
    Thyra::eval_f_W( *model_, *x_curr, &*f, &*J_ );
    if (showNewtonDetails)
      *out << "\nSolving the system J*dx = -f ...\n";
    Thyra::V_S(dx.ptr(),ST::zero()); // Initial guess is needed!
    Thyra::SolveCriteria<Scalar>
      linearSolveCriteria(
        Thyra::SolveMeasureType(
          Thyra::SOLVE_MEASURE_NORM_RESIDUAL, Thyra::SOLVE_MEASURE_NORM_RHS
          ),
        linearTolSafety*tol
        );
    VOTSLOWSB J_outputTempState(J_,out,incrVerbLevel(verbLevel,-1));
    Thyra::SolveStatus<Scalar> linearSolveStatus
      = J_->solve(Thyra::NOTRANS, *f, dx.ptr(), Teuchos::ptr(&linearSolveCriteria) );
    if (showNewtonDetails)
      *out << "\nLinear solve status:\n" << linearSolveStatus;
    Thyra::Vt_S(dx.ptr(),Scalar(-ST::one()));
    if (dumpAll)
      *out << "\ndx = " << Teuchos::describe(*dx,verbLevel);
    if (delta != NULL) {
      Thyra::Vp_V(ptr(delta),*dx);
      if (dumpAll)
        *out << "\ndelta = " << Teuchos::describe(*delta,verbLevel);
    }
    // Check the linear solve
    if(linearSolveStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) {
      sawFailedLinearSolve = true;
      failedLinearSolveStatus = linearSolveStatus;
      if (throwOnLinearSolveFailure_) {
        TEST_FOR_EXCEPTION(
          throwOnLinearSolveFailure_, Thyra::CatastrophicSolveFailure,
          "Error, the linear solver did not converge!"
          );
      }
      if (showNewtonDetails)
        *out << "\nWarning, linear solve did not converge!  Continuing anyway :-)\n";
    }
    // Update the solution: x_curr = x_curr + dx
    Vp_V( x_curr.ptr(), *dx );
    if (dumpAll)
      *out << "\nUpdated solution x = " << Teuchos::describe(*x_curr,verbLevel);
    // Convergence test
    nrm_dx = Thyra::norm(*dx);
    if ( R*nrm_dx <= nonlinearSafetyFactor_*tol )
      converged = true;
    if (showNewtonDetails)
      *out
        << "\nConvergence test:\n"
        << "  R*||dx|| = " << R << "*" << nrm_dx
        << " = " << (R*nrm_dx) << "\n"
        << "    <= nonlinearSafetyFactor*tol = " << nonlinearSafetyFactor_ << "*" << tol
        << " = " << (nonlinearSafetyFactor_*tol)
        << " : " << ( converged ? "converged!" : " unconverged" )
        << endl;
    if(converged)
      break; // We have converged!!!
    // Update convergence criteria for the next iteration ...
    if(iter > 1) {
      const Scalar
        MinR = RMinFraction_*R,
        nrm_dx_ratio = nrm_dx/nrm_dx_last;
      R = std::max(MinR,nrm_dx_ratio);
      if (showNewtonDetails)
      *out
        << "\nUpdated R\n"
        << "  = max(RMinFraction*R,||dx||/||dx_last||)\n"
        << "  = max("<<RMinFraction_<<"*"<<R<<","<<nrm_dx<<"/"<<nrm_dx_last<<")\n"
        << "  = max("<<MinR<<","<<nrm_dx_ratio<<")\n"
        << "  = " << R << endl;
    }
    // Save to old
    std::swap(dx_last,dx);
    nrm_dx_last = nrm_dx;
  }

  // Set the solution
  Thyra::assign(ptr(x),*x_curr);
  
  if (dumpAll)
    *out << "\nFinal solution x = " << Teuchos::describe(*x,verbLevel);

  // Check the status

  Thyra::SolveStatus<Scalar> solveStatus;

  std::ostringstream oss;
  Teuchos::FancyOStream omsg(rcp(&oss,false));

  omsg << "Solver: " << this->description() << endl;

  if(converged) {
    solveStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
    omsg << "CVODE status test converged!\n";
  }
  else {
    solveStatus.solveStatus = Thyra::SOLVE_STATUS_UNCONVERGED;
    omsg << "CVODE status test failed!\n";
  }

  if (sawFailedLinearSolve) {
    omsg << "Warning!  A failed linear solve was encountered with status:\n";
    OSTab tab(omsg);
    omsg << failedLinearSolveStatus;
  }

  omsg
    << "R*||dx|| = " << R << "*" << nrm_dx
    << " <= nonlinearSafetyFactor*tol = " << nonlinearSafetyFactor_ << "*" << tol << " : "
    << ( converged ? "converged!" : " unconverged" ) << endl;

  omsg
    << "Iterations = " << iter;
  // Above, we leave off the last newline since this is the convention for the
  // SolveStatus::message string!

  solveStatus.message = oss.str();

  // Update the solution state for external clients
  current_x_ = x->clone_v();
  J_is_current_ = false;
  // 2007/09/04: rabartl: Note, above the Jacobian J is always going to be out
  // of date since this algorithm computes x_curr = x_curr + dx for at least
  // one solve for dx = -inv(J)*f.  Therefore, J is never at the updated
  // x_curr, only the old x_curr!

  if (showNewtonDetails)
    *out << "\nLeaving TimeStepNonlinearSolver::solve(...) ...\n";

  return solveStatus;

}


template <class Scalar>
bool TimeStepNonlinearSolver<Scalar>::supportsCloning() const
{
  return true;
}


template <class Scalar>
RCP<Thyra::NonlinearSolverBase<Scalar> >
TimeStepNonlinearSolver<Scalar>::cloneNonlinearSolver() const
{
  RCP<TimeStepNonlinearSolver<Scalar> >
    nonlinearSolver = Teuchos::rcp(new TimeStepNonlinearSolver<Scalar>);
  nonlinearSolver->model_ = model_; // Shallow copy is okay, model is stateless
  nonlinearSolver->defaultTol_ = defaultTol_;
  nonlinearSolver->defaultMaxIters_ = defaultMaxIters_;
  nonlinearSolver->nonlinearSafetyFactor_ = nonlinearSafetyFactor_;
  nonlinearSolver->linearSafetyFactor_ = linearSafetyFactor_;
  nonlinearSolver->RMinFraction_ = RMinFraction_;
  nonlinearSolver->throwOnLinearSolveFailure_ = throwOnLinearSolveFailure_;
  // Note: The specification of this virtual function in the interface class
  // allows us to just copy the algorithm, not the entire state so we are
  // done!
  return nonlinearSolver;
}


template <class Scalar>
RCP<const Thyra::VectorBase<Scalar> >
TimeStepNonlinearSolver<Scalar>::get_current_x() const
{
  return current_x_;
}


template <class Scalar>
bool TimeStepNonlinearSolver<Scalar>::is_W_current() const
{
  return J_is_current_;
}


template <class Scalar>
RCP<Thyra::LinearOpWithSolveBase<Scalar> >
TimeStepNonlinearSolver<Scalar>::get_nonconst_W(const bool forceUpToDate)
{
  if (is_null(J_))
    return Teuchos::null;
  if (forceUpToDate) {
#ifdef RYTHMOS_DEBUG
    TEST_FOR_EXCEPT(is_null(current_x_));
#endif
    Thyra::eval_f_W<Scalar>( *model_, *current_x_, 0, &*J_ );
    J_is_current_ = true;
  }
  return J_;
}


template <class Scalar>
RCP<const Thyra::LinearOpWithSolveBase<Scalar> >
TimeStepNonlinearSolver<Scalar>::get_W() const
{
  return J_;
}


template <class Scalar>
void TimeStepNonlinearSolver<Scalar>::set_W_is_current(bool W_is_current)
{
  J_is_current_ = W_is_current;
}


} // namespace Rythmos


// Nonmember constructors


template <class Scalar>
Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<Scalar> >
Rythmos::timeStepNonlinearSolver()
{
  return Teuchos::rcp(new TimeStepNonlinearSolver<Scalar>);
}


template <class Scalar>
Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<Scalar> >
Rythmos::timeStepNonlinearSolver(const RCP<ParameterList> &pl)
{
  const RCP<Rythmos::TimeStepNonlinearSolver<Scalar> >
    solver = timeStepNonlinearSolver<Scalar>();
  solver->setParameterList(pl);
  return solver;
}


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_INSTANT(SCALAR) \
  \
  template class TimeStepNonlinearSolver< SCALAR >; \
  \
  template RCP<TimeStepNonlinearSolver< SCALAR > > timeStepNonlinearSolver(); \
  \
  template RCP<TimeStepNonlinearSolver<SCALAR > > \
  timeStepNonlinearSolver(const RCP<ParameterList> &pl);



#endif // RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_DEF_HPP
