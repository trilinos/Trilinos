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


#ifndef RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_HPP
#define RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_HPP


#include "Thyra_NonlinearSolverBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"


namespace Rythmos {


/** \brief Simple undampended Newton solver designed to solve time step
 * equations in accurate times-tepping methods.
 * 
 * ToDo: Finish documentation.
 *
 * 2007/05/18: rabartl: ToDo: Derive NonlinearSolverBase from
 * Teuchos::ParameterListAcceptor and accept options through a validated
 * parameter list!  Then remove these STANDARD_MEMBER_COMPOSITION_MEMBERS()
 * macros.
 */
template <class Scalar>
class TimeStepNonlinearSolver : public Thyra::NonlinearSolverBase<Scalar> {
public:

  /** \brief. */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief. */
  typedef typename ST::magnitudeType ScalarMag;
  /** \brief. */
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  /** \brief The default maximum number of iterations. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, defaultMaxIterations );

  /** \brief The default solution tolerance. */
   STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, defaultTol );

  /** \brief Stream that warnings are printed to. */
   STANDARD_NONCONST_COMPOSITION_MEMBERS( std::ostream, warningOut );

  /** \brief . */
  TimeStepNonlinearSolver(
    const int defaultMaxIterations = 3,
    const ScalarMag defaultTol = 1e-2,
    const warningOut_ptr_t &warningOut = Teuchos::rcp(&std::cerr,false)
    );

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(
    const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  Thyra::SolveStatus<Scalar> solve(
    Thyra::VectorBase<Scalar> *x,
    const Thyra::SolveCriteria<Scalar> *solveCriteria,
    Thyra::VectorBase<Scalar> *delta = NULL
    );
  /** \brief . */
  bool supportsCloning() const;
  /** \brief . */
  Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> >
  cloneNonlinearSolver() const;  
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_current_x() const;
  /** \brief . */
  bool is_W_current() const;
  /** \brief . */
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> > get_nonconst_W();
  /** \brief . */
  Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> > get_W() const;
  /** \brief . */
  void set_W_is_current(bool W_is_current);

  //@}

private:

  Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> > J_;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > current_x_;
  bool J_is_current_;

};


// ////////////////////////
// Defintions


template <class Scalar>
TimeStepNonlinearSolver<Scalar>::TimeStepNonlinearSolver(
  const int defaultMaxIterations,
  const ScalarMag defaultTol,
  const warningOut_ptr_t &warningOut
  )
  :defaultMaxIterations_(defaultMaxIterations),
   defaultTol_(defaultTol),
   warningOut_(warningOut),
   J_is_current_(false)
{}


// Overridden from NonlinearSolverBase


template <class Scalar>
void TimeStepNonlinearSolver<Scalar>::setModel(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
  current_x_ = Teuchos::null;
  J_is_current_ = false;
}


template <class Scalar>
Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >
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

  using std::endl;
  using Teuchos::rcp;
  using Teuchos::RefCountPtr;
  using Teuchos::OSTab;
  using Teuchos::getFancyOStream;

  TEST_FOR_EXCEPT(solveCriteria!=NULL); // ToDo: Pass to linear solver?

  // Initialize storage for algorithm
  if(!J_.get()) J_ = model_->create_W();
  RefCountPtr<Thyra::VectorBase<Scalar> > f = createMember(model_->get_f_space());
  RefCountPtr<Thyra::VectorBase<Scalar> > dx = createMember(model_->get_x_space());
  RefCountPtr<Thyra::VectorBase<Scalar> > dx_last = createMember(model_->get_x_space());
  RefCountPtr<Thyra::VectorBase<Scalar> > x_curr = createMember(model_->get_x_space());
  if (delta != NULL)
    Thyra::V_S(delta,ST::zero()); // delta stores the cumulative update to x over the whole Newton solve.
  Thyra::assign(&*x_curr,*x);

  // Initialize convergence criteria
  ScalarMag R = SMT::one();
  ScalarMag RMin = 0.3; // ToDo: Make this adjustable
  ScalarMag tolSafety = 0.1; // ToDo: Make this adjustable
  ScalarMag linearTolSafety = 0.05 * tolSafety; // ToDo: Make this adjustable
  int maxIters = this->defaultMaxIterations();
  ScalarMag tol = this->defaultTol();
  // ToDo: Get above from solveCriteria!

  // Do the undampened Newton iterations
  bool converged = false;
  ScalarMag nrm_dx_last;
  int iter = 1;
  for( ; iter <= maxIters; ++iter ) {
    // Evaluate f and W
    eval_f_W( *model_, *x_curr, ST::one(), &*f, &*J_ );
    // Zero out dx to deal with iterative linear solvers
    Thyra::V_S(&*dx,ST::zero());
    // Solve the system: J*dx = -f
    Thyra::SolveCriteria<Scalar>
      linearSolveCriteria(
        Thyra::SolveMeasureType(
          Thyra::SOLVE_MEASURE_NORM_RESIDUAL,Thyra::SOLVE_MEASURE_NORM_RHS
          )
        ,linearTolSafety*tol
        );
    Thyra::SolveStatus<Scalar> linearSolveStatus
      = Thyra::solve( *J_, Thyra::NOTRANS, *f, &*dx, &linearSolveCriteria );
    Thyra::Vt_S(&*dx,Scalar(-ST::one()));
    if (delta != NULL)
      Thyra::Vp_V(delta,*dx); 
    // Check the linear solve
    if(linearSolveStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) {
      warningOut()
        << "Rythmos::TimeStepNonlinearSolver<Scalar>::solve(...), Warning, linear solve did not converge with solve status:\n\n";
      OSTab(getFancyOStream(rcp(&warningOut(),false))).o() << linearSolveStatus;
      warningOut()
        <<endl<< "Continuing anyway :-)\n";
      // ToDo: Add option to throw exception failure
    }
    // Update the solution: x_curr = x_curr + dx
    Vp_V( &*x_curr, *dx );
    // Convergence test
    const ScalarMag nrm_dx = Thyra::norm(*dx);
    if(R*nrm_dx < tolSafety*tol) {
      converged = true;
      break;
    }
    // Update convergence criteria
    if(iter > 1) {
      R = std::max(RMin*R,nrm_dx/nrm_dx_last);
    }
    // Save to old
    std::swap(dx_last,dx);
    nrm_dx_last = nrm_dx;
  }

  // Set the solution
  Thyra::assign(x,*x_curr);

  // Check the status
  Thyra::SolveStatus<Scalar> solveStatus;
  std::ostringstream oss;
  if(converged) {
    solveStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
    oss << "CVODE status test converged!  Iterations = " << iter << ".";
  }
  else {
    solveStatus.solveStatus = Thyra::SOLVE_STATUS_UNCONVERGED;
    oss << "CVODE status test failed!  Iterations = " << iter << ".";
  }
  solveStatus.message = oss.str();

  // Update the solution state for external clients
  current_x_ = x->clone_v();
  J_is_current_ = true;

  return solveStatus;

}


template <class Scalar>
bool TimeStepNonlinearSolver<Scalar>::supportsCloning() const
{
  return true;
}


template <class Scalar>
Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> >
TimeStepNonlinearSolver<Scalar>::cloneNonlinearSolver() const
{
  Teuchos::RefCountPtr<TimeStepNonlinearSolver<Scalar> >
    nonlinearSolver = Teuchos::rcp(new TimeStepNonlinearSolver<Scalar>);
  nonlinearSolver->defaultMaxIterations_ = defaultMaxIterations_;
  nonlinearSolver->defaultTol_ = defaultTol_;
  nonlinearSolver->warningOut_ = warningOut_;
  nonlinearSolver->model_ = model_; // Shallow copy is okay, model is stateless
  // Note: The specification of this virtual function in the interface class
  // allows us to just copy the algorithm, not the entire state so we are
  // done!
  return nonlinearSolver;
}


template <class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
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
Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<Scalar> >
TimeStepNonlinearSolver<Scalar>::get_nonconst_W()
{
  return J_;
}


template <class Scalar>
Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
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


#endif // RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_HPP
