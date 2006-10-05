// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_NEWTON_NONLINEAR_SOLVER_HPP
#define THYRA_NEWTON_NONLINEAR_SOLVER_HPP

#include "Thyra_NonlinearSolverBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Thyra {

/** \brief Simple undampended Newton solver :-)
 * 
 * ToDo: Finish documentation.
 */
template <class Scalar>
class TimeStepNewtonNonlinearSolver : public NonlinearSolverBase<Scalar> {
public:

  /** \brief. */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief. */
  typedef typename ST::magnitudeType ScalarMag;
  /** \brief. */
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  /** \brief The default maximum number of iterations. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, defaultMaxIterations )

  /** \brief The default solution tolerance. */
   STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, defaultTol )

  /** \brief Stream that warnings are printed to. */
   STANDARD_NONCONST_COMPOSITION_MEMBERS( std::ostream, warningOut )

  TimeStepNewtonNonlinearSolver(
    const int                 defaultMaxIterations = 3
    ,const ScalarMag          defaultTol           = 1e-2
    ,const warningOut_ptr_t   &warningOut          = Teuchos::rcp(&std::cerr,false)
    );

  /** @name Overridden from NonlinearSolverBase */
  //@{


  /** \brief . */
  void setModel(
    const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  SolveStatus<Scalar> solve(
    VectorBase<Scalar>              *x
    ,const SolveCriteria<Scalar>    *solveCriteria
    ,VectorBase<Scalar>             *delta = NULL
    );
  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > get_current_x() const;
  /** \brief . */
  bool is_W_current() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > get_nonconst_W();
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > get_W() const;
  /** \brief . */
  void set_W_is_current(bool W_is_current);

  //@}

private:

  Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >    model_;
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >   J_;
  Teuchos::RefCountPtr<VectorBase<Scalar> >              current_x_;
  bool                                                   J_is_current_;

};

// ////////////////////////
// Defintions

template <class Scalar>
TimeStepNewtonNonlinearSolver<Scalar>::TimeStepNewtonNonlinearSolver(
  const int                 defaultMaxIterations
  ,const ScalarMag          defaultTol
  ,const warningOut_ptr_t   &warningOut
  )
  :defaultMaxIterations_(defaultMaxIterations)
  ,defaultTol_(defaultTol)
  ,warningOut_(warningOut)
  ,J_is_current_(false)
{}

// Overridden from NonlinearSolverBase

template <class Scalar>
void TimeStepNewtonNonlinearSolver<Scalar>::setModel(
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
  current_x_ = Teuchos::null;
  J_is_current_ = false;
}

template <class Scalar>
Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
TimeStepNewtonNonlinearSolver<Scalar>::getModel() const
{
  return model_;
}

template <class Scalar>
SolveStatus<Scalar> TimeStepNewtonNonlinearSolver<Scalar>::solve(
  VectorBase<Scalar>             *x
  ,const SolveCriteria<Scalar>   *solveCriteria
  ,VectorBase<Scalar>            *delta
  )
{
  using std::endl;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  using Teuchos::getFancyOStream;
  TEST_FOR_EXCEPT(solveCriteria!=NULL); // ToDo: Pass to linear solver?
  // Initialize storage for algorithm
  if(!J_.get())                                 J_ = model_->create_W();
  Teuchos::RefCountPtr<VectorBase<Scalar> >     f = createMember(model_->get_f_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >     dx = createMember(model_->get_x_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >     dx_last = createMember(model_->get_x_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >     x_curr = createMember(model_->get_x_space());
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
  // ToDo: Get these from solveCriteria
  bool converged = false;
  ScalarMag nrm_dx_last;
  int iter = 1;
  for( ; iter <= maxIters; ++iter ) {
    // Evaluate f and W
    eval_f_W( *model_, *x_curr, ST::one(), &*f, &*J_ );
    // Zero out dx to deal with iterative linear solvers
    Thyra::V_S(&*dx,ST::zero());
    // Solve the system: J*dx = -f
    SolveCriteria<Scalar>
      linearSolveCriteria(
        SolveMeasureType(
          SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS
          )
        ,linearTolSafety*tol
        );
    SolveStatus<Scalar>
      linearSolveStatus = Thyra::solve( *J_, NOTRANS, *f, &*dx, &linearSolveCriteria );
    Thyra::Vt_S(&*dx,Scalar(-ST::one()));
    if (delta != NULL)
      Thyra::Vp_V(delta,*dx); 
    // Check the linear solve
    if(linearSolveStatus.solveStatus != SOLVE_STATUS_CONVERGED) {
      warningOut()
        << "Thyra::TimeStepNewtonNonlinearSolver<Scalar>::solve(...), Warning, linear solve did not converge with solve status:\n\n";
      *OSTab(getFancyOStream(rcp(&warningOut(),false))).getOStream() << linearSolveStatus;
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
  SolveStatus<Scalar> solveStatus;
  std::ostringstream oss;
  if(converged) {
    solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
    oss << "CVODE status test converged!  Iterations = " << iter << ".";
  }
  else {
    solveStatus.solveStatus = SOLVE_STATUS_UNCONVERGED;
    oss << "CVODE status test failed!  Iterations = " << iter << ".";
  }
  solveStatus.message = oss.str();
  //
  current_x_ = x->clone_v();
  J_is_current_ = true;
  //
  return solveStatus;
}

template <class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
TimeStepNewtonNonlinearSolver<Scalar>::get_current_x() const
{
  return current_x_;
}

template <class Scalar>
bool TimeStepNewtonNonlinearSolver<Scalar>::is_W_current() const
{
  return J_is_current_;
}

template <class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
TimeStepNewtonNonlinearSolver<Scalar>::get_nonconst_W()
{
  return J_;
}

template <class Scalar>
Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> >
TimeStepNewtonNonlinearSolver<Scalar>::get_W() const
{
  return J_;
}

template <class Scalar>
void TimeStepNewtonNonlinearSolver<Scalar>::set_W_is_current(bool W_is_current)
{
  J_is_current_ = W_is_current;
}

} // namespace Thyra

#endif // THYRA_NEWTON_NONLINEAR_SOLVER_HPP
