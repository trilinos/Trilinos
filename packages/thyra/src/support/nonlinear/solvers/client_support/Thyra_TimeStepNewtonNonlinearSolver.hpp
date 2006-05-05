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
  SolveStatus<Scalar> solve(
    const ModelEvaluator<Scalar>          &model
    ,VectorBase<Scalar>                   *x
    ,const SolveCriteria<Scalar>          *solveCriteria
    ) const;

  //@}

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
{}

// Overridden from NonlinearSolverBase

template <class Scalar>
SolveStatus<Scalar> TimeStepNewtonNonlinearSolver<Scalar>::solve(
  const ModelEvaluator<Scalar>          &model
  ,VectorBase<Scalar>                   *x
  ,const SolveCriteria<Scalar>          *solveCriteria
  ) const
{
  using std::endl;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  using Teuchos::getFancyOStream;
  TEST_FOR_EXCEPT(solveCriteria!=NULL); // ToDo: Pass to linear solver?
  // Initialize storage for algorithm
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > J = model.create_W();
  Teuchos::RefCountPtr<VectorBase<Scalar> >            f = createMember(model.get_f_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >            dx = createMember(model.get_x_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >            dx_last = createMember(model.get_x_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >            x_curr = createMember(model.get_x_space());
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
    eval_f_W( model, *x_curr, ST::one(), &*f, &*J );
    // Solve the system: J*dx = -f
    SolveCriteria<Scalar> linearSolveCriteria(SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS),linearTolSafety*tol);
    SolveStatus<Scalar>
      linearSolveStatus = Thyra::solve( *J, NOTRANS, *f, &*dx, &linearSolveCriteria );
    Thyra::Vt_S(&*dx,Scalar(-ST::one()));
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
  return solveStatus;
}

} // namespace Thyra

#endif // THYRA_NEWTON_NONLINEAR_SOLVER_HPP
