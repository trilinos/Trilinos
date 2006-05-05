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

#ifndef THYRA_DAMPENED_NEWTON_NONLINEAR_SOLVER_HPP
#define THYRA_DAMPENED_NEWTON_NONLINEAR_SOLVER_HPP

#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Thyra {

/** \brief Simple dampended Newton solver using a Armijo line search :-)
 * 
 * This class derives from <tt>Teuchos::VerboseObject</tt> and therefore will
 * send output to <tt>*this->getOStream()</tt> if
 * <tt>!Teuchos::isNull(this->getOStream())</tt>. The amount of output sent to
 * <tt>*this->getOStream()</tt> depends on the verbosity level returned by
 * <tt>this->getVerbLevel()</tt>:
 * <ul>
 * <li><tt>Teuchos::VERB_DEFAULT</tt>: Same as <tt>Teuchos::VERB_LOW</tt>.
 * <li><tt>Teuchos::VERB_NONE</tt>: Output nothing
 * <li><tt>Teuchos::VERB_LOW</tt>: Ouput only two lines of output for each Newton iteration
 * <li><tt>Teuchos::VERB_MEDIUM</tt>: Output lines for each Newton iteration and line search iteration
 * <li><tt>Teuchos::VERB_HIGH</tt>: Output more details about the Newton and line search iterations (good for basic debugging) 
 * <li><tt>Teuchos::VERB_EXTREME</tt>: Dump all the matrices and vectors that are computed. 
 * </ul>
 *
 * ToDo: Finish documentation.
 */
template <class Scalar>
class DampenedNewtonNonlinearSolver : public NonlinearSolverBase<Scalar> {
public:

  /** \brief. */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief. */
  typedef typename ST::magnitudeType ScalarMag;
  /** \brief. */
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  /** \brief The default solution tolerance. */
   STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, defaultTol )

  /** \brief The default maximum number of iterations. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, defaultMaxNewtonIterations )
  
  /** \brief Set the armijo constant for the line search */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, armijoConstant )
  
  /** \brief Set the maximum number of backtracking line search iterations to take. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxLineSearchIterations )

  /** \brief . */
  DampenedNewtonNonlinearSolver(
    const ScalarMag           defaultTol                   = 1e-2
    ,const int                defaultMaxNewtonIterations   = 1000
    ,const Scalar             armijoConstant               = 1e-4
    ,const int                maxLineSearchIterations      = 20
    );

  /** \brief . */
  static Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidSolveCriteriaExtraParameters();

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
DampenedNewtonNonlinearSolver<Scalar>::DampenedNewtonNonlinearSolver(
  const ScalarMag           defaultTol
  ,const int                defaultMaxNewtonIterations
  ,const Scalar             armijoConstant
  ,const int                maxLineSearchIterations
  )
  :defaultTol_(defaultTol)
  ,defaultMaxNewtonIterations_(defaultMaxNewtonIterations)
  ,armijoConstant_(armijoConstant)
  ,maxLineSearchIterations_(maxLineSearchIterations)
{}

template <class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
DampenedNewtonNonlinearSolver<Scalar>::getValidSolveCriteriaExtraParameters()
{
  static Teuchos::RefCountPtr<const Teuchos::ParameterList> validSolveCriteriaExtraParameters;
  if(!validSolveCriteriaExtraParameters.get()) {
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      paramList = Teuchos::rcp(new Teuchos::ParameterList);
    paramList->set("Max Iters",int(1000));
    validSolveCriteriaExtraParameters = paramList;
  }
  return validSolveCriteriaExtraParameters;
}

// Overridden from NonlinearSolverBase

template <class Scalar>
SolveStatus<Scalar>
DampenedNewtonNonlinearSolver<Scalar>::solve(
  const ModelEvaluator<Scalar>          &model
  ,VectorBase<Scalar>                   *x_inout
  ,const SolveCriteria<Scalar>          *solveCriteria
  ) const
{
  using std::endl;
  // Validate input
  THYRA_ASSERT_VEC_SPACES(
    "DampenedNewtonNonlinearSolver<Scalar>::solve(...)",*x_inout->space(),*model.get_x_space());
  // Get the output stream and verbosity level
  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool showNewtonIters = (verbLevel==Teuchos::VERB_LOW);
  const bool showLineSearchIters = (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM));
  const bool showNewtonDetails = (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH));
  const bool dumpAll = (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_EXTREME)); 
  TEUCHOS_OSTAB;
  if(out.get() && showNewtonIters) *out
    << "\nBeginning dampended Newton solve of model = " << model.description() << "\n\n";
  // Initialize storage for algorithm
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > J     = model.create_W();
  Teuchos::RefCountPtr<VectorBase<Scalar> >            f     = createMember(model.get_f_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >            x     = Teuchos::rcp(x_inout,false);
  Teuchos::RefCountPtr<VectorBase<Scalar> >            dx    = createMember(model.get_x_space());
  Teuchos::RefCountPtr<VectorBase<Scalar> >            x_new = createMember(model.get_x_space());
  // Get convergence criteria
  ScalarMag tol = this->defaultTol();
  int maxIters = this->defaultMaxNewtonIterations();
  if(solveCriteria && !solveCriteria->solveMeasureType.useDefault()) {
    TEST_FOR_EXCEPTION(
      !solveCriteria->solveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS), CatastrophicSolveFailure
      ,"DampenedNewtonNonlinearSolver<Scalar>::solve(...): Error, can only support resudual-based"
      " convergence criteria!");
    tol = solveCriteria->requestedTol;
    if(solveCriteria->extraParameters.get()) {
      solveCriteria->extraParameters->validateParameters(*getValidSolveCriteriaExtraParameters());
      maxIters = solveCriteria->extraParameters->get("Max Iters",int(maxIters));
    }
  }
  if(out.get() && showNewtonDetails)
    *out << "\nCompute the initial starting point ...\n";
  eval_f_W( model, *x, &*f, &*J );
  if(out.get() && dumpAll) {
    *out << "\nInitial starting point:\n";
    *out << "\nx =\n" << *x;
    *out << "\nf =\n" << *f;
    *out << "\nJ =\n" << *J;
  }
  // Peform the Newton iterations
  int newtonIter, num_residual_evals = 1;
  SolveStatus<Scalar> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_UNCONVERGED;
  for( newtonIter = 1; newtonIter <= maxIters; ++newtonIter ) {
    if(out.get() && showNewtonDetails) *out << "\n*** newtonIter = " << newtonIter << endl;
    // Check convergence
    if(out.get() && showNewtonDetails) *out << "\nChecking for convergence ... : ";
    const Scalar phi = scalarProd(*f,*f), sqrt_phi = ST::squareroot(phi); // merit function: phi(f) = <f,f>
    solveStatus.achievedTol = sqrt_phi;
    const bool isConverged = sqrt_phi <= tol;
    if(out.get() && showNewtonDetails) *out
      << "sqrt(phi) = sqrt(<f,f>) = ||f|| = " << sqrt_phi << ( isConverged ? " <= " : " > " ) << "tol = " << tol << endl;
    if(out.get() && showNewtonIters) *out
      << "newton_iter="<<newtonIter<<": Check convergence: ||f|| = "
      << sqrt_phi << ( isConverged ? " <= " : " > " ) << "tol = " << tol << ( isConverged ? ", Converged!!!" : "" ) << endl;
    if(isConverged) {
      if(x_inout != x.get()) assign( x_inout, *x );  // Assign the solution if we have to
      if(out.get() && showNewtonDetails) {
        *out << "\nWe have converged :-)\n"
             << "\n||x||inf = " << norm_inf(*x) << endl;
        if(dumpAll) *out << "\nx =\n" << *x;
        *out << "\nExiting SimpleNewtonSolver::solve(...)\n";
      }
      std::ostringstream oss;
      oss << "Converged! ||f|| = " << sqrt_phi << ", num_newton_iters="<<newtonIter<<", num_residual_evals="<<num_residual_evals<<".";
      solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
      solveStatus.message = oss.str();
      break;
    }
    if(out.get() && showNewtonDetails) *out << "\nWe have to keep going :-(\n";
    // Compute the Jacobian if we have not already
    if(newtonIter > 1) {
      if(out.get() && showNewtonDetails) *out << "\nComputing the Jacobian J at current point ...\n";
      eval_f_W<Scalar>( model, *x, NULL, &*J );
      if(out.get() && dumpAll) *out << "\nJ =\n" << *J;
    }
    // Compute the newton step: dx = -inv(J)*f
    if(out.get() && showNewtonDetails) *out << "\nComputing the Newton step: dx = - inv(J)*f ...\n";
    if(out.get() && showNewtonIters) *out << "newton_iter="<<newtonIter<<": Computing Newton step ...\n";
    assign( &*dx, ST::zero() );       // Initial guess for the linear solve
    Thyra::solve(*J,NOTRANS,*f,&*dx); // Solve: J*dx = f
    Vt_S( &*dx, Scalar(-ST::one()) ); // dx *= -1.0
    if(out.get() && showNewtonDetails) *out << "\n||dx||inf = " << norm_inf(*dx) << endl;
    if(out.get() && dumpAll) *out << "\ndy =\n" << *dx;
    // Perform backtracking armijo line search
    if(out.get() && showNewtonDetails) *out << "\nStarting backtracking line search iterations ...\n";
    if(out.get() && showNewtonIters) *out << "newton_iter="<<newtonIter<<": Starting backtracking line search ...\n";
    const Scalar Dphi = -2.0*phi; // D(phi(x+alpha*dx))/D(alpha) at alpha=0.0 => -2.0*<f,c>: where dx = -inv(J)*f
    Scalar alpha = 1.0; // Try a full step initially since it will eventually be accepted near solution
    int lineSearchIter;
    ++num_residual_evals;
    for( lineSearchIter = 1; lineSearchIter <= maxLineSearchIterations(); ++lineSearchIter, ++num_residual_evals ) {
      TEUCHOS_OSTAB;
      if(out.get() && showNewtonDetails) *out << "\n*** lineSearchIter = " << lineSearchIter << endl;
      // x_new = x + alpha*dx
      assign( &*x_new, *x ); Vp_StV( &*x_new, alpha, *dx );
      if(out.get() && showNewtonDetails) *out << "\n||x_new||inf = " << norm_inf(*x_new) << endl;
      if(out.get() && dumpAll) *out << "\nx_new =\n" << *x_new;
      // Compute the residual at the updated point
      eval_f(model,*x_new,&*f);
      if(out.get() && dumpAll) *out << "\nf_new =\n" << *f;
      const Scalar phi_new = scalarProd(*f,*f), phi_frac = phi + alpha * armijoConstant() * Dphi;
      if(out.get() && showNewtonDetails) *out << "\nphi_new = <f_new,f_new> = " << phi_new << endl;
      if( Teuchos::ScalarTraits<Scalar>::isnaninf(phi_new) ) {
        if(out.get() && showNewtonDetails) *out << "\nphi_new is not a valid number, backtracking (alpha = 0.1*alpha) ...\n";
        alpha *= 0.1;
        continue;
      }
      const bool acceptPoint = (phi_new <= phi_frac);
      if(out.get() && showNewtonDetails) *out
        << "\nphi_new = " << phi_new << ( acceptPoint ? " <= " : " > " )
        << "phi + alpha * eta * Dphi = " << phi << " + " << alpha << " * " << armijoConstant() << " * " << Dphi
        << " = " << phi_frac << endl;
      if(out.get() && (showLineSearchIters || (showNewtonIters && acceptPoint))) *out
        << "newton_iter="<<newtonIter<<", ls_iter="<<lineSearchIter<<" : "
        << "phi(alpha="<<alpha<<") = "<<phi_new<<(acceptPoint ? " <=" : " >")<<" armijo_cord = " << phi_frac << endl;
      if( acceptPoint ) {
        if(out.get() && showNewtonDetails) *out << "\nAccepting the current step with step length alpha = " << alpha << "!\n";
        break;
      }
      if(out.get() && showNewtonDetails) *out << "\nBacktracking (alpha = 0.5*alpha) ...\n";
      alpha *= 0.5;
    }
    // Check for line search failure
    if( lineSearchIter > maxLineSearchIterations() ) {
      std::ostringstream oss;
      oss
        << "lineSearchIter = " << lineSearchIter << " > maxLineSearchIterations = " << maxLineSearchIterations()
        << ": Linear search failure!  Algorithm terminated!";
      solveStatus.message = oss.str();
      if(out.get() && (showNewtonIters || showNewtonDetails)) *out << endl << oss.str() << endl;
      goto exit;
    }
    // Take the Newton step
    std::swap<Teuchos::RefCountPtr<VectorBase<Scalar> > >( x_new, x ); // Now x is current point!
  }
exit:
  if(out.get() && showNewtonIters) *out
    << "\n[Final] newton_iters="<<newtonIter<<", num_residual_evals="<<num_residual_evals<<"\n";
  if(newtonIter > maxIters) {
    std::ostringstream oss;
    oss
      << "newton_iter = " << newtonIter << " > maxIters = " << maxIters
      << ": Newton algorithm terminated!";
    solveStatus.message = oss.str();
    if( out.get() && (showNewtonIters || showNewtonDetails)) *out << endl << oss.str() << endl;
  }
  if(x_inout != x.get()) assign( x_inout, *x ); // Assign the final point
  if(out.get() && showNewtonDetails) *out
    << "\n*** Ending dampended Newton solve." << endl; 
  return solveStatus;
}

} // namespace Thyra

#endif // THYRA_DAMPENED_NEWTON_NONLINEAR_SOLVER_HPP
