// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DAMPENED_NEWTON_NONLINEAR_SOLVER_HPP
#define THYRA_DAMPENED_NEWTON_NONLINEAR_SOLVER_HPP

#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_as.hpp"


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
 *
 * \ingroup Thyra_Nonlin_ME_solvers_grp
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
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, defaultTol );

  /** \brief The default maximum number of iterations. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, defaultMaxNewtonIterations );

  /** \brief The default maximum number of iterations. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, useDampenedLineSearch  );
  
  /** \brief Set the armijo constant for the line search */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, armijoConstant );
  
  /** \brief Set the maximum number of backtracking line search iterations to take. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxLineSearchIterations );

  /** \brief . */
  DampenedNewtonNonlinearSolver(
    const ScalarMag defaultTol = 1e-2
    ,const int defaultMaxNewtonIterations = 1000
    ,const bool useDampenedLineSearch = true
    ,const Scalar armijoConstant = 1e-4
    ,const int maxLineSearchIterations = 20
    );

  /** \brief . */
  static RCP<const Teuchos::ParameterList>
  getValidSolveCriteriaExtraParameters();

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(
    const RCP<const ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  RCP<const ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  SolveStatus<Scalar> solve(
    VectorBase<Scalar> *x,
    const SolveCriteria<Scalar> *solveCriteria,
    VectorBase<Scalar> *delta
    );
  /** \brief . */
  RCP<const VectorBase<Scalar> > get_current_x() const;
  /** \brief . */
  bool is_W_current() const;
  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> > get_nonconst_W(const bool forceUpToDate);
  /** \brief . */
  RCP<const LinearOpWithSolveBase<Scalar> > get_W() const;
  /** \brief . */
  void set_W_is_current(bool W_is_current);

  //@}

private:

  RCP<Teuchos::ParameterList> paramList_;
  RCP<const ModelEvaluator<Scalar> > model_;
  RCP<LinearOpWithSolveBase<Scalar> > J_;
  RCP<VectorBase<Scalar> > current_x_;
  bool J_is_current_;

};

// ////////////////////////
// Defintions

template <class Scalar>
DampenedNewtonNonlinearSolver<Scalar>::DampenedNewtonNonlinearSolver(
  const ScalarMag my_defaultTol
  ,const int my_defaultMaxNewtonIterations
  ,const bool my_useDampenedLineSearch
  ,const Scalar my_armijoConstant
  ,const int my_maxLineSearchIterations
  )
  :defaultTol_(my_defaultTol)
  ,defaultMaxNewtonIterations_(my_defaultMaxNewtonIterations)
  ,useDampenedLineSearch_(my_useDampenedLineSearch)
  ,armijoConstant_(my_armijoConstant)
  ,maxLineSearchIterations_(my_maxLineSearchIterations)
  ,J_is_current_(false)
{}

template <class Scalar>
RCP<const Teuchos::ParameterList>
DampenedNewtonNonlinearSolver<Scalar>::getValidSolveCriteriaExtraParameters()
{
  static RCP<const Teuchos::ParameterList> validSolveCriteriaExtraParameters;
  if(!validSolveCriteriaExtraParameters.get()) {
    RCP<Teuchos::ParameterList>
      paramList = Teuchos::rcp(new Teuchos::ParameterList);
    paramList->set("Max Iters",int(1000));
    validSolveCriteriaExtraParameters = paramList;
  }
  return validSolveCriteriaExtraParameters;
}

// Overridden from Teuchos::ParameterListAcceptor

template<class Scalar>
void DampenedNewtonNonlinearSolver<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  using Teuchos::get;
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*getValidParameters(),0);
  paramList_ = paramList;
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement!");
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG
}

template<class Scalar>
RCP<Teuchos::ParameterList>
DampenedNewtonNonlinearSolver<Scalar>::getNonconstParameterList()
{
  return paramList_;
}

template<class Scalar>
RCP<Teuchos::ParameterList>
DampenedNewtonNonlinearSolver<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
DampenedNewtonNonlinearSolver<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
DampenedNewtonNonlinearSolver<Scalar>::getValidParameters() const
{
  using Teuchos::setDoubleParameter; using Teuchos::setIntParameter;
  static RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();
    TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement!");
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}

// Overridden from NonlinearSolverBase

template <class Scalar>
void DampenedNewtonNonlinearSolver<Scalar>::setModel(
  const RCP<const ModelEvaluator<Scalar> > &model
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
  current_x_ = Teuchos::null;
  J_is_current_ = false;
}

template <class Scalar>
RCP<const ModelEvaluator<Scalar> >
DampenedNewtonNonlinearSolver<Scalar>::getModel() const
{
  return model_;
}

template <class Scalar>
SolveStatus<Scalar>
DampenedNewtonNonlinearSolver<Scalar>::solve(
  VectorBase<Scalar> *x_inout
  ,const SolveCriteria<Scalar> *solveCriteria
  ,VectorBase<Scalar> *delta
  ) 
{

  using std::endl;
  using Teuchos::as;

  // Validate input
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0==x_inout);
  THYRA_ASSERT_VEC_SPACES(
    "DampenedNewtonNonlinearSolver<Scalar>::solve(...)",
    *x_inout->space(), *model_->get_x_space() );
#endif

  // Get the output stream and verbosity level
  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool showNewtonIters = (verbLevel==Teuchos::VERB_LOW);
  const bool showLineSearchIters = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM));
  const bool showNewtonDetails = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH));
  const bool dumpAll = (as<int>(verbLevel) == as<int>(Teuchos::VERB_EXTREME)); 
  TEUCHOS_OSTAB;
  if(out.get() && showNewtonIters) {
    *out << "\nBeginning dampended Newton solve of model = " << model_->description() << "\n\n";
    if (!useDampenedLineSearch())
      *out << "\nDoing undampened newton ...\n\n";
  }

  // Initialize storage for algorithm
  if(!J_.get()) J_ = model_->create_W();
  RCP<VectorBase<Scalar> > f = createMember(model_->get_f_space());
  RCP<VectorBase<Scalar> > x = Teuchos::rcp(x_inout,false);
  RCP<VectorBase<Scalar> > dx = createMember(model_->get_x_space());
  RCP<VectorBase<Scalar> > x_new = createMember(model_->get_x_space());
  RCP<VectorBase<Scalar> > ee = createMember(model_->get_x_space());
  V_S(ee.ptr(),ST::zero());

  // Get convergence criteria
  ScalarMag tol = this->defaultTol();
  int maxIters = this->defaultMaxNewtonIterations();
  if(solveCriteria && !solveCriteria->solveMeasureType.useDefault()) {
    TEUCHOS_TEST_FOR_EXCEPTION(
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

  eval_f_W( *model_, *x, &*f, &*J_ );
  if(out.get() && dumpAll) {
    *out << "\nInitial starting point:\n";
    *out << "\nx =\n" << *x;
    *out << "\nf =\n" << *f;
    *out << "\nJ =\n" << *J_;
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
      if(x_inout != x.get()) assign( ptr(x_inout), *x ); // Assign the solution if we have to
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
      if(out.get() && showNewtonDetails) *out << "\nComputing the Jacobian J_ at current point ...\n";
      eval_f_W<Scalar>( *model_, *x, NULL, &*J_ );
      if(out.get() && dumpAll) *out << "\nJ =\n" << *J_;
    }

    // Compute the newton step: dx = -inv(J)*f
    if(out.get() && showNewtonDetails) *out << "\nComputing the Newton step: dx = - inv(J)*f ...\n";
    if(out.get() && showNewtonIters) *out << "newton_iter="<<newtonIter<<": Computing Newton step ...\n";
    assign( dx.ptr(), ST::zero() ); // Initial guess for the linear solve
    J_->solve(NOTRANS,*f,dx.ptr()); // Solve: J*dx = f
    Vt_S( dx.ptr(), Scalar(-ST::one()) ); // dx *= -1.0
    Vp_V( ee.ptr(), *dx); // ee += dx
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
      TEUCHOS_OSTAB_DIFF(lineSearchIter);
      if(out.get() && showNewtonDetails) *out << "\n*** lineSearchIter = " << lineSearchIter << endl;
      // x_new = x + alpha*dx
      assign( x_new.ptr(), *x ); Vp_StV( x_new.ptr(), alpha, *dx );
      if(out.get() && showNewtonDetails) *out << "\n||x_new||inf = " << norm_inf(*x_new) << endl;
      if(out.get() && dumpAll) *out << "\nx_new =\n" << *x_new;
      // Compute the residual at the updated point
      eval_f(*model_,*x_new,&*f);
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
      if (out.get() && showNewtonDetails && !useDampenedLineSearch())
        *out << "\nUndamped newton, always accpeting the point!\n";
      if( acceptPoint || !useDampenedLineSearch() ) {
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
        << ": Linear search failure! Algorithm terminated!";
      solveStatus.message = oss.str();
      if(out.get() && (showNewtonIters || showNewtonDetails)) *out << endl << oss.str() << endl;
      goto exit;
    }

    // Take the Newton step
    std::swap<RCP<VectorBase<Scalar> > >( x_new, x ); // Now x is current point!

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

  if(x_inout != x.get()) assign( ptr(x_inout), *x ); // Assign the final point
  if(delta != NULL) assign( ptr(delta), *ee );
  current_x_ = x_inout->clone_v(); // Remember the final point
  J_is_current_ = newtonIter==1; // J is only current with x if initial point was converged!

  if(out.get() && showNewtonDetails) *out
    << "\n*** Ending dampended Newton solve." << endl; 

  return solveStatus;

}

template <class Scalar>
RCP<const VectorBase<Scalar> >
DampenedNewtonNonlinearSolver<Scalar>::get_current_x() const
{
  return current_x_;
}

template <class Scalar>
bool DampenedNewtonNonlinearSolver<Scalar>::is_W_current() const
{
  return J_is_current_;
}

template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
DampenedNewtonNonlinearSolver<Scalar>::get_nonconst_W(const bool forceUpToDate)
{
  if (forceUpToDate) {
    TEUCHOS_TEST_FOR_EXCEPT(forceUpToDate);
  }
  return J_;
}

template <class Scalar>
RCP<const LinearOpWithSolveBase<Scalar> >
DampenedNewtonNonlinearSolver<Scalar>::get_W() const
{
  return J_;
}

template <class Scalar>
void DampenedNewtonNonlinearSolver<Scalar>::set_W_is_current(bool W_is_current)
{
  J_is_current_ = W_is_current;
}


} // namespace Thyra


#endif // THYRA_DAMPENED_NEWTON_NONLINEAR_SOLVER_HPP
