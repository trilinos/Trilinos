/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef OPTIPACK_NONLINEAR_CG_DEF_HPP
#define OPTIPACK_NONLINEAR_CG_DEF_HPP


#include "OptiPack_NonlinearCG_decl.hpp"
#include "OptiPack_DefaultPolyLineSearchPointEvaluator.hpp"
#include "OptiPack_UnconstrainedOptMeritFunc1D.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Tuple.hpp"


namespace OptiPack {


template<typename Scalar>
RCP<Teuchos::ParameterEntryValidator>
NonlinearCG<Scalar>::solverType_validator_ = Teuchos::null;


// Constructor/Initializers/Accessors


template<typename Scalar>
NonlinearCG<Scalar>::NonlinearCG()
  : paramIndex_(-1),
    responseIndex_(-1),
    solverType_(NonlinearCGUtils::solverType_default_integral_val),
    alpha_init_(NonlinearCGUtils::alpha_init_default),
    alpha_reinit_(NonlinearCGUtils::alpha_reinit_default),
    and_conv_tests_(NonlinearCGUtils::and_conv_tests_default),
    minIters_(NonlinearCGUtils::minIters_default),
    maxIters_(NonlinearCGUtils::maxIters_default),
    g_reduct_tol_(NonlinearCGUtils::g_reduct_tol_default),
    g_grad_tol_(NonlinearCGUtils::g_grad_tol_default),
    g_mag_(NonlinearCGUtils::g_mag_default),
    numIters_(0)
{}


template<typename Scalar>
void NonlinearCG<Scalar>::initialize(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const int paramIndex,
  const int responseIndex,
  const RCP<GlobiPack::LineSearchBase<Scalar> > &linesearch
  )
{
  // ToDo: Validate input objects!
  model_ = model.assert_not_null();
  paramIndex_ = paramIndex;
  responseIndex_ = responseIndex;
  linesearch_ = linesearch.assert_not_null();
}


template<typename Scalar>
NonlinearCGUtils::ESolverTypes NonlinearCG<Scalar>::get_solverType() const
{
  return solverType_;
}


template<typename Scalar>
typename NonlinearCG<Scalar>::ScalarMag
NonlinearCG<Scalar>::get_alpha_init() const
{
  return alpha_init_;
}


template<typename Scalar>
bool NonlinearCG<Scalar>::get_alpha_reinit() const
{
  return alpha_reinit_;
}


template<typename Scalar>
bool NonlinearCG<Scalar>::get_and_conv_tests() const
{
  return and_conv_tests_;
}


template<typename Scalar>
int NonlinearCG<Scalar>::get_minIters() const
{
  return minIters_;
}


template<typename Scalar>
int NonlinearCG<Scalar>::get_maxIters() const
{
  return maxIters_;
}


template<typename Scalar>
typename NonlinearCG<Scalar>::ScalarMag
NonlinearCG<Scalar>::get_g_reduct_tol() const
{
  return g_reduct_tol_;
}


template<typename Scalar>
typename NonlinearCG<Scalar>::ScalarMag
NonlinearCG<Scalar>::get_g_grad_tol() const
{
  return g_grad_tol_;
}


template<typename Scalar>
typename NonlinearCG<Scalar>::ScalarMag
NonlinearCG<Scalar>::get_g_mag() const
{
  return g_mag_;
}


// Overridden from ParameterListAcceptor (simple forwarding functions)


template<typename Scalar>
void NonlinearCG<Scalar>::setParameterList(RCP<ParameterList> const& paramList)
{
  typedef ScalarTraits<Scalar> ST;
  typedef ScalarTraits<ScalarMag> SMT;
  namespace NCGU = NonlinearCGUtils;
  using Teuchos::getParameter;
  using Teuchos::getIntegralValue;
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  solverType_ = getIntegralValue<NCGU::ESolverTypes>(*paramList, NCGU::solverType_name);
  alpha_init_ = getParameter<double>(*paramList, NCGU::alpha_init_name);
  alpha_reinit_ = getParameter<bool>(*paramList, NCGU::alpha_reinit_name);
  and_conv_tests_ = getParameter<bool>(*paramList, NCGU::and_conv_tests_name);
  minIters_ = getParameter<int>(*paramList, NCGU::minIters_name);
  maxIters_ = getParameter<int>(*paramList, NCGU::maxIters_name);
  g_reduct_tol_ = getParameter<double>(*paramList, NCGU::g_reduct_tol_name);
  g_grad_tol_ = getParameter<double>(*paramList, NCGU::g_grad_tol_name);
  g_mag_ = getParameter<double>(*paramList, NCGU::g_mag_name);
  TEUCHOS_ASSERT_INEQUALITY( alpha_init_, >, SMT::zero() );
  TEUCHOS_ASSERT_INEQUALITY( minIters_, >=, 0 );
  TEUCHOS_ASSERT_INEQUALITY( minIters_, <, maxIters_ );
  TEUCHOS_ASSERT_INEQUALITY( g_reduct_tol_, >=, SMT::zero() );
  TEUCHOS_ASSERT_INEQUALITY( g_grad_tol_, >=, SMT::zero() );
  TEUCHOS_ASSERT_INEQUALITY( g_mag_, >, SMT::zero() );
  Teuchos::readVerboseObjectSublist(&*paramList, this);
  setMyParamList(paramList);
}


template<typename Scalar>
RCP<const ParameterList>
NonlinearCG<Scalar>::getValidParameters() const
{
  using Teuchos::tuple;
  namespace NCGU = NonlinearCGUtils;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::rcp(new Teuchos::ParameterList());
    solverType_validator_ = 
      Teuchos::stringToIntegralParameterEntryValidator<NCGU::ESolverTypes>(
        tuple<std::string>(
          "FR",
          "PR+",
          "FR-PR",
          "HS"
          ),
        tuple<std::string>(
          "Fletcher-Reeves Method",
          "Polak-Ribiere Method",
          "Fletcher-Reeves Polak-Ribiere Hybrid Method",
          "Hestenes-Stiefel Method"
          ),
        tuple<NCGU::ESolverTypes>(
          NCGU::NONLINEAR_CG_FR,
          NCGU::NONLINEAR_CG_PR_PLUS,
          NCGU::NONLINEAR_CG_FR_PR,
          NCGU::NONLINEAR_CG_HS
          ),
        NCGU::solverType_name
        );
    pl->set( NCGU::solverType_name, NCGU::solverType_default,
      "Set the type of nonlinear CG solver algorithm to use.",
      solverType_validator_ );
    pl->set( NCGU::alpha_init_name, NCGU::alpha_init_default );
    pl->set( NCGU::alpha_reinit_name, NCGU::alpha_reinit_default );
    pl->set( NCGU::and_conv_tests_name, NCGU::and_conv_tests_default );
    pl->set( NCGU::minIters_name, NCGU::minIters_default );
    pl->set( NCGU::maxIters_name, NCGU::maxIters_default );
    pl->set( NCGU::g_reduct_tol_name, NCGU::g_reduct_tol_default );
    pl->set( NCGU::g_grad_tol_name, NCGU::g_grad_tol_default );
    pl->set( NCGU::g_mag_name, NCGU::g_mag_default );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
    // ToDo: Add documentation for these parameters
  }
  return validPL;
}


// Solve


template<typename Scalar>
NonlinearCGUtils::ESolveReturn
NonlinearCG<Scalar>::doSolve(
  const Ptr<Thyra::VectorBase<Scalar> > &p_inout,
  const Ptr<ScalarMag> &g_opt_out,
  const Ptr<const ScalarMag> &g_reduct_tol_in,
  const Ptr<const ScalarMag> &g_grad_tol_in,
  const Ptr<const ScalarMag> &alpha_init_in,
  const Ptr<int> &numIters_out
  )
{

  typedef ScalarTraits<Scalar> ST;
  typedef ScalarTraits<ScalarMag> SMT;
  
  using Teuchos::null;
  using Teuchos::as;
  using Teuchos::tuple;
  using Teuchos::rcpFromPtr;
  using Teuchos::optInArg;
  using Teuchos::inOutArg;
  using GlobiPack::computeValue;
  using GlobiPack::PointEval1D;
  using Thyra::VectorSpaceBase;
  using Thyra::VectorBase;
  using Thyra::MultiVectorBase;
  using Thyra::scalarProd;
  using Thyra::createMember;
  using Thyra::createMembers;
  using Thyra::get_ele;
  using Thyra::norm;
  using Thyra::V_StV;
  using Thyra::Vt_S;
  using Thyra::eval_g_DgDp;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  namespace NCGU = NonlinearCGUtils;
  using std::max;

  // Validate input

  g_opt_out.assert_not_null();

  // Set streams

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  linesearch_->setOStream(out);

  // Determine what step constants will be computed

  const bool compute_beta_PR =
    (
      solverType_ == NCGU::NONLINEAR_CG_PR_PLUS
      ||
      solverType_ == NCGU::NONLINEAR_CG_FR_PR
      );

  const bool compute_beta_HS = (solverType_ == NCGU::NONLINEAR_CG_HS);

  //
  // A) Set up the storage for the algorithm
  //
  
  const RCP<DefaultPolyLineSearchPointEvaluator<Scalar> >
    pointEvaluator = defaultPolyLineSearchPointEvaluator<Scalar>();

  const RCP<UnconstrainedOptMeritFunc1D<Scalar> >
    meritFunc = unconstrainedOptMeritFunc1D<Scalar>(
      model_, paramIndex_, responseIndex_ );

  const RCP<const VectorSpaceBase<Scalar> >
    p_space = model_->get_p_space(paramIndex_),
    g_space = model_->get_g_space(responseIndex_);

  // Stoarge for current iteration
  RCP<VectorBase<Scalar> >
    p_k = rcpFromPtr(p_inout),        // Current solution for p
    p_kp1 = createMember(p_space),    // Trial point for p (in line search)
    g_vec = createMember(g_space),    // Vector (size 1) form of objective g(p) 
    g_grad_k = createMember(p_space), // Gradient of g DgDp^T
    d_k = createMember(p_space),      // Search direction
    g_grad_k_diff_km1 = null;         // g_grad_k - g_grad_km1 (if needed)

  // Storage for previous iteration
  RCP<VectorBase<Scalar> >
    g_grad_km1 = null, // Will allocate if we need it!
    d_km1 = null; // Will allocate if we need it!
  ScalarMag
    alpha_km1 = SMT::zero(),
    g_km1 = SMT::zero(),
    g_grad_km1_inner_g_grad_km1 = SMT::zero(),
    g_grad_km1_inner_d_km1 = SMT::zero();
  
  if (compute_beta_PR || compute_beta_HS) {
    g_grad_km1 = createMember(p_space);
    g_grad_k_diff_km1 = createMember(p_space);
  }
  
  if (compute_beta_HS) {
    d_km1 = createMember(p_space);
  }

  //
  // B) Do the nonlinear CG iterations
  //

  *out << "\nStarting nonlinear CG iterations ...\n";

  if (and_conv_tests_) {
    *out << "\nNOTE: Using AND of convergence tests!\n";
  }
  else {
    *out << "\nNOTE: Using OR of convergence tests!\n";
  }

  const Scalar alpha_init =
    ( !is_null(alpha_init_in) ? *alpha_init_in : alpha_init_ );
  const Scalar g_reduct_tol =
    ( !is_null(g_reduct_tol_in) ? *g_reduct_tol_in : g_reduct_tol_ );
  const Scalar g_grad_tol =
    ( !is_null(g_grad_tol_in) ? *g_grad_tol_in : g_grad_tol_ );

  const Ordinal globalDim = p_space->dim();

  bool foundSolution = false;
  bool fatalLinesearchFailure = false;
  bool restart = true;
  int numConsecutiveLineSearchFailures = 0;

  int numConsecutiveIters = 0;

  for (numIters_ = 0; numIters_ < maxIters_; ++numIters_, ++numConsecutiveIters) {

    Teuchos::OSTab tab(out);

    *out << "\nNonlinear CG Iteration k = " << numIters_ << "\n";

    Teuchos::OSTab tab2(out);

    //
    // B.1) Evaluate the point (on first iteration)
    //
    
    eval_g_DgDp(
      *model_, paramIndex_, *p_k, responseIndex_,
      numIters_ == 0 ? g_vec.ptr() : null, // Only on first iteration
      MEB::Derivative<Scalar>(g_grad_k, MEB::DERIV_MV_GRADIENT_FORM) );

    const ScalarMag g_k = get_ele(*g_vec, 0);
    // Above: If numIters_ > 0, then g_vec was updated in meritFunc->eval(...).

    //
    // B.2) Check for convergence
    //

    // B.2.a) ||g_k - g_km1|| |g_k + g_mag| <= g_reduct_tol

    bool g_reduct_converged = false;

    if (numIters_ > 0) {

      const ScalarMag g_reduct = g_k - g_km1;
      
      *out << "\ng_k - g_km1 = "<<g_reduct<<"\n";
      
      const ScalarMag g_reduct_err =
        SMT::magnitude(g_reduct / SMT::magnitude(g_k + g_mag_));
      
      g_reduct_converged = (g_reduct_err <= g_reduct_tol);
      
      *out << "\nCheck convergence: |g_k - g_km1| / |g_k + g_mag| = "<<g_reduct_err
           << (g_reduct_converged ? " <= " : " > ")
           << "g_reduct_tol = "<<g_reduct_tol<<"\n";
      
    }

    // B.2.b) ||g_grad_k|| g_mag <= g_grad_tol

    const Scalar g_grad_k_inner_g_grad_k = scalarProd<Scalar>(*g_grad_k, *g_grad_k);
    const ScalarMag norm_g_grad_k = ST::magnitude(ST::squareroot(g_grad_k_inner_g_grad_k));

    *out << "\n||g_grad_k|| = "<<norm_g_grad_k << "\n";

    const ScalarMag g_grad_err = norm_g_grad_k / g_mag_;

    const bool g_grad_converged = (g_grad_err <= g_grad_tol);

    *out << "\nCheck convergence: ||g_grad_k|| / g_mag = "<<g_grad_err
         << (g_grad_converged ? " <= " : " > ")
         << "g_grad_tol = "<<g_grad_tol<<"\n";

    // B.2.c) Convergence status
    
    bool isConverged = false;
    if (and_conv_tests_) {
      isConverged = g_reduct_converged && g_grad_converged;
    }
    else {
      isConverged = g_reduct_converged || g_grad_converged;
    }

    if (isConverged) {
      if (numIters_ < minIters_) {
        *out << "\nnumIters="<<numIters_<<" < minIters="<<minIters_
             << ", continuing on!\n";
      }
      else {
        *out << "\nFound solution, existing algorithm!\n";
        foundSolution = true;
      }
    }
    else {
      *out << "\nNot converged!\n";
    }
    
    if (foundSolution) {
      break;
    }

    //
    // B.3) Compute the search direction d_k
    //

    if (numConsecutiveIters == globalDim) {

      *out << "\nThe number of consecutive iterations exceeds the"
           << " global dimension so restarting!\n";

      restart = true;

    }

    if (restart) {

      *out << "\nResetting search direction back to steppest descent!\n";

      // d_k = -g_grad_k
      V_StV( d_k.ptr(), as<Scalar>(-1.0), *g_grad_k );

      restart = false;

    }
    else {
      
      // g_grad_k - g_grad_km1
      if (!is_null(g_grad_k_diff_km1)) {
        V_VmV( g_grad_k_diff_km1.ptr(), *g_grad_k, *g_grad_km1 );
      }

      // beta_FR = inner(g_grad_k, g_grad_k) / inner(g_grad_km1, g_grad_km1)
      const Scalar beta_FR =
        g_grad_k_inner_g_grad_k / g_grad_km1_inner_g_grad_km1;
      *out << "\nbeta_FR = " << beta_FR << "\n";
      // NOTE: Computing beta_FR is free so we might as well just do it!

      // beta_PR = inner(g_grad_k, g_grad_k - g_grad_km1) /
      //    inner(g_grad_km1, g_grad_km1)
      Scalar beta_PR = ST::zero();
      if (compute_beta_PR) {
        beta_PR =
          inner(*g_grad_k, *g_grad_k_diff_km1) / g_grad_km1_inner_g_grad_km1;
        *out << "\nbeta_PR = " << beta_PR << "\n";
      }

      // beta_HS = inner(g_grad_k, g_grad_k - g_grad_km1) /
      //    inner(g_grad_k - g_grad_km1, d_km1)
      Scalar beta_HS = ST::zero();
      if (compute_beta_HS) {
        beta_HS =
          inner(*g_grad_k, *g_grad_k_diff_km1) / inner(*g_grad_k_diff_km1, *d_km1);
        *out << "\nbeta_HS = " << beta_HS << "\n";
      }
      
      Scalar beta_k = ST::zero();
      switch(solverType_) {
        case NCGU::NONLINEAR_CG_FR: {
          beta_k = beta_FR;
          break;
        }
        case NCGU::NONLINEAR_CG_PR_PLUS: {
          beta_k = max(beta_PR, ST::zero());
          break;
        }
        case NCGU::NONLINEAR_CG_FR_PR: {
          // NOTE: This does not seem to be working :-(
          if (numConsecutiveIters < 2) {
            beta_k = beta_PR;
          }
          else if (beta_PR < -beta_FR)
            beta_k = -beta_FR;
          else if (ST::magnitude(beta_PR) <= beta_FR)
            beta_k = beta_PR;
          else // beta_PR > beta_FR
            beta_k = beta_FR;
        }
        case NCGU::NONLINEAR_CG_HS: {
          beta_k = beta_HS;
          break;
        }
        default:
          TEST_FOR_EXCEPT(true);
      }
      *out << "\nbeta_k = " << beta_k << "\n";

      // d_k = beta_k * d_last + -g_grad_k
      if (!is_null(d_km1))
        V_StV( d_k.ptr(), beta_k, *d_km1 );
      else
        Vt_S( d_k.ptr(), beta_k );
      Vp_StV( d_k.ptr(), as<Scalar>(-1.0), *g_grad_k );

    }
    
    //
    // B.4) Perform the line search
    //

    // B.4.a) Compute the initial step length

    Scalar alpha_k = as<Scalar>(-1.0);

    if (numIters_ == 0) {
      alpha_k = alpha_init;
    }
    else {
      if (alpha_reinit_) {
        alpha_k = alpha_init;
      }
      else {
        alpha_k = alpha_km1;
        // ToDo: Implement better logic from Nocedal and Wright for selecting
        // this step length after first iteration!
      }
    }

    // B.4.b) Perform the linesearch (computing updated quantities in process)

    pointEvaluator->initialize(tuple<RCP<const VectorBase<Scalar> > >(p_k, d_k)());

    ScalarMag g_grad_k_inner_d_k = ST::zero();

    // Set up the merit function to only compute the value
    meritFunc->setEvaluationQuantities(pointEvaluator, p_kp1, g_vec, null);

    PointEval1D<ScalarMag> point_k(ST::zero(), g_k);
    if (linesearch_->requiresBaseDeriv()) {
      g_grad_k_inner_d_k = scalarProd(*g_grad_k, *d_k);
      point_k.Dphi = g_grad_k_inner_d_k;
    }

    ScalarMag g_kp1 = computeValue(*meritFunc, alpha_k);
    // NOTE: The above call updates p_kp1 and g_vec as well!

    PointEval1D<ScalarMag> point_kp1(alpha_k, g_kp1);

    const bool linesearchResult = linesearch_->doLineSearch(
      *meritFunc, point_k, inOutArg(point_kp1), null );

    alpha_k = point_kp1.alpha;
    g_kp1 = point_kp1.phi;

    if (linesearchResult) {
      numConsecutiveLineSearchFailures = 0;
    }
    else {
      if (numConsecutiveLineSearchFailures==0) {
        *out << "\nLine search failure, resetting the search direction!\n";
        restart = true;
      }
      if (numConsecutiveLineSearchFailures==1) {
        *out << "\nLine search failure on last iteration also, terminating algorithm!\n";
        fatalLinesearchFailure = true;
      }
      ++numConsecutiveLineSearchFailures;
    }

    if (fatalLinesearchFailure) {
      break;
    }

    //
    // B.5) Transition to the next iteration
    //
    
    alpha_km1 = alpha_k;
    g_km1 = g_k;
    g_grad_km1_inner_g_grad_km1 = g_grad_k_inner_g_grad_k;
    g_grad_km1_inner_d_km1 = g_grad_k_inner_d_k;
    std::swap(p_k, p_kp1);
    if (!is_null(g_grad_km1))
      std::swap(g_grad_km1, g_grad_k);
    if (!is_null(d_km1))
      std::swap(d_k, d_km1);
    
#ifdef TEUCHOS_DEBUG
    // Make sure we compute these correctly before they are used!
    V_S(g_grad_k.ptr(), ST::nan());
    V_S(p_kp1.ptr(), ST::nan());
#endif

  }

  //
  // C) Final clean up
  //
  
  // Get the most current value of g(p)
  *g_opt_out = get_ele(*g_vec, 0);

  // Make sure that the final value for p has been copied in!
  V_V( p_inout, *p_k );

  if (!is_null(numIters_out)) {
    *numIters_out = numIters_;
  }

  if (numIters_ == maxIters_) {
    *out << "\nMax nonlinear CG iterations exceeded!\n";
  }
  
  if (foundSolution) {
    return NonlinearCGUtils::SOLVE_SOLUTION_FOUND;
  }
  else if(fatalLinesearchFailure) {
    return NonlinearCGUtils::SOLVE_LINSEARCH_FAILURE;
  }

  // Else, the max number of iterations was exceeded
  return NonlinearCGUtils::SOLVE_MAX_ITERS_EXCEEDED;

}


} // namespace OptiPack


#endif // OPTIPACK_NONLINEAR_CG_DEF_HPP
