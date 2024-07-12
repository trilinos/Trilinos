// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>
#include "Piro_TransientDecorator.hpp"

template <typename Scalar>
Piro::TrapezoidRuleSolver<Scalar>::
TrapezoidRuleSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams_,
                          const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model_,
                          const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr_,
                          const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer_ ) :
  appParams(appParams_),
  model(Teuchos::rcp(new Piro::TrapezoidDecorator<Scalar>(model_))),
  observer(observer_),
  solMgr(solMgr_),
  out(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Only non-null when adaptation is used
  if(Teuchos::nonnull(solMgr)) {
    solMgr->initialize(rcp(new Thyra::TransAdaptiveState(model_)));
  }

  num_p = model->Np();
  num_g = model->Ng();

  TEUCHOS_TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::TrapezoidRuleSolver " <<
                     "Not Implemented for Np>1 : " << num_p << std::endl);
  TEUCHOS_TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::TrapezoidRuleSolver " <<
                     "Not Implemented for Ng>1 : " << num_g << std::endl);

  *out << "\nA) Get the base parameter list ...\n";
 
  RCP<Teuchos::ParameterList> trPL = sublist(appParams, "Trapezoid Rule", true);
  trPL->validateParameters(*getValidTrapezoidRuleParameters(),0);

  {
    const std::string verbosity = trPL->get("Verbosity Level", "VERB_DEFAULT");
    solnVerbLevel = Teuchos::VERB_DEFAULT;
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
  }

  numTimeSteps = trPL->get("Num Time Steps", 10);
  t_final = trPL->get("Final Time", 0.1);
  t_init  = trPL->get("Initial Time", 0.0);
  delta_t = (t_final - t_init) / numTimeSteps;

  *out << "\nB) Using Trapezoid Decorator and NOX Solver\n";

  // Construct NOX solver -- will look for NOX sublist -- this must be set!!
  trPL->sublist("NOX").set("Reset Initial Guess",true);
  noxSolver = Teuchos::rcp(new Piro::NOXSolver<Scalar>(appParams, model));

}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TrapezoidRuleSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::TrapezoidRuleSolver::get_p_map():  " <<
      "Invalid parameter index l = " <<
      l << std::endl);

  return model->get_p_space(l);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TrapezoidRuleSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::TrapezoidRuleSolver::get_g_map():  " <<
      "Invalid response index j = " <<
      j << std::endl);

  if (j < num_g) {
    return model->get_g_space(j);
  } else {
    // j == num_g
    return model->get_x_space();
  }
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TrapezoidRuleSolver<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TrapezoidRuleSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::TrapezoidRuleSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p, num_g + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();

  if (num_p > 0) {
    // Only one parameter supported
    const int l = 0;

    // Computing the DxDp sensitivity for a transient problem currently requires the evaluation of
    // the mutilivector-based, Jacobian-oriented DfDp derivatives of the underlying transient model.
    const Thyra::ModelEvaluatorBase::DerivativeSupport model_dfdp_support =
      modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
    if (!model_dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      // Ok to return early since only one parameter supported
      return outArgs;
    }

    // Solution sensitivity
    outArgs.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
        num_g,
        l,
        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

    if (num_g > 0) {
      // Only one response supported
      const int j = 0;

      const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdx_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
      if (!model_dgdx_support.none()) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdp_support =
          modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        // Response sensitivity
        Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support;
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        }
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        }
        outArgs.setSupports(
            Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
            j,
            l,
            dgdp_support);
      }
    }
  }

  return outArgs;

}

template <typename Scalar>
void
Piro::TrapezoidRuleSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(model->get_x_dotdot()),
                     std::logic_error,
                     std::endl << "Error in Piro::TrapezoidRuleSolver " <<
                     "Model must specify soln, soln_dot, and soln_dotdot." << std::endl);

  // noxSolver is the Piro::NOXSolver object
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nox_inargs = noxSolver->createInArgs();
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> nox_outargs = noxSolver->createOutArgs();

  // Parse InArgs
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(0);
    nox_inargs.set_p(0, p_in);
  }

  // Parse OutArgs: always 1 extra
  RCP<Thyra::VectorBase<Scalar> > g_out;
  if (num_g > 0) {
    g_out = outArgs.get_g(0);
    nox_outargs.set_g(0, g_out);
  }
  RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);
  if (Teuchos::is_null(gx_out)) {
    // Solution not requested by caller as a response, create local temporary instead
    // model is the Trapezoid decorator
    gx_out = Thyra::createMember<Scalar>(model->get_x_space());
  }
  nox_outargs.set_g(num_g, gx_out);

  // Build a multivector holding x (0th vector), v (1st vector), and a (2nd vector)
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > soln = createMembers(model->get_x_space(), 3);

// create a new vector and fill it with the contents of model->get_x()
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x = soln->col(0);
  assign(x.ptr(), *model->getNominalValues().get_x());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > v = soln->col(1);
  assign(v.ptr(), *model->getNominalValues().get_x_dot());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > a = soln->col(2);
  assign(a.ptr(), *model->get_x_dotdot());

// Note that Thyra doesn't have x_dotdot - go get it from the transient decorator around the Albany model
//  Teuchos::RCP<Thyra::VectorBase<Scalar> > a = model->get_x_dotdot()->clone_v();

  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_pred_a = Thyra::createMember<Scalar>(model->get_f_space());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_pred_v = Thyra::createMember<Scalar>(model->get_f_space());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > a_old = Thyra::createMember<Scalar>(model->get_f_space());

  TEUCHOS_TEST_FOR_EXCEPTION(v == Teuchos::null || x == Teuchos::null,
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::TrapezoidRuleSolver " <<
                     "Requires initial x and x_dot: " << std::endl);

  Scalar t = t_init;

  // Observe initial condition
  if (observer != Teuchos::null) observer->observeSolution(*soln, t);

  Scalar nrm = norm_2(*v);
  *out << "Initial Velocity = " << nrm << std::endl;

   //calculate intial acceleration using small time step (1.0e-3*delta_t)
   // AGS: Check this for inital velocity
  if (calc_init_accel_ == true) {
    Scalar pert= 1.0e6 * 4.0 / (delta_t * delta_t);
    assign(x_pred_a.ptr(), *x);
    assign(x_pred_v.ptr(), *x);

    Vp_StV(x_pred_v.ptr(), sqrt(pert), *v);
    model->injectData(x_pred_a, x_pred_a, pert, x_pred_v, sqrt(pert), t);

    noxSolver->evalModel(nox_inargs, nox_outargs);

    V_StVpStV(a.ptr(), pert, *gx_out,  -pert, *x_pred_a);
    nrm = norm_2(*a);
    *out << "Calculated a_init = " << nrm << std::endl;
   }

   // Start integration loop
   Scalar fdt2 = 4.0 / (delta_t * delta_t);
   Scalar dt2f =  delta_t * delta_t / 4.0;
   Scalar hdt  =  delta_t/ 2.0;
   Scalar tdt  =  2.0 / delta_t;

// GAH time step loop

  for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {

     t += delta_t;

     if(Teuchos::nonnull(solMgr)){
       solMgr->setIteration(timeStep);
       solMgr->setTime(t);
     }

     // Adapt the mesh if the user has turned on adaptation
     // and we have passed the criteria to adapt
     if(Teuchos::nonnull(solMgr) && 
         solMgr->isAdaptive() && solMgr->queryAdaptationCriteria()){

       // Adapt the mesh, send exception if adaptation fails
       TEUCHOS_TEST_FOR_EXCEPTION(
          !solMgr->adaptProblem(), 
          std::logic_error,
          "Error: Piro_TrapezoidRuleSolver, cannot adapt the mesh!" << std::endl);

       // Gets the Thyra MV directly from the updated discretization
       soln = solMgr->getCurrentSolution();

       // create a new vector and fill it with the contents of model->get_x()
       x = soln->col(0);
       v = soln->col(1);
       a = soln->col(2);
       nrm = norm_2(*a);

       x_pred_a = Thyra::createMember<Scalar>(model->get_f_space());
       x_pred_v = Thyra::createMember<Scalar>(model->get_f_space());
       a_old = Thyra::createMember<Scalar>(model->get_f_space());

       model->resize(x);
       gx_out = Thyra::createMember<Scalar>(model->get_x_space());

       nox_outargs.set_g(num_g, gx_out);

       noxSolver->reset();

     }

     assign(a_old.ptr(), *a);
     assign(x_pred_a.ptr(), *x);

     Vp_StV(x_pred_a.ptr(), dt2f, *a);
     Vp_StV(x_pred_a.ptr(), delta_t, *v);

     assign(x_pred_v.ptr(), *x);
     Vp_StV(x_pred_v.ptr(), hdt, *v);

     model->injectData(x, x_pred_a, fdt2, x_pred_v, tdt, t);

     noxSolver->evalModel(nox_inargs, nox_outargs);
     // Copy out final solution from nonlinear solver
//     *x =  *gx_out;
     assign(x.ptr(), *gx_out);

     // Compute a and v and new conditions
     V_StVpStV(a.ptr(), fdt2, *x, -fdt2, *x_pred_a);

     Vp_StV(v.ptr(), hdt, *a);
     Vp_StV(v.ptr(), hdt, *a_old);
     // Should be equivalent to: v->Update(tdt, *x, -tdt, *x_pred_v, 0.0);
     
     nrm = norm_2(*a);

     // Observe completed time step
     if (observer != Teuchos::null) observer->observeSolution(*soln, t);

     if (g_out != Teuchos::null)
       *out << "Responses at time step(time) = " << timeStep << "("<<t<<")\n" << g_out << std::endl;

   }

}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::TrapezoidRuleSolver<Scalar>::getValidTrapezoidRuleParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidTrapezoidRuleParams"));;
  validPL->set<int>("Num Time Steps", 0, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->set<double>("Initial Time", 0.0, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<bool>("Lump Mass Matrix", false, "Boolean to tell code to lump mass matrix");
  validPL->set<bool>("Constant Mass Matrix", false, "Boolean to tell code to if mass matrix is constant in time");
  validPL->sublist("Stratimikos", false, "");
  validPL->sublist("NOX", false, "");
  return validPL;
}

template <typename Scalar>
Teuchos::RCP<Piro::NOXSolver<Scalar> >
Piro::TrapezoidRuleSolver<Scalar>::getNOXSolver() const
{
  return noxSolver;
}

template <typename Scalar>
Teuchos::RCP<Piro::TrapezoidDecorator<Scalar> >
Piro::TrapezoidRuleSolver<Scalar>::getDecorator() const
{
  return model;
}

template <typename Scalar>
Teuchos::RCP<Thyra::AdaptiveSolutionManager>
Piro::TrapezoidRuleSolver<Scalar>::getSolutionManager() const
{
  return solMgr;
}

/****************************************************************************/

template <typename Scalar>
Piro::TrapezoidDecorator<Scalar>::TrapezoidDecorator(
                          const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& modelEvaluator) :
  Thyra::ModelEvaluatorDelegatorBase<Scalar>(modelEvaluator),
  DMEWSF(Teuchos::rcp_dynamic_cast<Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar> >(modelEvaluator)),
  fdt2(0.0),
  tdt(0.0),
  time(0.0),
  out(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(DMEWSF),
	std::logic_error,
    "Model passed into Trapezoid decorator does not cast to a Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>"
    << std::endl);

  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = this->getUnderlyingModel()->getNominalValues();

  xDotDot = Thyra::createMember<Scalar>(this->getUnderlyingModel()->get_x_space());
  xDot = Thyra::createMember<Scalar>(this->getUnderlyingModel()->get_x_space());

  x_pred_a = state_ic.get_x()->clone_v();
  x_pred_v = state_ic.get_x()->clone_v();
  x_save = state_ic.get_x()->clone_v();

}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Piro::TrapezoidDecorator<Scalar>::get_x() const
{

  return x_save;

}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Piro::TrapezoidDecorator<Scalar>::get_x_dot() const
{

  return xDot;

}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Piro::TrapezoidDecorator<Scalar>::get_x_dotdot() const
{
  Teuchos::RCP<const Piro::TransientDecorator<Scalar> > dec =
     Teuchos::rcp_dynamic_cast<const Piro::TransientDecorator<Scalar> >
         (DMEWSF->getUnderlyingModel());

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(dec),
    std::logic_error,
    "Underlying model in trapezoid decorator does not cast to a Piro::TransientDecorator<Scalar>\n"); 

  return dec->get_x_dotdot();

}

template <typename Scalar>
void
Piro::TrapezoidDecorator<Scalar>::reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint,
    const bool wasSolved)
{
  this->getNonconstUnderlyingModel()->reportFinalPoint(finalPoint,wasSolved);
}

template <typename Scalar>
void Piro::TrapezoidDecorator<Scalar>::injectData(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_pred_a_, Scalar fdt2_,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_pred_v_, Scalar tdt_,
    Scalar time_)
{

  assign(x_save.ptr(), *x_);
  assign(x_pred_a.ptr(), *x_pred_a_);
  assign(x_pred_v.ptr(), *x_pred_v_);

  fdt2 = fdt2_;
  tdt = tdt_;
  time = time_;

}

template <typename Scalar>
void Piro::TrapezoidDecorator<Scalar>::resize(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_)
{

  xDotDot = Thyra::createMember<Scalar>(x_->space());
  xDot = Thyra::createMember<Scalar>(x_->space());

  x_pred_a = x_->clone_v();
  x_pred_v = x_->clone_v();
  x_save   = x_->clone_v();

}


template <typename Scalar>
void Piro::TrapezoidDecorator<Scalar>::evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const {


  using Teuchos::RCP;
  using Teuchos::rcp;

  // Copy outArgs; add time term
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs(inArgs);
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs(outArgs);

  V_StVpStV(xDotDot.ptr(), fdt2, *inArgs.get_x(), -fdt2, *x_pred_a);

  V_StVpStV(xDot.ptr(), tdt, *inArgs.get_x(), -tdt, *x_pred_v);
  modelInArgs.set_x_dot(xDot);

  modelInArgs.set_alpha(tdt);    // tdt  = 2/dt
  modelInArgs.set_beta(1.0);
  modelInArgs.set_t(time);

  // No xdotdot support in Thyra, so set directly in the model
  // Need to set xdotdot and omega in the underlying model if the model is a DMEWSF

  Teuchos::RCP<const Piro::TransientDecorator<Scalar> > dec =
       Teuchos::rcp_dynamic_cast<const Piro::TransientDecorator<Scalar> >
           (DMEWSF->getUnderlyingModel());

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(dec),
      std::logic_error,
      "Underlying model in trapezoid decorator does not cast to a Piro::TransientDecorator<Scalar>\n");

  dec->set_omega(fdt2);  // fdt2 = 4/(dt)^2
  dec->set_x_dotdot_data(xDotDot); // no xdotdot in Thyra inArgs

  //Evaluate the underlying model
  this->getUnderlyingModel()->evalModel(modelInArgs, modelOutArgs);


}

