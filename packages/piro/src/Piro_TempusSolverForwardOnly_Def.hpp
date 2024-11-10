// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_TempusSolverForwardOnly.hpp"

#include "Piro_ObserverToTempusIntegrationObserverAdapter.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "Piro_MatrixFreeDecorator.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Piro_InvertMassMatrixDecorator.hpp"

#ifdef HAVE_PIRO_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#ifdef HAVE_PIRO_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

#include <string>
#include <stdexcept>
#include <iostream>

template <typename Scalar>
Piro::TempusSolverForwardOnly<Scalar>::TempusSolverForwardOnly() :
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized(false)
{
}

template <typename Scalar>
Piro::TempusSolverForwardOnly<Scalar>::TempusSolverForwardOnly(
  const Teuchos::RCP<Teuchos::ParameterList> &appParams,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
  const Teuchos::RCP<Tempus::IntegratorObserver<Scalar> > &observer)
  : out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    isInitialized(false)
{
  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");
  if (jacobianSource == "Matrix-Free") {
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > mf_model;
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      mf_model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model,
                           appParams->get<double>("Matrix-Free Perturbation")));
    }
    else mf_model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model));
    initialize(appParams, mf_model, observer);
  }
  else
    initialize(appParams, in_model, observer);
}

template <typename Scalar>
void Piro::TempusSolverForwardOnly<Scalar>::initialize(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Tempus::IntegratorObserver<Scalar> > &observer)
{

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  *out << "\nA) Get the base parameter list ...\n";
  //

  TEUCHOS_TEST_FOR_EXCEPTION(!appParams->isSublist("Tempus"),
      Teuchos::Exceptions::InvalidParameter, std::endl <<
      "Error! Piro::TempusSolverForwardOnly: must have Tempus sublist ");

  RCP<Teuchos::ParameterList> tempusPL = sublist(appParams, "Tempus", true);

  //
  *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
  //

  // Validation of Parameters in tempusPL is done through each objects
  // "create" function, e.g., createIntegratorBasic, createStepper,
  // and createTimeStepControl.
  fwdStateIntegrator = Tempus::createIntegratorBasic(tempusPL, in_model);


  // Tempus does not have a Verbosity parameter, so set to default.
  solnVerbLevel = Teuchos::VERB_DEFAULT;

  if (Teuchos::nonnull(observer))
    fwdStateIntegrator->setObserver(observer);

  fwdStateIntegrator->initialize();

  isInitialized = true;
}

template <typename Scalar>
Piro::TempusSolverForwardOnly<Scalar>::TempusSolverForwardOnly(
  const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
  const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
  Scalar finalTime,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel,
  Teuchos::EVerbosityLevel verbosityLevel)
  : fwdStateIntegrator(stateIntegrator),
    initialConditionModel(icModel),
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    solnVerbLevel(verbosityLevel),
    isInitialized(true)
{
  if (fwdStateIntegrator->getStepper() != stateStepper) {
    *out << "\nTempusSolverForwardOnly: Resetting the Stepper!\n";
    fwdStateIntegrator->setStepper(stateStepper);
  }

  if (getTimeStepSolver() != timeStepSolver) {
    *out << "\nTempusSolverForwardOnly: Resetting the solver!\n";
    fwdStateIntegrator->getStepper()->setSolver(timeStepSolver);
  }

  if (getUnderlyingModel() != underlyingModel) {
    *out << "\nTempusSolverForwardOnly: Resetting the model!\n";
    fwdStateIntegrator->getStepper()->setModel(underlyingModel);
  }

  fwdStateIntegrator->getNonConstTimeStepControl()->setInitTime(0.0);
  fwdStateIntegrator->getNonConstTimeStepControl()->setFinalTime(finalTime);
}

template <typename Scalar>
Piro::TempusSolverForwardOnly<Scalar>::TempusSolverForwardOnly(
  const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
  const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
  Scalar initialTime,
  Scalar finalTime,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel,
  Teuchos::EVerbosityLevel verbosityLevel)
  : fwdStateIntegrator(stateIntegrator),
    initialConditionModel(icModel),
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    solnVerbLevel(verbosityLevel),
    isInitialized(true)
{
  if (fwdStateIntegrator->getStepper() != stateStepper) {
    *out << "\nTempusSolverForwardOnly: Resetting the Stepper!\n";
    fwdStateIntegrator->setStepper(stateStepper);
  }

  if (fwdStateIntegrator->getStepper()->getSolver() != timeStepSolver) {
    *out << "\nTempusSolverForwardOnly: Resetting the solver!\n";
    fwdStateIntegrator->getStepper()->setSolver(timeStepSolver);
  }

  if (getUnderlyingModel() != underlyingModel) {
    *out << "\nTempusSolverForwardOnly: Resetting the model!\n";
    fwdStateIntegrator->getStepper()->setModel(underlyingModel);
  }

  fwdStateIntegrator->getNonConstTimeStepControl()->setInitTime(initialTime);
  fwdStateIntegrator->getNonConstTimeStepControl()->setFinalTime(finalTime);
}

template <typename Scalar>
Teuchos::RCP<const Tempus::IntegratorBasic<Scalar> >
Piro::TempusSolverForwardOnly<Scalar>::getTempusIntegrator() const
{
  return fwdStateIntegrator;
}


template <typename Scalar>
Teuchos::RCP<const Thyra::NonlinearSolverBase<Scalar> >
Piro::TempusSolverForwardOnly<Scalar>::getTimeStepSolver() const
{
  return fwdStateIntegrator->getStepper()->getSolver();
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TempusSolverForwardOnly<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= get_num_p() || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::TempusSolverForwardOnly::get_p_map():  " <<
      "Invalid parameter index l = " <<
      l << std::endl);

  return getUnderlyingModel()->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TempusSolverForwardOnly<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > get_num_g() || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::TempusSolverForwardOnly::get_g_map():  " <<
      "Invalid response index j = " <<
      j << std::endl);

  if (j < get_num_g()) {
    return getUnderlyingModel()->get_g_space(j);
  } else {
    // j == get_num_g()
    return getUnderlyingModel()->get_x_space();
  }
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TempusSolverForwardOnly<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = getUnderlyingModel()->getNominalValues();
  for (int l = 0; l < get_num_p(); ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TempusSolverForwardOnly<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(get_num_p());
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::TempusSolverForwardOnly<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(get_num_p(), get_num_g() + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = getUnderlyingModel()->createOutArgs();

  if (get_num_p() > 0) {
    // Only one parameter supported
    const int l = 0;

    if (Teuchos::nonnull(initialConditionModel)) {
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs =
        initialConditionModel->createOutArgs();
      const Thyra::ModelEvaluatorBase::DerivativeSupport init_dxdp_support =
        initCondOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, initCondOutArgs.Ng() - 1, l);
      if (!init_dxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
        // Ok to return early since only one parameter supported
        return outArgs;
      }
    }

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
        get_num_g(),
        l,
        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

    if (get_num_g() > 0) {
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
void Piro::TempusSolverForwardOnly<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // TODO: Support more than 1 parameter and 1 response
  const int j = 0;
  const int l = 0;

  // Parse InArgs
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (get_num_p() > 0) {
    p_in = inArgs.get_p(l);
  }
  RCP<const Thyra::VectorBase<Scalar> > p_in2;  //JF add for multipoint
  if (get_num_p() > 1) {
    p_in2 = inArgs.get_p(l+1);
  }

  // Parse OutArgs
  RCP<Thyra::VectorBase<Scalar> > g_out;
  if (get_num_g() > 0) {
    g_out = outArgs.get_g(j);
  }
  const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(get_num_g());

  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = getUnderlyingModel()->getNominalValues();

  // Set initial time in ME if needed

  if(get_t_initial() > 0.0 && state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_t))

    state_ic.set_t(get_t_initial());

  if (Teuchos::nonnull(initialConditionModel)) {
    // The initial condition depends on the parameter
    // It is found by querying the auxiliary model evaluator as the last response
    const RCP<Thyra::VectorBase<Scalar> > initialState =
      Thyra::createMember(getUnderlyingModel()->get_x_space());

    {
      Thyra::ModelEvaluatorBase::InArgs<Scalar> initCondInArgs = initialConditionModel->createInArgs();
      if (get_num_p() > 0) {
        initCondInArgs.set_p(l, inArgs.get_p(l));
      }

      Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs = initialConditionModel->createOutArgs();
      initCondOutArgs.set_g(initCondOutArgs.Ng() - 1, initialState);

      initialConditionModel->evalModel(initCondInArgs, initCondOutArgs);
    }

    state_ic.set_x(initialState);
  }

  // Set paramters p_in as part of initial conditions
  if (get_num_p() > 0) {
    if (Teuchos::nonnull(p_in)) {
      state_ic.set_p(l, p_in);
    }
  }
  if (get_num_p() > 1) { //JF added for multipoint
    if (Teuchos::nonnull(p_in2)) {
      state_ic.set_p(l+1, p_in2);
    }
  }

  *out << "\nstate_ic:\n" << Teuchos::describe(state_ic, solnVerbLevel);

  //JF  may need a version of the following for multipoint, i.e. get_num_p()>1, l+1, if we want sensitivities
  RCP<Thyra::MultiVectorBase<Scalar> > dgxdp_out;
  Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv_out;
  if (get_num_p() > 0) {
    const Thyra::ModelEvaluatorBase::DerivativeSupport dgxdp_support =
      outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, get_num_g(), l);
    if (dgxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgxdp_deriv =
        outArgs.get_DgDp(get_num_g(), l);
      dgxdp_out = dgxdp_deriv.getMultiVector();
    }

    if (get_num_g() > 0) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
      if (!dgdp_support.none()) {
        dgdp_deriv_out = outArgs.get_DgDp(j, l);
      }
    }
  }

  const bool requestedSensitivities =
    Teuchos::nonnull(dgxdp_out) || !dgdp_deriv_out.isEmpty();

  RCP<const Thyra::VectorBase<Scalar> > finalSolution;
  if (!requestedSensitivities) {
    //
    *out << "\nE) Solve the forward problem ...\n";
    //

    Teuchos::RCP<const Thyra::VectorBase<Scalar>> xinit, xdotinit, xdotdotinit;
    typedef Thyra::ModelEvaluatorBase MEB;
    if (state_ic.supports(MEB::IN_ARG_t)) state_ic.set_t(get_t_initial());
    if (state_ic.supports(MEB::IN_ARG_x))         xinit       = state_ic.get_x();
    if (state_ic.supports(MEB::IN_ARG_x_dot))     xdotinit    = state_ic.get_x_dot();
    if (state_ic.supports(MEB::IN_ARG_x_dot_dot)) xdotdotinit = state_ic.get_x_dot_dot();
    fwdStateIntegrator->initializeSolutionHistory(get_t_initial(), xinit, xdotinit, xdotdotinit);
    fwdStateIntegrator->initialize();

    Scalar t_final = get_t_final();
    fwdStateIntegrator->advanceTime(t_final);
    double time = fwdStateIntegrator->getTime();

    Scalar diff = 0.0;
    if (abs(t_final) == 0) diff = abs(time-t_final);
    else diff = abs(time-t_final)/abs(t_final);
    if (diff > 1.0e-10)
      *out << "\n WARNING: Piro::TempusSolverForwardOnly did not make it to final time.\n";

    finalSolution = fwdStateIntegrator->getX();

    if (Teuchos::VERB_MEDIUM <= solnVerbLevel)
      std::cout << "Final Solution\n" << *finalSolution << std::endl;

  }

  *out << "\nF) Check the solution to the forward problem ...\n";

  // As post-processing step, calculate responses at final solution
  {
    Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = getUnderlyingModel()->createInArgs();
    {
      modelInArgs.set_x(finalSolution);
      if (get_num_p() > 0) {
        modelInArgs.set_p(l, p_in);
      }
      if (get_num_p() > 1) {  //JF added for multipoint
        modelInArgs.set_p(l+1, p_in2);
      }
      //Set time to be final time at which the solve occurs (< t_final in the case we don't make it to t_final).
      modelInArgs.set_t(fwdStateIntegrator->getSolutionHistory()->getCurrentState()->getTime());
    }

    Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = getUnderlyingModel()->createOutArgs();
    if (Teuchos::nonnull(g_out)) {
      Thyra::put_scalar(Teuchos::ScalarTraits<Scalar>::zero(), g_out.ptr());
      modelOutArgs.set_g(j, g_out);
    }

    getUnderlyingModel()->evalModel(modelInArgs, modelOutArgs);
  }

  // Return the final solution as an additional g-vector, if requested
  if (Teuchos::nonnull(gx_out))
    Thyra::copy(*finalSolution, gx_out.ptr());
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::TempusSolverForwardOnly<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  TEUCHOS_ASSERT(j != get_num_g());
  const Teuchos::Array<Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > > dummy =
    Teuchos::tuple(Thyra::zero<Scalar>(this->get_g_space(j), this->get_p_space(l)));
  return Teuchos::rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::TempusSolverForwardOnly<Scalar>::getValidTempusParameters() const
{
  if (fwdStateIntegrator == Teuchos::null) {
    auto integrator = Teuchos::rcp(new Tempus::IntegratorBasic<Scalar>());
    return integrator->getValidParameters();
  }

  return fwdStateIntegrator->getValidParameters();
}

template <typename Scalar>
Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >
Piro::TempusSolverForwardOnly<Scalar>::getUnderlyingModel() const
{
  return Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar> > (
    fwdStateIntegrator->getStepper()->getModel());
}

template <typename Scalar>
Scalar
Piro::TempusSolverForwardOnly<Scalar>::get_t_initial() const
{
  return fwdStateIntegrator->getTimeStepControl()->getInitTime();
}

template <typename Scalar>
Scalar
Piro::TempusSolverForwardOnly<Scalar>::get_t_final() const
{
  return fwdStateIntegrator->getTimeStepControl()->getFinalTime();
}

template <typename Scalar>
int
Piro::TempusSolverForwardOnly<Scalar>::get_num_p() const
{
  return getUnderlyingModel()->Np();
}

template <typename Scalar>
int
Piro::TempusSolverForwardOnly<Scalar>::get_num_g() const
{
  return getUnderlyingModel()->Ng();
}

