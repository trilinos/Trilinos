// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_TempusSolver.hpp"

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

#ifdef HAVE_PIRO_IFPACK2
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

#ifdef HAVE_PIRO_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#ifdef HAVE_PIRO_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

//#define DEBUG_OUTPUT

#include <string>
#include <stdexcept>
#include <iostream>

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TempusSolver() :
#else
template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver() :
#endif
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized(false)
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TempusSolver(
#else
template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver(
#endif
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    bool computeSensitivities,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver):
  computeSensitivities_(computeSensitivities),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized(false),
  piroObserver_(piroObserver),
  supports_x_dotdot_(false)
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");
  if (jacobianSource == "Matrix-Free") {
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model,
                           appParams->get<double>("Matrix-Free Perturbation")));
    }
    else model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model));
    initialize(appParams, model);
  }
  else
    initialize(appParams, in_model);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initialize(
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::initialize(
#endif
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP< Thyra::ModelEvaluator<Scalar> > &in_model)
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // set some internals
  model = in_model;
  num_p = in_model->Np();
  num_g = in_model->Ng();

  //
  *out << "\nA) Get the base parameter list ...\n";
  //

  if (appParams->isSublist("Tempus")) {
    RCP<Teuchos::ParameterList> tempusPL = sublist(appParams, "Tempus", true);
    abort_on_failure_ = tempusPL->get<bool>("Abort on Failure", true);
    //*out << "tempusPL = " << *tempusPL << "\n";
    RCP<Teuchos::ParameterList> integratorPL = sublist(tempusPL, "Tempus Integrator", true);
    //*out << "integratorPL = " << *integratorPL << "\n";
    //IKT, 10/31/16, FIXME: currently there is no Verbosity Sublist in Tempus, but
    //Curt will add this at some point.  When this option is added, set Verbosity
    //based on that sublist, rather than hard-coding it here.
    solnVerbLevel = Teuchos::VERB_DEFAULT;

    RCP<Teuchos::ParameterList> timeStepControlPL = sublist(integratorPL, "Time Step Control", true);
    t_initial = timeStepControlPL->get<Scalar>("Initial Time", 0.0);
    t_final = timeStepControlPL->get<Scalar>("Final Time", 0.1);
    RCP<Teuchos::ParameterList> stepperPL = sublist(tempusPL, "Tempus Stepper", true);
    //*out << "stepperPL = " << *stepperPL << "\n";
    const std::string stepperType = stepperPL->get<std::string>("Stepper Type", "Backward Euler");
    //*out << "Stepper Type = " << stepperType << "\n";

    //
    // *out << "\nB) Create the Stratimikos linear solver factory ...\n";
    //
    // This is the linear solve strategy that will be used to solve for the
    // linear system with the W.
    //
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef HAVE_PIRO_IFPACK2
    typedef Thyra::PreconditionerFactoryBase<double> Base;
#ifdef ALBANY_BUILD
    typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double, LocalOrdinal, GlobalOrdinal, Node> > Impl;
#else
    typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double> > Impl;
#endif
    linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef HAVE_PIRO_MUELU
#ifdef ALBANY_BUILD
    Stratimikos::enableMueLu<LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);
#else
    Stratimikos::enableMueLu(linearSolverBuilder);
#endif
#endif

    linearSolverBuilder.setParameterList(sublist(tempusPL, "Stratimikos", true));
    tempusPL->validateParameters(*getValidTempusParameters(),0);
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    //
    *out << "\nC) Create and initalize the forward model ...\n";

    //
    // C.1) Create the underlying Thyra::ModelEvaluator
    // already constructed as "model". Decorate if needed.
    // IKT, 12/8/16: it may be necessary to expand the list of conditions
    // below, as more explicit schemes get added to Tempus
    // Explicit time-integrators for 1st order ODEs 
    if (
      stepperType == "RK Forward Euler" ||
      stepperType == "RK Explicit 4 Stage" ||
      stepperType == "RK Explicit 3/8 Rule" ||
      stepperType == "RK Explicit 4 Stage 3rd order by Runge" ||
      stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray"||
      stepperType == "RK Explicit 3 Stage 3rd order" ||
      stepperType == "RK Explicit 3 Stage 3rd order TVD" ||
      stepperType == "RK Explicit 3 Stage 3rd order by Heun" ||
      stepperType == "RK Explicit 2 Stage 2nd order by Runge" ||
      stepperType == "RK Explicit Trapezoidal" ||
      stepperType == "General ERK" ) {

      bool invertMassMatrix = tempusPL->get("Invert Mass Matrix", false); 
      if (!invertMassMatrix) {
        *out << "\n WARNING in Piro::TempusSolver!  You are attempting to run \n" 
             << "Explicit Stepper (" << stepperType << ") with 'Invert Mass Matrix' set to 'false'. \n" 
             << "This option should be set to 'true' unless your mass matrix is the identiy.\n"; 
      }
      else {
        Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > origModel = model;
#ifdef ALBANY_BUILD
        model = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
#else
        model = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar>(
#endif
        sublist(tempusPL,"Stratimikos", true), origModel, true, tempusPL->get("Lump Mass Matrix", false),false));
      }
    }

    //Explicit time-integrators for 2nd order ODEs
    //IKT, FIXME: fill this in as more explicit integrators for 2nd order ODEs are added to Tempus.
    else if (stepperType == "Newmark Explicit a-Form") {
      bool invertMassMatrix = tempusPL->get("Invert Mass Matrix", false); 
      if (!invertMassMatrix) {
        *out << "\n WARNING in Piro::TempusSolver!  You are attempting to run \n" 
             << "'Newmark Explicit a-Form' Stepper with 'Invert Mass Matrix' set to 'false'. \n" 
             << "This option should be set to 'true' unless your mass matrix is the identiy.\n"; 
      }
      else {
        Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > origModel = model;
#ifdef ALBANY_BUILD
        model = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
#else
        model = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar>(
#endif
        sublist(tempusPL,"Stratimikos", true), origModel, true, tempusPL->get("Lump Mass Matrix", false),true));
      }
    }
    // C.2) Create the Thyra-wrapped ModelEvaluator

    thyraModel = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>(model, lowsFactory));

    const RCP<const Thyra::VectorSpaceBase<double> > x_space = thyraModel->get_x_space();

    //
    *out << "\nD) Create the stepper and integrator for the forward problem ...\n";

    //Create Tempus integrator with observer using tempusPL and model.
    fwdStateIntegrator = Tempus::integratorBasic<Scalar>(tempusPL, model);

    //Get stepper from integrator
    fwdStateStepper = fwdStateIntegrator->getStepper();

    //Set observer
    supports_x_dotdot_ = model->createInArgs().supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot);
    setObserver();  

  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        appParams->isSublist("Tempus"),
        Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::TempusSolver: must have Tempus sublist ");

  }

  isInitialized = true;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TempusSolver(
#else
template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver(
#endif
    const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
    Scalar finalTime,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel,
    Teuchos::EVerbosityLevel verbosityLevel) :
  fwdStateIntegrator(stateIntegrator),
  fwdStateStepper(stateStepper),
  fwdTimeStepSolver(timeStepSolver),
  model(underlyingModel),
  initialConditionModel(icModel),
  t_initial(0.0),
  t_final(finalTime),
  num_p(model->Np()),
  num_g(model->Ng()),
  computeSensitivities_(false),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel(verbosityLevel),
  isInitialized(true)
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  if (fwdStateStepper->getModel() != underlyingModel) {
    fwdStateStepper->setNonConstModel(underlyingModel);
  }
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TempusSolver(
#else
template <typename Scalar>
Piro::TempusSolver<Scalar>::TempusSolver(
#endif
    const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
    Scalar initialTime,
    Scalar finalTime,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel,
    Teuchos::EVerbosityLevel verbosityLevel) :
  fwdStateIntegrator(stateIntegrator),
  fwdStateStepper(stateStepper),
  fwdTimeStepSolver(timeStepSolver),
  model(underlyingModel),
  initialConditionModel(icModel),
  t_initial(initialTime),
  t_final(finalTime),
  num_p(model->Np()),
  num_g(model->Ng()),
  computeSensitivities_(false),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel(verbosityLevel),
  isInitialized(true)
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //IKT, 12/5/16: the following exception check is needed until setInitialCondition method
  //is added to the Tempus::IntegratorBasic class.
  if (initialTime > 0.0) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TempusSolver: the constructor employed does not support initialTime > 0.0.  " <<
      "You have set initialTime = " << initialTime << "\n");
  }

  if (fwdStateStepper->getModel() != underlyingModel) {
    fwdStateStepper->setNonConstModel(underlyingModel);
  }
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_p_space(int l) const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TempusSolver<Scalar>::get_p_space(int l) const
#endif
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TempusSolver::get_p_map():  " <<
      "Invalid parameter index l = " <<
      l << "\n");

  return model->get_p_space(l);
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_g_space(int j) const
#else
template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::TempusSolver<Scalar>::get_g_space(int j) const
#endif
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      "\n Error in Piro::TempusSolver::get_g_map():  " <<
      "Invalid response index j = " <<
      j << "\n");

  if (j < num_g) {
    return model->get_g_space(j);
  } else {
    // j == num_g
    return model->get_x_space();
  }
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNominalValues() const
#else
template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TempusSolver<Scalar>::getNominalValues() const
#endif
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createInArgs() const
#else
template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::TempusSolver<Scalar>::createInArgs() const
#endif
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createOutArgsImpl() const
#else
template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::TempusSolver<Scalar>::createOutArgsImpl() const
#endif
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p, num_g + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();

  if (num_p > 0) {
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

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::evalModelImpl(
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::evalModelImpl(
#endif
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;
  using Teuchos::rcp;

  // TODO: Support more than 1 parameter and 1 response
  const int j = 0;
  const int l = 0;

  // Parse InArgs
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(l);
  }
  RCP<const Thyra::VectorBase<Scalar> > p_in2;  //JF add for multipoint
  if (num_p > 1) {
    p_in2 = inArgs.get_p(l+1);
  }

  // Parse OutArgs
  RCP<Thyra::VectorBase<Scalar> > g_out;
  if (num_g > 0) {
    g_out = outArgs.get_g(j);
  }
  const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);

  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = model->getNominalValues();

  // Set initial time in ME if needed

  if(t_initial > 0.0 && state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_t)) {
    state_ic.set_t(t_initial);
  }

  if (Teuchos::nonnull(initialConditionModel)) {
    // The initial condition depends on the parameter
    // It is found by querying the auxiliary model evaluator as the last response
    const RCP<Thyra::VectorBase<Scalar> > initialState =
      Thyra::createMember(model->get_x_space());

    {
      Thyra::ModelEvaluatorBase::InArgs<Scalar> initCondInArgs = initialConditionModel->createInArgs();
      if (num_p > 0) {
        initCondInArgs.set_p(l, inArgs.get_p(l));
      }

      Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs = initialConditionModel->createOutArgs();
      initCondOutArgs.set_g(initCondOutArgs.Ng() - 1, initialState);

      initialConditionModel->evalModel(initCondInArgs, initCondOutArgs);
    }

    state_ic.set_x(initialState);
  }

  // Set paramters p_in as part of initial conditions
  if (num_p > 0) {
    if (Teuchos::nonnull(p_in)) {
      state_ic.set_p(l, p_in);
    }
  }
  if (num_p > 1) { //JF added for multipoint
    if (Teuchos::nonnull(p_in2)) {
      state_ic.set_p(l+1, p_in2);
    }
  }

  //*out << "\nstate_ic:\n" << Teuchos::describe(state_ic, solnVerbLevel);

  //JF  may need a version of the following for multipoint, i.e. num_p>1, l+1, if we want sensitivities
  RCP<Thyra::MultiVectorBase<Scalar> > dgxdp_out;
  Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv_out;
  if (num_p > 0) {
    const Thyra::ModelEvaluatorBase::DerivativeSupport dgxdp_support =
      outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, num_g, l);
    if (dgxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgxdp_deriv =
        outArgs.get_DgDp(num_g, l);
      dgxdp_out = dgxdp_deriv.getMultiVector();
    }

    if (num_g > 0) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
      if (!dgdp_support.none()) {
        dgdp_deriv_out = outArgs.get_DgDp(j, l);
      }
    }
  }

  bool requestedSensitivities = true;
  if (computeSensitivities_ == true)
    requestedSensitivities = Teuchos::nonnull(dgxdp_out) || !dgdp_deriv_out.isEmpty();
  else
    requestedSensitivities = false;

  RCP<const Thyra::VectorBase<Scalar> > finalSolution;
  RCP<const Tempus::SolutionState<Scalar> > solutionState;
  RCP<const Tempus::SolutionHistory<Scalar> > solutionHistory;
  if (!requestedSensitivities)
  {
    //
    *out << "\nE) Solve the forward problem ...\n";
    //
    //
    *out << "T final requested: " << t_final << " \n";

    fwdStateIntegrator->advanceTime(t_final);

    double time = fwdStateIntegrator->getTime();

    *out << "T final actual: " << time << "\n";

    if (abs(time-t_final) > 1.0e-10) {
      if (abort_on_failure_ == true) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          Teuchos::Exceptions::InvalidParameter,
          "\n Error! Piro::TempusSolver: time-integrator did not make it to final time " <<
          "specified in Input File.  Final time in input file is " << t_final <<
          ", whereas actual final time is " << time << ".  If you'd like to " <<
          "suppress this exception, run with 'Abort on Failure' set to 'false' in " << 
          "Tempus sublist.\n" );
      }
      else {
         *out << "\n WARNING: Piro::TempusSolver did not make it to final time, but "
              << "solver will not abort since you have specified 'Abort on Failure' = 'false'.\n"; 
      }
    }

    finalSolution = fwdStateIntegrator->getX();

    solutionHistory = fwdStateIntegrator->getSolutionHistory();
    auto numStates = solutionHistory->getNumStates();
    solutionState = (*solutionHistory)[numStates-1];
    //Get final solution from solutionHistory.
    finalSolution = solutionState->getX();

    if (Teuchos::VERB_MEDIUM <= solnVerbLevel) {
      *out << "Final Solution\n" << *finalSolution << "\n";
    }

  }
  else {
    //
    *out << "\nE) Solve the forward problem with Sensitivities...\n";
    //
    TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        Teuchos::Exceptions::InvalidParameter,
        "\n Error! Piro::TempusSolver: sensitivities with Tempus are not yet supported!");
  }

  *out << "\nF) Check the solution to the forward problem ...\n";

  // As post-processing step, calculate responses at final solution
  {
    Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = model->createInArgs();
    {
      modelInArgs.set_x(finalSolution);
      if (num_p > 0) {
        modelInArgs.set_p(l, p_in);
      }
      if (num_p > 1) {  //JF added for multipoint
        modelInArgs.set_p(l+1, p_in2);
      }
      //Set time to be final time at which the solve occurs (< t_final in the case we don't make it to t_final).
      //IKT: get final time from solutionHistory workingSpace, which is different than how it is done in Piro::RythmosSolver class.
      //IKT, 11/1/16, FIXME? workingState pointer is null right now, so the following
      //code is commented out for now.  Use t_final and soln_dt in set_t instead for now.
      /*RCP<Tempus::SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();
      const Scalar time = workingState->getTime();
      const Scalar dt   = workingState->getTimeStep();
      const Scalar t = time + dt;
      modelInArgs.set_t(t);*/
      const Scalar soln_dt = solutionState->getTimeStep();
      modelInArgs.set_t(t_final - soln_dt);
    }

    Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
    if (Teuchos::nonnull(g_out)) {
      Thyra::put_scalar(Teuchos::ScalarTraits<Scalar>::zero(), g_out.ptr());
      modelOutArgs.set_g(j, g_out);
    }

    model->evalModel(modelInArgs, modelOutArgs);
  }

  // Return the final solution as an additional g-vector, if requested
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*finalSolution, gx_out.ptr());
  }
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::create_DgDp_op_impl(int j, int l) const
#else
template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::TempusSolver<Scalar>::create_DgDp_op_impl(int j, int l) const
#endif
{
  TEUCHOS_ASSERT(j != num_g);
  const Teuchos::Array<Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > > dummy =
    Teuchos::tuple(Thyra::zero<Scalar>(this->get_g_space(j), this->get_p_space(l)));
  return Teuchos::rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getValidTempusParameters() const
#else
template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::TempusSolver<Scalar>::getValidTempusParameters() const
#endif
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
    Teuchos::rcp(new Teuchos::ParameterList("ValidTempusSolverParams"));
  validPL->sublist("Tempus", false, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->sublist("NonLinear Solver", false, "");
  //validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");
  validPL->set<bool>("Lump Mass Matrix", false, "");
  validPL->set<bool>("Abort on Failure", true, "");
  validPL->set<std::string>("Integrator Name", "Tempus Integrator", "");
  validPL->sublist("Tempus Integrator", false, "");
  validPL->sublist("Tempus Stepper", false, "");
  validPL->sublist("Time Step Control", false, "");
  return validPL;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
addStepperFactory(const std::string & stepperName,const Teuchos::RCP<Piro::TempusStepperFactory<Scalar> > & factory)
{
  stepperFactories[stepperName] = factory;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
addStepControlFactory(const std::string & stepControlName,
                      const Teuchos::RCP<Piro::TempusStepControlFactory<Scalar>> & step_control_strategy)
{
  stepControlFactories[stepControlName] = step_control_strategy;
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
setStartTime(const Scalar start_time)
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc_const = fwdStateIntegrator->getTimeStepControl();
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc = Teuchos::rcp_const_cast<Tempus::TimeStepControl<Scalar> >(tsc_const); 
  tsc->setInitTime(start_time); 
} 

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Scalar Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
Scalar Piro::TempusSolver<Scalar>::
#endif
getStartTime() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc = fwdStateIntegrator->getTimeStepControl();
  Scalar start_time = tsc->getInitTime(); 
  return start_time; 
} 

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
setFinalTime(const Scalar final_time)
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc_const = fwdStateIntegrator->getTimeStepControl();
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc = Teuchos::rcp_const_cast<Tempus::TimeStepControl<Scalar> >(tsc_const); 
  t_final = final_time; 
  tsc->setFinalTime(final_time); 
} 

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Scalar Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
Scalar Piro::TempusSolver<Scalar>::
#endif
getFinalTime() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc = fwdStateIntegrator->getTimeStepControl();
  Scalar final_time = tsc->getFinalTime(); 
  return final_time; 
} 

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
setInitTimeStep(const Scalar init_time_step)
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc_const = fwdStateIntegrator->getTimeStepControl();
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc = Teuchos::rcp_const_cast<Tempus::TimeStepControl<Scalar> >(tsc_const); 
  tsc->setInitTimeStep(init_time_step); 
} 


#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Scalar Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
Scalar Piro::TempusSolver<Scalar>::
#endif
getInitTimeStep() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > tsc = fwdStateIntegrator->getTimeStepControl();
  auto init_time_step = tsc->getInitTimeStep(); 
  return init_time_step; 
} 
#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
setObserver()
{
  Teuchos::RCP<Tempus::IntegratorObserverBasic<Scalar> > observer = Teuchos::null;
  if (Teuchos::nonnull(piroObserver_)) {
    //Get solutionHistory from integrator
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > solutionHistory = fwdStateIntegrator->getSolutionHistory();
    const Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > timeStepControl = fwdStateIntegrator->getTimeStepControl();
    //Create Tempus::IntegratorObserverBasic object
    observer = Teuchos::rcp(new ObserverToTempusIntegrationObserverAdapter<Scalar>(solutionHistory,
                                timeStepControl, piroObserver_, supports_x_dotdot_));
  }
  if (Teuchos::nonnull(observer)) {
    //Set observer in integrator
    fwdStateIntegrator->getObserver()->clearObservers();
    fwdStateIntegrator->setObserver(observer);
    //Reinitialize everything in integrator class, since we have changed the observer.
    fwdStateIntegrator->initialize();
  }
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
setInitialState(Scalar t0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot0) 
{
   fwdStateIntegrator->setInitialState(t0, x0, xdot0, xdotdot0); 
   //Reset observer.  This is necessary for correct observation of solution
   //since setInitialState modifies the solutionHistory object.
   setObserver(); 
 
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
void Piro::TempusSolver<Scalar>::
#endif
setInitialGuess(Teuchos::RCP< const Thyra::VectorBase<Scalar> > initial_guess) 
{
   fwdStateStepper->setInitialGuess(initial_guess); 
}

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Tempus::SolutionHistory<Scalar> > Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
Teuchos::RCP<Tempus::SolutionHistory<Scalar> > Piro::TempusSolver<Scalar>::
#endif
getSolutionHistory() const
{
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > soln_history_const = fwdStateIntegrator->getSolutionHistory();
  Teuchos::RCP<Tempus::SolutionHistory<Scalar> > soln_history = Teuchos::rcp_const_cast<Tempus::SolutionHistory<Scalar> >(soln_history_const); 
  return soln_history;
}
  

#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > Piro::TempusSolver<Scalar>::
#endif
getSolver() const
{
  return fwdStateStepper->getSolver(); 
}


#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Tempus::Status Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
#else
template <typename Scalar>
Tempus::Status Piro::TempusSolver<Scalar>::
#endif
getTempusIntegratorStatus() const
{
  return fwdStateIntegrator->getStatus(); 
}


#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<Piro::TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#else
template <typename Scalar>
Teuchos::RCP<Piro::TempusSolver<Scalar> >
#endif
Piro::tempusSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out(Teuchos::VerboseObjectBase::getDefaultOStream());
#ifdef DEBUT_OUTPUT
  *out << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
   bool computeSensitivities = true;
   if (appParams->isSublist("Analysis")) {
     Teuchos::ParameterList& analysisPL = appParams->sublist("Analysis");
     if (analysisPL.isParameter("Compute Sensitivities"))
       computeSensitivities = analysisPL.get<bool>("Compute Sensitivities");
   }

#ifdef ALBANY_BUILD
  return Teuchos::rcp(new TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>(appParams, in_model, computeSensitivities, piroObserver));
#else
  return Teuchos::rcp(new TempusSolver<Scalar>(appParams, in_model, computeSensitivities, piroObserver));
#endif

}
