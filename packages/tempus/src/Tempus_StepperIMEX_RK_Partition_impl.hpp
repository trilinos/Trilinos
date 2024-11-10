//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperIMEX_RK_Partition_impl_hpp
#define Tempus_StepperIMEX_RK_Partition_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"

namespace Tempus {

template <class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
    std::string stepperType)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperType != "Partitioned IMEX RK 1st order" &&
          stepperType != "Partitioned IMEX RK SSP2" &&
          stepperType != "Partitioned IMEX RK ARS 233" &&
          stepperType != "General Partitioned IMEX RK",
      std::logic_error,
      "  'Stepper Type' (='"
          << stepperType
          << "')\n"
          << "  does not match one of the types for this Stepper:\n"
          << "    'Partitioned IMEX RK 1st order'\n"
          << "    'Partitioned IMEX RK SSP2'\n"
          << "    'Partitioned IMEX RK ARS 233'\n"
          << "    'General Partitioned IMEX RK'\n");

  this->setStepperName("Partitioned IMEX RK SSP2");
  this->setStepperType("Partitioned IMEX RK SSP2");
  this->setUseFSAL(false);
  this->setICConsistency("None");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);

  this->setStageNumber(-1);

  this->setTableaus(stepperType);
  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}

template <class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    std::string stepperType,
    Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau,
    Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau, Scalar order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperType != "Partitioned IMEX RK 1st order" &&
          stepperType != "Partitioned IMEX RK SSP2" &&
          stepperType != "Partitioned IMEX RK ARS 233" &&
          stepperType != "General Partitioned IMEX RK",
      std::logic_error,
      "  'Stepper Type' (='"
          << stepperType
          << "')\n"
          << "  does not match one of the types for this Stepper:\n"
          << "    'Partitioned IMEX RK 1st order'\n"
          << "    'Partitioned IMEX RK SSP2'\n"
          << "    'Partitioned IMEX RK ARS 233'\n"
          << "    'General Partitioned IMEX RK'\n");

  this->setStepperName(stepperType);
  this->setStepperType(stepperType);
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);
  this->setOrder(order);

  this->setStageNumber(-1);

  if (stepperType == "General Partitioned IMEX RK") {
    this->setExplicitTableau(explicitTableau);
    this->setImplicitTableau(implicitTableau);
  }
  else {
    this->setTableaus(stepperType);
  }
  this->setAppAction(stepperRKAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setTableaus(
    std::string stepperType,
    Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau,
    Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau)
{
  if (stepperType == "") stepperType = "Partitioned IMEX RK SSP2";

  if (stepperType == "Partitioned IMEX RK 1st order") {
    {
      // Explicit Tableau
      typedef Teuchos::ScalarTraits<Scalar> ST;
      int NumStages = 2;
      Teuchos::SerialDenseMatrix<int, Scalar> A(NumStages, NumStages);
      Teuchos::SerialDenseVector<int, Scalar> b(NumStages);
      Teuchos::SerialDenseVector<int, Scalar> c(NumStages);
      const Scalar one  = ST::one();
      const Scalar zero = ST::zero();

      // Fill A:
      A(0, 0) = zero;
      A(0, 1) = zero;
      A(1, 0) = one;
      A(1, 1) = zero;

      // Fill b:
      b(0) = one;
      b(1) = zero;

      // Fill c:
      c(0) = zero;
      c(1) = one;

      int order = 1;

      auto expTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
          "Explicit Tableau - Partitioned IMEX RK 1st order", A, b, c, order,
          order, order));
      expTableau->setTVD(true);
      expTableau->setTVDCoeff(2.0);

      this->setExplicitTableau(expTableau);
    }
    {
      // Implicit Tableau
      typedef Teuchos::ScalarTraits<Scalar> ST;
      int NumStages        = 2;
      const Scalar sspcoef = std::numeric_limits<Scalar>::max();
      Teuchos::SerialDenseMatrix<int, Scalar> A(NumStages, NumStages);
      Teuchos::SerialDenseVector<int, Scalar> b(NumStages);
      Teuchos::SerialDenseVector<int, Scalar> c(NumStages);
      const Scalar one  = ST::one();
      const Scalar zero = ST::zero();

      // Fill A:
      A(0, 0) = zero;
      A(0, 1) = zero;
      A(1, 0) = zero;
      A(1, 1) = one;

      // Fill b:
      b(0) = zero;
      b(1) = one;

      // Fill c:
      c(0) = zero;
      c(1) = one;

      int order = 1;

      auto impTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
          "Implicit Tableau - Partitioned IMEX RK 1st order", A, b, c, order,
          order, order));
      impTableau->setTVD(true);
      impTableau->setTVDCoeff(sspcoef);

      this->setImplicitTableau(impTableau);
    }
    this->setStepperName("Partitioned IMEX RK 1st order");
    this->setStepperType("Partitioned IMEX RK 1st order");
    this->setOrder(1);
  }
  else if (stepperType == "Partitioned IMEX RK SSP2") {
    // Explicit Tableau
    auto stepperERK = Teuchos::rcp(new StepperERK_Trapezoidal<Scalar>());
    this->setExplicitTableau(stepperERK->getTableau());

    // Implicit Tableau
    auto stepperSDIRK = Teuchos::rcp(new StepperSDIRK_2Stage3rdOrder<Scalar>());
    stepperSDIRK->setGammaType("2nd Order L-stable");
    this->setImplicitTableau(stepperSDIRK->getTableau());

    this->setStepperName("Partitioned IMEX RK SSP2");
    this->setStepperType("Partitioned IMEX RK SSP2");
    this->setOrder(2);
  }
  else if (stepperType == "Partitioned IMEX RK ARS 233") {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int, Scalar> A(NumStages, NumStages);
    Teuchos::SerialDenseVector<int, Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int, Scalar> c(NumStages);
    const Scalar one     = ST::one();
    const Scalar zero    = ST::zero();
    const Scalar onehalf = ST::one() / (2 * ST::one());
    const Scalar gamma   = (3 * one + ST::squareroot(3 * one)) / (6 * one);
    {
      // Explicit Tableau
      // Fill A:
      A(0, 0) = zero;
      A(0, 1) = zero;
      A(0, 2) = zero;
      A(1, 0) = gamma;
      A(1, 1) = zero;
      A(1, 2) = zero;
      A(2, 0) = (gamma - 1.0);
      A(2, 1) = (2.0 - 2.0 * gamma);
      A(2, 2) = zero;

      // Fill b:
      b(0) = zero;
      b(1) = onehalf;
      b(2) = onehalf;

      // Fill c:
      c(0) = zero;
      c(1) = gamma;
      c(2) = one - gamma;

      int order = 2;

      auto expTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
          "Explicit Tableau - Partitioned IMEX RK ARS 233", A, b, c, order,
          order, order));

      this->setExplicitTableau(expTableau);
    }
    {
      // Implicit Tableau
      // Fill A:
      A(0, 0) = zero;
      A(0, 1) = zero;
      A(0, 2) = zero;
      A(1, 0) = zero;
      A(1, 1) = gamma;
      A(1, 2) = zero;
      A(2, 0) = zero;
      A(2, 1) = (1.0 - 2.0 * gamma);
      A(2, 2) = gamma;

      // Fill b:
      b(0) = zero;
      b(1) = onehalf;
      b(2) = onehalf;

      // Fill c:
      c(0) = zero;
      c(1) = gamma;
      c(2) = one - gamma;

      int order = 3;

      auto impTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
          "Implicit Tableau - Partitioned IMEX RK ARS 233", A, b, c, order,
          order, order));

      this->setImplicitTableau(impTableau);
    }
    this->setStepperName("Partitioned IMEX RK ARS 233");
    this->setStepperType("Partitioned IMEX RK ARS 233");
    this->setOrder(3);
  }
  else if (stepperType == "General Partitioned IMEX RK") {
    if (explicitTableau == Teuchos::null) {
      // Default Explicit Tableau (i.e., Partitioned IMEX RK SSP2)
      auto stepperERK = Teuchos::rcp(new StepperERK_Trapezoidal<Scalar>());
      this->setExplicitTableau(stepperERK->getTableau());
    }
    else {
      this->setExplicitTableau(explicitTableau);
    }

    if (implicitTableau == Teuchos::null) {
      // Default Implicit Tablea (i.e., Partitioned IMEX RK SSP2)
      auto stepperSDIRK =
          Teuchos::rcp(new StepperSDIRK_2Stage3rdOrder<Scalar>());
      stepperSDIRK->setGammaType("2nd Order L-stable");
      this->setImplicitTableau(stepperSDIRK->getTableau());
    }
    else {
      this->setImplicitTableau(implicitTableau);
    }

    this->setStepperName("General Partitioned IMEX RK");
    this->setStepperType("General Partitioned IMEX RK");
    this->setOrder(1);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error - Not a valid StepperIMEX_RK_Partition type!  Stepper Type = "
            << stepperType << "\n"
            << "  Current valid types are: \n"
            << "      'Partitioned IMEX RK 1st order'\n"
            << "      'Partitioned IMEX RK SSP2'\n"
            << "      'Partitioned IMEX RK ARS 233'\n"
            << "      'General Partitioned IMEX RK'\n");
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      explicitTableau_ == Teuchos::null, std::runtime_error,
      "Error - StepperIMEX_RK_Partition - Explicit tableau is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(
      implicitTableau_ == Teuchos::null, std::runtime_error,
      "Error - StepperIMEX_RK_Partition - Implicit tableau is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(
      explicitTableau_->numStages() != implicitTableau_->numStages(),
      std::runtime_error,
      "Error - StepperIMEX_RK_Partition - Number of stages do not match!\n"
          << "  Explicit tableau = " << explicitTableau_->description() << "\n"
          << "    number of stages = " << explicitTableau_->numStages() << "\n"
          << "  Implicit tableau = " << implicitTableau_->description() << "\n"
          << "    number of stages = " << implicitTableau_->numStages()
          << "\n");

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setTableausPartition(
    Teuchos::RCP<Teuchos::ParameterList> pl, std::string stepperType)
{
  using Teuchos::RCP;
  if (stepperType == "") {
    if (pl == Teuchos::null)
      stepperType = "Partitioned IMEX RK SSP2";
    else
      stepperType =
          pl->get<std::string>("Stepper Type", "Partitioned IMEX RK SSP2");
  }

  if (stepperType != "General Partitioned IMEX RK") {
    this->setTableaus(stepperType);
  }
  else {
    if (pl != Teuchos::null) {
      Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau;
      Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau;
      if (pl->isSublist("IMEX-RK Explicit Stepper")) {
        RCP<Teuchos::ParameterList> explicitPL =
            Teuchos::rcp(new Teuchos::ParameterList(
                pl->sublist("IMEX-RK Explicit Stepper")));
        auto sf          = Teuchos::rcp(new StepperFactory<Scalar>());
        auto stepperTemp = sf->createStepper(explicitPL, Teuchos::null);
        auto stepperERK  = Teuchos::rcp_dynamic_cast<StepperExplicitRK<Scalar> >(
            stepperTemp, true);
        TEUCHOS_TEST_FOR_EXCEPTION(
            stepperERK == Teuchos::null, std::logic_error,
            "Error - The explicit component of a general partitioned IMEX RK "
            "stepper was not specified as an ExplicitRK stepper");
        explicitTableau = stepperERK->getTableau();
      }

      if (pl->isSublist("IMEX-RK Implicit Stepper")) {
        RCP<Teuchos::ParameterList> implicitPL =
            Teuchos::rcp(new Teuchos::ParameterList(
                pl->sublist("IMEX-RK Implicit Stepper")));
        auto sf          = Teuchos::rcp(new StepperFactory<Scalar>());
        auto stepperTemp = sf->createStepper(implicitPL, Teuchos::null);
        auto stepperDIRK =
            Teuchos::rcp_dynamic_cast<StepperDIRK<Scalar> >(stepperTemp, true);
        TEUCHOS_TEST_FOR_EXCEPTION(
            stepperDIRK == Teuchos::null, std::logic_error,
            "Error - The implicit component of a general partitioned IMEX RK "
            "stepper was not specified as an DIRK stepper");
        implicitTableau = stepperDIRK->getTableau();
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
          !(explicitTableau != Teuchos::null &&
            implicitTableau != Teuchos::null),
          std::logic_error,
          "Error - A parameter list was used to setup a general partitioned "
          "IMEX RK stepper, but did not "
          "specify both an explicit and an implicit tableau!\n");

      this->setTableaus(stepperType, explicitTableau, implicitTableau);

      this->setOrder(pl->get<int>("overall order", 1));
    }
  }
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setExplicitTableau(
    Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      explicitTableau->isImplicit() == true, std::logic_error,
      "Error - Received an implicit Tableau for setExplicitTableau()!\n"
          << "        Tableau = " << explicitTableau->description() << "\n");
  explicitTableau_ = explicitTableau;

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setImplicitTableau(
    Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      implicitTableau->isDIRK() != true, std::logic_error,
      "Error - Did not receive a DIRK Tableau for setImplicitTableau()!\n"
          << "        Tableau = " << implicitTableau->description() << "\n");
  implicitTableau_ = implicitTableau;

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  RCP<Thyra::ModelEvaluator<Scalar> > ncModel =
      rcp_const_cast<Thyra::ModelEvaluator<Scalar> >(appModel);
  RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> > modelPairIMEX =
      rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(
          ncModel);
  TEUCHOS_TEST_FOR_EXCEPTION(
      modelPairIMEX == Teuchos::null, std::logic_error,
      "Error - StepperIMEX_RK::setModel() was given a ModelEvaluator that\n"
          << "  could not be cast to a WrapperModelEvaluatorPairPartIMEX_Basic!\n"
          << "  From: " << appModel << "\n  To  : " << modelPairIMEX
          << "\n  Likely have given the wrong ModelEvaluator to this Stepper.\n");

  setModelPair(modelPairIMEX);

  this->isInitialized_ = false;
}

/** \brief Create WrapperModelPairIMEX from user-supplied ModelEvaluator pair.
 *
 *  The user-supplied ME pair can contain any user-specific IMEX interactions
 *  between explicit and implicit MEs.
 */
template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setModelPair(
    const Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >&
        mePairIMEX)
{
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      wrapperModelPairPartIMEX = Teuchos::rcp_dynamic_cast<
          WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(
          this->wrapperModel_);
  validExplicitODE(mePairIMEX->getExplicitModel());
  validImplicitODE_DAE(mePairIMEX->getImplicitModel());
  wrapperModelPairPartIMEX = mePairIMEX;
  wrapperModelPairPartIMEX->initialize();

  this->wrapperModel_ = wrapperModelPairPartIMEX;

  this->isInitialized_ = false;
}

/** \brief Create WrapperModelPairIMEX from explicit/implicit ModelEvaluators.
 *
 *  Use the supplied explicit/implicit MEs to create a WrapperModelPairIMEX
 *  with basic IMEX interactions between explicit and implicit MEs.
 */
template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setModelPair(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel)
{
  validExplicitODE(explicitModel);
  validImplicitODE_DAE(implicitModel);
  this->wrapperModel_ =
      Teuchos::rcp(new WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>(
          explicitModel, implicitModel));

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::initialize()
{
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      wrapperModelPairPartIMEX = Teuchos::rcp_dynamic_cast<
          WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(
          this->wrapperModel_);
  TEUCHOS_TEST_FOR_EXCEPTION(
      wrapperModelPairPartIMEX == Teuchos::null, std::logic_error,
      "Error - Can not cast the wrapper Model Evaluator to a IMEX Model Pair."
      "StepperIMEX_RK_Partition::initialize()\n");

  // Initialize the stage vectors
  const int numStages = explicitTableau_->numStages();
  stageF_.resize(numStages);
  stageGx_.resize(numStages);
  for (int i = 0; i < numStages; i++) {
    stageF_[i] = Thyra::createMember(
        wrapperModelPairPartIMEX->getExplicitModel()->get_f_space());
    stageGx_[i] = Thyra::createMember(
        wrapperModelPairPartIMEX->getImplicitModel()->get_f_space());
    assign(stageF_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
    assign(stageGx_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  xTilde_ = Thyra::createMember(
      wrapperModelPairPartIMEX->getImplicitModel()->get_x_space());
  assign(xTilde_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

  StepperImplicit<Scalar>::initialize();
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(
      numStates < 1, std::logic_error,
      "Error - setInitialConditions() needs at least one SolutionState\n"
      "        to set the initial condition.  Number of States = "
          << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out, 1, "StepperIMEX_RK::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"
         << std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x        = initialState->getX();

  // Use x from inArgs as ICs, if needed.
  auto inArgs = this->wrapperModel_->getNominalValues();
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        (x == Teuchos::null) && (inArgs.get_x() == Teuchos::null),
        std::logic_error,
        "Error - setInitialConditions() needs the ICs from the "
        "SolutionHistory\n        or getNominalValues()!\n");

    x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
    initialState->setX(x);
  }

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  TEUCHOS_TEST_FOR_EXCEPTION(
      icConsistency != "None", std::logic_error,
      "Error - setInitialConditions() requested a consistency of '"
          << icConsistency
          << "'.\n        But only 'None' is available for IMEX-RK!\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
      this->getUseFSAL(), std::logic_error,
      "Error - The First-Same-As-Last (FSAL) principle is not "
          << "available for IMEX-RK.  Set useFSAL=false.\n");
}

template <typename Scalar>
void StepperIMEX_RK_Partition<Scalar>::evalImplicitModelExplicitly(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& X,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& Y, Scalar time,
    Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& G) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      wrapperModelPairPartIMEX = Teuchos::rcp_dynamic_cast<
          WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(
          this->wrapperModel_);
  MEB::InArgs<Scalar> inArgs = wrapperModelPairPartIMEX->getInArgs();
  inArgs.set_x(X);
  inArgs.set_p(wrapperModelPairPartIMEX->getParameterIndex(), Y);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(stepSize);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(stageNumber);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);

  MEB::OutArgs<Scalar> outArgs = wrapperModelPairPartIMEX->getOutArgs();
  outArgs.set_f(G);

  wrapperModelPairPartIMEX->getImplicitModel()->evalModel(inArgs, outArgs);
  Thyra::Vt_S(G.ptr(), -1.0);
}

template <typename Scalar>
void StepperIMEX_RK_Partition<Scalar>::evalExplicitModel(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& Z, Scalar time,
    Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& F) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      wrapperModelPairPartIMEX = Teuchos::rcp_dynamic_cast<
          WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(
          this->wrapperModel_);
  MEB::InArgs<Scalar> inArgs =
      wrapperModelPairPartIMEX->getExplicitModel()->createInArgs();
  inArgs.set_x(Z);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(stepSize);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(stageNumber);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);

  MEB::OutArgs<Scalar> outArgs =
      wrapperModelPairPartIMEX->getExplicitModel()->createOutArgs();
  outArgs.set_f(F);

  wrapperModelPairPartIMEX->getExplicitModel()->evalModel(inArgs, outArgs);
  Thyra::Vt_S(F.ptr(), -1.0);
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::SerialDenseVector;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperIMEX_RK_Partition::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperIMEX_RK_Partition<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for IMEX_RK_Partition.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\nTry setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    const Scalar dt   = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages                        = explicitTableau_->numStages();
    const SerialDenseMatrix<int, Scalar>& AHat = explicitTableau_->A();
    const SerialDenseVector<int, Scalar>& bHat = explicitTableau_->b();
    const SerialDenseVector<int, Scalar>& cHat = explicitTableau_->c();
    const SerialDenseMatrix<int, Scalar>& A    = implicitTableau_->A();
    const SerialDenseVector<int, Scalar>& b    = implicitTableau_->b();
    const SerialDenseVector<int, Scalar>& c    = implicitTableau_->c();

    Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
        wrapperModelPairPartIMEX = Teuchos::rcp_dynamic_cast<
            WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(
            this->wrapperModel_);

    bool pass = true;
    Thyra::assign(workingState->getX().ptr(), *(currentState->getX()));
    RCP<Thyra::VectorBase<Scalar> > stageY =
        wrapperModelPairPartIMEX->getExplicitOnlyVector(workingState->getX());
    RCP<Thyra::VectorBase<Scalar> > stageX =
        wrapperModelPairPartIMEX->getIMEXVector(workingState->getX());

    RCP<StepperIMEX_RK_Partition<Scalar> > thisStepper =
        Teuchos::rcpFromRef(*this);
    this->stepperRKAppAction_->execute(
        solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    // Compute stage solutions
    for (int i = 0; i < numStages; ++i) {
      this->setStageNumber(i);

      Thyra::assign(stageY.ptr(),
                    *(wrapperModelPairPartIMEX->getExplicitOnlyVector(
                        currentState->getX())));
      Thyra::assign(
          xTilde_.ptr(),
          *(wrapperModelPairPartIMEX->getIMEXVector(currentState->getX())));
      for (int j = 0; j < i; ++j) {
        if (AHat(i, j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          RCP<Thyra::VectorBase<Scalar> > stageFy =
              wrapperModelPairPartIMEX->getExplicitOnlyVector(stageF_[j]);
          RCP<Thyra::VectorBase<Scalar> > stageFx =
              wrapperModelPairPartIMEX->getIMEXVector(stageF_[j]);
          Thyra::Vp_StV(stageY.ptr(), -dt * AHat(i, j), *stageFy);
          Thyra::Vp_StV(xTilde_.ptr(), -dt * AHat(i, j), *stageFx);
        }
        if (A(i, j) != Teuchos::ScalarTraits<Scalar>::zero())
          Thyra::Vp_StV(xTilde_.ptr(), -dt * A(i, j), *(stageGx_[j]));
      }

      this->stepperRKAppAction_->execute(
          solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STAGE);

      Scalar ts    = time + c(i) * dt;
      Scalar tHats = time + cHat(i) * dt;
      if (A(i, i) == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Explicit stage for the ImplicitODE_DAE
        bool isNeeded = false;
        for (int k = i + 1; k < numStages; ++k)
          if (A(k, i) != 0.0) isNeeded = true;
        if (b(i) != 0.0) isNeeded = true;
        if (isNeeded == false) {
          // stageGx_[i] is not needed.
          assign(stageGx_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
        }
        else {
          Thyra::assign(stageX.ptr(), *xTilde_);
          evalImplicitModelExplicitly(stageX, stageY, ts, dt, i, stageGx_[i]);
        }
      }
      else {
        // Implicit stage for the ImplicitODE_DAE
        const Scalar alpha = Scalar(1.0) / (dt * A(i, i));
        const Scalar beta  = Scalar(1.0);

        // Setup TimeDerivative
        RCP<TimeDerivative<Scalar> > timeDer =
            Teuchos::rcp(new StepperIMEX_RKPartTimeDerivative<Scalar>(
                alpha, xTilde_.getConst()));

        auto p = Teuchos::rcp(new ImplicitODEParameters<Scalar>(
            timeDer, dt, alpha, beta, SOLVE_FOR_X, i));

        this->stepperRKAppAction_->execute(
            solutionHistory, thisStepper,
            StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

        wrapperModelPairPartIMEX->setUseImplicitModel(true);
        this->solver_->setModel(wrapperModelPairPartIMEX);

        Thyra::SolveStatus<Scalar> sStatus;
        if (wrapperModelPairPartIMEX->getParameterIndex() >= 0)
          sStatus = this->solveImplicitODE(
              stageX, stageGx_[i], ts, p, stageY,
              wrapperModelPairPartIMEX->getParameterIndex());
        else
          sStatus = this->solveImplicitODE(stageX, stageGx_[i], ts, p);

        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass = false;

        wrapperModelPairPartIMEX->setUseImplicitModel(false);

        this->stepperRKAppAction_->execute(
            solutionHistory, thisStepper,
            StepperRKAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);

        // Update contributions to stage values
        Thyra::V_StVpStV(stageGx_[i].ptr(), -alpha, *stageX, alpha, *xTilde_);
      }

      this->stepperRKAppAction_->execute(
          solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);
      evalExplicitModel(workingState->getX(), tHats, dt, i, stageF_[i]);
      this->stepperRKAppAction_->execute(
          solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STAGE);
    }

    // Sum for solution: y_n = y_n-1 - dt*Sum{ bHat(i)*fy(i)            }
    // Sum for solution: x_n = x_n-1 - dt*Sum{ bHat(i)*fx(i) + b(i)*gx(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    RCP<Thyra::VectorBase<Scalar> > Z = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > X =
        wrapperModelPairPartIMEX->getIMEXVector(Z);
    for (int i = 0; i < numStages; ++i) {
      if (bHat(i) != Teuchos::ScalarTraits<Scalar>::zero())
        Thyra::Vp_StV(Z.ptr(), -dt * bHat(i), *(stageF_[i]));
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero())
        Thyra::Vp_StV(X.ptr(), -dt * b(i), *(stageGx_[i]));
    }

    if (pass == true)
      workingState->setSolutionStatus(Status::PASSED);
    else
      workingState->setSolutionStatus(Status::FAILED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    this->stepperRKAppAction_->execute(
        solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  // reset the stage number
  this->setStageNumber(-1);
  return;
}

/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template <class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperIMEX_RK_Partition<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperIMEX_RK_Partition ---\n";
  out << "  explicitTableau_   = " << explicitTableau_ << std::endl;
  if (verbLevel == Teuchos::VERB_HIGH)
    explicitTableau_->describe(out, verbLevel);
  out << "  implicitTableau_   = " << implicitTableau_ << std::endl;
  if (verbLevel == Teuchos::VERB_HIGH)
    implicitTableau_->describe(out, verbLevel);
  out << "  xTilde_            = " << xTilde_ << std::endl;
  out << "  stageF_.size()     = " << stageF_.size() << std::endl;
  int numStages = stageF_.size();
  for (int i = 0; i < numStages; ++i)
    out << "    stageF_[" << i << "] = " << stageF_[i] << std::endl;
  out << "  stageGx_.size()    = " << stageGx_.size() << std::endl;
  numStages = stageGx_.size();
  for (int i = 0; i < numStages; ++i)
    out << "    stageGx_[" << i << "] = " << stageGx_[i] << std::endl;
  out << "  stepperRKAppAction_= " << this->stepperRKAppAction_ << std::endl;
  out << "  order_             = " << order_ << std::endl;
  out << "--------------------------------" << std::endl;
}

template <class Scalar>
bool StepperIMEX_RK_Partition<Scalar>::isValidSetup(
    Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);

  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;

  Teuchos::RCP<WrapperModelEvaluatorPairIMEX<Scalar> > wrapperModelPairIMEX =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairIMEX<Scalar> >(
          this->wrapperModel_);

  if (wrapperModelPairIMEX->getExplicitModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The explicit ModelEvaluator is not set!\n";
  }

  if (wrapperModelPairIMEX->getImplicitModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The implicit ModelEvaluator is not set!\n";
  }

  if (this->wrapperModel_ == Teuchos::null) {
    isValidSetup = false;
    out << "The wrapper ModelEvaluator is not set!\n";
  }

  if (this->solver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The solver is not set!\n";
  }

  if (this->stepperRKAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The AppAction is not set!\n";
  }

  if (explicitTableau_ == Teuchos::null) {
    isValidSetup = false;
    out << "The explicit tableau is not set!\n";
  }

  if (implicitTableau_ == Teuchos::null) {
    isValidSetup = false;
    out << "The implicit tableau is not set!\n";
  }

  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperIMEX_RK_Partition<Scalar>::getValidParameters() const
{
  auto pl = this->getValidParametersBasicImplicit();
  pl->template set<int>("overall order", this->getOrder());

  auto explicitStepper = Teuchos::rcp(new StepperERK_General<Scalar>());
  explicitStepper->setTableau(
      explicitTableau_->A(), explicitTableau_->b(), explicitTableau_->c(),
      explicitTableau_->order(), explicitTableau_->orderMin(),
      explicitTableau_->orderMax(), explicitTableau_->bstar());
  pl->set("IMEX-RK Explicit Stepper", *explicitStepper->getValidParameters());

  auto implicitStepper = Teuchos::rcp(new StepperERK_General<Scalar>());
  implicitStepper->setTableau(
      implicitTableau_->A(), implicitTableau_->b(), implicitTableau_->c(),
      implicitTableau_->order(), implicitTableau_->orderMin(),
      implicitTableau_->orderMax(), implicitTableau_->bstar());
  pl->set("IMEX-RK Implicit Stepper", *implicitStepper->getValidParameters());

  return pl;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperIMEX_RK_Partition<Scalar> > createStepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper =
      Teuchos::rcp(new StepperIMEX_RK_Partition<Scalar>(stepperType));
  stepper->setStepperImplicitValues(pl);
  stepper->setTableausPartition(pl, stepperType);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperIMEX_RK_Partition_impl_hpp
