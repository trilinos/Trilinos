// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperIMEX_RK_Partition_impl_hpp
#define Tempus_StepperIMEX_RK_Partition_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_Basic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


template<class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition()
{
  this->setStepperType(        "Partitioned IMEX RK SSP2");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());
  this->setZeroInitialGuess(   false);

  this->setStageNumber(-1);

  this->setTableaus("Partitioned IMEX RK SSP2");
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  this->setObserver();
#endif
  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}


#ifndef TEMPUS_HIDE_DEPRECATED_CODE
template<class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  std::string stepperType,
  Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau,
  Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau,
  Scalar order)
{
  this->setStepperType(        stepperType);
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);

  this->setStageNumber(-1);

  this->setExplicitTableau(explicitTableau);
  this->setImplicitTableau(implicitTableau);
  this->setObserver(obs);
  this->setAppAction(Teuchos::null);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}
#endif
template<class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
  std::string stepperType,
  Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau,
  Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau,
  Scalar order)
{
  this->setStepperType(        stepperType);
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);
  this->setOrder(              order);

  this->setStageNumber(-1);

  this->setExplicitTableau(explicitTableau);
  this->setImplicitTableau(implicitTableau);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  this->setObserver(Teuchos::null);
#endif
  this->setAppAction(stepperRKAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}


template<class Scalar>
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
      const bool isTVD = true;
      const Scalar sspcoef = 2.0;
      Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
      Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
      Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
      const Scalar one = ST::one();
      const Scalar zero = ST::zero();

      // Fill A:
      A(0,0) = zero; A(0,1) = zero;
      A(1,0) =  one; A(1,1) = zero;

      // Fill b:
      b(0) = one; b(1) = zero;

      // Fill c:
      c(0) = zero; c(1) = one;

      int order = 1;

      auto expTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Explicit Tableau - Partitioned IMEX RK 1st order",
        A,b,c,order,order,order,isTVD,sspcoef));

      this->setExplicitTableau(expTableau);
    }
    {
      // Implicit Tableau
      typedef Teuchos::ScalarTraits<Scalar> ST;
      int NumStages = 2;
      const bool isTVD = true;
      const Scalar sspcoef =  std::numeric_limits<Scalar>::max();
      Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
      Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
      Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
      const Scalar one = ST::one();
      const Scalar zero = ST::zero();

      // Fill A:
      A(0,0) = zero; A(0,1) = zero;
      A(1,0) = zero; A(1,1) =  one;

      // Fill b:
      b(0) = zero; b(1) = one;

      // Fill c:
      c(0) = zero; c(1) = one;

      int order = 1;

      auto impTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Implicit Tableau - Partitioned IMEX RK 1st order",
        A,b,c,order,order,order,isTVD,sspcoef));

      this->setImplicitTableau(impTableau);
    }
    this->setStepperType("Partitioned IMEX RK 1st order");
    this->setOrder(1);

  } else if (stepperType == "Partitioned IMEX RK SSP2") {
    // Explicit Tableau
    auto stepperERK = Teuchos::rcp(new StepperERK_Trapezoidal<Scalar>());
    this->setExplicitTableau(stepperERK->getTableau());

    // Implicit Tableau
    auto stepperSDIRK = Teuchos::rcp(new StepperSDIRK_2Stage3rdOrder<Scalar>());
    stepperSDIRK->setGammaType("2nd Order L-stable");
    this->setImplicitTableau(stepperSDIRK->getTableau());

    this->setStepperType("Partitioned IMEX RK SSP2");
    this->setOrder(2);

  } else if (stepperType == "Partitioned IMEX RK ARS 233") {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 3;
    const bool isTVD = false;
    const Scalar sspcoef = -4.0;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = ST::one()/(2*ST::one());
    const Scalar gamma = (3*one+ST::squareroot(3*one))/(6*one);
    {
      // Explicit Tableau
      // Fill A:
      A(0,0) =        zero; A(0,1) =            zero; A(0,2) = zero;
      A(1,0) =       gamma; A(1,1) =            zero; A(1,2) = zero;
      A(2,0) = (gamma-1.0); A(2,1) = (2.0-2.0*gamma); A(2,2) = zero;

      // Fill b:
      b(0) = zero; b(1) = onehalf; b(2) = onehalf;

      // Fill c:
      c(0) = zero; c(1) = gamma; c(2) = one-gamma;

      int order = 2;

      auto expTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Explicit Tableau - Partitioned IMEX RK ARS 233",
        A,b,c,order,order,order,isTVD,sspcoef));

      this->setExplicitTableau(expTableau);
    }
    {
      // Implicit Tableau
      // Fill A:
      A(0,0) = zero; A(0,1) =            zero; A(0,2) =  zero;
      A(1,0) = zero; A(1,1) =           gamma; A(1,2) =  zero;
      A(2,0) = zero; A(2,1) = (1.0-2.0*gamma); A(2,2) = gamma;

      // Fill b:
      b(0) = zero; b(1) = onehalf; b(2) = onehalf;

      // Fill c:
      c(0) = zero; c(1) = gamma; c(2) = one-gamma;

      int order = 3;

      auto impTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Implicit Tableau - Partitioned IMEX RK ARS 233",
        A,b,c,order,order,order,isTVD,sspcoef));

      this->setImplicitTableau(impTableau);
    }
    this->setStepperType("Partitioned IMEX RK ARS 233");
    this->setOrder(3);

  } else if (stepperType == "General Partitioned IMEX RK") {
    this->setExplicitTableau(explicitTableau);
    this->setImplicitTableau(implicitTableau);
    this->setStepperType("General Partitioned IMEX RK");
    this->setOrder(1);

  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
       "Error - Not a valid StepperIMEX_RK_Partition type!  Stepper Type = "
       << stepperType <<  "\n"
       << "  Current valid types are: " << "\n"
       << "      'Partitioned IMEX RK 1st order'" << "\n"
       << "      'Partitioned IMEX RK SSP2'" << "\n"
       << "      'Partitioned IMEX RK ARS 233'" << "\n"
       << "      'General Partitioned IMEX RK'" << "\n");
  }

  TEUCHOS_TEST_FOR_EXCEPTION(explicitTableau_==Teuchos::null,
    std::runtime_error,
    "Error - StepperIMEX_RK_Partition - Explicit tableau is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(implicitTableau_==Teuchos::null,
    std::runtime_error,
    "Error - StepperIMEX_RK_Partition - Implicit tableau is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    explicitTableau_->numStages()!=implicitTableau_->numStages(),
    std::runtime_error,
       "Error - StepperIMEX_RK_Partition - Number of stages do not match!\n"
    << "  Explicit tableau = " << explicitTableau_->description() << "\n"
    << "    number of stages = " << explicitTableau_->numStages() << "\n"
    << "  Implicit tableau = " << implicitTableau_->description() << "\n"
    << "    number of stages = " << implicitTableau_->numStages() << "\n");

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setExplicitTableau(
  Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION(explicitTableau->isImplicit() == true,
    std::logic_error,
    "Error - Received an implicit Tableau for setExplicitTableau()!\n" <<
    "        Tableau = " << explicitTableau->description() << "\n");
  explicitTableau_ = explicitTableau;

  this->isInitialized_ = false;
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setImplicitTableau(
  Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION( implicitTableau->isDIRK() != true,
    std::logic_error,
    "Error - Did not receive a DIRK Tableau for setImplicitTableau()!\n" <<
    "        Tableau = " << implicitTableau->description() << "\n");
  implicitTableau_ = implicitTableau;

  this->isInitialized_ = false;
}

template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  RCP<Thyra::ModelEvaluator<Scalar> > ncModel =
    rcp_const_cast<Thyra::ModelEvaluator<Scalar> > (appModel);
  RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> > modelPairIMEX =
    rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(ncModel);
  TEUCHOS_TEST_FOR_EXCEPTION( modelPairIMEX == Teuchos::null, std::logic_error,
    "Error - StepperIMEX_RK::setModel() was given a ModelEvaluator that\n"
    "  could not be cast to a WrapperModelEvaluatorPairPartIMEX_Basic!\n"
    "  From: " << appModel << "\n"
    "  To  : " << modelPairIMEX << "\n"
    "  Likely have given the wrong ModelEvaluator to this Stepper.\n");

  setModelPair(modelPairIMEX);

  this->isInitialized_ = false;
}

/** \brief Create WrapperModelPairIMEX from user-supplied ModelEvaluator pair.
 *
 *  The user-supplied ME pair can contain any user-specific IMEX interactions
 *  between explicit and implicit MEs.
 */
template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setModelPair(
  const Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> > &
    mePairIMEX)
{
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
    wrapperModelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      (this->wrapperModel_);
  validExplicitODE    (mePairIMEX->getExplicitModel());
  validImplicitODE_DAE(mePairIMEX->getImplicitModel());
  wrapperModelPairIMEX = mePairIMEX;
  wrapperModelPairIMEX->initialize();

  this->wrapperModel_ = wrapperModelPairIMEX;

  this->isInitialized_ = false;
}

/** \brief Create WrapperModelPairIMEX from explicit/implicit ModelEvaluators.
 *
 *  Use the supplied explicit/implicit MEs to create a WrapperModelPairIMEX
 *  with basic IMEX interactions between explicit and implicit MEs.
 */
template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setModelPair(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel)
{
  validExplicitODE    (explicitModel);
  validImplicitODE_DAE(implicitModel);
  this->wrapperModel_ = Teuchos::rcp(
    new WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>(
                                             explicitModel, implicitModel));

  this->isInitialized_ = false;
}


#ifndef TEMPUS_HIDE_DEPRECATED_CODE
template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (this->stepperObserver_ == Teuchos::null)
     this->stepperObserver_  =
        Teuchos::rcp(new StepperRKObserverComposite<Scalar>());

  if (( obs == Teuchos::null ) and (this->stepperObserver_->getSize() >0 ) )
    return;

  if (( obs == Teuchos::null ) and (this->stepperObserver_->getSize() == 0) )
     obs = Teuchos::rcp(new StepperRKObserver<Scalar>());

    // Check that this casts to prevent a runtime error if it doesn't
  if (Teuchos::rcp_dynamic_cast<StepperRKObserver<Scalar> > (obs) != Teuchos::null) {
    this->stepperObserver_->addObserver(
         Teuchos::rcp_dynamic_cast<StepperRKObserver<Scalar> > (obs, true) );
  } else {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,0,"setObserver");
    *out << "Tempus::StepperIMEX_RK_Partition::setObserver: Warning: An observer has been provided that";
    *out << " does not support Tempus::StepperRKObserver. This observer WILL NOT be added.";
    *out << " In the future, this will result in a runtime error!" << std::endl;
  }

  this->isInitialized_ = false;
}
#endif


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::initialize()
{
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
    wrapperModelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      (this->wrapperModel_);
  TEUCHOS_TEST_FOR_EXCEPTION( wrapperModelPairIMEX == Teuchos::null,
    std::logic_error,
    "Error - Can not cast the wrapper Model Evaluator to a IMEX Model Pair."
    "StepperIMEX_RK_Partition::initialize()\n");

  // Initialize the stage vectors
  const int numStages = explicitTableau_->numStages();
  stageF_.resize(numStages);
  stageGx_.resize(numStages);
  for(int i=0; i < numStages; i++) {
    stageF_[i] = Thyra::createMember(wrapperModelPairIMEX->
                                     getExplicitModel()->get_f_space());
    stageGx_[i] = Thyra::createMember(wrapperModelPairIMEX->
                                     getImplicitModel()->get_f_space());
    assign(stageF_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
    assign(stageGx_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  xTilde_ = Thyra::createMember(wrapperModelPairIMEX->
                                getImplicitModel()->get_x_space());
  assign(xTilde_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

  StepperImplicit<Scalar>::initialize();
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperIMEX_RK::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();

  // Use x from inArgs as ICs, if needed.
  auto inArgs = this->wrapperModel_->getNominalValues();
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION( (x == Teuchos::null) &&
      (inArgs.get_x() == Teuchos::null), std::logic_error,
      "Error - setInitialConditions() needs the ICs from the SolutionHistory\n"
      "        or getNominalValues()!\n");

    x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
    initialState->setX(x);
  }


  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  TEUCHOS_TEST_FOR_EXCEPTION(icConsistency != "None", std::logic_error,
    "Error - setInitialConditions() requested a consistency of '"
             << icConsistency << "'.\n"
    "        But only 'None' is available for IMEX-RK!\n");

  TEUCHOS_TEST_FOR_EXCEPTION( this->getUseFSAL(), std::logic_error,
    "Error - The First-Step-As-Last (FSAL) principle is not "
         << "available for IMEX-RK.  Set useFSAL=false.\n");
}


template <typename Scalar>
void StepperIMEX_RK_Partition<Scalar>::evalImplicitModelExplicitly(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & Y,
  Scalar time, Scalar stepSize, Scalar stageNumber,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & G) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
    wrapperModelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      (this->wrapperModel_);
  MEB::InArgs<Scalar>  inArgs = wrapperModelPairIMEX->getInArgs();
  inArgs.set_x(X);
  inArgs.set_p(wrapperModelPairIMEX->getParameterIndex(), Y);
  if (inArgs.supports(MEB::IN_ARG_t))           inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size))   inArgs.set_step_size(stepSize);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(stageNumber);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);

  MEB::OutArgs<Scalar> outArgs = wrapperModelPairIMEX->getOutArgs();
  outArgs.set_f(G);

  wrapperModelPairIMEX->getImplicitModel()->evalModel(inArgs,outArgs);
  Thyra::Vt_S(G.ptr(), -1.0);
}


template <typename Scalar>
void StepperIMEX_RK_Partition<Scalar>::evalExplicitModel(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & Z,
  Scalar time, Scalar stepSize, Scalar stageNumber,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
    wrapperModelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      (this->wrapperModel_);
  MEB::InArgs<Scalar> inArgs =
    wrapperModelPairIMEX->getExplicitModel()->createInArgs();
  inArgs.set_x(Z);
  if (inArgs.supports(MEB::IN_ARG_t))           inArgs.set_t(time);
  if (inArgs.supports(MEB::IN_ARG_step_size))   inArgs.set_step_size(stepSize);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(stageNumber);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);

  MEB::OutArgs<Scalar> outArgs =
    wrapperModelPairIMEX->getExplicitModel()->createOutArgs();
  outArgs.set_f(F);

  wrapperModelPairIMEX->getExplicitModel()->evalModel(inArgs, outArgs);
  Thyra::Vt_S(F.ptr(), -1.0);
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::SerialDenseVector;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperIMEX_RK_Partition::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperIMEX_RK_Partition<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for IMEX_RK_Partition.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
#endif
    RCP<StepperIMEX_RK_Partition<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
      StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int numStages = explicitTableau_->numStages();
    const SerialDenseMatrix<int,Scalar> & AHat = explicitTableau_->A();
    const SerialDenseVector<int,Scalar> & bHat = explicitTableau_->b();
    const SerialDenseVector<int,Scalar> & cHat = explicitTableau_->c();
    const SerialDenseMatrix<int,Scalar> & A    = implicitTableau_->A();
    const SerialDenseVector<int,Scalar> & b    = implicitTableau_->b();
    const SerialDenseVector<int,Scalar> & c    = implicitTableau_->c();

    Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
     wrapperModelPairIMEX =
     Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
       (this->wrapperModel_);

    bool pass = true;
    Thyra::SolveStatus<Scalar> sStatus;
    stageZ_ = workingState->getX();
    Thyra::assign(stageZ_.ptr(), *(currentState->getX()));
    RCP<Thyra::VectorBase<Scalar> > stageY =
      wrapperModelPairIMEX->getExplicitOnlyVector(stageZ_);
    RCP<Thyra::VectorBase<Scalar> > stageX =
      wrapperModelPairIMEX->getIMEXVector(stageZ_);

    // Compute stage solutions
    for (int i = 0; i < numStages; ++i) {
      this->setStageNumber(i);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
      this->stepperObserver_->observeBeginStage(solutionHistory, *this);
#endif

      Thyra::assign(stageY.ptr(),
        *(wrapperModelPairIMEX->getExplicitOnlyVector(currentState->getX())));
      Thyra::assign(xTilde_.ptr(),
        *(wrapperModelPairIMEX->getIMEXVector(currentState->getX())));
      for (int j = 0; j < i; ++j) {
        if (AHat(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          RCP<Thyra::VectorBase<Scalar> > stageFy =
            wrapperModelPairIMEX->getExplicitOnlyVector(stageF_[j]);
          RCP<Thyra::VectorBase<Scalar> > stageFx =
            wrapperModelPairIMEX->getIMEXVector(stageF_[j]);
          Thyra::Vp_StV(stageY.ptr(),  -dt*AHat(i,j), *stageFy);
          Thyra::Vp_StV(xTilde_.ptr(), -dt*AHat(i,j), *stageFx);
        }
        if (A   (i,j) != Teuchos::ScalarTraits<Scalar>::zero())
          Thyra::Vp_StV(xTilde_.ptr(), -dt*A   (i,j), *(stageGx_[j]));
      }

      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEGIN_STAGE);

      Scalar ts    = time + c(i)*dt;
      Scalar tHats = time + cHat(i)*dt;
      if (A(i,i) == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Explicit stage for the ImplicitODE_DAE
        bool isNeeded = false;
        for (int k=i+1; k<numStages; ++k) if (A(k,i) != 0.0) isNeeded = true;
        if (b(i) != 0.0) isNeeded = true;
        if (isNeeded == false) {
          // stageGx_[i] is not needed.
          assign(stageGx_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
        } else {
          Thyra::assign(stageX.ptr(), *xTilde_);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
          this->stepperObserver_->observeBeforeImplicitExplicitly(solutionHistory, *this);
#endif
          evalImplicitModelExplicitly(stageX, stageY, ts, dt, i, stageGx_[i]);
        }
      } else {
        // Implicit stage for the ImplicitODE_DAE
        const Scalar alpha = Scalar(1.0)/(dt*A(i,i));
        const Scalar beta  = Scalar(1.0);

        // Setup TimeDerivative
        RCP<TimeDerivative<Scalar> > timeDer =
          Teuchos::rcp(new StepperIMEX_RKPartTimeDerivative<Scalar>(
            alpha, xTilde_.getConst()));

        // Setup InArgs and OutArgs
        typedef Thyra::ModelEvaluatorBase MEB;
        //MEB::InArgs<Scalar>  inArgs  = wrapperModelPairIMEX->getInArgs();
        //MEB::OutArgs<Scalar> outArgs = wrapperModelPairIMEX->getOutArgs();
        wrapperModelPairIMEX->setUseImplicitModel(true);
        MEB::InArgs<Scalar>  inArgs  = wrapperModelPairIMEX->createInArgs();
        MEB::OutArgs<Scalar> outArgs = wrapperModelPairIMEX->createOutArgs();
        inArgs.set_x(stageX);
        if (wrapperModelPairIMEX->getParameterIndex() >= 0)
          inArgs.set_p(wrapperModelPairIMEX->getParameterIndex(), stageY);
        if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(stageGx_[i]);
        if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (ts);
        if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
        if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (alpha);
        if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (beta);
        if (inArgs.supports(MEB::IN_ARG_stage_number))
          inArgs.set_stage_number(i);

        wrapperModelPairIMEX->setForSolve(timeDer, inArgs, outArgs);

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
        this->stepperObserver_->observeBeforeSolve(solutionHistory, *this);
#endif
        this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

        this->solver_->setModel(wrapperModelPairIMEX);
        sStatus = this->solveImplicitODE(stageX);
        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass = false;

        wrapperModelPairIMEX->setUseImplicitModel(false);

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
        this->stepperObserver_->observeAfterSolve(solutionHistory, *this);
#endif
        this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
          StepperRKAppAction<Scalar>::ACTION_LOCATION::AFTER_SOLVE);

        // Update contributions to stage values
        Thyra::V_StVpStV(stageGx_[i].ptr(), -alpha, *stageX, alpha, *xTilde_);
      }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
      this->stepperObserver_->observeBeforeExplicit(solutionHistory, *this);
#endif
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);
      evalExplicitModel(stageZ_, tHats, dt, i, stageF_[i]);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
      this->stepperObserver_->observeEndStage(solutionHistory, *this);
#endif
      this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
        StepperRKAppAction<Scalar>::ACTION_LOCATION::END_STAGE);
    }

    // Sum for solution: y_n = y_n-1 - dt*Sum{ bHat(i)*fy(i)            }
    // Sum for solution: x_n = x_n-1 - dt*Sum{ bHat(i)*fx(i) + b(i)*gx(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    RCP<Thyra::VectorBase<Scalar> > Z = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > X = wrapperModelPairIMEX->getIMEXVector(Z);
    for (int i=0; i < numStages; ++i) {
      if (bHat(i) != Teuchos::ScalarTraits<Scalar>::zero())
        Thyra::Vp_StV(Z.ptr(), -dt*bHat(i), *(stageF_[i]));
      if (b   (i) != Teuchos::ScalarTraits<Scalar>::zero())
        Thyra::Vp_StV(X.ptr(), -dt*b   (i), *(stageGx_[i]));
    }

    if (pass == true) workingState->setSolutionStatus(Status::PASSED);
    else              workingState->setSolutionStatus(Status::FAILED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
#endif
    this->stepperRKAppAction_->execute(solutionHistory, thisStepper,
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
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperIMEX_RK_Partition<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
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
  out << "  xTilde_            = " << xTilde_  << std::endl;
  out << "  stageZ_            = " << stageZ_  << std::endl;
  out << "  stageF_.size()     = " << stageF_.size() << std::endl;
  int numStages = stageF_.size();
  for (int i=0; i<numStages; ++i)
    out << "    stageF_["<<i<<"] = " << stageF_[i] << std::endl;
  out << "  stageGx_.size()    = " << stageGx_.size() << std::endl;
  numStages = stageGx_.size();
  for (int i=0; i<numStages; ++i)
    out << "    stageGx_["<<i<<"] = " << stageGx_[i] << std::endl;
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  out << "  stepperObserver_   = " << stepperObserver_ << std::endl;
#endif
  out << "  stepperRKAppAction_= " << this->stepperRKAppAction_ << std::endl;
  out << "  order_             = " << order_ << std::endl;
  out << "--------------------------------" << std::endl;
}


template<class Scalar>
bool StepperIMEX_RK_Partition<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;

  Teuchos::RCP<WrapperModelEvaluatorPairIMEX<Scalar> > wrapperModelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairIMEX<Scalar> >(
      this->wrapperModel_);

  if ( wrapperModelPairIMEX->getExplicitModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The explicit ModelEvaluator is not set!\n";
  }

  if ( wrapperModelPairIMEX->getImplicitModel() == Teuchos::null) {
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

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  if (stepperObserver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The stepper observer is not set!\n";
  }
#endif
  if (this->stepperRKAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The AppAction is not set!\n";
  }

  if ( explicitTableau_ == Teuchos::null ) {
    isValidSetup = false;
    out << "The explicit tableau is not set!\n";
  }

  if ( implicitTableau_ == Teuchos::null ) {
    isValidSetup = false;
    out << "The implicit tableau is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperIMEX_RK_Partition<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - Partitioned IMEX RK SSP2");
  pl->set<std::string>("Stepper Type", "Partitioned IMEX RK SSP2");
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<bool>("Initial Condition Consistency Check",
                this->getICConsistencyCheckDefault());
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>       ("Zero Initial Guess", false);
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperIMEX_RK_Partition_impl_hpp
