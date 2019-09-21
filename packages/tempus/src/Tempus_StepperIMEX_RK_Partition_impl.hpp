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

  this->setTableaus("Partitioned IMEX RK SSP2");
  this->setObserver();
}


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

  this->setExplicitTableau(explicitTableau);
  this->setImplicitTableau(implicitTableau);
  this->setObserver(obs);

  if (appModel != Teuchos::null) {

    this->setModel(appModel);
    this->setSolver(solver);
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

      auto explicitTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Explicit Tableau - Partitioned IMEX RK 1st order",
        A,b,c,order,order,order));

      this->setExplicitTableau(explicitTableau);
    }
    {
      // Implicit Tableau
      typedef Teuchos::ScalarTraits<Scalar> ST;
      int NumStages = 2;
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

      auto implicitTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Implicit Tableau - Partitioned IMEX RK 1st order",
        A,b,c,order,order,order));

      this->setImplicitTableau(implicitTableau);
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

      auto explicitTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Explicit Tableau - Partitioned IMEX RK ARS 233",
        A,b,c,order,order,order));

      this->setExplicitTableau(explicitTableau);
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

      auto implicitTableau = Teuchos::rcp(new RKButcherTableau<Scalar>(
        "Implicit Tableau - Partitioned IMEX RK ARS 233",
        A,b,c,order,order,order));

      this->setImplicitTableau(implicitTableau);
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
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (this->stepperObserver_ == Teuchos::null) {
      stepperIMEX_RKPartObserver_ =
        Teuchos::rcp(new StepperIMEX_RKPartObserver<Scalar>());
      this->stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >
          (stepperIMEX_RKPartObserver_);
     }
  } else {
    this->stepperObserver_ = obs;
    stepperIMEX_RKPartObserver_ =
      Teuchos::rcp_dynamic_cast<StepperIMEX_RKPartObserver<Scalar> >
        (this->stepperObserver_);
  }
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (explicitTableau_ == Teuchos::null) || (implicitTableau_ == Teuchos::null),
    std::logic_error,
    "Error - Need to set the Butcher Tableaus, setTableaus(), before calling "
    "StepperIMEX_RK_Partition::initialize()\n");

  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
    wrapperModelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      (this->wrapperModel_);
  TEUCHOS_TEST_FOR_EXCEPTION( wrapperModelPairIMEX == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
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

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
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
      if (!Teuchos::is_null(stepperIMEX_RKPartObserver_))
        stepperIMEX_RKPartObserver_->observeBeginStage(solutionHistory, *this);

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
          if (!Teuchos::is_null(stepperIMEX_RKPartObserver_))
            stepperIMEX_RKPartObserver_->
              observeBeforeImplicitExplicitly(solutionHistory, *this);
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

        if (!Teuchos::is_null(stepperIMEX_RKPartObserver_))
          stepperIMEX_RKPartObserver_->observeBeforeSolve(solutionHistory, *this);

        this->solver_->setModel(wrapperModelPairIMEX);
        sStatus = this->solveImplicitODE(stageX);
        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass = false;

        wrapperModelPairIMEX->setUseImplicitModel(false);

        if (!Teuchos::is_null(stepperIMEX_RKPartObserver_))
          stepperIMEX_RKPartObserver_->observeAfterSolve(solutionHistory, *this);

        // Update contributions to stage values
        Thyra::V_StVpStV(stageGx_[i].ptr(), -alpha, *stageX, alpha, *xTilde_);
      }

      if (!Teuchos::is_null(stepperIMEX_RKPartObserver_))
        stepperIMEX_RKPartObserver_->observeBeforeExplicit(solutionHistory,*this);
      evalExplicitModel(stageZ_, tHats, dt, i, stageF_[i]);
      if (!Teuchos::is_null(stepperIMEX_RKPartObserver_))
        stepperIMEX_RKPartObserver_->observeEndStage(solutionHistory, *this);
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
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
  }
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
   const Teuchos::EVerbosityLevel      /* verbLevel */) const
{
  out << this->getStepperType() << "::describe:" << std::endl
      << "wrapperModelPairIMEX = " << this->wrapperModel_->description()
      << std::endl;
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
