// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperIMEX_RK_Partition_impl_hpp
#define Tempus_StepperIMEX_RK_Partition_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
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
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType)
{
  this->setTableaus(Teuchos::null, stepperType);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  this->setTableaus(pList, "Partitioned IMEX RK SSP2");
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
StepperIMEX_RK_Partition<Scalar>::StepperIMEX_RK_Partition(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  this->setTableaus(pList, stepperType);
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setTableaus(
  Teuchos::RCP<Teuchos::ParameterList> pList,
  std::string stepperType)
{
  if (stepperType == "") {
    if (pList == Teuchos::null)
      stepperType = "Partitioned IMEX RK SSP2";
    else
      stepperType = pList->get<std::string>("Stepper Type",
                                            "Partitioned IMEX RK SSP2");
  }

  if (stepperType == "Partitioned IMEX RK 1st order") {
    {
      // Explicit Tableau
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
      pl->setName("IMEX-RK Explicit Stepper");
      pl->set<std::string>("Stepper Type", "General ERK");

      // Tableau ParameterList
      Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
      tableauPL->set<std::string>("A", "0.0 0.0; 1.0 0.0");
      tableauPL->set<std::string>("b", "1.0 0.0");
      tableauPL->set<std::string>("c", "0.0 1.0");
      tableauPL->set<int>("order", 1);
      pl->set("Tableau", *tableauPL);

      this->setExplicitTableau("General ERK", pl);
    }
    {
      // Implicit Tableau
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
      pl->setName("IMEX-RK Implicit Stepper");
      pl->set<std::string>("Stepper Type", "General DIRK");
      pl->set("Solver Name", "");

      // Tableau ParameterList
      Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
      tableauPL->set<std::string>("A", "0.0 0.0; 0.0 1.0");
      tableauPL->set<std::string>("b", "0.0 1.0");
      tableauPL->set<std::string>("c", "0.0 1.0");
      tableauPL->set<int>("order", 1);
      pl->set("Tableau", *tableauPL);

      this->setImplicitTableau("General DIRK", pl);
    }
    description_ = stepperType;
    order_ = 1;

  } else if (stepperType == "Partitioned IMEX RK SSP2") {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    // Explicit Tableau
    this->setExplicitTableau("RK Explicit Trapezoidal", Teuchos::null);

    // Implicit Tableau
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set<std::string>("Stepper Type", "SDIRK 2 Stage 3rd order");
    pl->set("Solver Name", "");
    const Scalar one = ST::one();
    Scalar gamma = one - one/ST::squareroot(2*one);
    pl->set<double>("gamma",gamma);
    this->setImplicitTableau("SDIRK 2 Stage 3rd order", pl);

    description_ = stepperType;
    order_ = 2;
  } else if (stepperType == "Partitioned IMEX RK ARS 233") {
    using std::to_string;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar gammaN = (3*one+ST::squareroot(3*one))/(6*one);
    std::string gamma      = to_string(        gammaN);
    std::string one_gamma  = to_string(1.0-    gammaN);
    std::string one_2gamma = to_string(1.0-2.0*gammaN);
    std::string two_2gamma = to_string(2.0-2.0*gammaN);
    std::string gamma_one  = to_string(        gammaN-1.0);
    {
      // Explicit Tableau
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
      pl->setName("IMEX-RK Explicit Stepper");
      pl->set<std::string>("Stepper Type", "General ERK");

      // Tableau ParameterList
      Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
      tableauPL->set<std::string>("A",
        "0.0 0.0 0.0; "+gamma+" 0.0 0.0; "+gamma_one+" "+two_2gamma+" 0.0");
      tableauPL->set<std::string>("b", "0.0 0.5 0.5");
      tableauPL->set<std::string>("c", "0.0 "+gamma+" "+one_gamma);
      tableauPL->set<int>("order", 2);
      pl->set("Tableau", *tableauPL);

      this->setExplicitTableau("General ERK", pl);
    }
    {
      // Implicit Tableau
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
      pl->setName("IMEX-RK Implicit Stepper");
      pl->set<std::string>("Stepper Type", "General DIRK");
      pl->set("Solver Name", "");

      // Tableau ParameterList
      Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
      tableauPL->set<std::string>("A",
        "0.0 0.0 0.0; 0.0 "+gamma+" 0.0; 0.0 "+one_2gamma+" "+gamma);
      tableauPL->set<std::string>("b", "0.0 0.5 0.5");
      tableauPL->set<std::string>("c", "0.0 "+gamma+" "+one_gamma);
      tableauPL->set<int>("order", 3);
      pl->set("Tableau", *tableauPL);

      this->setImplicitTableau("General DIRK", pl);
    }
    description_ = stepperType;
    order_ = 3;

  } else if (stepperType == "General Partitioned IMEX RK") {
    Teuchos::RCP<Teuchos::ParameterList> explicitPL = Teuchos::rcp(
      new Teuchos::ParameterList(pList->sublist("IMEX-RK Explicit Stepper")));

    Teuchos::RCP<Teuchos::ParameterList> implicitPL = Teuchos::rcp(
      new Teuchos::ParameterList(pList->sublist("IMEX-RK Implicit Stepper")));

    // TODO: should probably check the order of the tableau match
    this->setExplicitTableau("General ERK",  explicitPL);
    this->setImplicitTableau("General DIRK", implicitPL);
    description_ = stepperType;
    order_ = pList->get<int>("overall order", 0);

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
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  explicitTableau_ = createRKBT<Scalar>(stepperType,pList);
  TEUCHOS_TEST_FOR_EXCEPTION(explicitTableau_->isImplicit() == true,
    std::logic_error,
       "Error - Received an implicit Tableau for setExplicitTableau()!\n"
    << "  Stepper Type = " << stepperType << "\n");
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setExplicitTableau(
  Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION(explicitTableau->isImplicit() == true,
    std::logic_error,
       "Error - Received an implicit Tableau for setExplicitTableau()!\n"
    << "  explicitTableau = " << explicitTableau->description() << "\n");
  explicitTableau_ = explicitTableau;
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setImplicitTableau(
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  implicitTableau_ = createRKBT<Scalar>(stepperType,pList);
  //Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  //Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  //implicitTableau_->describe(*out,verbLevel);
  TEUCHOS_TEST_FOR_EXCEPTION( implicitTableau_->isDIRK() != true,
    std::logic_error,
       "Error - Did not receive a DIRK Tableau for setImplicitTableau()!\n"
    << "  Stepper Type = " << stepperType << "\n");
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setImplicitTableau(
  Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION( implicitTableau_->isDIRK() != true,
    std::logic_error,
       "Error - Did not receive a DIRK Tableau for setImplicitTableau()!\n"
    << "  implicitTableau = " << implicitTableau->description() << "\n");
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
  this->validExplicitODE    (mePairIMEX->getExplicitModel());
  this->validImplicitODE_DAE(mePairIMEX->getImplicitModel());
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
  this->validExplicitODE    (explicitModel);
  this->validImplicitODE_DAE(implicitModel);
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
    if (stepperObserver_ == Teuchos::null) {
      stepperIMEX_RKPartObserver_ =
        Teuchos::rcp(new StepperIMEX_RKPartObserver<Scalar>());
      stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >
          (stepperIMEX_RKPartObserver_);
     }
  } else {
    stepperObserver_ = obs;
    stepperIMEX_RKPartObserver_ =
      Teuchos::rcp_dynamic_cast<StepperIMEX_RKPartObserver<Scalar> >
        (stepperObserver_);
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

  this->setTableaus(this->stepperPL_);
  this->setParameterList(this->stepperPL_);
  this->setSolver();
  this->setObserver();

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

    stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
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
        Scalar alpha = 1.0/(dt*A(i,i));

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
        if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (1.0);
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
    stepperObserver_->observeEndTakeStep(solutionHistory, *this);
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
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperIMEX_RK_Partition<Scalar>::description() const
{
  return(description_);
}


template<class Scalar>
void StepperIMEX_RK_Partition<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModelPairIMEX = " << this->wrapperModel_->description()
      << std::endl;
}


template <class Scalar>
void StepperIMEX_RK_Partition<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null) this->stepperPL_ = this->getDefaultParameters();
  } else {
    this->stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperIMEX_RK_Partition<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - IMEX RK SSP2");
  pl->set("Stepper Type", "IMEX RK SSP2");
  pl->set("Zero Initial Guess", false);
  pl->set("Solver Name", "",
    "Name of ParameterList containing the solver specifications.");

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperIMEX_RK_Partition<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - IMEX RK SSP2");
  pl->set<std::string>("Stepper Type", "IMEX RK SSP2");
  pl->set<bool>       ("Zero Initial Guess", false);
  pl->set<std::string>("Solver Name", "Default Solver");
  Teuchos::RCP<Teuchos::ParameterList> solverPL=this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperIMEX_RK_Partition<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperIMEX_RK_Partition<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperIMEX_RK_Partition_impl_hpp
