// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperIMEX_RK_impl_hpp
#define Tempus_StepperIMEX_RK_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


// StepperIMEX_RK definitions:
template<class Scalar>
StepperIMEX_RK<Scalar>::StepperIMEX_RK(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  std::string stepperType,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  this->setTableaus(pList, stepperType);
  this->setParameterList(pList);
  this->setModel(appModel);
  this->setSolver();
  this->setObserver();
  this->initialize();
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::setTableaus(
  Teuchos::RCP<Teuchos::ParameterList> pList,
  std::string stepperType)
{
  if (stepperType == "") {
    if (pList == Teuchos::null)
      stepperType = "IMEX RK SSP2";
    else
      stepperType = pList->get<std::string>("Stepper Type");
  }

  if (stepperType == "IMEX RK 1st order") {
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

  } else if (stepperType == "IMEX RK SSP2") {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    // Explicit Tableau
    this->setExplicitTableau("RK Explicit Trapezoidal", Teuchos::null);

    // Implicit Tableau
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set<std::string>("Stepper Type", "SDIRK 2 Stage 3rd order");
    pl->set("Solver Name", "");
    Scalar gamma = 1.0 - 1.0/ST::squareroot(2.0);
    pl->set<double>("gamma",gamma);
    this->setImplicitTableau("SDIRK 2 Stage 3rd order", pl);

    description_ = stepperType;
    order_ = 2;
  } else if (stepperType == "IMEX RK ARS 233") {
    using std::to_string;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar gammaN = (3.0+ST::squareroot(3.0))/(6.0);
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

  } else if (stepperType == "General IMEX RK") {
    Teuchos::RCP<Teuchos::ParameterList> explicitPL = Teuchos::rcp(
      new Teuchos::ParameterList(pList->sublist("IMEX-RK Explicit Stepper")));

    Teuchos::RCP<Teuchos::ParameterList> implicitPL = Teuchos::rcp(
      new Teuchos::ParameterList(pList->sublist("IMEX-RK Implicit Stepper")));

    // TODO: should probably check the order of the tableau match
    this->setExplicitTableau("General ERK",  explicitPL);
    this->setImplicitTableau("General DIRK", implicitPL);
    description_ = stepperType;
    order_ = 0;  // TODO: Determine overall order

  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
       "Error - Not a valid StepperIMEX_RK type!  Stepper Type = "
       << stepperType <<  "\n"
       << "  Current valid types are: " << "\n"
       << "      'IMEX RK 1st order'" << "\n"
       << "      'IMEX RK SSP2'" << "\n"
       << "      'IMEX RK ARS 233'" << "\n"
       << "      'General IMEX RK'" << "\n");
  }

  TEUCHOS_TEST_FOR_EXCEPTION(explicitTableau_==Teuchos::null,
    std::runtime_error,
    "Error - StepperIMEX_RK - Explicit tableau is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(implicitTableau_==Teuchos::null,
    std::runtime_error,
    "Error - StepperIMEX_RK - Implicit tableau is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    explicitTableau_->numStages()!=implicitTableau_->numStages(),
    std::runtime_error,
       "Error - StepperIMEX_RK - Number of stages do not match!\n"
    << "  Explicit tableau = " << explicitTableau_->description() << "\n"
    << "    number of stages = " << explicitTableau_->numStages() << "\n"
    << "  Implicit tableau = " << implicitTableau_->description() << "\n"
    << "    number of stages = " << implicitTableau_->numStages() << "\n");
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::setExplicitTableau(
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
void StepperIMEX_RK<Scalar>::setExplicitTableau(
  Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION(explicitTableau->isImplicit() == true,
    std::logic_error,
       "Error - Received an implicit Tableau for setExplicitTableau()!\n"
    << "  explicitTableau = " << explicitTableau->description() << "\n");
  explicitTableau_ = explicitTableau;
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::setImplicitTableau(
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
void StepperIMEX_RK<Scalar>::setImplicitTableau(
  Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau)
{
  TEUCHOS_TEST_FOR_EXCEPTION( implicitTableau_->isDIRK() != true,
    std::logic_error,
       "Error - Did not receive a DIRK Tableau for setImplicitTableau()!\n"
    << "  implicitTableau = " << implicitTableau->description() << "\n");
  implicitTableau_ = implicitTableau;
}

template<class Scalar>
void StepperIMEX_RK<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  RCP<Thyra::ModelEvaluator<Scalar> > ncModel =
    rcp_const_cast<Thyra::ModelEvaluator<Scalar> > (appModel);
  RCP<WrapperModelEvaluatorPairIMEX_Basic<Scalar> > modelPairIMEX =
    rcp_dynamic_cast<WrapperModelEvaluatorPairIMEX_Basic<Scalar> > (ncModel);
  TEUCHOS_TEST_FOR_EXCEPTION( modelPairIMEX == Teuchos::null, std::logic_error,
    "Error - StepperIMEX_RK::setModel() was given a ModelEvaluator that\n"
    "  could not be cast to a WrapperModelEvaluatorPairIMEX_Basic!\n"
    "  From: " << appModel << "\n"
    "  To  : " << modelPairIMEX << "\n"
    "  Likely have given the wrong ModelEvaluator to this Stepper.\n");

  setModelPair(modelPairIMEX);
}

template<class Scalar>
void StepperIMEX_RK<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

/** \brief Create WrapperModelPairIMEX from user-supplied ModelEvaluator pair.
 *
 *  The user-supplied ME pair can contain any user-specific IMEX interactions
 *  between explicit and implicit MEs.
 */
template<class Scalar>
void StepperIMEX_RK<Scalar>::setModelPair(
  const Teuchos::RCP<WrapperModelEvaluatorPairIMEX_Basic<Scalar> > &
    modelPairIMEX)
{
  this->validExplicitODE    (modelPairIMEX->getExplicitModel());
  this->validImplicitODE_DAE(modelPairIMEX->getImplicitModel());
  wrapperModelPairIMEX_ = modelPairIMEX;
  wrapperModelPairIMEX_->initialize();
  int expXDim = wrapperModelPairIMEX_->getExplicitModel()->get_x_space()->dim();
  int impXDim = wrapperModelPairIMEX_->getImplicitModel()->get_x_space()->dim();
  TEUCHOS_TEST_FOR_EXCEPTION( expXDim != impXDim, std::logic_error,
    "Error - \n"
    "  Explicit and Implicit x vectors are incompatible!\n"
    "  Explicit vector dim = " << expXDim << "\n"
    "  Implicit vector dim = " << impXDim << "\n");

}

/** \brief Create WrapperModelPairIMEX from explicit/implicit ModelEvaluators.
 *
 *  Use the supplied explicit/implicit MEs to create a WrapperModelPairIMEX
 *  with basic IMEX interactions between explicit and implicit MEs.
 */
template<class Scalar>
void StepperIMEX_RK<Scalar>::setModelPair(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel)
{
  this->validExplicitODE    (explicitModel);
  this->validImplicitODE_DAE(implicitModel);
  wrapperModelPairIMEX_ = Teuchos::rcp(
    new WrapperModelEvaluatorPairIMEX_Basic<Scalar>(
                                              explicitModel, implicitModel));
}


/** \brief Set the solver to a pre-defined solver in the ParameterList.
 *  The solver is set to solverName sublist in the Stepper's ParameterList.
 *  The solverName sublist should already be defined in the Stepper's
 *  ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperIMEX_RK<Scalar>::setSolver(std::string solverName)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = Teuchos::sublist(stepperPL_, solverName, true);
  stepperPL_->set("Solver Name", solverName);
  solver_ = rcp(new Thyra::NOXNonlinearSolver());
  RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
  solver_->setParameterList(noxPL);
}


/** \brief Set the solver to the supplied Parameter sublist.
 *  This adds a new solver Parameter sublist to the Stepper's ParameterList.
 *  If the solver sublist is null, the solver is set to the solver name
 *  in the Stepper's ParameterList.
 */
template<class Scalar>
void StepperIMEX_RK<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string solverName = stepperPL_->get<std::string>("Solver Name");
  if (is_null(solverPL)) {
    // Create default solver, otherwise keep current solver.
    if (solver_ == Teuchos::null) {
      solverPL = Teuchos::sublist(stepperPL_, solverName, true);
      solver_ = rcp(new Thyra::NOXNonlinearSolver());
      RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
      solver_->setParameterList(noxPL);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( solverName == solverPL->name(),
      std::logic_error,
         "Error - Trying to add a solver that is already in ParameterList!\n"
      << "  Stepper Type = "<< stepperPL_->get<std::string>("Stepper Type")
      << "\n" << "  Solver Name  = "<<solverName<<"\n");
    solverName = solverPL->name();
    stepperPL_->set("Solver Name", solverName);
    stepperPL_->set(solverName, solverPL);      // Add sublist
    solver_ = rcp(new Thyra::NOXNonlinearSolver());
    RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
    solver_->setParameterList(noxPL);
  }
}


/** \brief Set the solver.
 *  This sets the solver to supplied solver and adds solver's ParameterList
 *  to the Stepper ParameterList.
 */
template<class Scalar>
void StepperIMEX_RK<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = solver->getNonconstParameterList();
  std::string solverName = solverPL->name();
  stepperPL_->set("Solver Name", solverName);
  stepperPL_->set(solverName, solverPL);      // Add sublist
  solver_ = solver;
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::setObserver(
  Teuchos::RCP<StepperIMEX_RKObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperIMEX_RKObserver_ == Teuchos::null) {
      stepperIMEX_RKObserver_ =
        Teuchos::rcp(new StepperIMEX_RKObserver<Scalar>());
    }
  } else {
    stepperIMEX_RKObserver_ = obs;
  }
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::initialize()
{
  // Initialize the stage vectors
  const int numStages = explicitTableau_->numStages();
  stageF_.resize(numStages);
  stageG_.resize(numStages);
  for(int i=0; i < numStages; i++) {
    stageF_[i] = Thyra::createMember(wrapperModelPairIMEX_->get_f_space());
    stageG_[i] = Thyra::createMember(wrapperModelPairIMEX_->get_f_space());
    assign(stageF_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
    assign(stageG_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  }

  xTilde_ = Thyra::createMember(wrapperModelPairIMEX_->get_x_space());
  assign(xTilde_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
}


template <typename Scalar>
void StepperIMEX_RK<Scalar>::evalImplicitModelExplicitly(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
  Scalar time, Scalar stepSize, Scalar stageNumber,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & G) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>  inArgs  = wrapperModelPairIMEX_->getInArgs();
  inArgs.set_x(X);
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

  MEB::OutArgs<Scalar> outArgs = wrapperModelPairIMEX_->getOutArgs();
  outArgs.set_f(G);

  wrapperModelPairIMEX_->getImplicitModel()->evalModel(inArgs,outArgs);
  Thyra::Vt_S(G.ptr(), -1.0);
}


template <typename Scalar>
void StepperIMEX_RK<Scalar>::evalExplicitModel(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
  Scalar time, Scalar stepSize, Scalar stageNumber,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inArgs =
    wrapperModelPairIMEX_->getExplicitModel()->createInArgs();
  inArgs.set_x(X);
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
    wrapperModelPairIMEX_->getExplicitModel()->createOutArgs();
  outArgs.set_f(F);

  wrapperModelPairIMEX_->getExplicitModel()->evalModel(inArgs, outArgs);
  Thyra::Vt_S(F.ptr(), -1.0);
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::SerialDenseVector;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperIMEX_RK::takeStep()");
  {
    stepperIMEX_RKObserver_->observeBeginTakeStep(solutionHistory, *this);
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

    bool pass = true;
    Thyra::SolveStatus<double> sStatus;
    stageX_ = workingState->getX();
    Thyra::assign(stageX_.ptr(), *(currentState->getX()));

    // Compute stage solutions
    for (int i = 0; i < numStages; ++i) {
      stepperIMEX_RKObserver_->observeBeginStage(solutionHistory, *this);
      Thyra::assign(xTilde_.ptr(), *(currentState->getX()));
      for (int j = 0; j < i; ++j) {
        if (AHat(i,j) != Teuchos::ScalarTraits<Scalar>::zero())
          Thyra::Vp_StV(xTilde_.ptr(), -dt*AHat(i,j), *(stageF_[j]));
        if (A   (i,j) != Teuchos::ScalarTraits<Scalar>::zero())
          Thyra::Vp_StV(xTilde_.ptr(), -dt*A   (i,j), *(stageG_[j]));
      }

      Scalar ts    = time + c(i)*dt;
      Scalar tHats = time + cHat(i)*dt;
      if (A(i,i) == Teuchos::ScalarTraits<Scalar>::zero()) {
        // Explicit stage for the ImplicitODE_DAE
        bool isNeeded = false;
        for (int k=i+1; k<numStages; ++k) if (A(k,i) != 0.0) isNeeded = true;
        if (b(i) != 0.0) isNeeded = true;
        if (isNeeded == false) {
          // stageG_[i] is not needed.
          assign(stageG_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
        } else {
          Thyra::assign(stageX_.ptr(), *xTilde_);
          stepperIMEX_RKObserver_->
            observeBeforeImplicitExplicitly(solutionHistory, *this);
          evalImplicitModelExplicitly(stageX_, ts, dt, i, stageG_[i]);
        }
      } else {
        // Implicit stage for the ImplicitODE_DAE
        Scalar alpha = 1.0/(dt*A(i,i));

        // Setup TimeDerivative
        Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
          Teuchos::rcp(new StepperIMEX_RKTimeDerivative<Scalar>(
            alpha, xTilde_.getConst()));

        // Setup InArgs and OutArgs
        typedef Thyra::ModelEvaluatorBase MEB;
        MEB::InArgs<Scalar>  inArgs  = wrapperModelPairIMEX_->getInArgs();
        MEB::OutArgs<Scalar> outArgs = wrapperModelPairIMEX_->getOutArgs();
        inArgs.set_x(stageX_);
        if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(stageG_[i]);
        if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (ts);
        if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
        if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (alpha);
        if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (1.0);
        if (inArgs.supports(MEB::IN_ARG_stage_number))
          inArgs.set_stage_number(i);

        wrapperModelPairIMEX_->setForSolve(timeDer, inArgs, outArgs);

        stepperIMEX_RKObserver_->observeBeforeSolve(solutionHistory, *this);

        sStatus = this->solveNonLinear(wrapperModelPairIMEX_,*solver_,stageX_);
        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass = false;

        stepperIMEX_RKObserver_->observeAfterSolve(solutionHistory, *this);

        // Update contributions to stage values
        Thyra::V_StVpStV(stageG_[i].ptr(), -alpha, *stageX_, alpha, *xTilde_);
      }

      stepperIMEX_RKObserver_->observeBeforeExplicit(solutionHistory, *this);
      evalExplicitModel(stageX_, tHats, dt, i, stageF_[i]);
      stepperIMEX_RKObserver_->observeEndStage(solutionHistory, *this);
    }

    // Sum for solution: x_n = x_n-1 - dt*Sum{ bHat(i)*f(i) + b(i)*g(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0; i < numStages; ++i) {
      if (bHat(i) != Teuchos::ScalarTraits<Scalar>::zero())
        Thyra::Vp_StV((workingState->getX()).ptr(), -dt*bHat(i), *(stageF_[i]));
      if (b   (i) != Teuchos::ScalarTraits<Scalar>::zero())
        Thyra::Vp_StV((workingState->getX()).ptr(), -dt*b   (i), *(stageG_[i]));
    }

    if (pass == true)
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;
    workingState->setOrder(this->getOrder());
    stepperIMEX_RKObserver_->observeEndTakeStep(solutionHistory, *this);
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
StepperIMEX_RK<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperIMEX_RK<Scalar>::description() const
{
  return(description_);
}


template<class Scalar>
void StepperIMEX_RK<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModelPairIMEX_ = " << wrapperModelPairIMEX_->description()
      << std::endl;
}


template <class Scalar>
void StepperIMEX_RK<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    stepperPL_ = this->getDefaultParameters();
  } else {
    stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperIMEX_RK<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - IMEX RK SSP2");
  pl->set("Stepper Type", "IMEX RK SSP2");
  pl->set("Solver Name", "",
    "Name of ParameterList containing the solver specifications.");

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperIMEX_RK<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - IMEX RK SSP2");
  pl->set<std::string>("Stepper Type", "IMEX RK SSP2");
  pl->set<std::string>("Solver Name", "Default Solver");
  Teuchos::RCP<Teuchos::ParameterList> solverPL=this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperIMEX_RK<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperIMEX_RK<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperIMEX_RK_impl_hpp
