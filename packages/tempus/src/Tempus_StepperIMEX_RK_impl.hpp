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
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
       "Error - 'General IMEX RK' is not implemented yet!\n");
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
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > nc_model =
    Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar> > (appModel);
  Teuchos::RCP<WrapperModelEvaluatorPairIMEX<Scalar> > modelPairIMEX =
    Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorPairIMEX<Scalar> > (nc_model);

  setModelPair(modelPairIMEX);
}

template<class Scalar>
void StepperIMEX_RK<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

/** \brief Create WrapperModelPairIMEX from user-supplied ModelEvaluator pair
 *  The user-supplied ME pair can contain any user-specific IMEX interactions
 *  between explicit and implicit MEs.
 */
template<class Scalar>
void StepperIMEX_RK<Scalar>::setModelPair(
  const Teuchos::RCP<WrapperModelEvaluatorPairIMEX<Scalar> > & modelPairIMEX)
{
  this->validExplicitODE    (modelPairIMEX->getExplicitModel());
  this->validImplicitODE_DAE(modelPairIMEX->getImplicitModel());
  wrapperModelPairIMEX_ = modelPairIMEX;

  inArgs_  = wrapperModelPairIMEX_->getImplicitModel()->getNominalValues();
  outArgs_ = wrapperModelPairIMEX_->getImplicitModel()->createOutArgs();
}

/** \brief Create WrapperModelPairIMEX from explicit/implicit ModelEvaluators.
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
void StepperIMEX_RK<Scalar>::initialize()
{
  // Initialize the stage vectors
  const int numStages = explicitTableau_->numStages();
  stageX_ = wrapperModelPairIMEX_->getNominalValues().get_x()->clone_v();
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
std::function<void (const Thyra::VectorBase<Scalar> &,
                          Thyra::VectorBase<Scalar> &)>
StepperIMEX_RK<Scalar>::xDotFunction(
  Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
{
  return [=](const Thyra::VectorBase<Scalar> & x,
                   Thyra::VectorBase<Scalar> & xDotTilde)
    {
      // ith stage
      // s = 1/(dt*a_ii)
      // xOld = solution at beginning of time step
      // xTilde = xOld + dt*(Sum_{j=1}^{i-1} a_ij x_dot_j)
      // xDotTilde = - (s*x_i - s*xTilde)
      Thyra::V_StVpStV(Teuchos::ptrFromRef(xDotTilde),s,x,-s,*xTilde);
    };
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

    // Compute stage solutions
    bool pass = true;
    Thyra::SolveStatus<double> sStatus;
    Thyra::assign(stageX_.ptr(), *(currentState->getX()));
    for (int i = 0; i < numStages; ++i) {
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
          typedef Thyra::ModelEvaluatorBase MEB;
          Thyra::assign(stageX_.ptr(), *xTilde_);
          inArgs_.set_x(stageX_);
          if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(ts);
          if (inArgs_.supports(MEB::IN_ARG_x_dot))
            inArgs_.set_x_dot(Teuchos::null);
          outArgs_.set_f(stageG_[i]);

          wrapperModelPairIMEX_->getImplicitModel()->evalModel(inArgs_,
                                                                outArgs_);
          Thyra::Vt_S(stageG_[i].ptr(), -1.0);
        }
      } else {
        // Implicit stage for the ImplicitODE_DAE
        Scalar alpha = 1.0/dt/A(i,i);
        Scalar beta  = 1.0;

        // function used to compute time derivative
        auto computeXDot = xDotFunction(alpha,xTilde_.getConst());

        wrapperModelPairIMEX_->initialize2(computeXDot, ts, tHats,
                                            alpha, beta, dt, i);

        sStatus = this->solveNonLinear(wrapperModelPairIMEX_,*solver_,stageX_);
        if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED) pass = false;

        // Update contributions to stage values
        Thyra::V_StVpStV(stageG_[i].ptr(), -alpha, *stageX_, alpha, *xTilde_);
      }

      // Evaluate the ExplicitODE
      wrapperModelPairIMEX_->evalExplicitModel(stageX_, tHats, stageF_[i]);
      Thyra::Vt_S(stageF_[i].ptr(), -1.0);
    }

    // Sum for solution: x_n = x_n-1 - dt*Sum{ bHat(i)*f(i) + b(i)*g(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0 ; i < numStages ; ++i) {
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
