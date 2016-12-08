#ifndef Tempus_StepperDIRK_impl_hpp
#define Tempus_StepperDIRK_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


// StepperDIRK definitions:
template<class Scalar>
StepperDIRK<Scalar>::StepperDIRK(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->setParameterList(pList);
  this->setModel(transientModel);
  this->setSolver();
  this->setTableau();
  this->initialize();
}

template<class Scalar>
void StepperDIRK<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->validImplicitODE_DAE(transientModel);
  residualModel_ =
    Teuchos::rcp(new ResidualModelEvaluator<Scalar>(transientModel));

  inArgs_  = residualModel_->getNominalValues();
  outArgs_ = residualModel_->createOutArgs();
}

template<class Scalar>
void StepperDIRK<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel)
{
  this->setModel(transientModel);
}

template<class Scalar>
void StepperDIRK<Scalar>::setSolver(
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &solver)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (solver != Teuchos::null) {
    solver_ = solver;
  } else {
    // Construct solver from ParameterList
    std::string solverName = pList_->get<std::string>("Solver Name");
    RCP<ParameterList> solverPL = Teuchos::sublist(pList_, solverName, true);
    RCP<ParameterList> noxPL    = Teuchos::sublist(solverPL, "NOX", true);
    solver_ = rcp(new Thyra::NOXNonlinearSolver());
    solver_->setParameterList(noxPL);
  }
}

template<class Scalar>
void StepperDIRK<Scalar>::setTableau(std::string stepperType)
{
  if (stepperType == "")
    stepperType = pList_->get<std::string>("Stepper Type");

  DIRK_ButcherTableau_ = createRKBT<Scalar>(stepperType,pList_);

  //Teuchos::SerialDenseMatrix<int,Scalar> A = DIRK_ButcherTableau_->A();
  //std::cout << " A = \n" << A << std::endl;

  TEUCHOS_TEST_FOR_EXCEPTION( DIRK_ButcherTableau_->isDIRK() != true,
    std::logic_error,
       "Error - StepperDIRK did not receive a DIRK Butcher Tableau!\n"
    << "  Stepper Type = "<< pList_->get<std::string>("Stepper Type") << "\n");
  description_ = DIRK_ButcherTableau_->description();

  // Initialize the stage vectors
  const int numStages = DIRK_ButcherTableau_->numStages();
  stageX_    = Thyra::createMember (residualModel_->get_x_space());
  stageXDot_ = Thyra::createMembers(residualModel_->get_x_space(), numStages);
  stageXPartial_ = Thyra::createMember(residualModel_->get_x_space());
}

template<class Scalar>
void StepperDIRK<Scalar>::initialize()
{
}

template <typename Scalar>
std::function<void (const Thyra::VectorBase<Scalar> &,
                          Thyra::VectorBase<Scalar> &)>
StepperDIRK<Scalar>::xDotFunction(
  Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > stageXPartial)
{
  return [=](const Thyra::VectorBase<Scalar> & x,
                   Thyra::VectorBase<Scalar> & x_dot)
    {
      // ith stage
      // s = 1/(dt*a_ii)
      // xOld = solution at beginning of time step
      // stageXPartial = xOld + dt*(Sum_{j=1}^{i-1} a_ij x_dot_j)
      // x_dot_i = s*x_i - s*stageXPartial
      Thyra::V_StVpStV(Teuchos::ptrFromRef(x_dot),s,x,-s,*stageXPartial);
    };
}

template<class Scalar>
void StepperDIRK<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperDIRK::takeStep()");
  {
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    const Scalar time = currentState->getTime();

    const int stages = DIRK_ButcherTableau_->numStages();
    Teuchos::SerialDenseMatrix<int,Scalar> A = DIRK_ButcherTableau_->A();
    Teuchos::SerialDenseVector<int,Scalar> b = DIRK_ButcherTableau_->b();
    Teuchos::SerialDenseVector<int,Scalar> c = DIRK_ButcherTableau_->c();

    // Compute stage solutions
    bool pass = true;
    Thyra::SolveStatus<double> sStatus;
    for (int i=0 ; i < stages ; ++i) {
      Thyra::assign(stageXPartial_.ptr(), *(currentState->getX()));
      for (int j=0 ; j < i ; ++j) {
        if (A(i,j) != Teuchos::ScalarTraits<Scalar>::zero()) {
          Thyra::Vp_StV(stageXPartial_.ptr(), dt*A(i,j), *(stageXDot_->col(j)));
        }
      }

      Scalar alpha = 1.0/dt/A(i,i);
      Scalar beta = 1.0;
      Scalar ts = time + c(i)*dt;

      // function used to compute time derivative
      auto computeXDot = xDotFunction(alpha,stageXPartial_.getConst());

      residualModel_->initialize(computeXDot, ts, alpha, beta);

      sStatus = this->solveNonLinear(residualModel_,*solver_,stageX_,inArgs_);
      if (sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED ) pass=false;

      computeXDot(*stageX_, *(stageXDot_->col(i)));
    }

    // Sum for solution: x_n = x_n-1 + Sum{ dt*b(i) * f(i) }
    Thyra::assign((workingState->getX()).ptr(), *(currentState->getX()));
    for (int i=0 ; i < stages ; ++i) {
      if (b(i) != Teuchos::ScalarTraits<Scalar>::zero()) {
        Thyra::Vp_StV((workingState->getX()).ptr(), dt*b(i),
                      *(stageXDot_->col(i)));
      }
    }

    if (DIRK_ButcherTableau_->isEmbedded() ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Implicit RK embedded methods not implemented yet!.\n");
    }

    if (pass == true)
      workingState->stepperState_->stepperStatus_ = Status::PASSED;
    else
      workingState->stepperState_->stepperStatus_ = Status::FAILED;
    workingState->setOrder(DIRK_ButcherTableau_->order());
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
StepperDIRK<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperDIRK<Scalar>::description() const
{
  return(description_);
}


template<class Scalar>
void StepperDIRK<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "residualModel_ = " << residualModel_->description() << std::endl;
}


template <class Scalar>
void StepperDIRK<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  //pList->validateParameters(*this->getValidParameters());
  pList_ = pList;

  Teuchos::readVerboseObjectSublist(&*pList_,this);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperDIRK<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperDIRK<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperDIRK<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperDIRK_impl_hpp
