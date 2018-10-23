// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperImplicit_impl_hpp
#define Tempus_StepperImplicit_impl_hpp

// Tempus
//#include "Tempus_Stepper.hpp"
//#include "Tempus_TimeDerivative.hpp"

// Thrya
//#include "Thyra_VectorBase.hpp"
//#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {


template<class Scalar>
void StepperImplicit<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validImplicitODE_DAE(appModel);
  wrapperModel_ =
    Teuchos::rcp(new WrapperModelEvaluatorBasic<Scalar>(appModel));
}


template<class Scalar>
void StepperImplicit<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}


/** \brief Set the solver to a pre-defined solver in the ParameterList.
 *
 *  The solver is set to solverName sublist in the Stepper's ParameterList.
 *  The solverName sublist should already be defined in the Stepper's
 *  ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperImplicit<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::ParameterList> solverPL =
    Teuchos::sublist(stepperPL_, solverName, true);
  this->setSolver(solverPL);
}


/** \brief Set the solver to the supplied Parameter sublist.
 *
 *  This adds a new solver Parameter sublist to the Stepper's ParameterList.
 *  If the solver sublist is null, the solver is set to the solver name
 *  in the Stepper's ParameterList.
 */
template<class Scalar>
void StepperImplicit<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string solverName = stepperPL_->get<std::string>("Solver Name");
  if (is_null(solverPL)) {
    if ( stepperPL_->isSublist(solverName) )
      solverPL = Teuchos::sublist(stepperPL_, solverName, true);
    else
      solverPL = this->defaultSolverParameters();
  }

  solverName = solverPL->name();
  stepperPL_->set("Solver Name", solverName);
  stepperPL_->set(solverName, *solverPL);      // Add sublist
  RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);

  solver_ = rcp(new Thyra::NOXNonlinearSolver());
  solver_->setParameterList(noxPL);

  TEUCHOS_TEST_FOR_EXCEPTION(wrapperModel_ == Teuchos::null, std::logic_error,
       "Error - StepperImplicit<Scalar>::setSolver() wrapperModel_ is unset!\n"
    << "  Should call setModel(...) first.\n");
  solver_->setModel(wrapperModel_);
}


/** \brief Set the solver.
 *
 *  This sets the solver to supplied solver and adds solver's ParameterList
 *  to the Stepper ParameterList.
 */
template<class Scalar>
void StepperImplicit<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  Teuchos::RCP<Teuchos::ParameterList> solverPL =
    solver->getNonconstParameterList();
  this->setSolver(solverPL);
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperImplicit<Scalar>::
getStepperXDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDot() != Teuchos::null) stepperXDot_ = state->getXDot();
  // Else use temporary storage stepperXDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDot_ == Teuchos::null, std::logic_error,
    "Error - stepperXDot_ has not been set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDot_;
}


template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperImplicit<Scalar>::
getStepperXDotDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDotDot() != Teuchos::null) stepperXDotDot_=state->getXDotDot();
  // Else use temporary storage stepperXDotDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDotDot_ == Teuchos::null,std::logic_error,
    "Error - stepperXDotDot_ has not been set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDotDot_;
}


template<class Scalar>
const Thyra::SolveStatus<Scalar>
StepperImplicit<Scalar>::solveImplicitODE(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x)
{
  if (getZeroInitialGuess())
    Thyra::assign(x.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

  const Thyra::SolveStatus<Scalar> sStatus = (*solver_).solve(&*x);

  return sStatus;
}


template<class Scalar>
const Thyra::SolveStatus<Scalar>
StepperImplicit<Scalar>::solveImplicitODE(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
  const Scalar time,
  const Teuchos::RCP<ImplicitODEParameters<Scalar> > & p )
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>  inArgs  = wrapperModel_->getInArgs();
  MEB::OutArgs<Scalar> outArgs = wrapperModel_->getOutArgs();
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
  if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
  if (inArgs.supports(MEB::IN_ARG_step_size))
    inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (p->alpha_);
  if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (p->beta_);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(p->stageNumber_);

  wrapperModel_->setForSolve(p->timeDer_, inArgs, outArgs, p->evaluationType_);

  Thyra::SolveStatus<Scalar> sStatus;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  switch (p->evaluationType_)
  {
    case SOLVE_FOR_X: {
      if (getZeroInitialGuess()) Thyra::assign(x.ptr(), ST::zero());
      sStatus = (*solver_).solve(&*x);
      break;
    }
    case SOLVE_FOR_XDOT_CONST_X: {
      //if (getZeroInitialGuess()) Thyra::assign(xDot.ptr(), ST::zero());
      sStatus = (*solver_).solve(&*xDot);
      break;
    }
    default: {
      TEUCHOS_TEST_FOR_EXCEPT("Invalid EVALUATION_TYPE!");
    }
  }

  return sStatus;
}


template<class Scalar>
void
StepperImplicit<Scalar>::evaluateImplicitODE(
        Teuchos::RCP<Thyra::VectorBase<Scalar> > & f,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & xDot,
  const Scalar time,
  const Teuchos::RCP<ImplicitODEParameters<Scalar> > & p )
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>  inArgs  = wrapperModel_->getInArgs();
  inArgs.set_x(x);
  if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
  if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
  if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(p->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (Scalar(0.0));
  if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (Scalar(0.0));

  MEB::OutArgs<Scalar> outArgs = wrapperModel_->getOutArgs();
  outArgs.set_f(f);

  wrapperModel_->setForSolve(Teuchos::null,inArgs,outArgs,p->evaluationType_);

  wrapperModel_->evalModel(inArgs, outArgs);
}


} // namespace Tempus
#endif // Tempus_StepperImplicit_impl_hpp
