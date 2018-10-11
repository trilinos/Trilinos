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
  } else {
    if (solverName == solverPL->name()) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperImplicit::setSolver()");
      *out << "Warning - Replacing the solver ParameterList.\n"
           << "  Stepper Type = "<< stepperPL_->get<std::string>("Stepper Type")
           << "\n  Solver PL = " << solverName << std::endl;
      stepperPL_->remove(solverName);
    }
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
const Thyra::SolveStatus<Scalar>
StepperImplicit<Scalar>::solveImplicitODE(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x)
{
  if (getZeroInitialGuess())
    Thyra::assign(x.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

  const Thyra::SolveStatus<Scalar> sStatus = (*solver_).solve(&*x);

  return sStatus;
}


} // namespace Tempus
#endif // Tempus_StepperImplicit_impl_hpp
