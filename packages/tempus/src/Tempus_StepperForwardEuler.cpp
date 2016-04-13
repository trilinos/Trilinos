#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperForwardEuler_impl.hpp"

namespace tempus {

template class StepperForwardEuler< double >;

template Teuchos::RCP< StepperForwardEuler< double > >
StepperForwardEuler(
  RCP<Teuchos::ParameterList> parameterList,
  RCP<Thyra::ModelEvaluator<double> > model,
  RCP<Thyra::NonlinearSolverBase<double> > solver = Teuchos::null);

template Teuchos::RCP< StepperForwardEuler< double > >
StepperForwardEuler(
  RCP<Teuchos::ParameterList>& parameterList,
  Array<RCP<Thyra::ModelEvaluator<double> > > models,
  Array<RCP<Thyra::NonlinearSolverBase<double> > >  solvers =Teuchos::null);

} // namespace tempus
