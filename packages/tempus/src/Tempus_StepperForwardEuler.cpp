#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperForwardEuler_impl.hpp"

namespace tempus {

template class StepperForwardEuler<double>;

template RCP<StepperForwardEuler<double> >
StepperForwardEuler(
  RCP<ParameterList> parameterList,
  RCP<Thyra::ModelEvaluator<double> > model,
  RCP<Thyra::NonlinearSolverBase<double> > solver = Teuchos::null);

template RCP<StepperForwardEuler<double> >
StepperForwardEuler(
  RCP<ParameterList>& parameterList,
  Array<RCP<Thyra::ModelEvaluator<double> > > models,
  Array<RCP<Thyra::NonlinearSolverBase<double> > >  solvers =Teuchos::null);

} // namespace tempus
