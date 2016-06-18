#include "Tempus_config.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_StepperBackwardEuler_impl.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class StepperBackwardEuler<double>;

template RCP<StepperBackwardEuler<double> >
StepperBackwardEuler(
  RCP<ParameterList> parameterList,
  RCP<Thyra::ModelEvaluator<double> > model,
  RCP<Thyra::NonlinearSolverBase<double> > solver = Teuchos::null);

template RCP<StepperBackwardEuler<double> >
StepperBackwardEuler(
  RCP<ParameterList>& parameterList,
  Array<RCP<Thyra::ModelEvaluator<double> > > models,
  Array<RCP<Thyra::NonlinearSolverBase<double> > >  solvers = Teuchos::null);

} // namespace Tempus

#endif
