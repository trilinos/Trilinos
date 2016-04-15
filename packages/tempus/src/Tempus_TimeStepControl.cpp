#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControl_impl.hpp"

namespace tempus {

template class TimeStepControl<double>;

template RCP<TimeStepControl<double> > TimeStepControl();

template RCP<TimeStepControl<double> >
TimeStepControl(RCP<ParameterList> pList_ = Teuchos::null);

} // namespace tempus
