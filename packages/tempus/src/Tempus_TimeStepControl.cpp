#include "Tempus_config.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControl_impl.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class TimeStepControl<double>;

template RCP<TimeStepControl<double> > TimeStepControl();

template RCP<TimeStepControl<double> >
TimeStepControl(RCP<ParameterList> pList_ = Teuchos::null);

} // namespace Tempus

#endif
