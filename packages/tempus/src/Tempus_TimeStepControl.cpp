#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControl_impl.hpp"

namespace tempus {

template class TimeStepControl< double >;

template Teuchos::RCP< TimeStepControl< double > > TimeStepControl();

template Teuchos::RCP< TimeStepControl< double > >
TimeStepControl( RCP<Teuchos::ParameterList> paramList_ = Teuchos::null );

} // namespace tempus
