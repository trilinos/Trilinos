#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SolutionHistory_impl.hpp"

namespace tempus {

template class SolutionHistory< double >;

template Teuchos::RCP< SolutionHistory< double > >
SolutionHistory( RCP<Teuchos::ParameterList> paramList_ = Teuchos::null );

} // namespace tempus
