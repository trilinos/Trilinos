#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SolutionHistory_impl.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class SolutionHistory<double>;

template RCP<SolutionHistory<double> >
SolutionHistory(RCP<ParameterList> pList_ = Teuchos::null);

} // namespace Tempus

#endif
