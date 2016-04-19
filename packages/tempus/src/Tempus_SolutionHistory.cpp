#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SolutionHistory_impl.hpp"

namespace Tempus {

template class SolutionHistory<double>;

template RCP<SolutionHistory<double> >
SolutionHistory(RCP<ParameterList> pList_ = Teuchos::null);

} // namespace Tempus
