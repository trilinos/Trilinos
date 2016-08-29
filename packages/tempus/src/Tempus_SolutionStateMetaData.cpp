#include "Tempus_config.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionStateMetaData_impl.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class SolutionStateMetaData<double>;

template Teuchos::RCP<SolutionStateMetaData<double> > SolutionStateMetaData();

template Teuchos::RCP<SolutionStateMetaData<double> >
SolutionStateMetaData(const double time,
                      const int    iStep,
                      const double dt,
                      const double errorAbs,
                      const double errorRel,
                      const int    order,
                      const int    nFailures,
                      const int    nConsecutiveFailures,
                      const Status solutionStatus,
                      const bool   output,
                      const bool   isRestartable,
                      const bool   isInterpolated,
                      const double accuracy);

} // namespace Tempus

#endif
