#include "Tempus_config.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionStateMetaData_impl.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class SolutionStateMetaData<double>;

template RCP<SolutionStateMetaData<double> > SolutionStateMetaData();

template RCP<SolutionStateMetaData<double> >
SolutionStateMetaData(const double time,
                      const double dt,
                      const int    iStep,
                      const double errorAbs,
                      const double errorRel,
                      const int    order,
                      const int    nFailures,
                      const int    nConsecutiveFailures,
                      const SolutionStatus status,
                      const bool   output,
                      const bool   isAccepted,
                      const bool   isRestartable,
                      const bool   isInterpolated,
                      const double accuracy);

} // namespace Tempus

#endif
