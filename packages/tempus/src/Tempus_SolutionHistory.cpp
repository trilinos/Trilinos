#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SolutionHistory_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(SolutionHistory)

  // Non-member ctor
  template Teuchos::RCP<SolutionHistory<double> >
  solutionHistory(Teuchos::RCP<Teuchos::ParameterList> pList_);

} // namespace Tempus

#endif
