#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Tempus_RKButcherTableauBuilder_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(RKButcherTableauBuilder)

  // Nonmember constructor
  template Teuchos::RCP<RKButcherTableauBuilder<double> >
  rKButcherTableauBuilder();

  // Nonmember helper function
  template Teuchos::RCP<RKButcherTableau<double> >
  createRKBT(const std::string& rkbt_name,
             Teuchos::RCP<Teuchos::ParameterList> pl);

} // namespace Tempus

#endif
