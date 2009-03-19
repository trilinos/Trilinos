#include "PB_LSCStrategy.hpp"

namespace PB {
namespace NS {

   // Staiblized constructor
StaticLSCStrategy::StaticLSCStrategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invBQBtmC,
                                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invD,
                                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invMass)
   : invF_(invF), invBQBtmC_(invBQBtmC), invD_(invD), invMass_(invMass)
{ }
 
   // Stable constructor
StaticLSCStrategy::StaticLSCStrategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invBQBtmC,
                                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invMass)
   : invF_(invF), invBQBtmC_(invBQBtmC), invD_(Teuchos::null), invMass_(invMass)
{ }

} // end namespace NS
} // end namespace PB
