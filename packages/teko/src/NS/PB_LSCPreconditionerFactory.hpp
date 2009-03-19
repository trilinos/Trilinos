#ifndef __PB_LSCPreconditionerFactory_hpp__
#define __PB_LSCPreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "PB_Block2x2PreconditionerFactory.hpp"
#include "PB_LSCStrategy.hpp"

namespace PB {
namespace NS { // Navier-Stokes specialization

class LSCPreconditionerFactory 
   : public Block2x2PreconditionerFactory {
   public:
      // constructors for a LSCPreconditionerFactory

      // Staiblized constructor
      LSCPreconditionerFactory(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invBQBtmC,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invD,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invMass);
 
      // Stable constructor
      LSCPreconditionerFactory(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invBQBtmC,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invMass);

      // fully generic constructor
      LSCPreconditionerFactory(const Teuchos::RCP<const LSCStrategy> & strategy);

      // for PreconditionerFactoryBase
      ///////////////////////////////////////////////////////////////////////

      // initialize a newly created preconditioner object
      void initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > & fwdOpSrc,
                          Thyra::PreconditionerBase<double> * precOp,
                          const Thyra::ESupportSolveUse supportSolveUse=Thyra::SUPPORT_SOLVE_UNSPECIFIED) const;
   protected:
      // main driver for code
      Teuchos::RCP<const LSCStrategy> invOpsStrategy_;
};

} // end namespace NS
} // end namespace PB

#endif
