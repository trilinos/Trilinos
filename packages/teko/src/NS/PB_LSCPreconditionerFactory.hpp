#ifndef __PB_LSCPreconditionerFactory_hpp__
#define __PB_LSCPreconditionerFactory_hpp__

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_LSCStrategy.hpp"

namespace PB {
namespace NS { // Navier-Stokes specialization

class LSCPreconditionerFactory 
   : public BlockPreconditionerFactory {
   public:
      // constructors for a LSCPreconditionerFactory

      // Staiblized constructor
      LSCPreconditionerFactory(const LinearOp & invF,const LinearOp & invBQBtmC,
                               const LinearOp & invD,const LinearOp & invMass);
 
      // Stable constructor
      LSCPreconditionerFactory(const LinearOp & invF,
                               const LinearOp & invBQBtmC,
                               const LinearOp & invMass);

      // fully generic constructor
      LSCPreconditionerFactory(const Teuchos::RCP<const LSCStrategy> & strategy);

      // for PreconditionerFactoryBase
      ///////////////////////////////////////////////////////////////////////

      virtual LinearOp buildPreconditionerOperator(BlockedLinearOp & blo) const;

   protected:
      // main driver for code
      Teuchos::RCP<const LSCStrategy> invOpsStrategy_;
};

} // end namespace NS
} // end namespace PB

#endif
