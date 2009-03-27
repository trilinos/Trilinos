#include "PB_LU2x2PreconditionerFactory.hpp"

// PB includes
#include "PB_SchurSolveLinearOp.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace PB {

// construct a PreconditionerFactory
LU2x2PreconditionerFactory::LU2x2PreconditionerFactory(LinearOp & invA00, LinearOp & invS)
      : invOpsStrategy_(rcp(new StaticLU2x2Strategy(invA00,invS)))
{ }

LU2x2PreconditionerFactory::LU2x2PreconditionerFactory(const RCP<const LU2x2Strategy> & strategy)
   : invOpsStrategy_(strategy)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// initialize a newly created preconditioner object
LinearOp LU2x2PreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & A) const
{
   LinearOp invA00 = invOpsStrategy_->getInvA00(A);
   LinearOp invS   = invOpsStrategy_->getInvS(A);

   // build the SchurSolve LinearOp
   return createNewSchurSolveLinearOp(A,invA00,invS);
}

} // end namespace PB
