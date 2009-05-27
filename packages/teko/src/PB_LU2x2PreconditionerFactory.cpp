#include "PB_LU2x2PreconditionerFactory.hpp"

// PB includes
#include "PB_LU2x2InverseOp.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace PB {

// construct a PreconditionerFactory
LU2x2PreconditionerFactory::LU2x2PreconditionerFactory(LinearOp & invA00, LinearOp & invS)
      : invOpsStrategy_(rcp(new StaticLU2x2Strategy(invA00,invA00,invS)))
{ }

/** @brief Build a simple static LU2x2 preconditioner */
LU2x2PreconditionerFactory::LU2x2PreconditionerFactory(LinearOp & hatInvA00,LinearOp & tildeInvA00,LinearOp & invS)
      : invOpsStrategy_(rcp(new StaticLU2x2Strategy(hatInvA00,tildeInvA00,invS)))
{ }

LU2x2PreconditionerFactory::LU2x2PreconditionerFactory(const RCP<const LU2x2Strategy> & strategy)
   : invOpsStrategy_(strategy)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// initialize a newly created preconditioner object
LinearOp LU2x2PreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LinearOp hatInvA00   = invOpsStrategy_->getHatInvA00(A,state);
   LinearOp tildeInvA00 = invOpsStrategy_->getTildeInvA00(A,state);
   LinearOp invS        = invOpsStrategy_->getInvS(A,state);

   // build the SchurSolve LinearOp
   return createLU2x2InverseOp(A,hatInvA00,tildeInvA00,invS);
}

} // end namespace PB
