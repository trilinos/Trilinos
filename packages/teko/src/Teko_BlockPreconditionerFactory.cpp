#include "Teko_BlockPreconditionerFactory.hpp"

#include "Teko_Preconditioner.hpp"
#include "Teko_InverseLibrary.hpp"

#include "Thyra_DefaultPreconditioner.hpp"

using namespace Thyra;

namespace Teko {

/////////////////////////////////////////////////////

LinearOp BlockPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const
{
   // get the blocked linear operator
   RCP<LinearOpBase<double> > loA = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(lo);
   BlockedLinearOp A = Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(loA);

   state.setInitialized(false);

   return buildPreconditionerOperator(A,dynamic_cast<BlockPreconditionerState &>(state));
}

//! is this operator compatiable with the preconditioner factory?
bool BlockPreconditionerFactory::isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const
{
   RCP<const Thyra::PhysicallyBlockedLinearOpBase<double> > A 
         = Teuchos::rcp_dynamic_cast<const Thyra::PhysicallyBlockedLinearOpBase<double> >(fwdOpSrc.getOp());
   return A!=Teuchos::null;
}

} // end namespace Teko
