#include "PB_EpetraBlockPreconditioner.hpp"

namespace PB {
namespace Epetra {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

/** \brief Constructor that takes the BlockPreconditionerFactory that will
  *        build the preconditioner.
  *
  * Constructor that takes the BlockPreconditionerFactory that will
  * build the preconditioner.
  */
EpetraBlockPreconditioner::EpetraBlockPreconditioner(const Teuchos::RCP<const BlockPreconditionerFactory> & bfp)
   : preconFactory_(bfp)
{ }

/** \brief Build this preconditioner from an Epetra_Operator 
  * passed in to this object. It is assume that this Epetra_Operator
  *
  * Build this preconditioner from an Epetra_Operator 
  * passed in to this object. It is assume that this Epetra_Operator
  * will be a EpetraOperatorWrapper object, so the block Thyra components
  * can be easily extracted.
  *
  * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
  */
void EpetraBlockPreconditioner::buildPreconditioner(const Epetra_Operator & A)
{
   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   RCP<const LinearOpBase<double> > thyraA = eow->getThyraOp(); 

   // set the mapping strategy
   SetMapStrategy(rcp(new InverseMappingStrategy(eow->getMapStrategy())));
   
   // actually build the preconditioner
   RCP<const LinearOpBase<double> > preconditioner = preconFactory_->buildPreconditionerOperator(thyraA);

   SetOperator(preconditioner,false);

   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
}

} // end namespace Epetra
} // end namespace PB
