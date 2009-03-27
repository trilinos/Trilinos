#include "PB_EpetraBlockJacobiPreconditioner.hpp"

#include "Teuchos_RCP.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "PB_EpetraHelpers.hpp"
#include "PB_JacobiPreconditionerFactory.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::dyn_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::null;

namespace PB {
namespace Epetra {

EpetraBlockJacobiPreconditioner::EpetraBlockJacobiPreconditioner(const Epetra_Operator * A, const Epetra_Operator * invD1, 
                                                             const Epetra_Operator * invD2)
   : EpetraInverseOpWrapper(rcp(new InverseMappingStrategy(dyn_cast<const EpetraOperatorWrapper>(*A).getMapStrategy())))
{
   TEUCHOS_ASSERT(A!=0);
   TEUCHOS_ASSERT(invD1!=0);
   TEUCHOS_ASSERT(invD2!=0);

   // Have to convert epetra inverses to something I can handle
   LinearOp thyraInvD1 = Thyra::epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(invD1)));
   LinearOp thyraInvD2 = Thyra::epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(invD2)));

   // get blocked version of the matrix 
   const RCP<const EpetraOperatorWrapper> wrapA = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcp(A,false));
   TEST_FOR_EXCEPTION(wrapA==null,std::runtime_error,"EpetraBlockJacobiPreconditioner: Argument \"Epetra_Operator * A\""
                                                     " could not be cast to PB::Epetra::EpetraOperatorWrapper");
  
   // build preconditioned operator
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory
         = rcp(new PB::JacobiPreconditionerFactory(thyraInvD1,thyraInvD2));
   const RCP<const Thyra::PreconditionerBase<double> > prec = Thyra::prec(*precFactory,wrapA->getThyraOp());
   const RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

   // set the preconditioner as the operator
   SetOperator(precOp);

   // build the label for correct printing
   label_ = "Epetra Block Jacobi : op0 = { "
          + std::string(invD1->Label()) + " }, op1 = { " 
          + std::string(invD2->Label()) +" }";
} 

} // end namespace Epetra
} // end namespace PB
