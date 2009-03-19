#include "PB_EpetraLSCStablePreconditioner.hpp"

#include "Teuchos_RCP.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "PB_EpetraHelpers.hpp"
#include "NS/PB_NavierStokesPreconditioners.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::dyn_cast;
using Teuchos::rcp_dynamic_cast;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::DefaultDiagonalLinearOp;

namespace PB {
namespace Epetra {

EpetraLSCStablePreconditioner::EpetraLSCStablePreconditioner(const Epetra_Operator * A, const Epetra_Operator * invF, 
                                                             const Epetra_Operator * invBQBt, const Epetra_Vector * invMass)
   : EpetraInverseOpWrapper(dyn_cast<const EpetraOperatorWrapper>(*A).getMapStrategy())
{
   TEUCHOS_ASSERT(A!=0);
   TEUCHOS_ASSERT(invF!=0);
   TEUCHOS_ASSERT(invBQBt!=0);

   // Have to convert epetra inverses to something I can handle
   const RCP<const Thyra::LinearOpBase<double> > thyraInvF = Thyra::epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(invF)));
   const RCP<const Thyra::LinearOpBase<double> > thyraInvS = Thyra::epetraLinearOp(rcp(PB::Epetra::mechanicalInverse(invBQBt)));
   RCP<const Thyra::LinearOpBase<double> > thyraInvM; 
   if(invMass!=0) {
      const RCP<const VectorBase<double> > thyraVec  // need a Thyra::VectorBase object
         = Thyra::create_Vector(rcp(invMass,false),Thyra::create_VectorSpace(rcp(&invF->OperatorDomainMap(),false)));

      // create a diagonal matrix
      thyraInvM = Teuchos::rcp(new DefaultDiagonalLinearOp<double>(thyraVec));

   }
   else 
      thyraInvM = Thyra::identity<double>(thyraInvF->domain());

   // get blocked version of the matrix 
   const RCP<const EpetraOperatorWrapper> wrapA = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcp(A,false));
   TEST_FOR_EXCEPTION(wrapA==Teuchos::null,std::runtime_error,"EpetraLSCStablePreconditioner: Argument \"Epetra_Operator * A\""
                                                     " could not be cast to PB::Epetra::EpetraOperatorWrapper");
  
   // build preconditioned operator
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory
         = rcp(new PB::NS::LSCPreconditionerFactory(thyraInvF,thyraInvS,thyraInvM));
         // = rcp(new PB::NS::LSCStablePreconditionerFactory(thyraInvF,thyraInvS,thyraInvM));
   const RCP<const Thyra::PreconditionerBase<double> > prec = Thyra::prec(*precFactory,wrapA->getThyraOp());
   const RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

   // set the preconditioner as the operator
   SetOperator(precOp);

   // build the label for correct printing
   label_ = "Epetra LSCStable : invF = "+std::string(invF->Label())+", invS = " + std::string(invBQBt->Label());
} 

} // end namespace Epetra
} // end namespace PB
