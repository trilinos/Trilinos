#include "Epetra/PB_StridedEpetraOperator.hpp"
#include "Epetra/PB_StridedMappingStrategy.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "EpetraExt_MultiVectorOut.h"

namespace PB {
namespace Epetra {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

StridedEpetraOperator::StridedEpetraOperator(int numVars,const Teuchos::RCP<Epetra_Operator> & content) 
      : PB::Epetra::EpetraOperatorWrapper()
{
   std::vector<int> vars;
   
   // build vector describing the sub maps
   for(int i=0;i<numVars;i++) vars.push_back(1);

   SetContent(vars,content);
}

StridedEpetraOperator::StridedEpetraOperator(const std::vector<int> & vars,const Teuchos::RCP<Epetra_Operator> & content) 
      : PB::Epetra::EpetraOperatorWrapper()
{
   SetContent(vars,content);
}

void StridedEpetraOperator::SetContent(const std::vector<int> & vars,const Teuchos::RCP<Epetra_Operator> & content)
{ 
   fullContent_ = content;
   mapStrategy_ = rcp(new StridedMappingStrategy(vars,Teuchos::rcpFromRef(fullContent_->OperatorDomainMap()),
                                                         fullContent_->Comm()));

   // build thyra operator
   BuildBlockedOperator(); 
}

int StridedEpetraOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
   if (!useTranspose_)
   {
       // allocate space for each vector
       const RCP<Thyra::MultiVectorBase<double> > tX = Thyra::createMembers(thyraOp_->domain(),X.NumVectors()); 
       const RCP<Thyra::MultiVectorBase<double> > tY = Thyra::createMembers(thyraOp_->range(),X.NumVectors());
       Thyra::assign(tX.ptr(),0.0);
       Thyra::assign(tY.ptr(),0.0);

       // copy epetra X into thyra X
       mapStrategy_->copyEpetraIntoThyra(X, tX.ptr(),*this);

       // perform matrix vector multiplication
       thyraOp_->apply(Thyra::NONCONJ_ELE,*tX,&*tY);

       // copy thyra Y into epetra Y
       mapStrategy_->copyThyraIntoEpetra(tY, Y,*this);
   }
   else
   {
       TEUCHOS_ASSERT(false);
   }
 
   return 0;
}

void StridedEpetraOperator::BuildBlockedOperator()
{
   TEUCHOS_ASSERT(mapStrategy_!=Teuchos::null);

   // get a CRS matrix
   const RCP<const Epetra_CrsMatrix> crsContent = rcp_dynamic_cast<const Epetra_CrsMatrix>(fullContent_);

   // ask the strategy to build the Thyra operator for you
   const RCP<const Thyra::LinearOpBase<double> > A 
         = rcp_dynamic_cast<const StridedMappingStrategy>(mapStrategy_)->buildBlockedThyraOp(crsContent);

   // set whatever is returned
   SetOperator(A,false);
}

const Teuchos::RCP<const Epetra_Operator> StridedEpetraOperator::GetBlock(int i,int j) const
{
   const RCP<const Thyra::BlockedLinearOpBase<double> > blkOp 
         = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(getThyraOp());

   return Thyra::get_Epetra_Operator(*blkOp->getBlock(i,j));
}

} // end namespace Epetra
} // end namespace PB
