#include "Teko_MLLinearOp.hpp"

#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_mlutils.hpp"

#include "ml_MultiLevelPreconditioner.h"

namespace Teko {

MLLinearOp::MLLinearOp(const Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> & mlPrecOp)
   : mlPrecOp_(mlPrecOp)
{
   extractConversionInformation(*mlPrecOp_);
}

void MLLinearOp::extractConversionInformation(ML_Epetra::MultiLevelPreconditioner & mlPrec)
{
   const ML * ml = mlPrec.GetML();
   const ML_Smoother * preSmoother = ml->pre_smoother;
   const ML_Smoother * postSmoother = ml->post_smoother;

   // grab data object
   const mlutils::SmootherData * smootherData;
   if(preSmoother!=0)
      smootherData = (const mlutils::SmootherData *) ML_Get_MySmootherData(preSmoother);
   else if(postSmoother!=0)
      smootherData = (const mlutils::SmootherData *) ML_Get_MySmootherData(preSmoother);
   else
      TEST_FOR_EXCEPTION(true,std::runtime_error,
                              "MLLinearOp::extractConversionInformation pre and post smoother " <<
                              "are both null, cannot build operator");

   Amat_ = Teuchos::rcp_dynamic_cast<Epetra::EpetraOperatorWrapper>(smootherData->Amat);

   // for doing Epetra_Vector -> Thyra vector conversion
   mappingStrategy_ = Amat_->getMapStrategy();

   // to construct vectorspace for this object
   BlockedLinearOp bloA = toBlockedLinearOp(Amat_->getThyraOp());
   productDomain_ = bloA->productRange(); // this operator is the inverse of bloA...hence swap range and domain
   productRange_ = bloA->productDomain();

}

void MLLinearOp::implicitApply(const BlockedMultiVector & x, BlockedMultiVector & y,
                               const double alpha, const double beta) const
{
   int columns = x->domain()->dim();
   TEUCHOS_ASSERT(columns==y->domain()->dim()); // check for equal columns

   // (re)allocate vectors conditinolly if required size changed
   if(eX_==Teuchos::null || columns!=eX_->NumVectors()) {
      eX_ = Teuchos::rcp(new Epetra_MultiVector(mlPrecOp_->OperatorDomainMap(),x->domain()->dim()));
      eY_ = Teuchos::rcp(new Epetra_MultiVector(mlPrecOp_->OperatorRangeMap(),y->domain()->dim()));
   }

   BlockedMultiVector yCopy;
   if(beta!=0)
      yCopy = deepcopy(y);
   else
      yCopy = y;

   // initialize Epetra vectors
   eY_->PutScalar(0.0);
   mappingStrategy_->copyThyraIntoEpetra(x, *eX_,*Amat_);

   // run multigrid
   mlPrecOp_->ApplyInverse(*eX_,*eY_);

   // fill thyra vectors
   mappingStrategy_->copyEpetraIntoThyra(*eY_, yCopy.ptr(),*Amat_);

   // scale result by alpha
   if(beta!=0)
      update(alpha,yCopy,beta,y); // y = alpha * yCopy + beta * y 
   else if(alpha!=1.0)
      scale(alpha,y); // y = alpha * y
}

void MLLinearOp::describe(Teuchos::FancyOStream & out_arg,
                          const Teuchos::EVerbosityLevel verbLevel) const
{
   out_arg << "MLLinearOp";
}

}
