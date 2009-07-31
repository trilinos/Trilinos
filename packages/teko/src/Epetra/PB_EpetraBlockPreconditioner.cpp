#include "PB_EpetraBlockPreconditioner.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"

// Teuchos includes
#include "Teuchos_Time.hpp"

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
   RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 

   // set the mapping strategy
   SetMapStrategy(rcp(new InverseMappingStrategy(eow->getMapStrategy())));
   
   // actually build the preconditioner
   preconObj_ = preconFactory_->createPrec();
   RCP<const Thyra::LinearOpSourceBase<double> > lOpSrc = Thyra::defaultLinearOpSource(thyraA);
   preconFactory_->initializePrec(lOpSrc,&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);

   // extract preconditioner operator
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();

   SetOperator(preconditioner,false);

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
}

/** \brief Build this preconditioner from an Epetra_Operator 
  * passed in to this object. It is assume that this Epetra_Operator
  *
  * Build this preconditioner from an Epetra_Operator 
  * passed in to this object. It is assume that this Epetra_Operator
  * will be a EpetraOperatorWrapper object, so the block Thyra components
  * can be easily extracted.
  *
  * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
  * \param[in] src A vector that was used to build the source operator.
  */
void EpetraBlockPreconditioner::buildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & epetra_mv)
{
   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 

   // set the mapping strategy
   SetMapStrategy(rcp(new InverseMappingStrategy(eow->getMapStrategy())));

   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
   
   // build the thyra version of the source multivector
   RCP<Thyra::MultiVectorBase<double> > thyra_mv = Thyra::createMembers(thyraA->range(),epetra_mv.NumVectors());
   getMapStrategy()->copyEpetraIntoThyra(epetra_mv,thyra_mv.ptr(),*eow);

   // actually build the preconditioner
   preconObj_ = preconFactory_->createPrec();
   preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA),thyra_mv,&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();

   SetOperator(preconditioner,false);

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
}

/** \brief Rebuild this preconditioner from an Epetra_Operator passed
  * in this to object. 
  *
  * Rebuild this preconditioner from an Epetra_Operator passed
  * in this to object.  If <code>buildPreconditioner</code> has not been called
  * the preconditioner will be built instead. Otherwise efforts are taken
  * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
  * will be an EpetraOperatorWrapper object, so the block Thyra components
  * can be easily extracted.
  *
  * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
  * \param[in] mv A vector that was used to build the source operator.
  */
void EpetraBlockPreconditioner::rebuildPreconditioner(const Epetra_Operator & A)
{
   // if the preconditioner hasn't been built yet, rebuild from scratch
   if(preconObj_==Teuchos::null) {
      buildPreconditioner(A);
      return;
   }
   PB_DEBUG_EXPR(Teuchos::Time timer(""));

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   PB_DEBUG_EXPR(timer.start(true));
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild get thyraop time =  " << timer.totalElapsedTime(),2);

   // reinitialize the preconditioner
   PB_DEBUG_EXPR(timer.start(true));
   preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA),&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild initialize prec time =  " << timer.totalElapsedTime(),2);

   PB_DEBUG_EXPR(timer.start(true));
   SetOperator(preconditioner,false);
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild set operator time =  " << timer.totalElapsedTime(),2);

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
}

/** \brief Rebuild this preconditioner from an Epetra_Operator passed
  * in this to object. 
  *
  * Rebuild this preconditioner from an Epetra_Operator passed
  * in this to object.  If <code>buildPreconditioner</code> has not been called
  * the preconditioner will be built instead. Otherwise efforts are taken
  * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
  * will be an EpetraOperatorWrapper object, so the block Thyra components
  * can be easily extracted.
  *
  * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
  * \param[in] mv A vector that was used to build the source operator.
  */
void EpetraBlockPreconditioner::rebuildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & epetra_mv)
{
   // if the preconditioner hasn't been built yet, rebuild from scratch
   if(preconObj_==Teuchos::null) {
      buildPreconditioner(A);
      return;
   }
   PB_DEBUG_EXPR(Teuchos::Time timer(""));

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   PB_DEBUG_EXPR(timer.start(true));
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild get thyraop time =  " << timer.totalElapsedTime(),2);

   // build the thyra version of the source multivector
   PB_DEBUG_EXPR(timer.start(true));
   RCP<Thyra::MultiVectorBase<double> > thyra_mv = Thyra::createMembers(thyraA->range(),epetra_mv.NumVectors());
   getMapStrategy()->copyEpetraIntoThyra(epetra_mv,thyra_mv.ptr(),*eow);
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild vector copy time =  " << timer.totalElapsedTime(),2);

   // reinitialize the preconditioner
   PB_DEBUG_EXPR(timer.start(true));
   preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA),thyra_mv,&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild initialize prec time =  " << timer.totalElapsedTime(),2);

   PB_DEBUG_EXPR(timer.start(true));
   SetOperator(preconditioner,false);
   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("EBP::rebuild set operator time =  " << timer.totalElapsedTime(),2);

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
}

} // end namespace Epetra
} // end namespace PB
