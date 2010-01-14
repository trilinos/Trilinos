#include "PB_EpetraBlockPreconditioner.hpp"

// Thyra includes
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_EpetraLinearOp.hpp"

// Teuchos includes
#include "Teuchos_Time.hpp"

// Teko includes
#include "Teko_BasicMappingStrategy.hpp"

namespace Teko {
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
   : preconFactory_(bfp), firstBuildComplete_(false)
{ 
}

void EpetraBlockPreconditioner::initPreconditioner(bool clearOld)
{
   if((not clearOld) && preconObj_!=Teuchos::null)
      return;
   preconObj_ = preconFactory_->createPrec();
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
  * 
  * \note This will clear any internal state stored by the state object
  */
void EpetraBlockPreconditioner::buildPreconditioner(const Epetra_Operator & A,bool clear)
{
   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   // const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   // RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 
   RCP<const Thyra::LinearOpBase<double> > thyraA = extractLinearOp(A);

   // set the mapping strategy
   // SetMapStrategy(rcp(new InverseMappingStrategy(eow->getMapStrategy())));
   SetMapStrategy(rcp(new InverseMappingStrategy(extractMappingStrategy(A))));

   // build preconObj_ 
   initPreconditioner(clear);
   
   // actually build the preconditioner
   RCP<const Thyra::LinearOpSourceBase<double> > lOpSrc = Thyra::defaultLinearOpSource(thyraA);
   preconFactory_->initializePrec(lOpSrc,&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);

   // extract preconditioner operator
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();

   SetOperator(preconditioner,false);

   firstBuildComplete_ = true;

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
   TEUCHOS_ASSERT(firstBuildComplete_==true);
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
  *
  * \note This will clear any internal state stored by the state object
  */
void EpetraBlockPreconditioner::buildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & epetra_mv,bool clear)
{
   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   // const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   // RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 
   RCP<const Thyra::LinearOpBase<double> > thyraA = extractLinearOp(A);

   // set the mapping strategy
   // SetMapStrategy(rcp(new InverseMappingStrategy(eow->getMapStrategy())));
   SetMapStrategy(rcp(new InverseMappingStrategy(extractMappingStrategy(A))));

   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
   
   // build the thyra version of the source multivector
   RCP<Thyra::MultiVectorBase<double> > thyra_mv = Thyra::createMembers(thyraA->range(),epetra_mv.NumVectors());
   // getMapStrategy()->copyEpetraIntoThyra(epetra_mv,thyra_mv.ptr(),*eow);
   getMapStrategy()->copyEpetraIntoThyra(epetra_mv,thyra_mv.ptr(),*this);

   // build preconObj_ 
   initPreconditioner(clear);

   // actually build the preconditioner
   preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA),thyra_mv,&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();

   SetOperator(preconditioner,false);

   firstBuildComplete_ = true;

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(getMapStrategy()!=Teuchos::null);
   TEUCHOS_ASSERT(firstBuildComplete_==true);
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
   if(not firstBuildComplete_) {
      buildPreconditioner(A,false);
      return;
   }
   Teko_DEBUG_EXPR(Teuchos::Time timer(""));

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   Teko_DEBUG_EXPR(timer.start(true));
   // const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   // RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 
   RCP<const Thyra::LinearOpBase<double> > thyraA = extractLinearOp(A);
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild get thyraop time =  " << timer.totalElapsedTime(),2);

   // reinitialize the preconditioner
   Teko_DEBUG_EXPR(timer.start(true));
   preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA),&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild initialize prec time =  " << timer.totalElapsedTime(),2);

   Teko_DEBUG_EXPR(timer.start(true));
   SetOperator(preconditioner,false);
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild set operator time =  " << timer.totalElapsedTime(),2);

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(firstBuildComplete_==true);
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
   if(not firstBuildComplete_) {
      buildPreconditioner(A,epetra_mv,false);
      return;
   }
   Teko_DEBUG_EXPR(Teuchos::Time timer(""));

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   Teko_DEBUG_EXPR(timer.start(true));
   // const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(rcpFromRef(A),true);
   // RCP<const Thyra::LinearOpBase<double> > thyraA = eow->getThyraOp(); 
   RCP<const Thyra::LinearOpBase<double> > thyraA = extractLinearOp(A);
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild get thyraop time =  " << timer.totalElapsedTime(),2);

   // build the thyra version of the source multivector
   Teko_DEBUG_EXPR(timer.start(true));
   RCP<Thyra::MultiVectorBase<double> > thyra_mv = Thyra::createMembers(thyraA->range(),epetra_mv.NumVectors());
   // getMapStrategy()->copyEpetraIntoThyra(epetra_mv,thyra_mv.ptr(),*eow);
   getMapStrategy()->copyEpetraIntoThyra(epetra_mv,thyra_mv.ptr(),*this);
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild vector copy time =  " << timer.totalElapsedTime(),2);

   // reinitialize the preconditioner
   Teko_DEBUG_EXPR(timer.start(true));
   preconFactory_->initializePrec(Thyra::defaultLinearOpSource(thyraA),thyra_mv,&*preconObj_,Thyra::SUPPORT_SOLVE_UNSPECIFIED);
   RCP<const Thyra::LinearOpBase<double> > preconditioner = preconObj_->getUnspecifiedPrecOp();
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild initialize prec time =  " << timer.totalElapsedTime(),2);

   Teko_DEBUG_EXPR(timer.start(true));
   SetOperator(preconditioner,false);
   Teko_DEBUG_EXPR(timer.stop());
   Teko_DEBUG_MSG("EBP::rebuild set operator time =  " << timer.totalElapsedTime(),2);

   TEUCHOS_ASSERT(preconObj_!=Teuchos::null);
   TEUCHOS_ASSERT(getThyraOp()!=Teuchos::null);
   TEUCHOS_ASSERT(firstBuildComplete_==true);
}

/** Try to get a <code>Teko::BlockPreconditionerState</code> object. This method
  * attempts to cast its internal representation of a preconditioner 
  * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a 
  * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
  *
  * \returns Get the state object associated with this preconditioner.
  *          If it doesn't exist for this type of preconditioner factory
  *          this method returns null.
  */
Teuchos::RCP<BlockPreconditionerState> EpetraBlockPreconditioner::getPreconditionerState()
{
   Teuchos::RCP<BlockPreconditioner> bp = rcp_dynamic_cast<BlockPreconditioner>(preconObj_);

   if(bp!=Teuchos::null)
      return bp->getStateObject();
 
   return Teuchos::null;
}

/** Try to get a <code>Teko::BlockPreconditionerState</code> object. This method
  * attempts to cast its internal representation of a preconditioner 
  * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a 
  * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
  *
  * \returns Get the state object associated with this preconditioner.
  *          If it doesn't exist for this type of preconditioner factory
  *          this method returns null.
  */
Teuchos::RCP<const BlockPreconditionerState> EpetraBlockPreconditioner::getPreconditionerState() const
{
   Teuchos::RCP<const BlockPreconditioner> bp = rcp_dynamic_cast<const BlockPreconditioner>(preconObj_);

   if(bp!=Teuchos::null)
      return bp->getStateObject();
 
   return Teuchos::null;
}

Teuchos::RCP<const Thyra::LinearOpBase<double> > EpetraBlockPreconditioner::extractLinearOp(const Epetra_Operator & A) const
{
   Teuchos::RCP<const Epetra_Operator> ptrA = rcpFromRef(A);

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(ptrA);
  
   // if it is an EpetraOperatorWrapper, then get the Thyra operator
   if(eow!=Teuchos::null)
      return eow->getThyraOp(); 

   // otherwise wrap it up as a thyra operator 
   return Thyra::epetraLinearOp(ptrA);
}

Teuchos::RCP<const MappingStrategy> EpetraBlockPreconditioner::extractMappingStrategy(const Epetra_Operator & A) const
{
   Teuchos::RCP<const Epetra_Operator> ptrA = rcpFromRef(A);

   // extract EpetraOperatorWrapper (throw on failure) and corresponding thyra operator
   const RCP<const EpetraOperatorWrapper> & eow = rcp_dynamic_cast<const EpetraOperatorWrapper>(ptrA);
  
   // if it is an EpetraOperatorWrapper, then get the Thyra operator
   if(eow!=Teuchos::null)
      return eow->getMapStrategy(); 

   // otherwise wrap it up as a thyra operator 
   RCP<const Epetra_Map> range = rcpFromRef(A.OperatorRangeMap());
   RCP<const Epetra_Map> domain = rcpFromRef(A.OperatorDomainMap());
   return rcp(new BasicMappingStrategy(range,domain,A.Comm()));
}

} // end namespace Epetra
} // end namespace Teko
