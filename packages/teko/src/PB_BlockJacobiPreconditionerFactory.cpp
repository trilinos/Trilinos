#include "PB_BlockJacobiPreconditionerFactory.hpp"

// Thyra includes
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

// Teuchos includes
#include "Teuchos_TestForException.hpp"

using namespace PB;
using namespace Teuchos;
using namespace Thyra;

// construct a PreconditionerFactory
BlockJacobiPreconditionerFactory::BlockJacobiPreconditionerFactory(const RCP<const LinearOpBase<double> > & invD0,
                                                                   const RCP<const LinearOpBase<double> > & invD1)
         : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0,invD1)))
{ }

BlockJacobiPreconditionerFactory::BlockJacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy)
         : invOpsStrategy_(strategy)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// is this operator compatiable with the preconditioner factory?
bool BlockJacobiPreconditionerFactory::isCompatible(const LinearOpSourceBase<double> &opSrc) const 
{ 
   // probably this should be moved entirely into the strategy object.
   // however I think that will just complicate the strategy, probably
   // unnecessarily 

   const RCP<const BlockedLinearOpBase<double> > A 
         = rcp_dynamic_cast<const BlockedLinearOpBase<double> >(opSrc.getOp(),false);

   // make sure cast succeeded
   bool compatible = (A!=Teuchos::null);

   // make sure the block counts are correct
   if(compatible) {
      compatible &= (A->productRange()->numBlocks()==invOpsStrategy_->numDiagonalBlocks());
      compatible &= (A->productDomain()->numBlocks()==invOpsStrategy_->numDiagonalBlocks());
   }

   return compatible; 
}

// create an instance of the preconditioner
RCP<PreconditionerBase<double> > BlockJacobiPreconditionerFactory::createPrec() const 
{
   return rcp(new DefaultPreconditioner<double>());
}

// initialize a newly created preconditioner object
void BlockJacobiPreconditionerFactory::initializePrec(const RCP<const LinearOpSourceBase<double> > & opSrc,
      PreconditionerBase<double> * prec, const ESupportSolveUse supportSolveUse) const
{
   // get diagonal blocks
   const std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> >  > & invDiag
         = invOpsStrategy_->getInvD(opSrc->getOp());

   TEUCHOS_ASSERT(invOpsStrategy_->numDiagonalBlocks()==invDiag.size()); // sanity check

   // allocate new linear operator
   const RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blkPrecOp 
         = Thyra::defaultBlockedLinearOp<double>();

   // set digaonals on newly block operator and give it a label
   blkPrecOp->beginBlockFill(invDiag.size(),invDiag.size());
   for(int i=0;i<invDiag.size();i++)
      blkPrecOp->setBlock(i,i,invDiag[i]);
   blkPrecOp->endBlockFill();
   blkPrecOp->setObjectLabel("inv(diag("+opSrc->getOp()->getObjectLabel()+"))");

   // set preconditioner: this cast is needed to force it 
   //    to use the const version of initializeUnspecified
   const RCP<const LinearOpBase<double> > precOp = blkPrecOp;

   // must first cast that to be initialized
   DefaultPreconditioner<double> & dPrec = dyn_cast<DefaultPreconditioner<double> >(*prec);
   dPrec.initializeUnspecified(precOp);
}

// wipe clean a already initialized preconditioner object
void BlockJacobiPreconditionerFactory::uninitializePrec(PreconditionerBase<double> * prec, 
      RCP<const LinearOpSourceBase<double> > * fwdOpSrc,ESupportSolveUse *supportSolveUse) const
{
   TEST_FOR_EXCEPT_MSG(true,"\"BlockJacobiPreconditionerFactory::uninitializePrec not implemented\"");
}

// for ParameterListAcceptor
///////////////////////////////////////////////////////////////////////

// Set parameters from a parameter list and return with default values.
void BlockJacobiPreconditionerFactory::setParameterList(const RCP<ParameterList> & paramList)
{
   paramList_ = paramList;
}

// Get the parameter list that was set using setParameterList().
RCP<ParameterList> BlockJacobiPreconditionerFactory::getNonconstParameterList() 
{ 
   return paramList_;
}

// Unset the parameter list that was set using setParameterList(). 
RCP<ParameterList> BlockJacobiPreconditionerFactory::unsetParameterList()
{
   RCP<ParameterList> _paramList = paramList_;
   paramList_ = Teuchos::null;
   return _paramList;
}
