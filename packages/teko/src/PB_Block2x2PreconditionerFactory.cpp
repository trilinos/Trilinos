#include "PB_Block2x2PreconditionerFactory.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

// Teuchos includes
#include "Teuchos_TestForException.hpp"

// PB includes
#include "PB_SchurSolveLinearOp.hpp"

using namespace Teuchos;
using namespace Thyra;

namespace PB {

// construct a PreconditionerFactory
Block2x2PreconditionerFactory::Block2x2PreconditionerFactory(const RCP<const LinearOpBase<double> > & invA00,
      const RCP<const LinearOpBase<double> > & invS)
         : invOpsStrategy_(Teuchos::null)
{ 
   // create a static strategy for this object
   invOpsStrategy_ = rcp(new StaticBlock2x2Strategy(invA00,invS));
}

Block2x2PreconditionerFactory::Block2x2PreconditionerFactory(const RCP<const Block2x2Strategy> & strategy)
   : invOpsStrategy_(strategy)
{ }

// for PreconditionerFactoryBase
///////////////////////////////////////////////////////////////////////

// is this operator compatiable with the preconditioner factory?
bool Block2x2PreconditionerFactory::isCompatible(const LinearOpSourceBase<double> &fwdOpSrc) const 
{ 
   TEST_FOR_EXCEPT_MSG(true,"\"Block2x2PreconditionerFactory::isCompatible not implemented\"");
   return true; 
}

// create an instance of the preconditioner
RCP<PreconditionerBase<double> > Block2x2PreconditionerFactory::createPrec() const 
{
   return rcp(new DefaultPreconditioner<double>());
}

// initialize a newly created preconditioner object
void Block2x2PreconditionerFactory::initializePrec(const RCP<const LinearOpSourceBase<double> > & opSrc,
      PreconditionerBase<double> * prec, const ESupportSolveUse supportSolveUse) const
{
   const RCP<const LinearOpBase<double> > A = opSrc->getOp();
   const RCP<const LinearOpBase<double> > invA00 = invOpsStrategy_->getInvA00(A);
   const RCP<const LinearOpBase<double> > invS   = invOpsStrategy_->getInvS(A);

   const RCP<const LinearOpBase<double> > precOp = build2x2InverseOperator(A,invA00,invS);
 
   // must first cast that to be initialized
   DefaultPreconditioner<double> & dPrec = dyn_cast<DefaultPreconditioner<double> >(*prec);
   dPrec.initializeUnspecified(precOp); // should this be initializeUnspecified?
}


// wipe clean a already initialized preconditioner object
void Block2x2PreconditionerFactory::uninitializePrec(PreconditionerBase<double> * prec, 
      RCP<const LinearOpSourceBase<double> > * fwdOpSrc,ESupportSolveUse *supportSolveUse) const
{
   TEST_FOR_EXCEPT_MSG(true,"\"Block2x2PreconditionerFactory::uninitializePrec not implemented\"");
}

// class specific functionality
///////////////////////////////////////////////////////////////////////

// build a 2x2 inverse operator
const RCP<const LinearOpBase<double> > 
Block2x2PreconditionerFactory::build2x2InverseOperator(const RCP<const LinearOpBase<double> > & A,
                                                       const RCP<const LinearOpBase<double> > & invA00, 
                                                       const RCP<const LinearOpBase<double> > & invS) const
{
   // cast to a blocked operator
   const RCP<const BlockedLinearOpBase<double> > blkA 
         = rcp_dynamic_cast<const BlockedLinearOpBase<double> >(A);

   // return a Schur complement solve operator
   return rcp(new SchurSolveLinearOp(blkA,invA00,invS));
}

// for ParameterListAcceptor
///////////////////////////////////////////////////////////////////////

// Set parameters from a parameter list and return with default values.
void Block2x2PreconditionerFactory::setParameterList(const RCP<ParameterList> & paramList)
{
   paramList_ = paramList;
}

// Get the parameter list that was set using setParameterList().
RCP<ParameterList> Block2x2PreconditionerFactory::getNonconstParameterList() 
{ 
   return paramList_;
}

// Unset the parameter list that was set using setParameterList(). 
RCP<ParameterList> Block2x2PreconditionerFactory::unsetParameterList()
{
   RCP<ParameterList> _paramList = paramList_;
   paramList_ = Teuchos::null;
   return _paramList;
}

} // end namespace PB
