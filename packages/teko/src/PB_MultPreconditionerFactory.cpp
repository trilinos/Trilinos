#include "PB_MultPreconditionerFactory.hpp"

namespace PB {

using Teuchos::RCP;

void MultPrecsLinearOp::implicitApply(const PB::BlockedMultiVector & r, PB::BlockedMultiVector & y,
                     const double alpha, const double beta) const 
{ 
   // Casting is a bit delicate. We basically use 
   //
   //  1) deepcopy      to copy & cast BlockedMultiVectors to MultiVectors.
   //
   //  2) toMultiVector to cast BlockedMultiVectors to MultiVectors.
   //
   PB::MultiVector MOne_r = PB::deepcopy(r);
   PB::MultiVector t      = PB::deepcopy(r);
   PB::MultiVector w      = PB::toMultiVector(y);

   PB::applyOp(M1_, r, MOne_r);
   PB::applyOp(A_, MOne_r,  t);
   PB::update(1.,r,-1.,t);
   PB::applyOp(M2_, t,  w);
   PB::update(1.,MOne_r, 1.,  w);
}

//! Constructor
MultPreconditionerFactory 
   ::MultPreconditionerFactory(const RCP<const PB::BlockPreconditionerFactory> & FirstFactory,
                                  const RCP<const PB::BlockPreconditionerFactory> & SecondFactory)
   : FirstFactory_(FirstFactory), SecondFactory_(SecondFactory)
{ } 

//! Build the MultPrecondState object
RCP<PB::BlockPreconditionerState> MultPreconditionerFactory::buildPreconditionerState() const
{ 
   MultPrecondState*   mystate = new MultPrecondState(); 
   mystate->StateOne_ = FirstFactory_->buildPreconditionerState();
   mystate->StateTwo_ = SecondFactory_->buildPreconditionerState();
   return rcp(mystate);
}


//! Use the factory to build the preconditioner (this is where the work goes)
PB::LinearOp MultPreconditionerFactory
   ::buildPreconditionerOperator(PB::BlockedLinearOp & blockOp,
                                 PB::BlockPreconditionerState & state) const
{
   
   MultPrecondState *MyState = dynamic_cast<MultPrecondState *> (&state);

   TEUCHOS_ASSERT(MyState != 0);

   PB::LinearOp M1 = FirstFactory_->buildPreconditionerOperator(blockOp, *MyState->StateOne_);
   PB::LinearOp M2 = SecondFactory_->buildPreconditionerOperator(blockOp, *MyState->StateTwo_);


   /*************************************************************************
      A different way to create the same preconditioner using the funky 
      matrix representation discussed above. At the present time, there
      appears to be some kind of bug in Thrya so this doesn't work.

   const RCP<const Thyra::LinearOpBase<double>> Mat1= Thyra::block2x1(PB::identity(PB::rangeSpace(M1)) ,M1);
   const RCP<const Thyra::LinearOpBase<double>> Mat3= Thyra::block1x2(M2,PB::identity(PB::rangeSpace(M1)));
   const RCP<const Thyra::LinearOpBase<double>> Mat2= Thyra::block2x2(
                     PB::identity(PB::rangeSpace(M1)),                            PB::scale(-1.,PB::toLinearOp(blockOp)),
                     Thyra::zero<double>(PB::rangeSpace(M1),PB::domainSpace(M1)), PB::identity(PB::rangeSpace(M1)));
   PB::LinearOp invA = PB::multiply(Mat3,Mat2,Mat1);

   return invA; 
    *************************************************************************/

   // construct an implicit operator corresponding to multiplicative 
   // preconditioning, wrap it in an rcp pointer and return.

   return Teuchos::rcp(new MultPrecsLinearOp(blockOp,M1,M2));
}

} // end namespace PB
