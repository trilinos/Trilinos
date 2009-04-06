#include "NS/PB_LSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

// PB includes
#include "PB_Helpers.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"
#include "Epetra/PB_EpetraHelpers.hpp"

namespace PB {
namespace NS {

   // Staiblized constructor
StaticLSCStrategy::StaticLSCStrategy(const LinearOp & invF,
                                     const LinearOp & invBQBtmC,
                                     const LinearOp & invD,
                                     const LinearOp & invMass)
   : invF_(invF), invBQBtmC_(invBQBtmC), invD_(invD), invMass_(invMass)
{ }
 
   // Stable constructor
StaticLSCStrategy::StaticLSCStrategy(const LinearOp & invF,
                                     const LinearOp & invBQBtmC,
                                     const LinearOp & invMass)
   : invF_(invF), invBQBtmC_(invBQBtmC), invD_(Teuchos::null), invMass_(invMass)
{ }

//////////////////////////////////////////////

// constructors
InvLSCStrategy::InvLSCStrategy(InverseFactory & factory)
   : massMatrix_(Teuchos::null), invFactory_(factory)
{ }

InvLSCStrategy::InvLSCStrategy(InverseFactory & factory,LinearOp & mass)
   : massMatrix_(mass), invFactory_(factory)
{ }

// functions inherited from LSCStrategy
LinearOp InvLSCStrategy::getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);

   // if neccessary save state information
   if(not lscState->isInitialized())
      initializeState(A,lscState);
   else 
      reinitializeState(A,lscState);

   return buildInverse(invFactory_,lscState->BQBtmC_);
}

LinearOp InvLSCStrategy::getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   const LinearOp F  = getBlock(0,0,A);
 
   return buildInverse(invFactory_,F);
}

LinearOp InvLSCStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);

   // if neccessary save state information
   if(not lscState->isInitialized())
      initializeState(A,lscState);
   else 
      reinitializeState(A,lscState);
 
   return lscState->aiD_;
}

LinearOp InvLSCStrategy::getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);

   return lscState->invMass_;
}

//! Initialize the state object using this blocked linear operator
void InvLSCStrategy::initializeState(const BlockedLinearOp & A,LSCPrecondState * state) const
{
   const LinearOp B  = getBlock(1,0,A);
   const LinearOp Bt = getBlock(0,1,A);

   if(massMatrix_!=Teuchos::null) {
      state->invMass_ = getInvDiagonalOp(massMatrix_);
      state->BQBt_ = explicitMultiply(B,state->invMass_,Bt);  
   }
   else {
      state->BQBt_ = explicitMultiply(B,Bt);  
   }

   state->setInitialized(true);

   // do some real work
   reinitializeState(A,state);
}

void InvLSCStrategy::reinitializeState(const BlockedLinearOp & A,LSCPrecondState * state) const
{
   const LinearOp F  = getBlock(0,0,A);
   const LinearOp Bt = getBlock(0,1,A);
   const LinearOp B  = getBlock(1,0,A);
   const LinearOp C  = getBlock(1,1,A);

   bool isStabilized = (not isZeroOp(C));

   // compute gamma
   LinearOp iQuF;
   if(state->invMass_!=Teuchos::null) 
      iQuF = multiply(state->invMass_,F);
   else
      iQuF = F;

   // if this is a stable discretization...we are done!
   if(not isStabilized) {
      state->BQBtmC_ = state->BQBt_;
      state->gamma_ = 0.0;
      state->alpha_ = 0.0;
      state->aiD_ = Teuchos::null;

      return;
   }

   // do 6 power iterations to compute spectral radius: EHSST2007 Eq. 4.28
   state->gamma_ = std::abs(PB::computeSpectralRad(iQuF,5e-2,false,5))/3.0; 

   // compute alpha scaled inv(D): EHSST2007 Eq. 4.29
   const LinearOp invDiagF = getInvDiagonalOp(F);
   const LinearOp B_idF_Bt = explicitMultiply(B,invDiagF,Bt);

   MultiVector vec_D = getDiagonal(B_idF_Bt);
   update(1.0,getDiagonal(C),1.0,vec_D); // vec_D = diag(B*inv(diag(F))*Bt)+diag(C)
   const LinearOp invD = buildInvDiagonal(vec_D);

   const LinearOp BidFBtidD = multiply<double>(B_idF_Bt,invD);
   state->alpha_ = 1.0/std::abs(PB::computeSpectralRad(BidFBtidD,5e-2,false,5));
   state->aiD_ = Thyra::scale(state->alpha_,invD);

   // now build B*Q*Bt-gamma*C
   state->BQBtmC_ = explicitMultiply(state->BQBt_,scale(-state->gamma_,C));
}

} // end namespace NS
} // end namespace PB
