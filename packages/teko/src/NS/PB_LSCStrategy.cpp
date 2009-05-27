#include "NS/PB_LSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// PB includes
#include "PB_Helpers.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"
#include "Epetra/PB_EpetraHelpers.hpp"
#include "Epetra/PB_EpetraLSCHelpers.hpp"

using Teuchos::rcp_dynamic_cast;

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
InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & factory,bool rzn)
   : massMatrix_(Teuchos::null), invFactoryF_(factory), invFactoryS_(factory), eigSolveParam_(5), rowZeroingNeeded_(rzn), useFullLDU_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & invFactF,
                               const Teuchos::RCP<const InverseFactory> & invFactS,
                               bool rzn)
   : massMatrix_(Teuchos::null), invFactoryF_(invFactF), invFactoryS_(invFactS), eigSolveParam_(5), rowZeroingNeeded_(rzn), useFullLDU_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & factory,LinearOp & mass,bool rzn)
   : massMatrix_(mass), invFactoryF_(factory), invFactoryS_(factory), eigSolveParam_(5), rowZeroingNeeded_(rzn), useFullLDU_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & invFactF,
                               const Teuchos::RCP<const InverseFactory> & invFactS,
                               LinearOp & mass,bool rzn)
   : massMatrix_(mass), invFactoryF_(invFactF), invFactoryS_(invFactS), eigSolveParam_(5), rowZeroingNeeded_(rzn), useFullLDU_(false)
{ }

void InvLSCStrategy::buildState(BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);

   // if neccessary save state information
   if(not lscState->isInitialized())
      initializeState(A,lscState);
   else 
      reinitializeState(A,lscState);
}

// functions inherited from LSCStrategy
LinearOp InvLSCStrategy::getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   return buildInverse(*invFactoryS_,lscState->BQBtmC_);
}

LinearOp InvLSCStrategy::getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   const LinearOp F  = getBlock(0,0,A);
 
   return buildInverse(*invFactoryF_,F);
}

LinearOp InvLSCStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   return lscState->aiD_;
}

LinearOp InvLSCStrategy::getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

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

   // now we can just reintialize
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

   if(massMatrix_==Teuchos::null) {
      state->invMass_ = getInvDiagonalOp(F);
      state->BQBt_ = explicitMultiply(B,state->invMass_,Bt);  
   }

   // if this is a stable discretization...we are done!
   if(not isStabilized) {
      state->BQBtmC_ = state->BQBt_;
      state->gamma_ = 0.0;
      state->alpha_ = 0.0;
      state->aiD_ = Teuchos::null;

      return;
   }

   // for Epetra_CrsMatrix...zero out certain rows: this ensures spectral radius is correct
   LinearOp modF = F;
   const RCP<const Epetra_Operator> epF = Thyra::get_Epetra_Operator(*F);
   if(epF!=Teuchos::null && rowZeroingNeeded_) {
      // try to get a CRS matrix
      const RCP<const Epetra_CrsMatrix> crsF = rcp_dynamic_cast<const Epetra_CrsMatrix>(epF);

      // if it is a CRS matrix get rows that need to be zeroed
      if(crsF!=Teuchos::null) {
         std::vector<int> zeroIndicies;
          
         // get rows in need of zeroing
         PB::Epetra::identityRowIndicies(crsF->RowMap(), *crsF,zeroIndicies);

         // build an operator that zeros those rows
         modF = Thyra::epetraLinearOp(rcp(new PB::Epetra::ZeroedOperator(zeroIndicies,crsF)));
      }
   }

   // compute gamma
   LinearOp iQuF;
   if(state->invMass_!=Teuchos::null) 
      iQuF = multiply(state->invMass_,modF);
   else
      iQuF = modF;

   // do 6 power iterations to compute spectral radius: EHSST2007 Eq. 4.28
   state->gamma_ = std::fabs(PB::computeSpectralRad(iQuF,5e-2,false,eigSolveParam_))/3.0; 

   // compute alpha scaled inv(D): EHSST2007 Eq. 4.29
   const LinearOp invDiagF = getInvDiagonalOp(F);
   const LinearOp B_idF_Bt = explicitMultiply(B,invDiagF,Bt);

   MultiVector vec_D = getDiagonal(B_idF_Bt);
   update(-1.0,getDiagonal(C),1.0,vec_D); // vec_D = diag(B*inv(diag(F))*Bt)-diag(C)
   const LinearOp invD = buildInvDiagonal(vec_D,"inv(D)");

   const LinearOp BidFBtidD = multiply<double>(B_idF_Bt,invD);
   double num = std::fabs(PB::computeSpectralRad(BidFBtidD,5e-2,false,eigSolveParam_));
   state->alpha_ = 1.0/num;
   state->aiD_ = Thyra::scale(state->alpha_,invD);

   // now build B*Q*Bt-gamma*C
   state->BQBtmC_ = explicitAdd(state->BQBt_,scale(-state->gamma_,C));
}

} // end namespace NS
} // end namespace PB
