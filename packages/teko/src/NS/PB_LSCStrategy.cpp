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
#include "PB_Utilities.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"
#include "Epetra/PB_EpetraHelpers.hpp"
#include "Epetra/PB_EpetraOperatorWrapper.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;

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
// InvLSCStrategy Implementation
//////////////////////////////////////////////

// helper function for a lid driven cavity-like problem
// This function _WILL_ change the operator
PB::LinearOp reduceCrsOperator(PB::LinearOp & op,const std::vector<int> & zeroIndicies)
{
   // Extract a non-const version of the operator
   RCP<Epetra_Operator> eOp = get_Epetra_Operator(*rcp_const_cast<Thyra::LinearOpBase<double> >(op));
   RCP<Epetra_CrsMatrix> eCrsOp = rcp_dynamic_cast<Epetra_CrsMatrix>(eOp);

   // build vector for scaling
   Epetra_Vector scaling(eCrsOp->OperatorRangeMap());
   scaling.PutScalar(1.0);

   // figure out which zero index (if any) this processor owns
   std::vector<std::pair<int,int> > localIndicies; // local, global indicies
   for(int i=0;i<zeroIndicies.size();i++) {
      // check if this zero index is owned by this processor
      int local = eCrsOp->LRID(zeroIndicies[i]);
      if(local>=0) 
         localIndicies.push_back(std::make_pair(local,zeroIndicies[i]));
   }

   // set a number of rows to zero
   for(int i=0;i<localIndicies.size();i++) 
      TEST_FOR_EXCEPT(scaling.ReplaceGlobalValue(localIndicies[i].second,0,0.0));

   // wipe out all but the desired rows and columns
   eCrsOp->RightScale(scaling);
   eCrsOp->LeftScale(scaling);

   #if 0
   // so the matrix is still invertable...set the digaonals of the
   // wiped rows to 1
   for(int i=0;i<localIndicies.size();i++) {
      double value = 1.0;
      int index = localIndicies[i].second;

      // set diagonal to one
      TEST_FOR_EXCEPT(eCrsOp->ReplaceGlobalValues(index,1,&value,&index));
   } 
   #endif

   // now wrap the crs matrix in a ZeroedEpetraOperator
   return Thyra::epetraLinearOp(rcp(new PB::Epetra::ZeroedOperator(zeroIndicies,eCrsOp)));
}

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
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::buildState",10);

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);

   // if neccessary save state information
   if(not lscState->isInitialized())
      initializeState(A,lscState);
   else 
      reinitializeState(A,lscState);

   PB_DEBUG_MSG("END InvLSCStrategy::buildState",10);
}

// functions inherited from LSCStrategy
LinearOp InvLSCStrategy::getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::getInvBQBt",10);

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   // (re)build the inverse of the Schur complement
   PB::ModifiableLinearOp BQBtmC = state.getInverse("BQBtmC");
   if(lscState->invBQBtmC_==Teuchos::null)
      lscState->invBQBtmC_ = buildInverse(*invFactoryS_,BQBtmC);
   else
      rebuildInverse(*invFactoryS_,BQBtmC,lscState->invBQBtmC_);

   // if(nullPresIndicies_.size()>0) {
   //    PB::LinearOp inv = lscState->invBQBtmC_;
   //    RCP<Epetra_Operator> eOp = rcp(new PB::Epetra::EpetraOperatorWrapper(inv));
   //    return Thyra::epetraLinearOp(rcp(new PB::Epetra::ZeroedOperator(nullPresIndicies_,eOp)));
   // } 

   PB_DEBUG_MSG("END InvLSCStrategy::getInvBQBt",10);

   return lscState->invBQBtmC_;
}

LinearOp InvLSCStrategy::getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::getInvF",10);

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   const LinearOp F  = getBlock(0,0,A);

   // (re)build the inverse of F
   InverseLinearOp invF = state.getInverse("invF");
   if(invF==Teuchos::null) {
      invF = buildInverse(*invFactoryF_,F);
      state.addInverse("invF",invF); 
   } else {
      rebuildInverse(*invFactoryF_,F,invF);
   }

   PB_DEBUG_MSG("END InvLSCStrategy::getInvF",10);
   return invF;
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
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::initiailzeState",10);

   const LinearOp B  = getBlock(1,0,A);
   const LinearOp Bt = getBlock(0,1,A);

   if(massMatrix_!=Teuchos::null) {
      state->invMass_ = getInvDiagonalOp(massMatrix_);
      state->BQBt_ = explicitMultiply(B,state->invMass_,Bt,state->BQBt_);
   }

   // now we can just reintialize
   state->setInitialized(true);

   // do some real work
   reinitializeState(A,state);

   PB_DEBUG_MSG("END InvLSCStrategy::initiailzeState",10);
}

void InvLSCStrategy::reinitializeState(const BlockedLinearOp & A,LSCPrecondState * state) const
{
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::reinitiailzeState",10);

   const LinearOp F  = getBlock(0,0,A);
   const LinearOp Bt = getBlock(0,1,A);
   const LinearOp B  = getBlock(1,0,A);
   const LinearOp C  = getBlock(1,1,A);

   bool isStabilized = (not isZeroOp(C));

   if(massMatrix_==Teuchos::null) {
      state->invMass_ = getInvDiagonalOp(F);
      state->BQBt_ = explicitMultiply(B,state->invMass_,Bt,state->BQBt_);  
   }

   // if this is a stable discretization...we are done!
   if(not isStabilized) {
      state->addInverse("BQBtmC",state->BQBt_);
      state->gamma_ = 0.0;
      state->alpha_ = 0.0;
      state->aiD_ = Teuchos::null;

      PB_DEBUG_MSG("END InvLSCStrategy::reinitiailzeState",10);
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

   if(graphLaplacian_!=Teuchos::null) {
      PB::LinearOp invDGl = PB::getInvDiagonalOp(graphLaplacian_);
      PB::LinearOp gammaOp = multiply(invDGl,C);
      state->gamma_ *= std::fabs(PB::computeSpectralRad(gammaOp,5e-2,false,eigSolveParam_));
   }

   // compute alpha scaled inv(D): EHSST2007 Eq. 4.29
   const LinearOp invDiagF = getInvDiagonalOp(F);
   
   // construct B_idF_Bt and save it for refilling later
   PB::ModifiableLinearOp modB_idF_Bt = state->getInverse("BidFBt");
   modB_idF_Bt = explicitMultiply(B,invDiagF,Bt,modB_idF_Bt);
   state->addInverse("BidFBt",modB_idF_Bt);
   const LinearOp B_idF_Bt = modB_idF_Bt;

   MultiVector vec_D = getDiagonal(B_idF_Bt);
   update(-1.0,getDiagonal(C),1.0,vec_D); // vec_D = diag(B*inv(diag(F))*Bt)-diag(C)
   const LinearOp invD = buildInvDiagonal(vec_D,"inv(D)");

   const LinearOp BidFBtidD = multiply<double>(B_idF_Bt,invD);
   double num = std::fabs(PB::computeSpectralRad(BidFBtidD,5e-2,false,eigSolveParam_));
   state->alpha_ = 1.0/num;
   state->aiD_ = Thyra::scale(state->alpha_,invD);

   // now build B*Q*Bt-gamma*C
   PB::ModifiableLinearOp BQBtmC = state->getInverse("BQBtmC");
   if(graphLaplacian_==Teuchos::null)
      BQBtmC = explicitAdd(state->BQBt_,scale(-state->gamma_,C),BQBtmC);
   else
      BQBtmC = explicitAdd(state->BQBt_,scale(state->gamma_,graphLaplacian_),BQBtmC);
   state->addInverse("BQBtmC",BQBtmC);

   PB_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LSC Gamma Parameter = " << state->gamma_ << std::endl;
      DEBUG_STREAM << "LSC Alpha Parameter = " << state->alpha_ << std::endl;
   PB_DEBUG_MSG_END()
   PB_DEBUG_MSG("END InvLSCStrategy::reinitiailzeState",10);
}

} // end namespace NS
} // end namespace PB
