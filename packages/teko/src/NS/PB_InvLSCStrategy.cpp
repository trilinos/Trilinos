#include "NS/PB_InvLSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "Teuchos_Time.hpp"

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

//////////////////////////////////////////////
// InvLSCStrategy Implementation
//////////////////////////////////////////////

// helper function for a lid driven cavity-like problem
// This function _WILL_ change the operator
PB::ModifiableLinearOp reduceCrsOperator(PB::ModifiableLinearOp & op,const std::vector<int> & zeroIndicies)
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

   #if 1
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
   // return Thyra::epetraLinearOp(rcp(new PB::Epetra::ZeroedOperator(zeroIndicies,eCrsOp)));
   return Thyra::nonconstEpetraLinearOp(eCrsOp);
}

// constructors
InvLSCStrategy::InvLSCStrategy()
   : massMatrix_(Teuchos::null), invFactoryF_(Teuchos::null), invFactoryS_(Teuchos::null), eigSolveParam_(5)
   , rowZeroingNeeded_(false), useFullLDU_(false), useMass_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & factory,bool rzn)
   : massMatrix_(Teuchos::null), invFactoryF_(factory), invFactoryS_(factory), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & invFactF,
                               const Teuchos::RCP<const InverseFactory> & invFactS,
                               bool rzn)
   : massMatrix_(Teuchos::null), invFactoryF_(invFactF), invFactoryS_(invFactS), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & factory,LinearOp & mass,bool rzn)
   : massMatrix_(mass), invFactoryF_(factory), invFactoryS_(factory), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false)
{ }

InvLSCStrategy::InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & invFactF,
                               const Teuchos::RCP<const InverseFactory> & invFactS,
                               LinearOp & mass,bool rzn)
   : massMatrix_(mass), invFactoryF_(invFactF), invFactoryS_(invFactS), eigSolveParam_(5), rowZeroingNeeded_(rzn)
   , useFullLDU_(false), useMass_(false)
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
   PB_DEBUG_EXPR(Teuchos::Time timer(""));
   PB_DEBUG_EXPR(timer.start());

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   // (re)build the inverse of the Schur complement
   PB::ModifiableLinearOp BQBtmC = state.getInverse("BQBtmC");
   if(lscState->invBQBtmC_==Teuchos::null)
      lscState->invBQBtmC_ = buildInverse(*invFactoryS_,BQBtmC);
   else
      rebuildInverse(*invFactoryS_,BQBtmC,lscState->invBQBtmC_);

   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("LSC::getInvBQBt RunTime = " << timer.totalElapsedTime(),1);
   PB_DEBUG_MSG("END InvLSCStrategy::getInvBQBt",10);

   return lscState->invBQBtmC_;
}

LinearOp InvLSCStrategy::getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::getInvF",10);
   PB_DEBUG_EXPR(Teuchos::Time timer(""));
   PB_DEBUG_EXPR(timer.start());

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

   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("LSC::getInvF RunTime = " << timer.totalElapsedTime(),1);
   PB_DEBUG_MSG("END InvLSCStrategy::getInvF",10);
   return invF;
}

LinearOp InvLSCStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   PB_DEBUG_EXPR(Teuchos::Time timer(""));
   PB_DEBUG_EXPR(timer.start());

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("LSC::getInvD RunTime = " << timer.totalElapsedTime(),1);

   return lscState->aiD_;
}

LinearOp InvLSCStrategy::getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const
{
   PB_DEBUG_EXPR(Teuchos::Time timer(""));
   PB_DEBUG_EXPR(timer.start());

   LSCPrecondState * lscState = dynamic_cast<LSCPrecondState*>(&state);
   TEUCHOS_ASSERT(lscState!=0);
   TEUCHOS_ASSERT(lscState->isInitialized())

   PB_DEBUG_EXPR(timer.stop());
   PB_DEBUG_MSG("LSC::getInvMass RunTime = " << timer.totalElapsedTime(),1);

   return lscState->invMass_;
}

//! Initialize the state object using this blocked linear operator
void InvLSCStrategy::initializeState(const BlockedLinearOp & A,LSCPrecondState * state) const
{
   PB_DEBUG_MSG("BEGIN InvLSCStrategy::initiailzeState",10);

   const LinearOp B  = getBlock(1,0,A);
   const LinearOp Bt = getBlock(0,1,A);

   if(massMatrix_!=Teuchos::null) {
      PB_DEBUG_MSG("LSC::init using passed in mass matrix for scaling",1); 
      state->invMass_ = getInvDiagonalOp(massMatrix_);

      PB_DEBUG_EXPR(Teuchos::Time matTimer(""));
      PB_DEBUG_EXPR(matTimer.start(true));
      state->BQBt_ = explicitMultiply(B,state->invMass_,Bt,state->BQBt_);
      PB_DEBUG_EXPR(matTimer.stop());
      PB_DEBUG_MSG("LSC::init BQBt MultTime = " << matTimer.totalElapsedTime(),1); 
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
   PB_DEBUG_EXPR(Teuchos::Time matTimer(""));
   PB_DEBUG_EXPR(Teuchos::Time totalTimer(""));
   PB_DEBUG_EXPR(totalTimer.start());

   const LinearOp F  = getBlock(0,0,A);
   const LinearOp Bt = getBlock(0,1,A);
   const LinearOp B  = getBlock(1,0,A);
   const LinearOp C  = getBlock(1,1,A);

   bool isStabilized = (not isZeroOp(C));

   if(massMatrix_==Teuchos::null) {
      PB_DEBUG_MSG("LSC::reinit using diag(F) for scaling",1); 
      state->invMass_ = getInvDiagonalOp(F);

      PB_DEBUG_EXPR(matTimer.start(true));
      state->BQBt_ = explicitMultiply(B,state->invMass_,Bt,state->BQBt_);  
      PB_DEBUG_EXPR(matTimer.stop());
      PB_DEBUG_MSG("LSC::reinit BQBt MultTime = " << matTimer.totalElapsedTime(),1); 
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
   PB::LinearOp stabMatrix; // this is the pressure stabilization matrix to use
   state->gamma_ = std::fabs(PB::computeSpectralRad(iQuF,5e-2,false,eigSolveParam_))/3.0; 
   if(userPresStabMat_!=Teuchos::null) {
      PB::LinearOp invDGl = PB::getInvDiagonalOp(userPresStabMat_);
      PB::LinearOp gammaOp = multiply(invDGl,C);
      state->gamma_ *= std::fabs(PB::computeSpectralRad(gammaOp,5e-2,false,eigSolveParam_));
      stabMatrix = userPresStabMat_;
   } else 
      stabMatrix = C;

   // compute alpha scaled inv(D): EHSST2007 Eq. 4.29
   const LinearOp invDiagF = getInvDiagonalOp(F);
   
   // construct B_idF_Bt and save it for refilling later: This could reuse BQBt graph
   PB::ModifiableLinearOp modB_idF_Bt = state->getInverse("BidFBt");

   PB_DEBUG_EXPR(matTimer.start(true));
   modB_idF_Bt = explicitMultiply(B,invDiagF,Bt,modB_idF_Bt);
   PB_DEBUG_EXPR(matTimer.stop());
   PB_DEBUG_MSG("LSC::reinit BFBt MultTime = " << matTimer.totalElapsedTime(),1); 
   state->addInverse("BidFBt",modB_idF_Bt);
   const LinearOp B_idF_Bt = modB_idF_Bt;

   MultiVector vec_D = getDiagonal(B_idF_Bt); // this memory could be reused
   update(-1.0,getDiagonal(C),1.0,vec_D); // vec_D = diag(B*inv(diag(F))*Bt)-diag(C)
   const LinearOp invD = buildInvDiagonal(vec_D,"inv(D)");

   const LinearOp BidFBtidD = multiply<double>(B_idF_Bt,invD);
   double num = std::fabs(PB::computeSpectralRad(BidFBtidD,5e-2,false,eigSolveParam_));
   state->alpha_ = 1.0/num;
   state->aiD_ = Thyra::scale(state->alpha_,invD);

   // now build B*Q*Bt-gamma*C
   PB::ModifiableLinearOp BQBtmC = state->getInverse("BQBtmC");
   PB_DEBUG_EXPR(matTimer.start(true));
   BQBtmC = explicitAdd(state->BQBt_,scale(-state->gamma_,stabMatrix),BQBtmC);
   PB_DEBUG_EXPR(matTimer.stop());
   PB_DEBUG_MSG("LSC::reinit BFBtgC AddTime = " << matTimer.totalElapsedTime(),1); 
   state->addInverse("BQBtmC",BQBtmC);

   PB_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LSC Gamma Parameter = " << state->gamma_ << std::endl;
      DEBUG_STREAM << "LSC Alpha Parameter = " << state->alpha_ << std::endl;
   PB_DEBUG_MSG_END()
   PB_DEBUG_EXPR(totalTimer.stop());
   PB_DEBUG_MSG("LSC reinit TotalTime = " << totalTimer.totalElapsedTime(),1);
   PB_DEBUG_MSG("END InvLSCStrategy::reinitializeState: " << totalTimer.totalElapsedTime(),10);
}

//! Initialize from a parameter list
void InvLSCStrategy::initializeFromParameterList(const Teuchos::ParameterList & pl,const InverseLibrary & invLib) 
{
   // get string specifying inverse
   std::string invStr="", invVStr="", invPStr="";
   bool rowZeroing = true;
   bool useLDU = false;

   // "parse" the parameter list
   if(pl.isParameter("Inverse Type"))
      invStr = pl.get<std::string>("Inverse Type");
   if(pl.isParameter("Inverse Velocity Type"))
      invVStr = pl.get<std::string>("Inverse Velocity Type");
   if(pl.isParameter("Inverse Pressure Type")) 
      invPStr = pl.get<std::string>("Inverse Pressure Type");
   if(pl.isParameter("Ignore Boundary Rows"))
      rowZeroing = pl.get<bool>("Ignore Boundary Rows");
   if(pl.isParameter("Use LDU"))
      useLDU = pl.get<bool>("Use LDU");
   if(pl.isParameter("Use Mass Scaling"))
      useMass_ = pl.get<bool>("Use Mass Scaling");

   PB_DEBUG_MSG_BEGIN(5)
      DEBUG_STREAM << "LSC Inverse Strategy Parameters: " << std::endl;
      DEBUG_STREAM << "   inv type   = \"" << invStr  << "\"" << std::endl;
      DEBUG_STREAM << "   inv v type = \"" << invVStr << "\"" << std::endl;
      DEBUG_STREAM << "   inv p type = \"" << invPStr << "\"" << std::endl;
      DEBUG_STREAM << "   bndry rows = " << rowZeroing << std::endl;
      DEBUG_STREAM << "   use ldu    = " << useLDU << std::endl;
      DEBUG_STREAM << "   use mass    = " << useMass_ << std::endl;
      DEBUG_STREAM << "LSC  Inverse Strategy Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   PB_DEBUG_MSG_END()

   // set defaults as needed
   if(invStr=="") invVStr = "Amesos";
   if(invVStr=="") invVStr = invStr;
   if(invPStr=="") invPStr = invStr;

   //  two inverse factory objects
   RCP<const InverseFactory> invVFact, invPFact;

   // build velocity inverse factory
   invFactoryF_ = invLib.getInverseFactory(invVStr);
   invFactoryS_ = invFactoryF_; // by default these are the same
   if(invVStr!=invPStr) // if different, build pressure inverse factory
      invFactoryS_ = invLib.getInverseFactory(invPStr);

   // set other parameters
   setUseFullLDU(useLDU);
   setRowZeroing(rowZeroing);
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> InvLSCStrategy::getRequestedParameters() const 
{
   Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

   // use the mass matrix
   if(useMass_)
      pl->set<PB::LinearOp>("Velocity Mass Matrix", Teuchos::null,"Velocity mass matrix");

   return Teuchos::null;
}

//! For assiting in construction of the preconditioner
bool InvLSCStrategy::updateRequestedParameters(const Teuchos::ParameterList & pl) 
{
   // set the mass matrix: throw if the strategy is not the right type
   if(useMass_) {
      PB::LinearOp mass = pl.get<PB::LinearOp>("Velocity Mass Matrix");

      // we must have a mass matrix
      if(mass==Teuchos::null) return false;

      // set the mass matrix
      setMassMatrix(mass);
   }

   return true;
}

} // end namespace NS
} // end namespace PB
