#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockLowerTriInverseOp.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

using Teuchos::RCP;

// Declaration of the example preconditioner state
class ExamplePreconditionerState /*@ \label{lne2:begin-decl-state} @*/
   : public PB::BlockPreconditionerState {
public:
   // Default constructor
   ExamplePreconditionerState() 
      : PB::BlockPreconditionerState() {}

   // class members
   PB::LinearOp P_; // store P = A_00+\alpha A_01
   
}; /*@ \label{lne2:end-decl-state} @*/

// Declaration of the preconditioner factory
class ExamplePreconditionerFactory /*@ \label{lne2:begin-decl} @*/
   : public PB::BlockPreconditionerFactory {
public:
   // Constructor
   ExamplePreconditionerFactory(const RCP<const PB::InverseFactory> & inverse,
                                double alpha);

   // Function inherited from PB::BlockPreconditionerFactory
   PB::LinearOp buildPreconditionerOperator(PB::BlockedLinearOp & blo,
                                            PB::BlockPreconditionerState & state) const;

   // Function that returns the correct type of state object
   virtual RCP<PB::BlockPreconditionerState> buildPreconditionerState() const;
    
protected:
   // class members
   RCP<const PB::InverseFactory> inverse_;
   double alpha_;
   std::string invP_str_;
   
};/*@ \label{lne2:end-decl} @*/

// Constructor definition
ExamplePreconditionerFactory /*@ \label{lne2:begin-constructor} @*/
   ::ExamplePreconditionerFactory(const RCP<const PB::InverseFactory> & inverse,
                                 double alpha)
   : inverse_(inverse), alpha_(alpha)
{ 
   // store the string name for retrieving the inverse
   invP_str_ = "invP";
} /*@ \label{lne2:end-constructor} @*/

// Use the factory to build the preconditioner (this is where the work goes)
PB::LinearOp ExamplePreconditionerFactory /*@ \label{lne2:begin-bpo} @*/
   ::buildPreconditionerOperator(PB::BlockedLinearOp & blockOp,
                                 PB::BlockPreconditionerState & state) const
{
   int rows = PB::blockRowCount(blockOp); /*@ \label{lne2:begin-extraction} @*/
   int cols = PB::blockColCount(blockOp);
 
   TEUCHOS_ASSERT(rows==2); // sanity checks
   TEUCHOS_ASSERT(cols==2);
 
   // get the casted version of the state object
   ExamplePreconditionerState & exampleState  /*@ \label{lne2:state-cast} @*/
         = dynamic_cast<ExamplePreconditionerState &>(state);

   // extract subblocks
   const PB::LinearOp A_00 = PB::getBlock(0,0,blockOp);
   const PB::LinearOp A_01 = PB::getBlock(0,1,blockOp);
   const PB::LinearOp A_10 = PB::getBlock(1,0,blockOp);
   const PB::LinearOp A_11 = PB::getBlock(1,1,blockOp); /*@ \label{lne2:end-extraction} @*/

   // get inverse of diag(A11)
   const PB::LinearOp invH = PB::getInvDiagonalOp(A_11); /*@ \label{lne2:invH} @*/

   // build or rebuild inverse P /*@ \label{lne2:invP} @*/
   PB::InverseLinearOp invP = exampleState.getInverse(invP_str_);
   if(invP==Teuchos::null) {
      // build 0,0 block in the preconditioner
      exampleState.P_ = PB::explicitAdd(A_00, PB::scale(alpha_,A_01)); /*@ \label{lne2:P} @*/

      invP = PB::buildInverse(*inverse_,exampleState.P_); // build inverse P
      exampleState.addInverse(invP_str_,invP);  // add inverse operator to state
   } 

   // build lower triangular inverse matrix
   PB::BlockedLinearOp L = PB::zeroBlockedOp(blockOp); /*@ \label{lne2:begin-trisolve} @*/
   PB::setBlock(1,0,L,A_10);
   PB::endBlockFill(L);

   std::vector<PB::LinearOp> invDiag(2); // vector storing inverses /*@ \label{lne2:begin-invdiags} @*/
   invDiag[0] = invP;
   invDiag[1] = invH; /*@ \label{lne2:end-invdiags} @*/

   PB::LinearOp invTildeA = PB::createBlockLowerTriInverseOp(L,invDiag); /*@ \label{lne2:invLower} @*/

   // tell the state object it has been initialized for this operator
   exampleState.setInitialized(true);

   // return fully constructed preconditioner
   return invTildeA; /*@ \label{lne2:end-trisolve} @*/
} /*@ \label{lne2:end-bpo} @*/

// Function that returns the correct type of state object
RCP<PB::BlockPreconditionerState> ExamplePreconditionerFactory /*@ \label{lne2:begin-bps} @*/
   ::buildPreconditionerState() const
{
   // build the state object
   return Teuchos::rcp(new ExamplePreconditionerState());

} /*@ \label{lne2:end-bps} @*/
