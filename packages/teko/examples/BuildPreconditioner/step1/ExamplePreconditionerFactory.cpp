#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockLowerTriInverseOp.hpp"
#include "PB_BlockUpperTriInverseOp.hpp"

using Teuchos::RCP;

// Declaration of the preconditioner factory
class ExamplePreconditionerFactory /*@ \label{lne1:begin-decl} @*/
   : public PB::BlockPreconditionerFactory {
public:
   // Constructor
   ExamplePreconditionerFactory(const RCP<const PB::InverseFactory> & inverse,
                                double alpha);

   // Function inherited from PB::BlockPreconditionerFactory
   PB::LinearOp buildPreconditionerOperator(PB::BlockedLinearOp & blo,
                                            PB::BlockPreconditionerState & state) const;
    
protected:
   // class members
   RCP<const PB::InverseFactory> inverse_;
   double alpha_;
   
};/*@ \label{lne1:end-decl} @*/

// Constructor definition
ExamplePreconditionerFactory /*@ \label{lne1:begin-constructor} @*/
   ::ExamplePreconditionerFactory(const RCP<const PB::InverseFactory> & inverse,
                                 double alpha)
   : inverse_(inverse), alpha_(alpha)
{ } /*@ \label{lne1:end-constructor} @*/

// Use the factory to build the preconditioner (this is where the work goes)
PB::LinearOp ExamplePreconditionerFactory /*@ \label{lne1:begin-bpo} @*/
   ::buildPreconditionerOperator(PB::BlockedLinearOp & blockOp,
                                 PB::BlockPreconditionerState & state) const
{
   int rows = PB::blockRowCount(blockOp); /*@ \label{lne1:begin-extraction} @*/
   int cols = PB::blockColCount(blockOp);
 
   TEUCHOS_ASSERT(rows==2); // sanity checks
   TEUCHOS_ASSERT(cols==2);

   // extract subblocks
   const PB::LinearOp A_00 = PB::getBlock(0,0,blockOp);
   const PB::LinearOp A_01 = PB::getBlock(0,1,blockOp);
   const PB::LinearOp A_10 = PB::getBlock(1,0,blockOp);
   const PB::LinearOp A_11 = PB::getBlock(1,1,blockOp); /*@ \label{lne1:end-extraction} @*/

   // get inverse of diag(A11)
   const PB::LinearOp invH = PB::getInvDiagonalOp(A_11); /*@ \label{lne1:invH} @*/

   // build 0,0 block in the preconditioner
   const PB::LinearOp P = PB::explicitAdd(A_00, PB::scale(alpha_,A_01)); /*@ \label{lne1:P} @*/
   const PB::LinearOp invP = PB::buildInverse(*inverse_,P); // build inverse P /*@ \label{lne1:invP} @*/

   // build lower triangular inverse matrix
   PB::BlockedLinearOp L = PB::zeroBlockedOp(blockOp); /*@ \label{lne1:begin-trisolve} @*/
   PB::setBlock(1,0,L,A_10);
   PB::endBlockFill(L);

   std::vector<PB::LinearOp> invDiag(2); // vector storing inverses /*@ \label{lne1:begin-invdiags} @*/
   invDiag[0] = invP;
   invDiag[1] = invH; /*@ \label{lne1:end-invdiags} @*/

   PB::LinearOp invTildeA = PB::createBlockLowerTriInverseOp(L,invDiag); /*@ \label{lne1:invLower} @*/

   // return fully constructed preconditioner
   return invTildeA; /*@ \label{lne1:end-trisolve} @*/
} /*@ \label{lne1:end-bpo} @*/
