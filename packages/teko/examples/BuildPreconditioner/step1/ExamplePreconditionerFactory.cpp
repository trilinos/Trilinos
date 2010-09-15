/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//  
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//  
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//  
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission. 
//  
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

using Teuchos::RCP;

// Declaration of the preconditioner factory
class ExamplePreconditionerFactory /*@ \label{lne1:begin-decl} @*/
   : public Teko::BlockPreconditionerFactory {
public:
   // Constructor
   ExamplePreconditionerFactory(const RCP<const Teko::InverseFactory> & inverse,
                                double alpha);

   // Function inherited from Teko::BlockPreconditionerFactory
   Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp & blo,
                                            Teko::BlockPreconditionerState & state) const;
    
protected:
   // class members
   RCP<const Teko::InverseFactory> inverse_;
   double alpha_;
   
};/*@ \label{lne1:end-decl} @*/

// Constructor definition
ExamplePreconditionerFactory /*@ \label{lne1:begin-constructor} @*/
   ::ExamplePreconditionerFactory(const RCP<const Teko::InverseFactory> & inverse,
                                 double alpha)
   : inverse_(inverse), alpha_(alpha)
{ } /*@ \label{lne1:end-constructor} @*/

// Use the factory to build the preconditioner (this is where the work goes)
Teko::LinearOp ExamplePreconditionerFactory /*@ \label{lne1:begin-bpo} @*/
   ::buildPreconditionerOperator(Teko::BlockedLinearOp & blockOp,
                                 Teko::BlockPreconditionerState & state) const
{
   int rows = Teko::blockRowCount(blockOp); /*@ \label{lne1:begin-extraction} @*/
   int cols = Teko::blockColCount(blockOp);
 
   TEUCHOS_ASSERT(rows==2); // sanity checks
   TEUCHOS_ASSERT(cols==2);

   // extract subblocks
   const Teko::LinearOp A_00 = Teko::getBlock(0,0,blockOp);
   const Teko::LinearOp A_01 = Teko::getBlock(0,1,blockOp);
   const Teko::LinearOp A_10 = Teko::getBlock(1,0,blockOp);
   const Teko::LinearOp A_11 = Teko::getBlock(1,1,blockOp); /*@ \label{lne1:end-extraction} @*/

   // get inverse of diag(A11)
   const Teko::LinearOp invH = Teko::getInvDiagonalOp(A_11); /*@ \label{lne1:invH} @*/

   // build 0,0 block in the preconditioner
   const Teko::LinearOp P = Teko::explicitAdd(A_00, Teko::scale(alpha_,A_01)); /*@ \label{lne1:P} @*/
   const Teko::LinearOp invP = Teko::buildInverse(*inverse_,P); // build inverse P /*@ \label{lne1:invP} @*/

   // build lower triangular inverse matrix
   Teko::BlockedLinearOp L = Teko::zeroBlockedOp(blockOp); /*@ \label{lne1:begin-trisolve} @*/
   Teko::setBlock(1,0,L,A_10);
   Teko::endBlockFill(L);

   std::vector<Teko::LinearOp> invDiag(2); // vector storing inverses /*@ \label{lne1:begin-invdiags} @*/
   invDiag[0] = invP;
   invDiag[1] = invH; /*@ \label{lne1:end-invdiags} @*/

   Teko::LinearOp invTildeA = Teko::createBlockLowerTriInverseOp(L,invDiag); /*@ \label{lne1:invLower} @*/

   // return fully constructed preconditioner
   return invTildeA; /*@ \label{lne1:end-trisolve} @*/
} /*@ \label{lne1:end-bpo} @*/
