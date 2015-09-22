/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include "Teko_BlockInvDiagonalStrategy.hpp"

namespace Teko {

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const Teuchos::RCP<InverseFactory> & factory)
{
   // only one factory to use!
   invDiagFact_.resize(1,factory);
   defaultInvFact_ = factory;
}

InvFactoryDiagStrategy::InvFactoryDiagStrategy(const std::vector<Teuchos::RCP<InverseFactory> > & factories,
                                               const Teuchos::RCP<InverseFactory> & defaultFact)
{
   invDiagFact_ = factories;

   if(defaultFact==Teuchos::null)
      defaultInvFact_ = invDiagFact_[0];
   else
      defaultInvFact_ = defaultFact;
}

/** returns an (approximate) inverse of the diagonal blocks of A
  * where A is closely related to the original source for invD0 and invD1
  */
void InvFactoryDiagStrategy::getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state,std::vector<LinearOp> & invDiag) const
{ 
   Teko_DEBUG_SCOPE("InvFactoryDiagStrategy::getInvD",10);

   // loop over diagonals, build an inverse operator for each
   int diagCnt = A->productRange()->numBlocks();
   int invCnt = invDiagFact_.size();

   Teko_DEBUG_MSG("# diags = " << diagCnt << ", # inverses = " << invCnt,6);

   const std::string opPrefix = "JacobiDiagOp";
   if(diagCnt<=invCnt) {
      for(int i=0;i<diagCnt;i++) 
         invDiag.push_back(buildInverse(*invDiagFact_[i],getBlock(i,i,A),state,opPrefix,i));
   }
   else {
      for(int i=0;i<invCnt;i++) 
         invDiag.push_back(buildInverse(*invDiagFact_[i],getBlock(i,i,A),state,opPrefix,i));

      for(int i=invCnt;i<diagCnt;i++) 
         invDiag.push_back(buildInverse(*defaultInvFact_,getBlock(i,i,A),state,opPrefix,i));
   }
}

LinearOp InvFactoryDiagStrategy::buildInverse(const InverseFactory & invFact,const LinearOp & matrix,
                                              BlockPreconditionerState & state,
                                              const std::string & opPrefix,int i) const
{
   std::stringstream ss;
   ss << opPrefix << "_" << i;

   ModifiableLinearOp & invOp = state.getModifiableOp(ss.str());
   if(invOp==Teuchos::null) 
      invOp = Teko::buildInverse(invFact,matrix);
   else 
      rebuildInverse(invFact,matrix,invOp);
 
   return invOp;
}

} // end namespace Teko
