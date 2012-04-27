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

#include "Teko_JacobiPreconditionerFactory.hpp"

using Teuchos::rcp;

namespace Teko {

JacobiPreconditionerFactory::JacobiPreconditionerFactory(const LinearOp & invD0,const LinearOp & invD1)
      : invOpsStrategy_(rcp(new StaticInvDiagStrategy(invD0,invD1)))
{ }

JacobiPreconditionerFactory::JacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy)
         : invOpsStrategy_(strategy)
{ }

/** Build a Jacobi preconditioner factory from a parameter list 
  */
JacobiPreconditionerFactory::JacobiPreconditionerFactory()
{ }

LinearOp JacobiPreconditionerFactory::buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const
{
   int rows = blo->productRange()->numBlocks();
   int cols = blo->productDomain()->numBlocks();
 
   TEUCHOS_ASSERT(rows==cols);

   // get diagonal blocks
   std::vector<LinearOp> invDiag;
   invOpsStrategy_->getInvD(blo,state,invDiag);
   TEUCHOS_ASSERT(rows==(int) invDiag.size());

   // create a blocked linear operator
   BlockedLinearOp precond = createBlockedOp();
   std::stringstream ss;
   ss << "Jacobi Preconditioner ( ";

   // start filling the blocked operator
   beginBlockFill(precond,rows,rows); // this is assuming the matrix is square

   // build blocked diagonal matrix
   for(int i=0;i<rows;i++) {
      ss << " op" << i << " = " << invDiag[i]->description() << ", ";
      precond->setBlock(i,i,invDiag[i]);
   }
   ss << " )";
   
   endBlockFill(precond);
   // done filling the blocked operator

   // precond->setObjectLabel(ss.str());
   precond->setObjectLabel("Jacobi");
   
   return precond; 
}

//! Initialize from a parameter list
void JacobiPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
#if 0
{
   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // get string specifying inverse
   std::string invStr = pl.get<std::string>("Inverse Type");
   if(invStr=="") invStr = "Amesos";

   // based on parameter type build a strategy
   invOpsStrategy_ = rcp(new InvFactoryDiagStrategy(invLib->getInverseFactory(invStr)));
}
#endif 
{
   Teko_DEBUG_SCOPE("JacobiPreconditionerFactory::initializeFromParameterList",10);
   Teko_DEBUG_MSG_BEGIN(9);
      DEBUG_STREAM << "Parameter list: " << std::endl;
      pl.print(DEBUG_STREAM);
   Teko_DEBUG_MSG_END();

   const std::string inverse_type = "Inverse Type";
   std::vector<RCP<InverseFactory> > inverses;

   RCP<const InverseLibrary> invLib = getInverseLibrary();

   // get string specifying default inverse
   std::string invStr ="Amesos"; 
   if(pl.isParameter(inverse_type))
      invStr = pl.get<std::string>(inverse_type);

   Teko_DEBUG_MSG("JacobiPrecFact: Building default inverse \"" << invStr << "\"",5);
   RCP<InverseFactory> defaultInverse = invLib->getInverseFactory(invStr);

   // now check individual solvers
   Teuchos::ParameterList::ConstIterator itr;
   for(itr=pl.begin();itr!=pl.end();++itr) {
      std::string fieldName = itr->first;
      Teko_DEBUG_MSG("JacobiPrecFact: checking fieldName = \"" << fieldName << "\"",9);

      // figure out what the integer is
      if(fieldName.compare(0,inverse_type.length(),inverse_type)==0 && fieldName!=inverse_type) {
         int position = -1;
         std::string inverse,type;

         // figure out position
         std::stringstream ss(fieldName);
         ss >> inverse >> type >> position;

         if(position<=0) {
            Teko_DEBUG_MSG("Jacobi \"Inverse Type\" must be a (strictly) positive integer",1);
         }

         // inserting inverse factory into vector
         std::string invStr = pl.get<std::string>(fieldName);
         Teko_DEBUG_MSG("JacobiPrecFact: Building inverse " << position << " \"" << invStr << "\"",5);
         if(position>(int) inverses.size()) {
            inverses.resize(position,defaultInverse);
            inverses[position-1] = invLib->getInverseFactory(invStr);
         }
         else
            inverses[position-1] = invLib->getInverseFactory(invStr);
      }
   }

   // use default inverse
   if(inverses.size()==0) 
      inverses.push_back(defaultInverse);

   // based on parameter type build a strategy
   invOpsStrategy_ = rcp(new InvFactoryDiagStrategy(inverses,defaultInverse));
}

} // end namspace Teko
