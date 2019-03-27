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

#include "Teko_PreconditionerState.hpp"

#include "Thyra_DefaultPreconditioner.hpp"

using namespace Thyra;

namespace Teko {
/////////////////////////////////////////////////////

//! Set parameters from a parameter list and return with default values.
void PreconditionerState::setParameterList(const RCP<Teuchos::ParameterList> & paramList)
{
   paramList_ = paramList;
}

//! Get the parameter list that was set using setParameterList().
RCP<Teuchos::ParameterList> PreconditionerState::getNonconstParameterList()
{
   if(paramList_==Teuchos::null)
      paramList_ = Teuchos::rcp(new Teuchos::ParameterList());

   return paramList_;
}

//! Unset the parameter list that was set using setParameterList(). 
RCP<Teuchos::ParameterList> PreconditionerState::unsetParameterList()
{
   RCP<Teuchos::ParameterList> paramList = paramList_;
   paramList_ = Teuchos::null;
   return paramList;
}

//! Merge internal storage of another PreconditionerState object into this one
void PreconditionerState::merge(const PreconditionerState & ps,int /* position */)
{
   // merge the two linearOps lists
   linearOps_.insert(ps.linearOps_.begin(),ps.linearOps_.end());
 
   // merge two parameter lists
   Teuchos::ParameterList::ConstIterator itr;
   if(ps.paramList_!=Teuchos::null) {
      Teuchos::RCP<Teuchos::ParameterList> paramList = getNonconstParameterList();
      for(itr=ps.paramList_->begin();itr!=ps.paramList_->end();++itr)
         paramList->setEntry(itr->first,itr->second);
   }
}

//! Get the tag for this operator
unsigned int PreconditionerState::getTag() const
{
   return tag_;
}

//! Set the tag for this operator
void PreconditionerState::setTag(unsigned int tag)
{
   tag_ = tag;
}

/////////////////////////////////////////////////////

} // end namespace Teko
