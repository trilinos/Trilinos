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

#include "Teko_DiagonalPreconditionerOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "EpetraExt_PointToBlockDiagPermute.h"
#include "Epetra_MultiVector.h"

using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;
using Teuchos::RCP;

using Thyra::MultiVectorBase;

namespace Teko {

DiagonalPreconditionerOp::DiagonalPreconditionerOp(Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP, const VectorSpace range, const VectorSpace domain):
  BDP_(BDP),
  range_(range),
  domain_(domain)
{}

void DiagonalPreconditionerOp::implicitApply(const MultiVector & x, MultiVector & y,
					      const double alpha, const double beta) const
{
  // Get the Multivectors into Epetra land
  // NTS: Thyra inexplicably wants maps, even when they are completely unecessary.
  const Epetra_Map & rangemap_=BDP_->OperatorRangeMap();
  const Epetra_Map & domainmap_=BDP_->OperatorDomainMap();

  RCP<const Epetra_MultiVector> x_=Thyra::get_Epetra_MultiVector(domainmap_,x);
  RCP<Epetra_MultiVector> y_=Thyra::get_Epetra_MultiVector(rangemap_,y); 
  TEUCHOS_ASSERT(x_!=Teuchos::null);
  TEUCHOS_ASSERT(y_!=Teuchos::null);


  // y = \alpha M x + \beta y $
  if(beta==0.0){
    BDP_->ApplyInverse(*x_,*y_);
    scale(alpha,y);
  }
  else{
    MultiVector y0=deepcopy(y);
    BDP_->ApplyInverse(*x_,*y_);
    update(alpha,y,beta,y0);  
  }
}

void DiagonalPreconditionerOp::describe(Teuchos::FancyOStream & out_arg,
                                      const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;

  switch(verbLevel) {
     case Teuchos::VERB_DEFAULT:
     case Teuchos::VERB_LOW:
        out_arg << this->description() << std::endl;
        break;
     case Teuchos::VERB_MEDIUM:
     case Teuchos::VERB_HIGH:
     case Teuchos::VERB_EXTREME:
       if(BDP_!=Teuchos::null)
	 BDP_->Print(out_arg);
       break;
     default:
        TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

} // end namespace Teko
