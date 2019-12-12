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

#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_BlockedMappingStrategy.hpp"
#include "Teko_ReorderedMappingStrategy.hpp"

#include "Teuchos_VerboseObject.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"

#include "Teko_Utilities.hpp"

namespace Teko {
namespace Epetra {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

BlockedEpetraOperator::BlockedEpetraOperator(const std::vector<std::vector<int> > & vars,
                                             const Teuchos::RCP<const Epetra_Operator> & content,
                                             const std::string & label) 
      : Teko::Epetra::EpetraOperatorWrapper(), label_(label)
{
   SetContent(vars,content);
}

void BlockedEpetraOperator::SetContent(const std::vector<std::vector<int> > & vars,
                                       const Teuchos::RCP<const Epetra_Operator> & content)
{ 
   fullContent_ = content;
   blockedMapping_ = rcp(new BlockedMappingStrategy(vars,Teuchos::rcpFromRef(fullContent_->OperatorDomainMap()),
                                                         fullContent_->Comm()));
   SetMapStrategy(blockedMapping_);

   // build thyra operator
   BuildBlockedOperator(); 
}

void BlockedEpetraOperator::BuildBlockedOperator()
{
   TEUCHOS_ASSERT(blockedMapping_!=Teuchos::null);

   // get a CRS matrix
   const RCP<const Epetra_CrsMatrix> crsContent 
         = rcp_dynamic_cast<const Epetra_CrsMatrix>(fullContent_);

   // ask the strategy to build the Thyra operator for you
   if(blockedOperator_==Teuchos::null) {
      blockedOperator_ = blockedMapping_->buildBlockedThyraOp(crsContent,label_);
   }
   else {
      const RCP<Thyra::BlockedLinearOpBase<double> > blkOp 
            = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockedOperator_,true);
      blockedMapping_->rebuildBlockedThyraOp(crsContent,blkOp);
   }

   // set whatever is returned
   SetOperator(blockedOperator_,false);

   // reorder if neccessary
   if(reorderManager_!=Teuchos::null) 
      Reorder(*reorderManager_);
}

const Teuchos::RCP<const Epetra_Operator> BlockedEpetraOperator::GetBlock(int i,int j) const
{
   const RCP<const Thyra::BlockedLinearOpBase<double> > blkOp 
         = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(getThyraOp());

   return Thyra::get_Epetra_Operator(*blkOp->getBlock(i,j));
}

/** Use a reorder manager to block this operator as desired.
  * Multiple calls to the function reorder only the underlying object. 
  */
void BlockedEpetraOperator::Reorder(const BlockReorderManager & brm)
{
   reorderManager_ = rcp(new BlockReorderManager(brm));

   // build reordered objects
   RCP<const MappingStrategy> reorderMapping = rcp(new ReorderedMappingStrategy(*reorderManager_,blockedMapping_));
   RCP<const Thyra::BlockedLinearOpBase<double> > blockOp
         = rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(blockedOperator_);

   RCP<const Thyra::LinearOpBase<double> > A = buildReorderedLinearOp(*reorderManager_,blockOp);

   // set them as working values
   SetMapStrategy(reorderMapping);
   SetOperator(A,false);
}

//! Remove any reordering on this object
void BlockedEpetraOperator::RemoveReording()
{
   SetMapStrategy(blockedMapping_);
   SetOperator(blockedOperator_,false);
   reorderManager_ = Teuchos::null;
}

/** Write out this operator to matrix market files
  */
void BlockedEpetraOperator::WriteBlocks(const std::string & prefix) const
{
   RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blockOp
         = rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(blockedOperator_);

   // get size of blocked block operator
   int rows = Teko::blockRowCount(blockOp);

   for(int i=0;i<rows;i++) {
      for(int j=0;j<rows;j++) {
         // build the file name
         std::stringstream ss;
         ss << prefix << "_" << i << j << ".mm";

         // get the row matrix object
         RCP<const Epetra_RowMatrix> mat
               = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(Thyra::get_Epetra_Operator(*blockOp->getBlock(i,j)));

         // write to file
         EpetraExt::RowMatrixToMatrixMarketFile(ss.str().c_str(),*mat);
      }
   }
}

bool BlockedEpetraOperator::testAgainstFullOperator(int count,double tol) const
{
   Epetra_Vector xf(OperatorRangeMap());
   Epetra_Vector xs(OperatorRangeMap());
   Epetra_Vector y(OperatorDomainMap());

   // test operator many times
   bool result = true;
   double diffNorm=0.0,trueNorm=0.0;
   for(int i=0;i<count;i++) {
      xf.PutScalar(0.0);
      xs.PutScalar(0.0);
      y.Random();

      // apply operator
      Apply(y,xs); // xs = A*y
      fullContent_->Apply(y,xf); // xf = A*y

      // compute norms
      xs.Update(-1.0,xf,1.0);
      xs.Norm2(&diffNorm);
      xf.Norm2(&trueNorm);

      // check result
      result &= (diffNorm/trueNorm < tol);
   }
   return result;
}

} // end namespace Epetra
} // end namespace Teko
