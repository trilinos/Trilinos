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

#ifndef __Teko_BlockedTpetraOperator_hpp__
#define __Teko_BlockedTpetraOperator_hpp__

// Tpetra includes
#include "Tpetra_Operator.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

// Teko includes
#include "Teko_BlockedReordering.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraBlockedMappingStrategy.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace TpetraHelpers {

/** \brief Tear about a user specified Tpetra::Operator<ST,LO,GO,NT>  (CrsMatrix)
  *        using a vector of vectors of GIDs for each block. 
  */
class BlockedTpetraOperator : public TpetraOperatorWrapper {
public:

   /** Build a blocked operator based on a vector of vector
     * of global IDs.
     *
     * \param[in] vars Vector of vectors of global ids specifying
     *                 how the operator is to be blocked.
     * \param[in] content Operator to be blocked
     * \param[in] label Label for name the operator
     */
   BlockedTpetraOperator(const std::vector<std::vector<GO> > & vars,
                         const Teuchos::RCP<const Tpetra::Operator<ST,LO,GO,NT> > & content,
                         const std::string & label="<ANYM>");

   /** Build a blocked operator based on a vector of vector
     * of global IDs. This function basically sets up the mapping
     * strategy used by this operator.
     *
     * \param[in] vars Vector of vectors of global ids specifying
     *                 how the operator is to be blocked.
     * \param[in] content Operator to be blocked
     */
   virtual void SetContent(const std::vector<std::vector<GO> > & vars,
                           const Teuchos::RCP<const Tpetra::Operator<ST,LO,GO,NT> > & content);

   /** Force a rebuild of the blocked operator from the stored
     * content operator.
     */
   virtual void RebuildOps()
   { BuildBlockedOperator(); }

   virtual const Teuchos::RCP<const Tpetra::Operator<ST,LO,GO,NT> > GetContent() const
   { return fullContent_; }

   virtual const Teuchos::RCP<const Tpetra::Operator<ST,LO,GO,NT> > GetContent()
   { return fullContent_; }

   const Teuchos::RCP<const Tpetra::Operator<ST,LO,GO,NT> > GetBlock(int i,int j) const;

   /** Use a reorder manager to block this operator as desired.
     * Multiple calls to the function reorder only the underlying object. 
     */
   void Reorder(const BlockReorderManager & brm);

   //! Remove any reordering on this object
   void RemoveReording();

   /** Write out this operator to matrix market files
     */
   virtual void WriteBlocks(const std::string & prefix) const;

   // functions overloading Tpetra::Operator<ST,LO,GO,NT> 
   ////////////////////////////////////////////////

   // destructor
   virtual ~BlockedTpetraOperator() {}

   // attribute set methods
   
   // don't use transpose...ever!
   virtual int SetUseTranspose(bool /* useTranspose */)
   { return -1; }

   virtual int ApplyInverse(const Tpetra::MultiVector<ST,LO,GO,NT> &/* X */, Tpetra::MultiVector<ST,LO,GO,NT> &/* Y */) const
   { TEUCHOS_ASSERT(false); return -1; }

   virtual double NormInf() const
   { TEUCHOS_ASSERT(false); return 0.0; }

   // attribute access functions
   virtual bool UseTranspose() const { return false; }
   virtual bool HasNormInf() const { return false; }
   virtual const Teuchos::Comm<int> & Comm() const { return *fullContent_->getRangeMap()->getComm(); }

   #ifndef Teko_DEBUG_OFF
   //! Helps perform sanity checks
   bool testAgainstFullOperator(int count,ST tol) const;
   #endif

protected:
   // gooey center of this shell
   Teuchos::RCP<const Tpetra::Operator<ST,LO,GO,NT> > fullContent_;
   Teuchos::RCP<TpetraBlockedMappingStrategy> blockedMapping_;
   Teuchos::RCP<Thyra::LinearOpBase<ST> > blockedOperator_;
   Teuchos::RCP<const BlockReorderManager> reorderManager_;

   std::string label_;

   void BuildBlockedOperator();
};

} // end namespace TpetraHelpers
} // end namespace Teko

#endif
