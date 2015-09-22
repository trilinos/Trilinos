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

#ifndef __Teko_InverseFactoryOperator_hpp__
#define __Teko_InverseFactoryOperator_hpp__

#include "Teuchos_ConstNonconstObjectContainer.hpp"

#include "Teko_InverseFactory.hpp"
#include "Teko_EpetraInverseOpWrapper.hpp"

namespace Teko {
namespace Epetra {

/** \brief A single Epetra wrapper for all operators constructed
  *        from an inverse operator.
  *
  * This class uses the Teko inverse factories to
  * build an Epetra_Operator that behaves like the inverse operatotr.
  * This is done by using the InverseFactory, and letting
  * it build whatever operator is neccessary. Thus the Epetra
  * "layer" is just a single class that handles any generic
  * InverseFactory.
  */
class InverseFactoryOperator : public EpetraInverseOpWrapper {
public:
   /** \brief Constructor that takes the InverseFactory that will
     *        build the operator.
     *
     * Constructor that takes the InverseFactory that will
     * build the operator.
     */
   InverseFactoryOperator(const Teuchos::RCP<const InverseFactory> & bfp); 

   /** \brief Build the underlying data structure for the inverse operator.
     *        
     * Build the underlying data structure for the inverse operator. This
     * permits the manipulation of the state object for an inverse operator.
     *
     * \param[in] clearOld If true any previously constructed
     *                     operator will be wiped out and
     *                     a new one created. If false, anoperator 
     *                     will be created only if the current one is
     *                     empty (i.e. <code>initPreconditioner</code>
     *                     had not been called).
     */
   virtual void initInverse(bool clearOld=false);

   /** \brief Build this inverse operator from an Epetra_Operator 
     * passed in to this object.
     *
     * Build this inverse opeerator from an Epetra_Operator 
     * passed in to this object. If this Epetra_Operator
     * is an EpetraOperatorWrapper object then the block Thyra components
     * are extracted.
     *
     * \param[in] A The Epetra source operator.
     * \param[in] clear If true, than any previous state saved by the operator 
     *                  is discarded.
     */
   virtual void buildInverseOperator(const Teuchos::RCP<const Epetra_Operator> & A,bool clear=true);

   /** \brief Build this inverse operator from an Epetra_Operator 
     * passed in to this object.
     *
     * Build this inverse opeerator from an Epetra_Operator 
     * passed in to this object. If this Epetra_Operator
     * is an EpetraOperatorWrapper object then the block Thyra components
     * are extracted.
     *
     * \param[in] A The Epetra source operator.
     * \param[in] clear If true, than any previous state saved by the operator 
     *                  is discarded.
     */
   virtual void buildInverseOperator(const Teuchos::RCP<Epetra_Operator> & A,bool clear=true);

   /** \brief Rebuild this inverse from an Epetra_Operator passed
     * in this to object. 
     *
     * Rebuild this inverse from an Epetra_Operator passed
     * in this to object.  If <code>buildInverseOperator</code> has not been called
     * the inverse operator will be built instead. Otherwise efforts are taken
     * to only rebuild what is neccessary. Also, that this Epetra_Operator
     * may be an EpetraOperatorWrapper object, so the block Thyra components
     * can be extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     */
   virtual void rebuildInverseOperator(const Teuchos::RCP<const Epetra_Operator> & A);

   /** \brief Rebuild this inverse from an Epetra_Operator passed
     * in this to object. 
     *
     * Rebuild this inverse from an Epetra_Operator passed
     * in this to object.  If <code>buildInverseOperator</code> has not been called
     * the inverse operator will be built instead. Otherwise efforts are taken
     * to only rebuild what is neccessary. Also, that this Epetra_Operator
     * may be an EpetraOperatorWrapper object, so the block Thyra components
     * can be extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     */
   virtual void rebuildInverseOperator(const Teuchos::RCP<Epetra_Operator> & A);

   /** Extract the forward op used by <code>buildInverseOperator</code>
     * or <code>rebuildInverseOperator</code>.
     */
   Teuchos::RCP<const Epetra_Operator> getForwardOp() const { return fwdOp_.getConstObj(); }

   /** Extract the forward op used by <code>buildInverseOperator</code>
     * or <code>rebuildInverseOperator</code>.
     */
   Teuchos::RCP<Epetra_Operator> getNonconstForwardOp() const { return fwdOp_.getNonconstObj(); }

protected:
   Teuchos::RCP<const Thyra::LinearOpBase<double> > extractLinearOp(const Teuchos::RCP<const Epetra_Operator> & A) const;
   Teuchos::RCP<const MappingStrategy> extractMappingStrategy(const Teuchos::RCP<const Epetra_Operator> & A) const;

   InverseFactoryOperator(); 
   InverseFactoryOperator(const InverseFactoryOperator &); 

   Teuchos::RCP<const Teko::InverseFactory> inverseFactory_;
   Teko::ModifiableLinearOp invOperator_;
   bool firstBuildComplete_;

   Teuchos::ConstNonconstObjectContainer<Epetra_Operator> fwdOp_;
   bool setConstFwdOp_;
};

} // end namespace Epetra
} // end namespace Teko

#endif
