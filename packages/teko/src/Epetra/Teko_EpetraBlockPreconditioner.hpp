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

#ifndef __Teko_EpetraBlockPreconditioner_hpp__
#define __Teko_EpetraBlockPreconditioner_hpp__

// #include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_EpetraInverseOpWrapper.hpp"

namespace Teko {
namespace Epetra {

/** \brief A single Epetra wrapper for all the BlockPreconditioners.
  *
  * This class uses the Thyra based preconditioner factories to
  * build an Epetra_Operator that behaves like a preconditioner.
  * This is done by using the BlockPreconditionerFactory, and letting
  * it build whatever preconditioner is neccessary. Thus the Epetra
  * "layer" is just a single class that handles any generic
  * BlockPreconditionerFactory.
  */
class EpetraBlockPreconditioner : public EpetraInverseOpWrapper {
public:
   /** \brief Constructor that takes the BlockPreconditionerFactory that will
     *        build the preconditioner.
     *
     * Constructor that takes the BlockPreconditionerFactory that will
     * build the preconditioner.
     */
   EpetraBlockPreconditioner(const Teuchos::RCP<const PreconditionerFactory> & bfp); 

   /** \brief Build the underlying data structure for the preconditioner.
     *        
     * Build the underlying data structure for the preconditioner. This
     * permits the manipulation of the state object for a preconditioner.
     * and is useful in that case some extra data needs to fill the
     * preconditioner state.
     *
     * \param[in] clearOld If true any previously constructed
     *                     preconditioner will be wiped out and
     *                     a new one created. If false, a preconditioner
     *                     will be create only if the current one is
     *                     empty (i.e. <code>initPreconditioner</code>
     *                     had not been called).
     */
   virtual void initPreconditioner(bool clearOld=false);

   /** \brief Build this preconditioner from an Epetra_Operator 
     * passed in to this object.
     *
     * Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assumed that this Epetra_Operator
     * will be a EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     * \param[in] clear If true, than any previous state saved by the preconditioner
     *                  is discarded.
     */
   virtual void buildPreconditioner(const Epetra_Operator & A,bool clear=true);

   /** \brief Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assumed that this Epetra_Operator
     *
     * Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assumed that this Epetra_Operator
     * will be a EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     * \param[in] mv A vector that was used to build the source operator.
     * \param[in] clear If true, than any previous state saved by the preconditioner
     *                  is discarded.
     */
   virtual void buildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & mv,bool clear=true);

   /** \brief Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object. 
     *
     * Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object.  If <code>buildPreconditioner</code> has not been called
     * the preconditioner will be built instead. Otherwise efforts are taken
     * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
     * will be an EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     */
   virtual void rebuildPreconditioner(const Epetra_Operator & A);

   /** \brief Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object. 
     *
     * Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object.  If <code>buildPreconditioner</code> has not been called
     * the preconditioner will be built instead. Otherwise efforts are taken
     * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
     * will be an EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     * \param[in] mv A vector that was used to build the source operator.
     */
   virtual void rebuildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & mv);

   /** Try to get a <code>Teko::PreconditionerState</code> object. This method
     * attempts to cast its internal representation of a preconditioner 
     * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a 
     * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
     *
     * \returns Get the state object associated with this preconditioner.
     *          If it doesn't exist for this type of preconditioner factory
     *          this method returns null.
     */
   virtual Teuchos::RCP<PreconditionerState> getPreconditionerState();

   /** Try to get a <code>Teko::PreconditionerState</code> object. This method
     * attempts to cast its internal representation of a preconditioner 
     * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a 
     * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
     *
     * \returns Get the state object associated with this preconditioner.
     *          If it doesn't exist for this type of preconditioner factory
     *          this method returns null.
     */
   virtual Teuchos::RCP<const PreconditionerState> getPreconditionerState() const;

protected:
   Teuchos::RCP<const Thyra::LinearOpBase<double> > extractLinearOp(const Epetra_Operator & A) const;
   Teuchos::RCP<const MappingStrategy> extractMappingStrategy(const Epetra_Operator & A) const;

   EpetraBlockPreconditioner(); 
   EpetraBlockPreconditioner(const EpetraBlockPreconditioner &); 

   // Teuchos::RCP<const BlockPreconditionerFactory> preconFactory_;
   Teuchos::RCP<const PreconditionerFactory> preconFactory_;
   Teuchos::RCP<Thyra::PreconditionerBase<double> > preconObj_;
   bool firstBuildComplete_;
};

} // end namespace Epetra
} // end namespace Teko

#endif
