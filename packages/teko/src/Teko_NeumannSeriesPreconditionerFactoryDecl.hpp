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

#ifndef __Teko_NeumannSeriesPreconditionerFactoryDecl_hpp__
#define __Teko_NeumannSeriesPreconditionerFactoryDecl_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
#include "Thyra_PreconditionerFactoryBase.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

namespace Teko {

using Teuchos::RCP;

template <typename ScalarT>
class NeumannSeriesPreconditionerFactory
   : public virtual Thyra::PreconditionerFactoryBase<ScalarT> {
public:

   NeumannSeriesPreconditionerFactory();

   //! is this operator compatiable with the preconditioner factory?
   bool isCompatible(const Thyra::LinearOpSourceBase<ScalarT> &fwdOpSrc) const;

   //! create an instance of the preconditioner
   RCP<Thyra::PreconditionerBase<ScalarT> > createPrec() const;

   /** \brief initialize a newly created preconditioner object
     *
     * Initialize a newly created preconditioner object. For use with
     * nonlinear solvers.
     *
     * \param[in] fwdOpSrc Forward operator to be preconditioned
     * \param[in] solnVec Vector associated with this linear operator.
     * \param[in,out] precOp Return location for the preconditioner
     * \param[in] supportSolveUse Thyra information (?)
     */
   void initializePrec(const RCP<const Thyra::LinearOpSourceBase<ScalarT> > & fwdOpSrc,
                       const RCP<const Thyra::MultiVectorBase<ScalarT> > & solnVec,
                       Thyra::PreconditionerBase<ScalarT> * precOp,
                       const Thyra::ESupportSolveUse supportSolveUse) const;

   /** \brief initialize a newly created preconditioner object
     *
     * Initialize a newly created preconditioner object. 
     *
     * \param[in] fwdOpSrc Forward operator to be preconditioned
     * \param[in,out] precOp Return location for the preconditioner
     * \param[in] supportSolveUse Thyra information (?)
     */
   void initializePrec(const RCP<const Thyra::LinearOpSourceBase<ScalarT> > & fwdOpSrc,
                       Thyra::PreconditionerBase<ScalarT> * precOp,
                       const Thyra::ESupportSolveUse supportSolveUse) const;

   //! wipe clean a already initialized preconditioner object
   void uninitializePrec(Thyra::PreconditionerBase<ScalarT> * prec, 
                         RCP<const Thyra::LinearOpSourceBase<ScalarT> > * fwdOpSrc,
                         Thyra::ESupportSolveUse *supportSolveUse) const;

   /** @name Overridden from Teuchos::ParameterListAcceptor */
   //@{
 
   //! \brief Set parameters from a parameter list
   void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

   //! \brief Get the parameter list that was set using setParameterList().
   Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

   //! \brief Unset the parameter list that was set using setParameterList(). 
   Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

   //! \brief Get the parameter list that was set using setParameterList().
   Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

   /** \brief Get the valid parameters */
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
   //@}
 
   /** \name Public functions overridden from Describable. */
   //@{
 
   /** \brief . */
   std::string description() const;
 
   //@}
  
protected:
   //! for ParameterListAcceptor
   Teuchos::RCP<Teuchos::ParameterList> paramList_;

   int numberOfTerms_;
   Teko::DiagonalType scalingType_;
};

} // end namespace Teko

#endif
