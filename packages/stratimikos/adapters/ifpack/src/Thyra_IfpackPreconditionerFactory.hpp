// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_IFPACK_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_IFPACK_PRECONDITIONER_FACTORY_DECL_HPP

#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Ifpack.h"

namespace Thyra {

/** \brief Concrete preconditioner factory subclass based on Ifpack.
 *
 * ToDo: Finish documentation!
 */
class IfpackPreconditionerFactory : public PreconditionerFactoryBase<double> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  IfpackPreconditionerFactory();
    
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS( EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor );

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;
  /** \brief . */
  bool applySupportsConj(EConj conj) const;
  /** \brief . */
  bool applyTransposeSupportsConj(EConj conj) const;
  /** \brief . */
  Teuchos::RCP<PreconditionerBase<double> > createPrec() const;
  /** \brief . */
  void initializePrec(
    const Teuchos::RCP<const LinearOpSourceBase<double> >    &fwdOpSrc
    ,PreconditionerBase<double>                                      *prec
    ,const ESupportSolveUse                                           supportSolveUse
    ) const;
  /** \brief . */
  void uninitializePrec(
    PreconditionerBase<double>                                *prec
    ,Teuchos::RCP<const LinearOpSourceBase<double> >  *fwdOpSrc
    ,ESupportSolveUse                                         *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // ////////////////////////////////
  // Private data members

  Teuchos::RCP<Teuchos::ParameterList>       paramList_;
  ::Ifpack::EPrecType                                precType_;
  int                                                overlap_;

  // ////////////////////////////////
  // Private member functions

  static void initializeTimers();

};

} // namespace Thyra

#endif // THYRA_IFPACK_PRECONDITIONER_FACTORY_DECL_HPP
