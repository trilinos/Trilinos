// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
