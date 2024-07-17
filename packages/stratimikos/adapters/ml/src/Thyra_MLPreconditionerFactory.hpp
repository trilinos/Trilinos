// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_ML_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_ML_PRECONDITIONER_FACTORY_DECL_HPP


#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Thyra {


/** \brief Concrete preconditioner factory subclass based on ML.
 *
 * For information on ML and its available options, please see the ML <A HREF=http://trilinos.sandia.gov/packages/ml/index.html>home page</A> or the ML User <A HREF="http://trilinos.sandia.gov/packages/ml/manuals.html">manuals</A>.
 */
class MLPreconditionerFactory : public PreconditionerFactoryBase<double> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  MLPreconditionerFactory();
    
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS(
    EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor );

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOp ) const;
  /** \brief . */
  bool applySupportsConj(EConj conj) const;
  /** \brief . */
  bool applyTransposeSupportsConj(EConj conj) const;
  /** \brief . */
  Teuchos::RCP<PreconditionerBase<double> > createPrec() const;
  /** \brief . */
  void initializePrec(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOp,
    PreconditionerBase<double> *prec,
    const ESupportSolveUse supportSolveUse
    ) const;
  /** \brief . */
  void uninitializePrec(
    PreconditionerBase<double> *prec
    ,Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOp
    ,ESupportSolveUse *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  // ToDo: Add an override of describe(...) to give more detail!

  //@}

private:

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

};


} // namespace Thyra


#endif // THYRA_ML_PRECONDITIONER_FACTORY_DECL_HPP
