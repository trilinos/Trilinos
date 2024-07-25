// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_LinearOpWithSolveFactoryBase.hpp"


namespace Thyra {


/** \brief Create a DefaultDiagonalLinearOpWithSolve out of a diagonal
 * Epetra_RowMatrix object.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class DiagonalEpetraLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<double> {
public:

  /** @name Overridden from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;

  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<double> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  void uninitializeOp(
    LinearOpWithSolveBase<double> *Op,
    Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOpSrc,
    Teuchos::RCP<const PreconditionerBase<double> > *prec,
    Teuchos::RCP<const LinearOpSourceBase<double> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
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

};


} // namespace Thyra


#endif // THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

