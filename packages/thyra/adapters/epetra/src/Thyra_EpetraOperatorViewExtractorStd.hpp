// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_STD_HPP
#define THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_STD_HPP

#include "Thyra_EpetraOperatorViewExtractorBase.hpp"


namespace Thyra {


/** \brief Standard strategy subclass for extracting an
 * <tt>Epetra_Operator</tt> view out of a <tt>Thyra::LinearOpBase<double></tt>
 * object by dynamic casting to the <tt>EpetraLinearOpBase</tt> interface.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraOperatorViewExtractorStd : virtual public EpetraOperatorViewExtractorBase
{
public:

  /** \name Overridden from EpetraOperatorViewExtractorBase. */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpBase<double> &fwdOp ) const;
  /** \brief . */
  void getNonconstEpetraOpView(
    const RCP<LinearOpBase<double> > &fwdOp,
    const Ptr<RCP<Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport,
    const Ptr<double> &epetraOpScalar
    ) const;
  /** \brief . */
  void getEpetraOpView(
    const RCP<const LinearOpBase<double> > &fwdOp,
    const Ptr<RCP<const Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport,
    const Ptr<double> &epetraOpScalar
    ) const;

  //@}

};


} // namespace Thyra


#endif // THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_STD_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

