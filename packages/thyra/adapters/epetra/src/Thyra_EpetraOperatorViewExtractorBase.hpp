// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_BASE_HPP
#define THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_BASE_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraTypes.hpp"


namespace Thyra {


/** \brief Strategy interface for extracting an <tt>Epetra_Operator</tt> view
 * out of a <tt>Thyra::LinearOpBase<double></tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraOperatorViewExtractorBase : virtual public Teuchos::Describable
{
public:

  /** \name Pure virtual functions that must be overridden in subclasses. */
  //@{

  /** \brief Check that a <tt>LinearOpBase</tt> object is compatible with
   * <tt>*this</tt> factory object.
   */
  virtual bool isCompatible( const LinearOpBase<double> &fwdOp ) const = 0;

  /** \brief Gat a smart pointer to a non-<tt>const</tt>
   * <tt>Epetra_Operator</tt> view of a <tt>Thyra::LinearOpBase</tt> object
   * and how the object is applied to implement the forward linear operator.
   *
   * \param fwdOp [in] The forward linear operator that the view will be
   * extracted from.  This object may be a wrapped scaled/adjoint operator on
   * input.  On output this object is "remembered" in the returned RCP for
   * <tt>epetraOp</tt>.
   *
   * \param epetraOp [out] The non-<tt>const</tt> epetra operator view of
   * <tt>*this</tt>.
   *
   * \param epetraOpTransp [out] Determines if the operator is applied as its
   * transpose or its non-transpose.  The Client should use this value and
   * ignore the value in <tt>(*epetraOp)->UseTranspose()</tt> since it has
   * been shown to be problematic and error prone.
   *
   * \param epetraOpApplyAs [out] Determines if the operator should be applied
   * using <tt>(*epetraOp)->Apply(...)</tt> or using
   * <tt>(*epetraOp)->ApplyInverse(...)</tt>.
   *
   * \param epetraOpAdjointSupport [out] Determines if the operator supports
   * transposes or not.
   *
   * \param epetraOpScalar [out] The scalar from the wrapped scaled/adjoint
   * linear operator
   *
   * <b>Preconditions:</b></ul>
   * <li><tt>epetraOp!=NULL</tt>
   * <li><tt>epetraOpOpTransp!=NULL</tt>
   * <li><tt>epetraOpApplyAs!=NULL</tt>
   * <li><tt>epetraOpAdjointSupport!=NULL</tt>
   * </ul>
   *
   * <b>Posconditions:</b></ul>
   * <li><tt>epetraOp->get() != NULL</tt>
   * <li><tt>fwdOp.count()</tt> is greater on output than on input and hense
   *   <tt>fwdOp</tt> is "remembered"
   * </ul>
   *
   * The object accessed from <tt>*epetraOp</tt> is only guaranteed to be
   * valid while the returned <tt>Teuchos::RCP</tt> object exits.
   * This allows for some very specialized implementations where a
   * <tt>Epetra_Operator</tt> view of <tt>*this</tt> can be acquired and
   * released according to the lifetime of the returned
   * <tt>Teuchos::RCP</tt> object.
   *
   * The <tt>Epetra_Operator</tt> object may be dynamic casted to more
   * specialized interfaces and therefore modified.  Then, when the last
   * <tt>RCP</tt> object ancestor returned from this function goes
   * away, then <tt>*this</tt> will be updated to relect the change.
   */
  virtual void getNonconstEpetraOpView(
    const RCP<LinearOpBase<double> > &fwdOp,
    const Ptr<RCP<Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport,
    const Ptr<double> &epetraOpScalar
    ) const = 0;

  /** \brief Gat a smart pointer to a <tt>const</tt>
   * <tt>Epetra_Operator</tt> view of a <tt>Thyra::LinearOpBase</tt> object
   * and how the object is applied to implement the forward linear operator.
   *
   * \param fwdOp [in] The forward linear operator that the view will be
   * extracted from.  This object may be a wrapped scaled/adjoint operator on
   * input.  On output this object is "remembered" in the returned RCP for
   * <tt>epetraOp</tt>.
   *
   * \param epetraOp [out] The <tt>const</tt> epetra operator view of
   * <tt>*this</tt>.
   *
   * \param epetraOpTransp [out] Determines if the operator is applied as its
   * transpose or its non-transpose.  The Client should use this value and
   * ignore the value in <tt>(*epetraOp)->UseTranspose()</tt> since it has
   * been shown to be problematic and error prone.
   *
   * \param epetraOpApplyAs [out] Determines if the operator should be applied
   * using <tt>(*epetraOp)->Apply(...)</tt> or using
   * <tt>(*epetraOp)->ApplyInverse(...)</tt>.
   *
   * \param epetraOpAdjointSupport [out] Determines if the operator supports
   * transposes or not.
   *
   * \param epetraOpScalar [out] The scalar from the wrapped scaled/adjoint
   * linear operator.
   *
   * <b>Preconditions:</b></ul>
   * <li><tt>epetraOp!=NULL</tt>
   * <li><tt>epetraOpOpTransp!=NULL</tt>
   * <li><tt>epetraOpApplyAs!=NULL</tt>
   * <li><tt>epetraOpAdjointSupport!=NULL</tt>
   * </ul>
   *
   * <b>Posconditions:</b></ul>
   * <li><tt>epetraOp->get() != NULL</tt>
   * <li><tt>fwdOp.count()</tt> is greater on output than on input and hense <tt>fwdOp</tt> is "remembered"
   * </ul>
   *
   * The object accessed from <tt>*return</tt> is only guaranteed to be valid
   * while the returned <tt>Teuchos::RCP</tt> object exits.  This
   * allows for some very specialized implementations where a
   * <tt>Epetra_Operator</tt> view of <tt>*this</tt> can be acquired and
   * released according to the lifetime of the returned
   * <tt>Teuchos::RCP</tt> object.
   *
   * Note that if the client tries to constant cast the returned object and
   * modify it that this returned view is not guaranteed to update
   * <tt>*this</tt>.  If the goal is to modify <tt>*this</tt> then the client
   * should call the non-<tt>const</tt> version of this function.
   */
  virtual void getEpetraOpView(
    const RCP<const LinearOpBase<double> > &fwdOp,
    const Ptr<RCP<const Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport,
    const Ptr<double> &epetraOpScalar
    ) const = 0;

  //@}

};


} // namespace Thyra


#endif // THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_BASE_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

