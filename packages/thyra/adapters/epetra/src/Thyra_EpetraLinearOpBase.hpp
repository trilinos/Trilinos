// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRA_LINEAR_OP_BASE_HPP
#define THYRA_EPETRA_LINEAR_OP_BASE_HPP

#include "Thyra_EpetraTypes.hpp"
#include "Teuchos_Describable.hpp"


class Epetra_Operator;


namespace Thyra {


/** \brief Abstract base class for all <tt>LinearOpBase</tt> objects that can
 * return an <tt>Epetra_Operator</tt> view of themselves and details about how
 * to apply the view.
 *
 * This interface defines a key interoperability interface that allows a
 * general client to extract an <tt>Epetra_Operator</tt> view of a linear
 * operator object.  In some cases, <tt>*this</tt> linear operator can be
 * modifed through the <tt>Epetra_Operator</tt> view and in some cases, it can
 * not.
 *
 * ToDo: Finish documentatation!
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraLinearOpBase : virtual public Teuchos::Describable {
public:

  /** \name Pure virtual functions that must be overridden in subclasses. */
  //@{

  /** \brief Return a smart pointer to a non-<tt>const</tt>
   * <tt>Epetra_Operator</tt> view of this object and how the object is
   * applied to implement <tt>*this</tt> linear operator.
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
   * <b>Preconditions:</b></ul>
   * <li><tt>epetraOp!=NULL</tt>
   * <li><tt>epetraOpOpTransp!=NULL</tt>
   * <li><tt>epetraOpApplyAs!=NULL</tt>
   * <li><tt>epetraOpAdjointSupport!=NULL</tt>
   * </ul>
   *
   * <b>Posconditions:</b></ul>
   * <li><tt>epetraOp->get() != NULL</tt>
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
   *
   * <b>Warning!</b> The client can not assume that the view in
   * <tt>*(*epetraOp)</tt> will be valid past the lifetime of <tt>*this</tt>
   * object which is providing the view!  The client must take special care in
   * the case!
   */
  virtual void getNonconstEpetraOpView(
    const Ptr<RCP<Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    ) = 0;

  /** \brief Return a smart pointer to a <tt>const</tt>
   * <tt>Epetra_Operator</tt> view of this object and how the object is
   * applied to implement <tt>*this</tt> linear operator.
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
   * <b>Preconditions:</b></ul>
   * <li><tt>epetraOp!=NULL</tt>
   * <li><tt>epetraOpOpTransp!=NULL</tt>
   * <li><tt>epetraOpApplyAs!=NULL</tt>
   * <li><tt>epetraOpAdjointSupport!=NULL</tt>
   * </ul>
   *
   * <b>Posconditions:</b></ul>
   * <li><tt>epetraOp->get() != NULL</tt>
   * </ul>
   *
   * The object accessed from <tt>*return</tt> is only guaranteed to be valid
   * while the returned <tt>Teuchos::RCP</tt> object exits.  This allows for
   * some very specialized implementations where a <tt>Epetra_Operator</tt>
   * view of <tt>*this</tt> can be acquired and released according to the
   * lifetime of the returned <tt>Teuchos::RCP</tt> object.
   *
   * Note that if the client tries to constant cast the returned object and
   * modify it that this returned view is not guaranteed to update
   * <tt>*this</tt>.  If the goal is to modify <tt>*this</tt> then the client
   * should call the non-<tt>const</tt> version of this function.
   *
   * <b>Warning!</b> The client can not assume that the view in
   * <tt>*(*epetraOp)</tt> will be valid past the lifetime of <tt>*this</tt>
   * object which is providing the view!  The client must take special care in
   * the case!
   */
  virtual void getEpetraOpView(
    const Ptr<RCP<const Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    ) const = 0;

  //@}

};	// end class EpetraLinearOpBase


}	// end namespace Thyra


#endif	// THYRA_EPETRA_LINEAR_OP_BASE_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

