// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
   * \param  fwdOp     [in] The forward linear operator that the view will be extracted from.
   *                   This object may be a wrapped scaled/adjoint operator on input.  On output
   *                   this object is "remembered" in the returned RCP for <tt>epetraOp</tt>.
   * \param  epetraOp  [out] The non-<tt>const</tt> epetra operator view of <tt>*this</tt>.
   * \param  epetraOpTransp
   *                   [out] Determines if the operator is applied
   *                   as its transpose or its non-transpose.
   *                   The Client should use this value and ignore the value in
   *                   <tt>(*epetraOp)->UseTranspose()</tt> since it has been shown to be
   *                   problematic and error prone.
   * \param  epetraOpApplyAs
   *                  [out] Determines if the operator should be applied using
   *                  <tt>(*epetraOp)->Apply(...)</tt> or using <tt>(*epetraOp)->ApplyInverse(...)</tt>.
   * \param  epetraOpAdjointSupport
   *                  [out] Determines if the operator supports transposes or not.
   * \param  epetraOpScalar
   *                  [out] The scalar from the wrapped scaled/adjoint linear operator
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
   * The object accessed from <tt>*epetraOp</tt> is only guaranteed to be
   * valid while the returned <tt>Teuchos::RefCountPtr</tt> object exits.
   * This allows for some very specialized implementations where a
   * <tt>Epetra_Operator</tt> view of <tt>*this</tt> can be acquired and
   * released according to the lifetime of the returned
   * <tt>Teuchos::RefCountPtr</tt> object.
   *
   * The <tt>Epetra_Operator</tt> object may be dynamic casted to more
   * specialized interfaces and therefore modified.  Then, when the last
   * <tt>RefCountPtr</tt> object ancestor returned from this function goes
   * away, then <tt>*this</tt> will be updated to relect the change.
   */
  virtual void getEpetraOpView(
    const Teuchos::RefCountPtr<LinearOpBase<double> >   &fwdOp
    ,Teuchos::RefCountPtr<Epetra_Operator>              *epetraOp
    ,ETransp                                            *epetraOpTransp
    ,EApplyEpetraOpAs                                   *epetraOpApplyAs
    ,EAdjointEpetraOp                                   *epetraOpAdjointSupport
    ,double                                             *epetraOpScalar
    ) const = 0;

  /** \brief Gat a smart pointer to a <tt>const</tt>
   * <tt>Epetra_Operator</tt> view of a <tt>Thyra::LinearOpBase</tt> object
   * and how the object is applied to implement the forward linear operator.
   *
   * \param  fwdOp     [in] The forward linear operator that the view will be extracted from.
   *                   This object may be a wrapped scaled/adjoint operator on input.  On output
   *                   this object is "remembered" in the returned RCP for <tt>epetraOp</tt>.
   * \param  epetraOp  [out] The <tt>const</tt> epetra operator view of <tt>*this</tt>.
   * \param  epetraOpTransp
   *                   [out] Determines if the operator is applied
   *                   as its transpose or its non-transpose.
   *                   The Client should use this value and ignore the value in
   *                   <tt>(*epetraOp)->UseTranspose()</tt> since it has been shown to be
   *                   problematic and error prone.
   * \param  epetraOpApplyAs
   *                  [out] Determines if the operator should be applied using
   *                  <tt>(*epetraOp)->Apply(...)</tt> or using <tt>(*epetraOp)->ApplyInverse(...)</tt>.
   * \param  epetraOpAdjointSupport
   *                  [out] Determines if the operator supports transposes or not.
   * \param  epetraOpScalar
   *                  [out] The scalar from the wrapped scaled/adjoint linear operator
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
   * while the returned <tt>Teuchos::RefCountPtr</tt> object exits.  This
   * allows for some very specialized implementations where a
   * <tt>Epetra_Operator</tt> view of <tt>*this</tt> can be acquired and
   * released according to the lifetime of the returned
   * <tt>Teuchos::RefCountPtr</tt> object.
   *
   * Note that if the client tries to constant cast the returned object and
   * modify it that this returned view is not guaranteed to update
   * <tt>*this</tt>.  If the goal is to modify <tt>*this</tt> then the client
   * should call the non-<tt>const</tt> version of this function.
   */
  virtual void getEpetraOpView(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >   &fwdOp
    ,Teuchos::RefCountPtr<const Epetra_Operator>              *epetraOp
    ,ETransp                                                  *epetraOpTransp
    ,EApplyEpetraOpAs                                         *epetraOpApplyAs
    ,EAdjointEpetraOp                                         *epetraOpAdjointSupport
    ,double                                                   *epetraOpScalar
    ) const = 0;

  //@}

};

} // namespace Thyra

#endif // THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_BASE_HPP
