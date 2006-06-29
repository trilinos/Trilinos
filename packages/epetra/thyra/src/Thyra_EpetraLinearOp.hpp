// @HEADER
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

#ifndef THYRA_EPETRA_LINEAR_OP_HPP
#define THYRA_EPETRA_LINEAR_OP_HPP

#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpBase</tt> adapter subclass for
 * <tt>Epetra_Operator</tt> object.
 *
 * This subclass can be used to represent the non-transposed operator or
 * transposed operator defined by an <tt>Epetra_Operator</tt> object.  This
 * class can implement <tt>apply()</tt> using either
 * <tt>Epetra_Operator::Apply()</tt> or
 * <tt>Epetra_Operator::ApplyInverse()</tt>.  In addition, the user can
 * specify whether adjoints are supported or not.
 *
 * <b>Partial Automatic Change Propagation:</b> This class shall maintain no
 * state with respect to the <em>values</em> of the internally stored
 * <tt>Epetra_Operator</tt> object.  Therefore, as long as the domain and
 * range maps do not change, the the <tt>Epetra_Operator</tt> can be changed
 * and this will automatically update <tt>*this</tt> object.  This simplifies
 * some types of update operations.  Since this is a simple concrete class,
 * this is harmless.  However, if the range and domain maps change, then one
 * must call the <tt>this->initialize()</tt> function.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraLinearOp : virtual public EpetraLinearOpBase {
public:

  /** \brief . */
  using SingleScalarEuclideanLinearOpBase<double>::euclideanApply;

  /** @name Public types */
  //@{

  /** \brief . */
  typedef double Scalar;

  //@}

  /** @name Constructors / initializers / accessors */
  //@{

  /** \brief Construct to uninitialized.
   *
   * See the postconditions for <tt>uninitialize()</tt>
   */
  EpetraLinearOp();

  /// Calls <tt>initialize()</tt>.
  EpetraLinearOp(
    const Teuchos::RefCountPtr<Epetra_Operator>                        &op
    ,ETransp                                                           opTrans         = NOTRANS
    ,EApplyEpetraOpAs                                                  applyAs         = EPETRA_OP_APPLY_APPLY
    ,EAdjointEpetraOp                                                  adjointSupport  = EPETRA_OP_ADJOINT_SUPPORTED
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdRange       = Teuchos::null
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdDomain      = Teuchos::null
    );

  /** \brief Initialize
   *
   * @param  op       [in] The <tt>Epetra_Operator</tt> this <tt>*this</tt> will wrap.
   * @param  opTrans  [in] If <tt>opTrans==NOTRANS</tt> then <tt>op</tt> will be viewed as <tt>op</tt>
   *                  and if <tt>opTrans==TRANS</tt> then <tt>op</tt> will be viewed as its transpose
   *                  <tt>op'</tt> for the behavior of <tt>apply()</tt>.
   * @param  applyAs  [in] If <tt>applyAs==APPLY_APPLY</tt> then <tt>op->Apply()</tt> will be used
   *                  and if <tt>applyAs==APPLY_APPLY_INVERSE</tt> then <tt>op->ApplyInverse()</tt>
   *                  is used instead.
   * @param  adjointSupport
   *                  [in] Determines if it is to be assumed that adjoints are supported on the
   *                  underlying <tt>Epetra_Operator</tt> object <tt>op</tt>.  If
   *                  <tt>adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt> then <tt>this->opSupported(TRANS)</tt>
   *                  will return <tt>true</tt>.  If <tt>adjointSupport==EPETRA_OP_ADJOINT_UNSUPPORTED</tt> then
   *                  <tt>this->opSupported(TRANS)</tt> will return <tt>false</tt>.
   * @param  spmdRange
   *                  [in] Smart pointer to the range space for the <tt>Epetra_Operator</tt>.  The default
   *                  value is <tt>Teuchos::null</tt> in which case <tt>*this</tt> will allocate
   *                  a new <tt>SpmdVectorSpace</tt> given range map from <tt>op</tt>.  A client may only bother
   *                  to specify this space if one wants to override the defintion of the scalar product.
   * @param  spmdDomain
   *                  [in] Smart pointer to the domain space for the <tt>Epetra_Operator</tt>.  The default
   *                  value is <tt>Teuchos::null</tt> in which case <tt>*this</tt> will allocate
   *                  a new <tt>DefaultSpmdVectorSpace</tt> given map from <tt>op</tt>.  A client may only bother
   *                  to specify this space if one wants to override the defintion of the scalar product.
   *
   * Preconditions:<ul>
   * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->epetra_op().get() == op.get()</tt>
   * <li> [<tt>spmdRange.get() != NULL</tt>] <tt>this->spmdRange().get() == spmdRange.get()</tt>
   * <li> [<tt>spmdDomain.get() != NULL</tt>] <tt>this->spmdDomain().get() == spmdDomain.get()</tt>
   * <li> [<tt>spmdRange.get() == NULL</tt>] <tt>this->spmdRange().get() != NULL</tt>
   * <li> [<tt>spmdDomain.get() == NULL</tt>] <tt>this->spmdDomain().get() != NULL</tt>
   * <li> <tt>this->opSupported(NOTRANS) == true</tt>
   * <li> <tt>this->opSupported(TRNAS) == adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<Epetra_Operator>                        &op
    ,ETransp                                                           opTrans         = NOTRANS
    ,EApplyEpetraOpAs                                                  applyAs         = EPETRA_OP_APPLY_APPLY
    ,EAdjointEpetraOp                                                  adjointSupport  = EPETRA_OP_ADJOINT_SUPPORTED
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdRange       = Teuchos::null
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdDomain      = Teuchos::null
    );
  
  /** \brief Set to uninitialized and optionally return the current state.
   *
   * Postconditions:<ul>
   * <li> <tt>this->domain().get() == NULL</tt>
   * <li> <tt>this->range().get() == NULL</tt>
   * </ul>
   */
  void uninitialize(
    Teuchos::RefCountPtr<Epetra_Operator>                       *op             = NULL
    ,ETransp                                                    *opTrans        = NULL
    ,EApplyEpetraOpAs                                           *applyAs        = NULL
    ,EAdjointEpetraOp                                           *adjointSupport = NULL
    ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >  *spmdRange       = NULL
    ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >  *spmdDomain      = NULL
    );

  /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the range.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->range().get() != NULL</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->range().get() == NULL</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > spmdRange() const;

  /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the domain.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->domain().get() != NULL</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->domain().get() == NULL</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > spmdDomain() const;

  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> epetra_op();

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Operator> epetra_op() const;

  //@}

  /** @name Overridden from EpetraLinearOpBase */
  //@{

  /** \brief . */
  void getEpetraOpView(
    Teuchos::RefCountPtr<Epetra_Operator>   *epetraOp
    ,ETransp                                *epetraOpTransp
    ,EApplyEpetraOpAs                       *epetraOpApplyAs
    ,EAdjointEpetraOp                       *epetraOpAdjointSupport
    );
  /** \brief . */
  void getEpetraOpView(
    Teuchos::RefCountPtr<const Epetra_Operator>   *epetraOp
    ,ETransp                                      *epetraOpTransp
    ,EApplyEpetraOpAs                             *epetraOpApplyAs
    ,EAdjointEpetraOp                             *epetraOpAdjointSupport
    ) const;

  //@}

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  
  //@}
  
  /** @name Overridden from EuclideanLinearOpBase */
  //@{

  /// Returns <tt>this->spmdRange()</tt>
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;
  /// Returns <tt>this->spmdDomain()</tt>
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
  /** \brief . */
  void euclideanApply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}
  
  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    Teuchos::FancyOStream                &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const;
  
  //@}
  
protected:

  /** \name Allocators for domain and range spaces */
  //@{

  /** \brief Allocate the domain space of the operator.
   *
   * Purpose: In TSFExtended, both EpetraLinearOp and
   * EpetraVectorSpace are extended from the Thyra versions by
   * inheritance, and the TSFExtended operator subclasses expect to
   * work with an extended vector space subclass. Thus, it is
   * necessary for the base operator class to never directly allocate
   * vector space objects, and allocation is delegated to a virtual
   * allocator function.
   */
  virtual Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > 
  allocateDomain(
    const Teuchos::RefCountPtr<Epetra_Operator>  &op 
    ,ETransp                                     op_trans 
    ) const; 
  
  /** \brief Allocate the range space of the operator.
   *
   * Purpose: In TSFExtended, both EpetraLinearOp and
   * EpetraVectorSpace are extended from the Thyra versions by
   * inheritance, and the TSFExtended operator subclasses expect to
   * work with an extended vector space subclass. Thus, it is
   * necessary for the base operator class to never directly allocate
   * vector space objects, and allocation is delegated to a virtual
   * allocator function.
   */
  virtual Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >
  allocateRange( 
    const Teuchos::RefCountPtr<Epetra_Operator>  &op 
    ,ETransp                                     op_trans 
    ) const; 

  //@}

private:

  // ////////////////////////////////////
  // Private data members

  Teuchos::RefCountPtr<Epetra_Operator>                     op_;
  ETransp                                                   opTrans_;
  EApplyEpetraOpAs                                          applyAs_;
  EAdjointEpetraOp                                          adjointSupport_;
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >  range_;
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >  domain_;
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >  sp_range_;
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >  sp_domain_;
  // ////////////////////////////////////
  // Private member functions

  const Epetra_Map& getRangeMap() const;
  const Epetra_Map& getDomainMap() const;

};	// end class EpetraLinearOp

}	// end namespace Thyra

#endif	// THYRA_EPETRA_LINEAR_OP_HPP
