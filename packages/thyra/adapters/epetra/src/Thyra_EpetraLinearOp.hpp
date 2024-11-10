// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRA_LINEAR_OP_HPP
#define THYRA_EPETRA_LINEAR_OP_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

#include "Epetra_RowMatrix.h"


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
class EpetraLinearOp
  : virtual public LinearOpBase<double>,
    virtual public ScaledLinearOpBase<double>,
    virtual public RowStatLinearOpBase<double>,
    virtual public EpetraLinearOpBase
{
public:

  /** \name Constructors / initializers / accessors */
  //@{

  /** \brief Construct to uninitialized.
   *
   * See the postconditions for <tt>uninitialize()</tt>
   */
  EpetraLinearOp();

  /** \brief Fully initialize.
   *
   * \param op [in] The <tt>Epetra_Operator</tt> this <tt>*this</tt> will
   * wrap.
   *
   * \param opTrans [in] If <tt>opTrans==NOTRANS</tt> then <tt>op</tt> will be
   * viewed as <tt>op</tt> and if <tt>opTrans==TRANS</tt> then <tt>op</tt>
   * will be viewed as its transpose <tt>op'</tt> for the behavior of
   * <tt>apply()</tt>.
   *
   * \param applyAs [in] If <tt>applyAs==APPLY_APPLY</tt> then
   * <tt>op->Apply()</tt> will be used and if
   * <tt>applyAs==APPLY_APPLY_INVERSE</tt> then <tt>op->ApplyInverse()</tt> is
   * used instead.
   *
   * \param adjointSupport [in] Determines if it is to be assumed that
   * adjoints are supported on the underlying <tt>Epetra_Operator</tt> object
   * <tt>op</tt>.  If <tt>adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt>
   * then <tt>this->opSupported(TRANS)</tt> will return <tt>true</tt>.  If
   * <tt>adjointSupport==EPETRA_OP_ADJOINT_UNSUPPORTED</tt> then
   * <tt>this->opSupported(TRANS)</tt> will return <tt>false</tt>.
   *
   * \param range [in] Smart pointer to the range space for the
   * <tt>Epetra_Operator</tt>.  The default value is <tt>Teuchos::null</tt> in
   * which case <tt>*this</tt> will allocate a new <tt>SpmdVectorSpace</tt>
   * given range map from <tt>op</tt>.  A client may only bother to specify
   * this space if one wants to override the defintion of the scalar product.
   *
   * \param domain [in] Smart pointer to the domain space for the
   * <tt>Epetra_Operator</tt>.  The default value is <tt>Teuchos::null</tt> in
   * which case <tt>*this</tt> will allocate a new
   * <tt>DefaultSpmdVectorSpace</tt> given map from <tt>op</tt>.  A client may
   * only bother to specify this space if one wants to override the defintion
   * of the scalar product.
   *
   * Preconditions:<ul>
   * <li> <tt>!is_null(op)</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->epetra_op().get() == op.get()</tt>
   * <li> [<tt>range.get() != NULL</tt>] <tt>this->range().get() == range.get()</tt>
   * <li> [<tt>domain.get() != NULL</tt>] <tt>this->domain().get() == domain.get()</tt>
   * <li> [<tt>range.get() == NULL</tt>] <tt>this->range().get() != NULL</tt>
   * <li> [<tt>domain.get() == NULL</tt>] <tt>this->domain().get() != NULL</tt>
   * <li> <tt>this->opSupported(NOTRANS) == true</tt>
   * <li> <tt>this->opSupported(TRNAS) == adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt>
   * </ul>
   *
   * After this function is called, <tt>this</tt> will be fully initialized
   * and ready to go.
   */
  void initialize(
    const RCP<Epetra_Operator> &op,
    EOpTransp opTrans = NOTRANS,
    EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
    EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED,
    const RCP<const VectorSpaceBase<double> > &range = Teuchos::null,
    const RCP<const VectorSpaceBase<double> > &domain = Teuchos::null
    );

  /** \brief Partially initialize.
   *
   * \param range [in] Smart pointer to the range space for the
   * <tt>Epetra_Operator</tt>.
   *
   * \param domain [in] Smart pointer to the domain space for the
   * <tt>Epetra_Operator</tt>.
   *
   * \param op [in] The <tt>Epetra_Operator</tt> this <tt>*this</tt> will
   * wrap.  This object is assumed to not be fully unitialized.
   *
   * \param opTrans [in] If <tt>opTrans==NOTRANS</tt> then <tt>op</tt> will be
   * viewed as <tt>op</tt> and if <tt>opTrans==TRANS</tt> then <tt>op</tt>
   * will be viewed as its transpose <tt>op'</tt> for the behavior of
   * <tt>apply()</tt>.
   *
   * \param applyAs [in] If <tt>applyAs==APPLY_APPLY</tt> then
   * <tt>op->Apply()</tt> will be used and if
   * <tt>applyAs==APPLY_APPLY_INVERSE</tt> then <tt>op->ApplyInverse()</tt> is
   * used instead.
   *
   * \param adjointSupport [in] Determines if it is to be assumed that
   * adjoints are supported on the underlying <tt>Epetra_Operator</tt> object
   * <tt>op</tt>.  If <tt>adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt>
   * then <tt>this->opSupported(TRANS)</tt> will return <tt>true</tt>.  If
   * <tt>adjointSupport==EPETRA_OP_ADJOINT_UNSUPPORTED</tt> then
   * <tt>this->opSupported(TRANS)</tt> will return <tt>false</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>!is_null(range)</tt>
   * <li> <tt>!is_null(domain)</tt>
   * <li> <tt>!is_null(op)</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->epetra_op().get() == op.get()</tt>
   * <li> [<tt>range.get() != NULL</tt>] <tt>this->range().get() == range.get()</tt>
   * <li> [<tt>domain.get() != NULL</tt>] <tt>this->domain().get() == domain.get()</tt>
   * <li> [<tt>range.get() == NULL</tt>] <tt>this->range().get() != NULL</tt>
   * <li> [<tt>domain.get() == NULL</tt>] <tt>this->domain().get() != NULL</tt>
   * <li> <tt>this->opSupported(NOTRANS) == true</tt>
   * <li> <tt>this->opSupported(TRNAS) == adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt>
   * </ul>
   *
   * After this function is called, only the range and domain spaces will be
   * supported and this must be followed up by a call to
   * <tt>setFullyInitialized()</tt>.
   */
  void partiallyInitialize(
    const RCP<const VectorSpaceBase<double> > &range,
    const RCP<const VectorSpaceBase<double> > &domain,
    const RCP<Epetra_Operator> &op,
    EOpTransp opTrans = NOTRANS,
    EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
    EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED
    );

  /** \brief Set to fully initialized.
   *
   * In debug mode, asserts will be performed to ensure that everything
   * matches up as it should.
   *
   * The functions <tt>initialize()</tt> ore <tt>partiallyInitialize()</tt>
   * must have been called prior to calling this function.
   */
  void setFullyInitialized(bool isFullyInitialized = true);
  
  /** \brief Set to uninitialized and optionally return the current state.
   *
   * Postconditions:<ul>
   * <li> <tt>this->domain().get() == NULL</tt>
   * <li> <tt>this->range().get() == NULL</tt>
   * </ul>
   */
  void uninitialize(
    RCP<Epetra_Operator> *op= NULL,
    EOpTransp *opTrans = NULL,
    EApplyEpetraOpAs *applyAs = NULL,
    EAdjointEpetraOp *adjointSupport = NULL,
    RCP<const VectorSpaceBase<double> > *range = NULL,
    RCP<const VectorSpaceBase<double> > *domain = NULL
    );

  /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the
   * range.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->range().get() != NULL</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->range().get() == NULL</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  RCP<const SpmdVectorSpaceBase<double> > spmdRange() const;

  /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the
   * domain.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->domain().get() != NULL</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->domain().get() == NULL</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  RCP<const SpmdVectorSpaceBase<double> > spmdDomain() const;

  /** \brief . */
  RCP<Epetra_Operator> epetra_op();

  /** \brief . */
  RCP<const Epetra_Operator> epetra_op() const;

  //@}

  /** \name Overridden from EpetraLinearOpBase */
  //@{

  /** \brief . */
  void getNonconstEpetraOpView(
    const Ptr<RCP<Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    );
  /** \brief . */
  void getEpetraOpView(
    const Ptr<RCP<const Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    ) const;

  //@}

  /** \name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<double> > range() const;

  /** \brief . */
  RCP<const VectorSpaceBase<double> > domain() const;

  /** \brief . */
  RCP<const LinearOpBase<double> > clone() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;
  
  //@}
  
protected:

  /** \name Protected member functions overridden from LinearOpBase. */
  //@{

  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &X,
    const Ptr<MultiVectorBase<double> > &Y,
    const double alpha,
    const double beta
    ) const;

  //@}

  /** \name Protected member functions overridden from ScaledLinearOpBase. */
  //@{

  /** \brief . */
  virtual bool supportsScaleLeftImpl() const;

  /** \brief . */
  virtual bool supportsScaleRightImpl() const;

  /** \brief . */
  virtual void scaleLeftImpl(const VectorBase<double> &row_scaling);

  /** \brief . */
  virtual void scaleRightImpl(const VectorBase<double> &col_scaling);
  
  //@}

  /** \name Protected member functions overridden from RowStatLinearOpBase. */
  //@{

  /** \brief . */
  virtual bool rowStatIsSupportedImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat) const;

  /** \brief . */
  virtual void getRowStatImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat,
    const Ptr<VectorBase<double> > &rowStatVec) const;

  //@}

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
  virtual RCP< const SpmdVectorSpaceBase<double> > 
  allocateDomain(
    const RCP<Epetra_Operator> &op, 
    EOpTransp op_trans 
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
  virtual RCP< const SpmdVectorSpaceBase<double> >
  allocateRange( 
    const RCP<Epetra_Operator> &op, 
    EOpTransp op_trans 
    ) const; 

  //@}

private:

  // ////////////////////////////////////
  // Private data members

  bool isFullyInitialized_;
  RCP<Epetra_Operator> op_;
  RCP<Epetra_RowMatrix> rowMatrix_;
  EOpTransp opTrans_;
  EApplyEpetraOpAs applyAs_;
  EAdjointEpetraOp adjointSupport_;
  RCP<const SpmdVectorSpaceBase<double> > range_;
  RCP<const SpmdVectorSpaceBase<double> > domain_;

  // ////////////////////////////////////
  // Private member functions

  const Epetra_Map& getRangeMap() const;
  const Epetra_Map& getDomainMap() const;

  /** \brief Compute the absolute row sum for this matrix.
    * 
    * A concrete implementation is required because Epetra
    * does not support absolute row sums.
    * 
    * \note This only works for Epetra_CrsMatrix objects.
    */
  void computeAbsRowSum(Epetra_Vector & rowStatVec_in) const;

};	// end class EpetraLinearOp


/** \brief Default nonmember constructor.
 *
 * \relates EpetraLinearOp
 */
RCP<EpetraLinearOp> nonconstEpetraLinearOp();


/** \brief Partially initialized EpetraLinearOp
 *
 * \relates EpetraLinearOp
 */
RCP<EpetraLinearOp>
partialNonconstEpetraLinearOp(
  const RCP<const VectorSpaceBase<double> > &range,
  const RCP<const VectorSpaceBase<double> > &domain,
  const RCP<Epetra_Operator> &op,
  EOpTransp opTrans = NOTRANS,
  EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED
  );


/** \brief Dynamically allocate an const EpetraLinearOp to wrap a const
 * Epetra_Operator object.
 *
 * \relates EpetraLinearOp
 */
RCP<EpetraLinearOp>
nonconstEpetraLinearOp(
  const RCP<Epetra_Operator> &op,
  EOpTransp opTrans = NOTRANS,
  EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED,
  const RCP< const VectorSpaceBase<double> > &range = Teuchos::null,
  const RCP< const VectorSpaceBase<double> > &domain = Teuchos::null
  );


/** \brief Dynamically allocate a nonconst EpetraLinearOp to wrap a const
 * Epetra_Operator object.
 *
 * \relates EpetraLinearOp
 */
RCP<const EpetraLinearOp>
epetraLinearOp(
  const RCP<const Epetra_Operator> &op,
  EOpTransp opTrans = NOTRANS,
  EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED,
  const RCP<const VectorSpaceBase<double> > &range = Teuchos::null,
  const RCP<const VectorSpaceBase<double> > &domain = Teuchos::null
  );


/** \brief Dynamically allocate an const EpetraLinearOp to wrap a const
 * Epetra_Operator object and give it a string label.
 *
 * \relates EpetraLinearOp
 */
RCP<EpetraLinearOp>
nonconstEpetraLinearOp(
  const RCP<Epetra_Operator> &op,
  const std::string &label,
  EOpTransp opTrans = NOTRANS,
  EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED,
  const RCP<const VectorSpaceBase<double> > &range = Teuchos::null,
  const RCP<const VectorSpaceBase<double> > &domain = Teuchos::null
  );


/** \brief Dynamically allocate a nonconst EpetraLinearOp to wrap a const
 * Epetra_Operator object.
 *
 * \relates EpetraLinearOp
 */
RCP<const EpetraLinearOp>
epetraLinearOp(
  const RCP<const Epetra_Operator> &op,
  const std::string &label,
  EOpTransp opTrans = NOTRANS,
  EApplyEpetraOpAs applyAs = EPETRA_OP_APPLY_APPLY,
  EAdjointEpetraOp adjointSupport = EPETRA_OP_ADJOINT_SUPPORTED,
  const RCP< const SpmdVectorSpaceBase<double> > &range = Teuchos::null,
  const RCP< const SpmdVectorSpaceBase<double> > &domain = Teuchos::null
  );


}	// end namespace Thyra


#endif	// THYRA_EPETRA_LINEAR_OP_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

