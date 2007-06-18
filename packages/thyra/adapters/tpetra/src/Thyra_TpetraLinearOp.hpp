
#ifndef THYRA_TPETRA_LINEAR_OP_HPP
#define THYRA_TPETRA_LINEAR_OP_HPP

#include "Thyra_TpetraLinearOpBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_SingleRhsEuclideanLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_getConst.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpBase</tt> adapter subclass for
 * <tt>Tpetra::Operator</tt> objects.
 *
 * This subclass can be used to represent the non-transposed operator or
 * transposed operator defined by an <tt>Tpetra::Operator</tt> object. In
 * addition, the user can specify whether adjoints are supported or not.
 *
 * <b>Partial Automatic Change Propagation:</b> This class shall maintain no
 * state with respect to the <em>values</em> of the internally stored
 * <tt>Tpetra::Operator</tt> object.  Therefore, as long as the domain and
 * range spaces do not change, the the <tt>Tpetra::Operator</tt> can be
 * changed and this will automatically update <tt>*this</tt> object.  This
 * simplifies some types of update operations.  Since this is a simple
 * concrete class, this behavior is harmless.  However, if the range and
 * domain maps change, then one must re-call the <tt>this->initialize()</tt>
 * function to update <tt>*this</tt>.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template<class Ordinal, class Scalar>
class TpetraLinearOp
  : virtual public TpetraLinearOpBase<Ordinal,Scalar>
  , virtual public SingleRhsEuclideanLinearOpBase<Scalar>
{
public:

  /** \brief . */
  using SingleRhsEuclideanLinearOpBase<Scalar>::euclideanApply;

  /** @name Constructors / initializers / accessors */
  //@{

  /** \brief Construct to uninitialized.
   *
   * See the postconditions for <tt>uninitialize()</tt>
   */
  TpetraLinearOp();

  /** \brief Calls <tt>initialize()</tt>. */
  TpetraLinearOp(
    const Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >       &op
    ,EAdjointTpetraOp                                                   adjointSupport  = TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     &mpiRange       = Teuchos::null
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     &mpiDomain      = Teuchos::null
    );

  /** \brief Calls <tt>initialize()</tt>. */
  TpetraLinearOp(
    const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> > &op
    ,EAdjointTpetraOp                                                   adjointSupport  = TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     &mpiRange       = Teuchos::null
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     &mpiDomain      = Teuchos::null
    );

  /** \brief Initialize given a non-constant <tt>Tpetra::Operator</tt> object.
   *
   * @param  op       [in] The <tt>Tpetra::Operator</tt> this <tt>*this</tt> will wrap.
   * @param  adjointSupport
   *                  [in] Determines if it is to be assumed that adjoints are supported on the
   *                  underlying <tt>Tpetra::Operator</tt> object <tt>op</tt>.  If
   *                  <tt>adjointSupport==TPETRA_OP_ADJOINT_SUPPORTED</tt> then <tt>this->opSupported(TRANS)</tt>
   *                  will return <tt>true</tt>.  If <tt>adjointSupport==TPETRA_OP_ADJOINT_UNSUPPORTED</tt> then
   *                  <tt>this->opSupported(TRANS)</tt> will return <tt>false</tt>.
   * @param  mpiRange
   *                  [in] Smart pointer to the range space for the <tt>Tpetra::Operator</tt>.  The default
   *                  value is <tt>Teuchos::null</tt> in which case <tt>*this</tt> will allocate
   *                  a new <tt>MPIVectorSpace</tt> given range map from <tt>op</tt>.  A client may only bother
   *                  to specify this space if one wants to override the defintion of the scalar product.
   * @param  mpiDomain
   *                  [in] Smart pointer to the domain space for the <tt>Tpetra::Operator</tt>.  The default
   *                  value is <tt>Teuchos::null</tt> in which case <tt>*this</tt> will allocate
   *                  a new <tt>DefaultSpmdVectorSpace</tt> given map from <tt>op</tt>.  A client may only bother
   *                  to specify this space if one wants to override the defintion of the scalar product.
   *
   * Preconditions:<ul>
   * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->tpetra_op().get() == op.get()</tt>
   * <li> [<tt>mpiRange.get() != NULL</tt>] <tt>this->mpiRange().get() == mpiRange.get()</tt>
   * <li> [<tt>mpiDomain.get() != NULL</tt>] <tt>this->mpiDomain().get() == mpiDomain.get()</tt>
   * <li> [<tt>mpiRange.get() == NULL</tt>] <tt>this->mpiRange().get() != NULL</tt>
   * <li> [<tt>mpiDomain.get() == NULL</tt>] <tt>this->mpiDomain().get() != NULL</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >      &op
    ,EAdjointTpetraOp                                                  adjointSupport  = TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >    &mpiRange       = Teuchos::null
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >    &mpiDomain      = Teuchos::null
    );

  /** \brief Initialize given a constant <tt>Tpetra::Operator</tt> object.
   *
   * @param  op       [in] The <tt>Tpetra::Operator</tt> this <tt>*this</tt> will wrap.
   * @param  adjointSupport
   *                  [in] Determines if it is to be assumed that adjoints are supported on the
   *                  underlying <tt>Tpetra::Operator</tt> object <tt>op</tt>.  If
   *                  <tt>adjointSupport==TPETRA_OP_ADJOINT_SUPPORTED</tt> then <tt>this->opSupported(TRANS)</tt>
   *                  will return <tt>true</tt>.  If <tt>adjointSupport==TPETRA_OP_ADJOINT_UNSUPPORTED</tt> then
   *                  <tt>this->opSupported(TRANS)</tt> will return <tt>false</tt>.
   * @param  mpiRange
   *                  [in] Smart pointer to the range space for the <tt>Tpetra::Operator</tt>.  The default
   *                  value is <tt>Teuchos::null</tt> in which case <tt>*this</tt> will allocate
   *                  a new <tt>MPIVectorSpace</tt> given range map from <tt>op</tt>.  A client may only bother
   *                  to specify this space if one wants to override the defintion of the scalar product.
   * @param  mpiDomain
   *                  [in] Smart pointer to the domain space for the <tt>Tpetra::Operator</tt>.  The default
   *                  value is <tt>Teuchos::null</tt> in which case <tt>*this</tt> will allocate
   *                  a new <tt>DefaultSpmdVectorSpace</tt> given map from <tt>op</tt>.  A client may only bother
   *                  to specify this space if one wants to override the defintion of the scalar product.
   *
   * Preconditions:<ul>
   * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->tpetra_op().get() == op.get()</tt>
   * <li> [<tt>mpiRange.get() != NULL</tt>] <tt>this->mpiRange().get() == mpiRange.get()</tt>
   * <li> [<tt>mpiDomain.get() != NULL</tt>] <tt>this->mpiDomain().get() == mpiDomain.get()</tt>
   * <li> [<tt>mpiRange.get() == NULL</tt>] <tt>this->mpiRange().get() != NULL</tt>
   * <li> [<tt>mpiDomain.get() == NULL</tt>] <tt>this->mpiDomain().get() != NULL</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> > &op
    ,EAdjointTpetraOp                                                   adjointSupport  = TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     &mpiRange       = Teuchos::null
    ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     &mpiDomain      = Teuchos::null
    );
  
  /** \brief Set to uninitialized and optionally return the current state.
   *
   * Postconditions:<ul>
   * <li> <tt>this->domain().get() == NULL</tt>
   * <li> <tt>this->range().get() == NULL</tt>
   * </ul>
   */
  void uninitialize(
    Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >     *op             = NULL
    ,EAdjointTpetraOp                                           *adjointSupport = NULL
    ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   *mpiRange       = NULL
    ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   *mpiDomain      = NULL
    );
  
  /** \brief Set to uninitialized and optionally return the current state.
   *
   * Postconditions:<ul>
   * <li> <tt>this->domain().get() == NULL</tt>
   * <li> <tt>this->range().get() == NULL</tt>
   * </ul>
   */
  void uninitialize(
    Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> > *op             = NULL
    ,EAdjointTpetraOp                                             *adjointSupport = NULL
    ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     *mpiRange       = NULL
    ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     *mpiDomain      = NULL
    );

  /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the range.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->range().get() != NULL</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->range().get() == NULL</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > mpiRange() const;

  /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the domain.
   *
   * Postconditions:<ul>
   * <li> [<tt>this->domain().get() != NULL</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->domain().get() == NULL</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > mpiDomain() const;

  /** \brief . */
  bool isTpetraOpConst() const;

  /** \brief . */
  Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> > getNonconstTpetraOp();

  /** \brief . */
  Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> > getTpetraOp() const;

  //@}

  /** @name Overridden from TpetraLinearOpBase */
  //@{

  /** \brief . */
  void getTpetraOpView(
    Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >   *tpetraOp
    ,EAdjointTpetraOp                                         *tpetraOpAdjointSupport
    );
  /** \brief . */
  void getTpetraOpView(
    Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >   *tpetraOp
    ,EAdjointTpetraOp                                               *tpetraOpAdjointSupport
    ) const;

  //@}

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  
  //@}
  
  /** @name Overridden from EuclideanLinearOpBase */
  //@{

  /// Returns <tt>this->mpiRange()</tt>
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;
  /// Returns <tt>this->mpiDomain()</tt>
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;

  //@}
  
  /** @name Overridden from SingleRhsEuclideanLinearOpBase */
  //@{

  /** \brief . */
  void euclideanApply(
    const ETransp                     M_trans
    ,const VectorBase<Scalar>         &x
    ,VectorBase<Scalar>               *y
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

  /** \brief Allocate the domain space of the operator. */
  virtual Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > 
  allocateDomain(
    const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >  &op 
    ) const; 
  
  /** \brief Allocate the range space of the operator. */
  virtual Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >
  allocateRange( 
    const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >  &op 
    ) const; 
  
  //@}

private:

  // ////////////////////////////////////
  // Private data members

  Teuchos::ConstNonconstObjectContainer<Tpetra::Operator<Ordinal,Scalar> >  op_;
  EAdjointTpetraOp                                                          adjointSupport_;
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >                  range_;
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >                  domain_;
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >           sp_range_;
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >           sp_domain_;

};	// end class TpetraLinearOp

} // namespace Thyra

// /////////////////////////
// Implementations

#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

namespace Thyra {

// Constructors / initializers / accessors

template<class Ordinal, class Scalar>
TpetraLinearOp<Ordinal,Scalar>::TpetraLinearOp()
  :adjointSupport_(TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED)
{}

template<class Ordinal, class Scalar>
TpetraLinearOp<Ordinal,Scalar>::TpetraLinearOp(
  const Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >      &op
  ,EAdjointTpetraOp                                                  adjointSupport
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >    &mpiRange
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >    &mpiDomain
  )
{
  this->initialize(op,adjointSupport,mpiRange,mpiDomain);
}

template<class Ordinal, class Scalar>
TpetraLinearOp<Ordinal,Scalar>::TpetraLinearOp(
  const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >  &op
  ,EAdjointTpetraOp                                                    adjointSupport
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >      &mpiRange
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >      &mpiDomain
  )
{
  this->initialize(op,adjointSupport,mpiRange,mpiDomain);
}

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::initialize(
  const Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >      &op
  ,EAdjointTpetraOp                                                  adjointSupport
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >    &mpiRange
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >    &mpiDomain
  )
{
  TEST_FOR_EXCEPT(adjointSupport==TPETRA_OP_ADJOINT_SUPPORTED||adjointSupport==TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED);
  op_.initialize(op);
  adjointSupport_ = adjointSupport;
  range_  = ( mpiRange.get()  ? mpiRange  : allocateRange(op)  );
  domain_ = ( mpiDomain.get() ? mpiDomain : allocateDomain(op) );
  sp_range_ = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(range_);
  sp_domain_ = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(domain_);
}

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::initialize(
  const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >  &op
  ,EAdjointTpetraOp                                                    adjointSupport
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >      &mpiRange
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >      &mpiDomain
  )
{
  TEST_FOR_EXCEPT(adjointSupport==TPETRA_OP_ADJOINT_SUPPORTED||adjointSupport==TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED);
  op_.initialize(op);
  adjointSupport_ = adjointSupport;
  range_  = ( mpiRange.get()  ? mpiRange  : allocateRange(op)  );
  domain_ = ( mpiDomain.get() ? mpiDomain : allocateDomain(op) );
  sp_range_ = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(range_);
  sp_domain_ = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(domain_);
}

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::uninitialize(
  Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >     *op
  ,EAdjointTpetraOp                                           *adjointSupport
  ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   *mpiRange
  ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   *mpiDomain
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::uninitialize(
  Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> > *op
  ,EAdjointTpetraOp                                             *adjointSupport
  ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     *mpiRange
  ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >     *mpiDomain
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >
TpetraLinearOp<Ordinal,Scalar>::mpiRange() const
{
  return range_;
}

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >
TpetraLinearOp<Ordinal,Scalar>::mpiDomain() const
{
  return domain_;
}

template<class Ordinal, class Scalar>
bool TpetraLinearOp<Ordinal,Scalar>::isTpetraOpConst() const
{
  return op_.isConst();
}

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >
TpetraLinearOp<Ordinal,Scalar>::getNonconstTpetraOp()
{
  return op_.getNonconstObj();
}

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >
TpetraLinearOp<Ordinal,Scalar>::getTpetraOp() const
{
  return op_.getConstObj();
}

// Overridden from TpetraLinearOpBase

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::getTpetraOpView(
  Teuchos::RefCountPtr<Tpetra::Operator<Ordinal,Scalar> >   *tpetraOp
  ,EAdjointTpetraOp                                         *tpetraOpAdjointSupport
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(tpetraOp==NULL);
  TEST_FOR_EXCEPT(tpetraOpAdjointSupport==NULL);
#endif
  *tpetraOp = op_.getNonconstObj();
  *tpetraOpAdjointSupport = adjointSupport_;
}

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::getTpetraOpView(
  Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >   *tpetraOp
  ,EAdjointTpetraOp                                               *tpetraOpAdjointSupport
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(tpetraOp==NULL);
  TEST_FOR_EXCEPT(tpetraOpAdjointSupport==NULL);
#endif
  *tpetraOp = op_.getConstObj();
  *tpetraOpAdjointSupport = adjointSupport_;
}

// Overridden from SingleScalarLinearOpBase

template<class Ordinal, class Scalar>
bool TpetraLinearOp<Ordinal,Scalar>::opSupported(ETransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  switch(M_trans) {
    case NOTRANS:
      return true;
    case CONJ:
      return ( ST::isComplex ? false : true );
    case TRANS:
      return ( ST::isComplex
               ? ( adjointSupport_==TPETRA_OP_TRANSPOSE_SUPPORTED
                   || adjointSupport_==TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED )
               : adjointSupport_!=TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED );
    case CONJTRANS:
      return ( ST::isComplex
               ? ( adjointSupport_==TPETRA_OP_ADJOINT_SUPPORTED
                   || adjointSupport_==TPETRA_OP_TRANSPOSE_ADJOINT_SUPPORTED )
               : adjointSupport_!=TPETRA_OP_TRANSPOSE_ADJOINT_UNSUPPORTED );
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  return false; // Will never be executed!
}

// Overridden from EuclideanLinearOpBase

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
TpetraLinearOp<Ordinal,Scalar>::rangeScalarProdVecSpc() const
{
  return sp_range_;
}

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
TpetraLinearOp<Ordinal,Scalar>::domainScalarProdVecSpc() const
{
  return sp_domain_;
}

// Overridden from SingleRhsEuclideanLinearOpBase

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::euclideanApply(
  const ETransp                     M_trans
  ,const VectorBase<Scalar>         &x_in
  ,VectorBase<Scalar>               *y_inout
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  const ETransp real_M_trans = real_trans(M_trans);
#ifdef TEUCHOS_DEBUG
  // ToDo: Assert vector spaces!
  TEST_FOR_EXCEPTION(
    !opSupported(M_trans), Exceptions::OpNotSupported
    ,"TpetraLinearOp::apply(...): *this was informed that adjoints are not supported when initialized." 
    );
#endif
  //
  // Get Tpetra::Vector objects for the arguments
  //
  Teuchos::RefCountPtr<const Tpetra::Vector<Ordinal,Scalar> >
    x = get_Tpetra_Vector<Ordinal,Scalar>(
      real_M_trans==NOTRANS ? op_.getConstObj()->getDomainDist() : op_.getConstObj()->getRangeDist()
      ,Teuchos::rcp(&x_in,false)
      );
  Teuchos::RefCountPtr<Tpetra::Vector<Ordinal,Scalar> >
    y;
  if( beta == ST::zero() ) {
    y = get_Tpetra_Vector<Ordinal,Scalar>(
      real_M_trans==NOTRANS ? op_.getConstObj()->getRangeDist() : op_.getConstObj()->getDomainDist()
      ,Teuchos::rcp(y_inout,false)
      );
  }
  //
  // Perform the operation
  //
  if( beta == ST::zero() ) {
    // y = M * x
    op_.getConstObj()->apply( *x, *y, real_M_trans==TRANS );
    // y = alpha * y
    if( alpha != ST::one() ) y->scale(alpha);
  }
  else {
    // y_inout = beta * y_inout
    if( beta != ST::zero() ) scale( beta, y_inout );
    else assign( y_inout, ST::zero() );
    // t = M * x
    Tpetra::Vector<Ordinal,Scalar>
      t(real_M_trans==NOTRANS ? op_.getConstObj()->getRangeDist() : op_.getConstObj()->getDomainDist());
    op_.getConstObj()->apply( *x, t, real_M_trans==TRANS );
    // y_inout += alpha * t
    Vp_StV(
      y_inout
      ,alpha
      ,*create_Vector(
        Teuchos::rcp(&Teuchos::getConst(t),false)
        ,Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(y_inout->range(),true)
        )
      );
  }
}

// Overridden from LinearOpBase

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
TpetraLinearOp<Ordinal,Scalar>::clone() const
{
  return Teuchos::null; // We can not support this at this time!
}

// Overridden from Teuchos::Describable

template<class Ordinal, class Scalar>
std::string TpetraLinearOp<Ordinal,Scalar>::description() const
{
  using Teuchos::ScalarTraits;
  std::ostringstream oss;
  oss
    << "Thyra::TpetraLinearOp<"
    << Teuchos::ScalarTraits<Ordinal>::name()
    <<","
    << Teuchos::ScalarTraits<Scalar>::name()
    << ">";
  oss << "{";
  if(op_.getConstObj().get()) {
    oss << "op=\'"<<typeName(*op_.getConstObj())<<"\'";
  }
  else {
    oss << "op=NULL";
  }
  oss << "}";
  return oss.str();
}

template<class Ordinal, class Scalar>
void TpetraLinearOp<Ordinal,Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << "Thyra::TpetraLinearOp<"
        << Teuchos::ScalarTraits<Ordinal>::name()
        <<","
        << Teuchos::ScalarTraits<Scalar>::name()
        << ">,"
        << "rangeDim = " << this->range()->dim() << ", domainDim = " << this->domain()->dim() << std::endl;
      OSTab tab(out);
      if(op_.getConstObj().get()) {
        *out << "op=\'"<<typeName(*op_.getConstObj())<<"\'\n";
        *out << "adjointSupport="<<toString(adjointSupport_)<<"\n";
      }
      else {
        *out << "op=NULL"<<"\n";
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

// protected

// Allocators for domain and range spaces

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > 
TpetraLinearOp<Ordinal,Scalar>::allocateDomain(
  const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >  &op 
  ) const
{
  return create_VectorSpace<Ordinal,Scalar>(
    Teuchos::rcp(new Tpetra::VectorSpace<Ordinal,Scalar>(op->getDomainDist()))
    );
} 

template<class Ordinal, class Scalar>
Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >
TpetraLinearOp<Ordinal,Scalar>::allocateRange( 
  const Teuchos::RefCountPtr<const Tpetra::Operator<Ordinal,Scalar> >  &op 
  ) const
{
  return create_VectorSpace<Ordinal,Scalar>(
    Teuchos::rcp(new Tpetra::VectorSpace<Ordinal,Scalar>(op->getRangeDist()))
    );
}

}	// end namespace Thyra

#endif	// THYRA_TPETRA_LINEAR_OP_HPP
