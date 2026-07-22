// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SCALAR_PROD_BASE_DECL_HPP
#define THYRA_SCALAR_PROD_BASE_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Describable.hpp"


namespace Thyra {


/** \brief Abstract interface for scalar products.
 * 
 * This interface is not considered a user-level interface.  Instead, this
 * interface is designed to be sub-classed off of and used with
 * <tt>ScalarProdVectorSpaceBase</tt> objects to define their scalar products.
 * Applications should create subclasses of this interface to define
 * application-specific scalar products (i.e. such as PDE finite-element codes
 * often do).
 *
 * This interface requires subclasses to override a multi-vector version of
 * the scalar product function <tt>scalarProds()</tt>.  This version yields
 * the most efficient implementation in a distributed memory environment by
 * requiring only a single global reduction operation and a single
 * communication.
 *
 * Note that one of the preconditions on the vector and multi-vector arguments
 * in <tt>scalarProds()</tt> is a little vague in stating that the vector or
 * multi-vector objects must be "compatible" with the underlying
 * implementation of <tt>*this</tt>.  The reason that this precondition must
 * be vague is that we can not expose a method to return a
 * <tt>VectorSpaceBase</tt> object that could be checked for compatibility
 * since <tt>%ScalarProdBase</tt> is used to define a <tt>VectorSpaceBase</tt>
 * object (through the <tt>ScalarProdVectorSpaceBase</tt> node subclass).
 * Also, some definitions of <tt>%ScalarProdBase</tt>
 * (i.e. <tt>EuclideanScalarProd</tt>) will work for any vector space
 * implementation since they only rely on <tt>RTOp</tt> operators.  In other
 * cases, however, an application-specific scalar product may a have
 * dependency on the data-structure of vector and multi-vector objects in
 * which case one can not just use this with any vector or multi-vector
 * implementation.
 *
 * This interface class also defines functions to modify the application of a
 * Euclidean linear operator to insert the definition of the application
 * specific scalar product.
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class ScalarProdBase : virtual public Teuchos::Describable {
public:
  
  /** @name Non-virtual public interface */
  //@{

  /** \brief Return if this is a Euclidean (identity) scalar product is the
   * same as the dot product.
   *
   * The default implementation returns <tt>false</tt> (evenn though on average
   * the truth is most likely <tt>true</tt>).
   */
  bool isEuclidean() const
    { return isEuclideanImpl(); }

  /** \brief Return the scalar product of two vectors in the vector space.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li>The vectors <tt>x</tt> and <tt>y</tt> are <em>compatible</em> with
   * <tt>*this</tt> implementation or an exception will be thrown.
   *
   * <li><tt>x.space()->isCompatible(*y.space())</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>The scalar product is returned.
   *
   * </ul>
   *
   * The default implementation calls on the multi-vector version
   * <tt>scalarProds()</tt>.
   */
  Scalar scalarProd(
    const VectorBase<Scalar>& x, const VectorBase<Scalar>& y
    ) const
    { return scalarProdImpl(x, y); }

  /** \brief Return the scalar product of each column in two multi-vectors in
   * the vector space.
   *
   * \param X [in] Multi-vector.
   *
   * \param Y [in] Multi-vector.
   *
   * \param scalar_prod [out] Array (length <tt>X.domain()->dim()</tt>)
   * containing the scalar products <tt>scalar_prod[j] =
   * this->scalarProd(*X.col(j),*Y.col(j))</tt>, for <tt>j = 0
   * ... X.domain()->dim()-1</tt>.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>X.domain()->isCompatible(*Y.domain())</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li><tt>X.range()->isCompatible(*Y.range())</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li>The MultiVectorBase objects <tt>X</tt> and <tt>Y</tt> are
   * <em>compatible</em> with this implementation or an exception will be
   * thrown.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li><tt>scalar_prod[j] = this->scalarProd(*X.col(j),*Y.col(j))</tt>, for
   * <tt>j = 0 ... X.domain()->dim()-1</tt>
   *
   * </ul>
   */
  void scalarProds(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out
    ) const
    { scalarProdsImpl(X, Y, scalarProds_out); }

  /** \brief Return a linear operator representing the scalar product
   * <tt>Q</tt>.
   *
   * All scalar products are not required to return this operator so a return
   * value of <tt>null</tt> is allowed.  Note that if <tt>this->isEuclidean()
   * == true</tt> then there is no reason to return an identity operator.
   */
  RCP<const LinearOpBase<Scalar> > getLinearOp() const
    { return getLinearOpImpl(); }

  //@}

protected:

  /** \name Protected virtual functions. */
  //@{

  /** \brief . */
  virtual bool isEuclideanImpl() const = 0;
  
  /** \brief Default implementation calls scalarProdsImpl(). */
  virtual Scalar scalarProdImpl(
    const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const;

  /** \brief . */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out
    ) const = 0;

  /** \brief . */
  virtual RCP<const LinearOpBase<Scalar> > getLinearOpImpl() const
    {
      return Teuchos::null;
    }

  //@}

};


} // end namespace Thyra


#endif  // THYRA_SCALAR_PROD_BASE_DECL_HPP
