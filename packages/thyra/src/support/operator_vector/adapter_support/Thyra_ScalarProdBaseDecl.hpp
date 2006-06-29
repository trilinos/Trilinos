// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_SCALAR_PROD_DECL_HPP
#define THYRA_SCALAR_PROD_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_EuclideanLinearOpBaseDecl.hpp"

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
 * object (through the <tt>VectorSpaceStdBase</tt> node subclass).  Also, some
 * definitions of <tt>%ScalarProdBase</tt> (i.e. <tt>EuclideanScalarProd</tt>)
 * will work for any vector space implementation since they only rely on
 * <tt>RTOp</tt> operators.  In other cases, however, an application-specific
 * scalar product may a have dependency on the data-structure of vector and
 * multi-vector objects in which case one can not just use this with any
 * vector or multi-vector implementation.
 *
 * This interface class also defines functions to modify the application of a
 * Euclidean linear operator to insert the definition of the application
 * specific scalar product.
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class ScalarProdBase {
public:
  
  /** @name Destructor */
  //@{
  
  /** \brief . */
  virtual ~ScalarProdBase() {}

  //@}
  
  /** @name Public pure virtual functions that must be overridden */
  //@{

  /** \brief Return the scalar product of each column in two multi-vectors in the vector space.
   *
   * @param  X            [in] Multi-vector.
   * @param  Y            [in] Multi-vector.
   * @param  scalar_prod  [out] Array (length <tt>X.domain()->dim()</tt>) containing the
   *                      scalar products <tt>scalar_prod[j] = this->scalarProd(*X.col(j),*Y.col(j))</tt>,
   *                      for <tt>j = 0 ... X.domain()->dim()-1</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>X.domain()->isCompatible(*Y.domain())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li><tt>X.range()->isCompatible(*Y.range())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li>The MultiVectorBase objects <tt>X</tt> and <tt>Y</tt> are <em>compatible</em> with this implementation or
   *     an exception will be thrown.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>scalar_prod[j] = this->scalarProd(*X.col(j),*Y.col(j))</tt>, for <tt>j = 0 ... X.domain()->dim()-1</tt>
   * </ul>
   */
  virtual void scalarProds( const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y, Scalar scalar_prods[] ) const = 0;

  /** \brief Modify the application of a Euclidean linear operator by
   * inserting the vector space's scalar product.
   *
   * Note that one responsibility of an implementation of this function is to
   * provide the block scalar product implementation of
   * <tt>MultiVectorBase</tt> objects that derive from
   * <tt>EuclideanLinearOpBase</tt>.  For example, let <tt>M</tt> be a
   * <tt>%MultiVectorBase</tt> object and consider the operation
   
   <tt>Y = adjoint(M)*X</tt>

   * where <tt>M_trans==CONJTRANS</tt>.  This function may, or many not, call
   * the <tt>EuclideanLinearOpBase::euclideanApplyTranspose()</tt> function in
   * order to implement this block Scalar product.
   *
   * Note that the special case of <tt>M==X</tt> should also be supported
   * which provides the symmetric operation
   
   <tt>Y = adjoint(X)*X</tt>

   * that can be performed in half the flops as the general case.
   *
   * ToDo: Finish documentation!
   */
  virtual void apply(
    const EuclideanLinearOpBase<Scalar>   &M
    ,const ETransp                        M_trans
    ,const MultiVectorBase<Scalar>        &X
    ,MultiVectorBase<Scalar>              *Y
    ,const Scalar                         alpha
    ,const Scalar                         beta
    ) const = 0;

  //@}

  /** @name Public virtual functions with default implementations */
  //@{

  /** \brief Return if this is a Euclidean (identity) scalar product is the
   * same as the dot product.
   *
   * The default implementation returns <tt>false</tt> (evenn though on average
   * the truth is most likely <tt>true</tt>).
   */
  virtual bool isEuclidean() const;

  /** \brief Return the scalar product of two vectors in the vector space.
   *
   * <b>Preconditions:</b><ul>
   * <li>The vectors <tt>x</tt> and <tt>y</tt> are <em>compatible</em> with <tt>*this</tt>
   *     implementation or an exception will be thrown.
   * <li><tt>x.space()->isCompatible(*y.space())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>The scalar product is returned.
   * </ul>
   *
   * The default implementation calls on the multi-vector version
   * <tt>scalarProds()</tt>.
   */
  virtual Scalar scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const;

  //@}

}; // end class ScalarProdBase

} // end namespace Thyra

#endif  // THYRA_SCALAR_PROD_DECL_HPP
