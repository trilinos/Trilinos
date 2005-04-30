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

namespace Thyra {

/** \brief Abstract interface for scalar products.
 * 
 * This interface is not considered a user-level interface.  Instead,
 * this interface is designed to subclassed off of and used with
 * <tt>VectorSpaceStdBase</tt> subclasses.  Applications should create
 * subclasses of this interface to define application-specific scalar
 * products (i.e. such as PDE code often do).
 *
 * This interface requires subclasses to override the multi-vector
 * version of scalar product function <tt>scalarProds()</tt> since
 * this will yield the most efficient implementation in a distributed
 * memory environment by requiring only a single global reduction
 * operation and a single communication.
 *
 * Note that one of the preconditions on the vector and multi-vector
 * arguments in <tt>scalarProd()</tt> and <tt>scalarProds()</tt> is a
 * little vague in stating that the vector or multi-vector objects
 * must be "compatible" with the underlying implementation of
 * <tt>*this</tt>.  The reason that this precondition must be vague is
 * that we can not expose a method to return a <tt>VectorSpaceBase</tt>
 * object that could be checked for compatibility since
 * <tt>%ScalarProdBase</tt> is used to define a <tt>VectorSpaceBase</tt>
 * object (through the <tt>VectorSpaceStdBase</tt> node subclass).
 * Also, some definitions of <tt>%ScalarProdBase</tt>
 * (i.e. <tt>DotProd</tt>) can work for any vector space
 * implementation since they only rely on <tt>RTOp</tt> operators.  In
 * other cases, however, an application-specific scalar product may
 * a have dependancy of the data-structure of vector and multi-vector
 * objects.
 *
 * This interface class also defines functions to modify the
 * application of a Euclidean linear operator to insert the definition
 * of the application specific scalar product.
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
  
  /** @name Pure virtual functions that must be overridden */
  //@{

  /** \brief Return the scalar product of each column in two multi-vectors in the vector space.
   *
   * @param  X            [in] Multi-vector.
   * @param  Y            [in] Multi-vector.
   * @param  scalar_prod  [out] Array (length <tt>X.domain()->dim()</tt>) containing the
   *                      scalar products <tt>scalar_prod[j-1] = this->scalarProd(*X.col(j),*Y.col(j))</tt>,
   *                      for <tt>j = 1 ... X.domain()->dim()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>X.domain()->isCompatible(*Y.domain())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li><tt>X.range()->isCompatible(*Y.range())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li>The MultiVectors <tt>X</tt> and <tt>Y</tt> are compatible with this implementation or
   *     an exception will be thrown.
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>scalar_prod[j-1] = this->scalarProd(*X.col(j),*Y.col(j))</tt>, for <tt>j = 1 ... X.domain()->dim()</tt>
   * </ul>
   */
  virtual void scalarProds( const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y, Scalar scalar_prods[] ) const = 0;

  /** \brief Modify the application of a Euclidean linear operator by
   * inserting the vector spaces scalar product.
   *
   *
   * The default implementation calls the single-vector version of
   * this function.
   */
  virtual void apply(
    const EuclideanLinearOpBase<Scalar>   &M
    ,const ETransp                        M_trans
    ,const VectorBase<Scalar>             &x
    ,VectorBase<Scalar>                   *y
    ,const Scalar                         alpha
    ,const Scalar                         beta
    ) const = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief Return the scalar product of two vectors in the vector space.
   *
   * Preconditions:<ul>
   * <li>The vectors <tt>X</tt> and <tt>Y</tt> are compatible with this implementation or
   *     an exception will be thrown.
   * <li><tt>x.space()->isCompatible(*y.space())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li>The scalar product is returned.
   * </ul>
   *
   * The default implementation calls on the multi-vector version
   * <tt>scalarProds()</tt>.
   */
  virtual Scalar scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const;

  /** \brief Modify the application of a Euclidean linear operator by
   * inserting the vector spaces scalar product.
   *
   *
   * The default implementation calls the single-vector version of
   * this function.
   */
  virtual void apply(
    const EuclideanLinearOpBase<Scalar>   &M
    ,const ETransp                        M_trans
    ,const MultiVectorBase<Scalar>        &X
    ,MultiVectorBase<Scalar>              *Y
    ,const Scalar                         alpha
    ,const Scalar                         beta
    ) const;

  //@}

}; // end class ScalarProdBase

} // end namespace Thyra

#endif  // THYRA_SCALAR_PROD_DECL_HPP
