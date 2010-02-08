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

#ifndef THYRA_SCALAR_PROD_VECTOR_SPACE_BASE_DECL_HPP
#define THYRA_SCALAR_PROD_VECTOR_SPACE_BASE_DECL_HPP

#include "Thyra_OperatorVectorAdapterSupportTypes.hpp"
#include "Thyra_VectorSpaceDefaultBase.hpp"


namespace Thyra {


/** \brief Base subclass for <tt>VectorSpaceBase</tt> that allows the
 * definition of an application-specific scalar product to be swapped
 * in and out.
 *
 * This subclass defines machinery for extracting out the definition of a
 * scalar product as an object that can be replaced.  The default
 * implementation of scalar product is the Euclidean scalar product (i.e. dot
 * product).  The idea is that, in most cases, the definition of a scalar
 * product may be more general than a specific concrete vector implementation
 * (i.e. a single scalar product may work with all serial and all MPI-based
 * vectors if, for example, it is implemented through an
 * <tt>RTOpPack::RTOpT</tt> object).  Or, a scalar product way work with any
 * MPI SPMD vector or multi-vector.  This subclass allows an application code
 * to set a specialized scalar product without having to depend on a
 * particular concrete vector (and vector space) implementation.
 *
 * Almost every data-structure centric concrete <tt>VectorSpaceBase</tt>
 * subclass should inherit from this subclass since it makes it easy for
 * application developers to redefine the scalar product without having to
 * create a new <tt>VectorSpaceBase</tt> subclass which can have many
 * repercussions.
 *
 * The reason that this machinery in this base subclass is separated out from
 * the <tt>VectorSpaceDefaultBase</tt> interface class is that, first it would
 * clutter the base interface since this machinery is an implementation
 * artifact and, second, every <tt>VectorSpaceBase</tt> subclass will not
 * utilize this machinery.  For example, composite (see
 * <tt>ProductVectorSpaceBase</tt>) and decorator subclasses should not derive
 * from this subclass.
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class ScalarProdVectorSpaceBase : virtual public VectorSpaceDefaultBase<Scalar> {
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Construct to use dot product as the default.
   *
   * <b>Postconditions:</b><ul>
   *
   * <li><tt>dynamic_cast<const
   * EuclideanScalarProd<Scalar>*>(&*this->getScalarProd()) != NULL</tt>
   *
   * </ul>
   */
  ScalarProdVectorSpaceBase();

  /** \brief Construct with a different scalar product.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>scalarProd.get()!=NULL</tt> (throw
   * <tt>std::invalid_argument</tt>)
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li><tt>this->getScalarProd().get() == scalarProd.get()</tt>
   *
   * </ul>
   */
  ScalarProdVectorSpaceBase(
    const RCP<const ScalarProdBase<Scalar> > &scalarProd );

  /** \brief Set a different scalar product.
   *
   * This function is made virtual so that subclasses can override it
   * and take control of what happens.  However, any override should
   * call back on this base implementation to set the actual scalar
   * product object.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>scalarProd.get()!=NULL</tt> (throw
   * <tt>std::invalid_argument</tt>)
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li><tt>this->getScalarProd().get() == scalarProd.get()</tt>
   *
   * </ul>
   */
  virtual void setScalarProd(
    const RCP<const ScalarProdBase<Scalar> > &scalarProd );

  /** \brief Return the current scalar product.
   */
  RCP<const ScalarProdBase<Scalar> > getScalarProd() const;

  //@}

  /** @name Overridden from VectorSpaceBase */
  //@{
  
  /// Returns <tt>getScalarProd()->isEuclidean()</tt>
  bool isEuclidean() const;
  /// Returns <tt>getScalarProd()->scalarProd(x,y)</tt>
  Scalar scalarProd(
    const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const;
  /// Calls <tt>getScalarProd()->scalarProds(X,Y,scalar_prods)</tt>
  void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out ) const;
  
  //@}

private:

  RCP<const ScalarProdBase<Scalar> > scalarProd_;

};


/** \brief Create a small vector space casted to ScalarProdVectorSpaceBase.
 *
 * \relates ScalarProdVectorSpaceBase
 */
template<class Scalar>
RCP<const ScalarProdVectorSpaceBase<Scalar> >
createSmallScalarProdVectorSpaceBase(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const Ordinal dim
  )
{
  return Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
    vs->smallVecSpcFcty()->createVecSpc(dim), true);
}


} // end namespace Thyra


#endif  // THYRA_SCALAR_PROD_VECTOR_SPACE_BASE_DECL_HPP
