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

#ifndef THYRA_LINEAR_OP_SCALAR_PROD_DECL_HPP
#define THYRA_LINEAR_OP_SCALAR_PROD_DECL_HPP

#include "Thyra_ScalarProdBase.hpp"

namespace Thyra {

/** \brief Concrete implementation of a scalar product using a
 * symmetric positive definite linear operator..
 *
 * This subclass will work with any <tt>VectorBase</tt> or
 * <tt>MultiVectorBase</tt> implementation who's vector spaces are
 * compatible with the underlying linear operator object..
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class LinearOpScalarProd : public ScalarProdBase<Scalar> {
public:

  /** @name Constructors, initializers, accessors */
  //@{

  /** \brief . */
  LinearOpScalarProd();

  /** \brief . */
  LinearOpScalarProd( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &op );

  /** \brief . */
  void initialize( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &op );

  /** \brief . */
  const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >& op() const;

  /** \brief . */
  void uninitialize( Teuchos::RefCountPtr<const LinearOpBase<Scalar> > *op = NULL );

  //@}
  
  /** @name Overridden from ScalarProdBase */
  //@{

  /** \brief . */
  void scalarProds( const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y, Scalar scalar_prods[] ) const;

  /** \brief . */
  void apply(
    const EuclideanLinearOpBase<Scalar>   &M
    ,const ETransp                        M_trans
    ,const VectorBase<Scalar>             &x
    ,VectorBase<Scalar>                   *y
    ,const Scalar                         alpha
    ,const Scalar                         beta
    ) const;

  /** \brief . */
  void apply(
    const EuclideanLinearOpBase<Scalar>   &M
    ,const ETransp                        M_trans
    ,const MultiVectorBase<Scalar>        &X
    ,MultiVectorBase<Scalar>              *Y
    ,const Scalar                         alpha
    ,const Scalar                         beta
    ) const;

  //@}

private:

  Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  op_;

}; // end class LinearOpScalarProd

// //////////////////////////////////
// Inline members

template<class Scalar>
inline
const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >& LinearOpScalarProd<Scalar>::op() const
{
  return op_;
}

} // end namespace Thyra

#endif  // THYRA_LINEAR_OP_SCALAR_PROD_DECL_HPP
