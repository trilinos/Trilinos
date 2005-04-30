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

#ifndef THYRA_EUCLIDEAN_SCALAR_PROD_HPP
#define THYRA_EUCLIDEAN_SCALAR_PROD_HPP

#include "Thyra_EuclideanScalarProdDecl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_EuclideanLinearOpBase.hpp"

namespace Thyra {

template<class Scalar>
void EuclideanScalarProd<Scalar>::scalarProds( const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y, Scalar scalar_prods[] ) const
{
  dots(X,Y,scalar_prods);
}

template<class Scalar>
void EuclideanScalarProd<Scalar>::apply(
  const EuclideanLinearOpBase<Scalar>   &M
  ,const ETransp                        M_trans
  ,const VectorBase<Scalar>             &x
  ,VectorBase<Scalar>                   *y
  ,const Scalar                         alpha
  ,const Scalar                         beta
  ) const
{
  M.euclideanApply(M_trans,x,y,alpha,beta);
}

template<class Scalar>
void EuclideanScalarProd<Scalar>::apply(
  const EuclideanLinearOpBase<Scalar>   &M
  ,const ETransp                        M_trans
  ,const MultiVectorBase<Scalar>        &X
  ,MultiVectorBase<Scalar>              *Y
  ,const Scalar                         alpha
  ,const Scalar                         beta
  ) const
{
  M.euclideanApply(M_trans,X,Y,alpha,beta);
}

} // end namespace Thyra

#endif  // THYRA_EUCLIDEAN_SCALAR_PROD_HPP
