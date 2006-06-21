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

#ifndef THYRA_LINEAR_OP_SCALAR_PROD_HPP
#define THYRA_LINEAR_OP_SCALAR_PROD_HPP

#include "Thyra_LinearOpScalarProdDecl.hpp"
#include "Thyra_ScalarProdBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_EuclideanLinearOpBase.hpp"

namespace Thyra {

// Constructors, initializers, accessors

template<class Scalar>
LinearOpScalarProd<Scalar>::LinearOpScalarProd()
{}

template<class Scalar>
LinearOpScalarProd<Scalar>::LinearOpScalarProd( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &op )
{
  this->initialize(op);
}

template<class Scalar>
void LinearOpScalarProd<Scalar>::initialize( const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &op )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(op.get()==NULL);
#endif
  op_ = op;
}

template<class Scalar>
void LinearOpScalarProd<Scalar>::uninitialize( Teuchos::RefCountPtr<const LinearOpBase<Scalar> > *op )
{
  if(op) *op = op_;
  op_ = Teuchos::null;
}

// Overridden from ScalarProdBase

template<class Scalar>
void LinearOpScalarProd<Scalar>::scalarProds( const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y, Scalar scalar_prods[] ) const
{
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    T = createMembers(Y.range(),Y.domain()->dim());
  Thyra::apply(*op_,NOTRANS,Y,&*T);
  dots(X,*T,scalar_prods);
}

template<class Scalar>
void LinearOpScalarProd<Scalar>::apply(
  const EuclideanLinearOpBase<Scalar>   &M
  ,const ETransp                        M_trans
  ,const MultiVectorBase<Scalar>        &X
  ,MultiVectorBase<Scalar>              *Y
  ,const Scalar                         alpha
  ,const Scalar                         beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(Y==NULL);
#endif
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    T = createMembers(X.range(),X.domain()->dim());
  Thyra::apply(*op_,NOTRANS,X,&*T);
  Thyra::euclideanApply(M,M_trans,*T,Y,alpha,beta);
}

} // end namespace Thyra

#endif  // THYRA_LINEAR_OP_SCALAR_PROD_HPP
