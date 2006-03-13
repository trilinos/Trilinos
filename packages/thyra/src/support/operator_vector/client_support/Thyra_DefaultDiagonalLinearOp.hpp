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

#ifndef THYRA_DIAGONAL_LINEAR_OP_HPP
#define THYRA_DIAGONAL_LINEAR_OP_HPP

#include "Thyra_DefaultDiagonalLinearOpDecl.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_VectorBase.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp()
{}

template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >   &space
  )
{
  initialize(space);
}

template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp(
  const Teuchos::RefCountPtr<VectorBase<Scalar> >   &diag
  )
{
  initialize(diag);
}

template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp(
  const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
  )
{
  initialize(diag);
}

template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::initialize(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &space
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(space.get()==NULL);
#endif
  initialize(createMember(space));
}

template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::initialize(
  const Teuchos::RefCountPtr<VectorBase<Scalar> >   &diag
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(diag.get()==NULL);
#endif
  nonconstDiag_ = diag;
  diag_         = diag;
}

template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::initialize(
  const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(diag.get()==NULL);
#endif
  nonconstDiag_ = Teuchos::null;
  diag_         = diag;
}

template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::uninitialize(
  Teuchos::RefCountPtr<VectorBase<Scalar> >   *diag
  )
{
  if(diag) *diag = nonconstDiag_;
  nonconstDiag_ = Teuchos::null;
  diag_         = Teuchos::null;
}

template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const VectorBase<Scalar> >   *diag
  )
{
  if(diag) *diag = diag_;
  nonconstDiag_ = Teuchos::null;
  diag_         = Teuchos::null;
}

// Overridden from DiagonalLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >  
DefaultDiagonalLinearOp<Scalar>::getNonconstDiag()
{
  if(diag_.get())
    return nonconstDiag_.assert_not_null();
  return nonconstDiag_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >  
DefaultDiagonalLinearOp<Scalar>::getDiag() const
{
  return diag_;
}

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultDiagonalLinearOp<Scalar>::range() const
{
  return diag_->space();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultDiagonalLinearOp<Scalar>::domain() const
{
  return diag_->space();
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultDiagonalLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool DefaultDiagonalLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
  return true; // ToDo: Update this!
}

// Overridden from SingleRhsLinearOpBase

template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::apply(
  const ETransp                M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  ,const Scalar                beta
  ) const
{
  // RAB: 4/16/2005: Warning! this does not work if Scalar is a complex type
  // and M_trans==CONJTRANS!
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if( beta != ST::one() ) Vt_S( y, beta );
  ele_wise_prod( alpha, x, *diag_, y );
}

}	// end namespace Thyra

#endif	// THYRA_DIAGONAL_LINEAR_OP_HPP
