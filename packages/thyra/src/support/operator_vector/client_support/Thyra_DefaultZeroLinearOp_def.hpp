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

#ifndef THYRA_DEFAULT_ZERO_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_ZERO_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultZeroLinearOp_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultZeroLinearOp<Scalar>::DefaultZeroLinearOp()
{}


template<class Scalar>
DefaultZeroLinearOp<Scalar>::DefaultZeroLinearOp(
  const RCP<const VectorSpaceBase<Scalar> > &range_in,
  const RCP<const VectorSpaceBase<Scalar> > &domain_in
  )
{
  initialize(range_in,domain_in);
}


template<class Scalar>
void DefaultZeroLinearOp<Scalar>::initialize(
  const RCP<const VectorSpaceBase<Scalar> > &range_in,
  const RCP<const VectorSpaceBase<Scalar> > &domain_in
  )
{
  range_ = range_in.assert_not_null();
  domain_ = domain_in.assert_not_null();
}


template<class Scalar>
void DefaultZeroLinearOp<Scalar>::uninitialize()
{
  range_ = Teuchos::null;
  domain_ = Teuchos::null;
}


// Overridden from LinearOpBase

  
template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultZeroLinearOp<Scalar>::range() const
{
  return range_;
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultZeroLinearOp<Scalar>::domain() const
{
  return domain_;
}

  
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultZeroLinearOp<Scalar>::clone() const
{
  typedef DefaultZeroLinearOp<Scalar> this_t;
  if(range_.get())
    return Teuchos::rcp(new this_t(range_,domain_));
  return Teuchos::rcp(new this_t());
}

  
// Overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultZeroLinearOp<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  std::ostringstream oss;
  oss
    << "Thyra::DefaultZeroLinearOp<" << ST::name() << ">{"
    << "range="<<(range_.get()?range_->description():"NULL")
    << ",domain="<<(domain_.get()?domain_->description():"NULL")
    << "}";
  return oss.str();
}


// protected


// Overridden from LinearOpBase


template<class Scalar>
bool DefaultZeroLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  return true;
}


template<class Scalar>
void DefaultZeroLinearOp<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultZeroLinearOp<Scalar>::apply(...)", *this, M_trans, X, &*Y
    );
#endif // TEUCHOS_DEBUG  
  scale(beta, Y);
}


}	// end namespace Thyra


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::zero(
  const RCP<const VectorSpaceBase<Scalar> > &range_in,
  const RCP<const VectorSpaceBase<Scalar> > &domain_in
  )
{
  return Teuchos::rcp(new DefaultZeroLinearOp<Scalar>(range_in, domain_in));
}


//
// Explicit instantaition
//


#define THYRA_DEFAULT_ZERO_LINEAR_OP_INSTANT(SCALAR) \
  \
  template class DefaultZeroLinearOp<SCALAR >; \
  \
  template RCP<const Thyra::LinearOpBase<SCALAR > >  \
  zero(  \
    const RCP<const VectorSpaceBase<SCALAR > > &range,  \
    const RCP<const VectorSpaceBase<SCALAR > > &domain  \
    );  \


#endif	// THYRA_DEFAULT_ZERO_LINEAR_OP_DEF_HPP
