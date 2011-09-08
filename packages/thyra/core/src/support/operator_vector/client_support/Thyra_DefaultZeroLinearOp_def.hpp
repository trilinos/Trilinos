// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
