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

#ifndef THYRA_DEFAULT_IDENTITY_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_IDENTITY_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultIdentityLinearOp_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultIdentityLinearOp<Scalar>::DefaultIdentityLinearOp()
{}


template<class Scalar>
DefaultIdentityLinearOp<Scalar>::DefaultIdentityLinearOp(
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >   &space
  )
{
  initialize(space);
}


template<class Scalar>
void DefaultIdentityLinearOp<Scalar>::initialize(
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >   &space
  )
{
  space_ = space.assert_not_null();
}


template<class Scalar>
void DefaultIdentityLinearOp<Scalar>::uninitialize()
{
  space_ = Teuchos::null;
}


// Overridden from LinearOpBase

  
template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultIdentityLinearOp<Scalar>::range() const
{
  return space_;
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultIdentityLinearOp<Scalar>::domain() const
{
  return space_;
}

  
template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultIdentityLinearOp<Scalar>::clone() const
{
  typedef DefaultIdentityLinearOp<Scalar> this_t;
  if(space_.get())
    return Teuchos::rcp(new this_t(space_));
  return Teuchos::rcp(new this_t());
}

  
// Overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultIdentityLinearOp<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  std::ostringstream oss;
  oss
    << "Thyra::DefaultIdentityLinearOp<" << ST::name() << ">{"
    << "space="<<(space_.get()?space_->description():"NULL")
    << "}";
  return oss.str();
}


// protected


// Overridden from LinearOpBase


template<class Scalar>
bool DefaultIdentityLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  return true;
}


template<class Scalar>
void DefaultIdentityLinearOp<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  using Teuchos::tuple;
  using Teuchos::ptrFromRef;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultIdentityLinearOp<Scalar>::apply(...)", *this, M_trans, X, &*Y
    );
#endif // TEUCHOS_DEBUG  
  Thyra::linear_combination<Scalar>(
    tuple<Scalar>(alpha)(),
    tuple<Ptr<const MultiVectorBase<Scalar> > >(ptrFromRef(X))(),
    beta, Y
    );
}


}	// end namespace Thyra


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::identity(
  const Teuchos::RCP<const VectorSpaceBase<Scalar> > &space,
  const std::string &label
  )
{
  RCP<Thyra::LinearOpBase<Scalar> > ilo =
    Teuchos::rcp(new DefaultIdentityLinearOp<Scalar>(space));
  if (label.length())
    ilo->setObjectLabel(label);
  return ilo;
}


//
// Explicit instantaition
//


#define THYRA_DEFAULT_IDENTITY_LINEAR_OP_INSTANT(SCALAR) \
  \
  template class DefaultIdentityLinearOp<SCALAR >; \
  \
  template RCP<const LinearOpBase<SCALAR > > \
  identity( \
    const RCP<const VectorSpaceBase<SCALAR > > &space, \
    const std::string &label \
    ); \


#endif	// THYRA_DEFAULT_IDENTITY_LINEAR_OP_DEF_HPP
