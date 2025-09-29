// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_IDENTITY_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_IDENTITY_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultIdentityLinearOp_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template <class Scalar>
DefaultIdentityLinearOp<Scalar>::DefaultIdentityLinearOp() {}

template <class Scalar>
DefaultIdentityLinearOp<Scalar>::DefaultIdentityLinearOp(
    const Teuchos::RCP<const VectorSpaceBase<Scalar> > &space) {
  initialize(space);
}

template <class Scalar>
void DefaultIdentityLinearOp<Scalar>::initialize(
    const Teuchos::RCP<const VectorSpaceBase<Scalar> > &space) {
  space_ = space.assert_not_null();
}

template <class Scalar>
void DefaultIdentityLinearOp<Scalar>::uninitialize() {
  space_ = Teuchos::null;
}

// Overridden from LinearOpBase

template <class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultIdentityLinearOp<Scalar>::range() const {
  return space_;
}

template <class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultIdentityLinearOp<Scalar>::domain() const {
  return space_;
}

template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultIdentityLinearOp<Scalar>::clone() const {
  typedef DefaultIdentityLinearOp<Scalar> this_t;
  if (space_.get())
    return Teuchos::rcp(new this_t(space_));
  return Teuchos::rcp(new this_t());
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string DefaultIdentityLinearOp<Scalar>::description() const {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::ostringstream oss;
  oss
      << "Thyra::DefaultIdentityLinearOp<" << ST::name() << ">{"
      << "space=" << (space_.get() ? space_->description() : "NULL")
      << "}";
  return oss.str();
}

// protected

// Overridden from LinearOpBase

template <class Scalar>
bool DefaultIdentityLinearOp<Scalar>::opSupportedImpl(EOpTransp /* M_trans */) const {
  return true;
}

template <class Scalar>
void DefaultIdentityLinearOp<Scalar>::applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta) const {
  using Teuchos::ptrFromRef;
  using Teuchos::tuple;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
      "DefaultIdentityLinearOp<Scalar>::apply(...)", *this, M_trans, X, &*Y);
#else
  (void)M_trans;
#endif  // TEUCHOS_DEBUG
  Thyra::linear_combination<Scalar>(
      tuple<Scalar>(alpha)(),
      tuple<Ptr<const MultiVectorBase<Scalar> > >(ptrFromRef(X))(),
      beta, Y);
}

}  // end namespace Thyra

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::identity(
    const Teuchos::RCP<const VectorSpaceBase<Scalar> > &space,
    const std::string &label) {
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
  template class DefaultIdentityLinearOp<SCALAR>;        \
                                                         \
  template RCP<const LinearOpBase<SCALAR> >              \
  identity(                                              \
      const RCP<const VectorSpaceBase<SCALAR> > &space,  \
      const std::string &label);

#endif  // THYRA_DEFAULT_IDENTITY_LINEAR_OP_DEF_HPP
