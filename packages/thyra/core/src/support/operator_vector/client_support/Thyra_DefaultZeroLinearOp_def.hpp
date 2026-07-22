// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_ZERO_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_ZERO_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultZeroLinearOp_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
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
bool DefaultZeroLinearOp<Scalar>::opSupportedImpl(EOpTransp /* M_trans */) const
{
  return true;
}


template<class Scalar>
void DefaultZeroLinearOp<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar /* alpha */,
  const Scalar beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultZeroLinearOp<Scalar>::apply(...)", *this, M_trans, X, &*Y
    );
#else
  (void)M_trans;
  (void)X;
#endif // TEUCHOS_DEBUG  
  scale(beta, Y);
}

template<class Scalar>
bool DefaultZeroLinearOp<Scalar>::
rowStatIsSupportedImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat) const
{
  if(   rowStat==RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM 
     || rowStat==RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM)
    return true;

  // else(   rowStat==RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM 
  //      || rowStat==RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM)
  
  // inverse of a zero diagonal is bad news, we won't allow it, return false
  return false;
}

template<class Scalar>
void DefaultZeroLinearOp<Scalar>::
getRowStatImpl(
    const RowStatLinearOpBaseUtils::ERowStat /* rowStat */, 
    const Teuchos::Ptr<VectorBase< Scalar> > &rowStatVec) const
{ 
  Thyra::put_scalar(Teuchos::ScalarTraits<Scalar>::zero(),rowStatVec); 
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

template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstZero(
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
  \
  template RCP<Thyra::LinearOpBase<SCALAR > >  \
  nonconstZero(  \
    const RCP<const VectorSpaceBase<SCALAR > > &range,  \
    const RCP<const VectorSpaceBase<SCALAR > > &domain  \
    );  \


#endif	// THYRA_DEFAULT_ZERO_LINEAR_OP_DEF_HPP
