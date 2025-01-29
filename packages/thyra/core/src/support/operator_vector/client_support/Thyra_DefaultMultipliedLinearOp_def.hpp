// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_MULTIPLIED_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_MULTIPLIED_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultMultipliedLinearOp_decl.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// Inline members only used in implementation


template<class Scalar>
inline
void DefaultMultipliedLinearOp<Scalar>::assertInitialized() const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( numOps() > 0 ) );
#endif
}


template<class Scalar>
inline
std::string DefaultMultipliedLinearOp<Scalar>::getClassName() const
{
  return Teuchos::Describable::description();
}


template<class Scalar>
inline
Ordinal DefaultMultipliedLinearOp<Scalar>::getRangeDim() const
{
  return (numOps() > 0 ? this->range()->dim() : 0);
}


template<class Scalar>
inline
Ordinal DefaultMultipliedLinearOp<Scalar>::getDomainDim() const
{
  return (numOps() > 0 ? this->domain()->dim() : 0);
}


// Constructors/initializers/accessors


template<class Scalar>
DefaultMultipliedLinearOp<Scalar>::DefaultMultipliedLinearOp()
{}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::initialize(
  const ArrayView<const RCP<LinearOpBase<Scalar> > > &Ops )
{
  const int nOps = Ops.size();
  Ops_.resize(nOps);
  for( int k = 0; k < nOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::initialize(
  const ArrayView<const RCP<const LinearOpBase<Scalar> > > &Ops )
{
  const int nOps = Ops.size();
  Ops_.resize(nOps);
  for( int k = 0; k < nOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::uninitialize()
{
  Ops_.resize(0);
  setupDefaultObjectLabel();
}


// Overridden form MultipliedLinearOpBase


template<class Scalar>
int DefaultMultipliedLinearOp<Scalar>::numOps() const
{
  return Ops_.size();
}


template<class Scalar>
bool DefaultMultipliedLinearOp<Scalar>::opIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].isConst();
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::getNonconstOp(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].getNonconstObj();
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::getOp(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k];
}


// Overridden from LinearOpBase


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::range() const
{
  if (numOps()) {
    return getOp(0)->range();
  }
  return Teuchos::null;
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::domain() const
{
  if (numOps()) {
    return getOp(numOps()-1)->domain();
  }
  return Teuchos::null;
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultMultipliedLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be!
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultMultipliedLinearOp<Scalar>::description() const
{
  std::ostringstream oss;
  oss << getClassName() << "{numOps="<<numOps()
      <<",rangeDim=" << getRangeDim()
      << ",domainDim="<< getDomainDim() <<"}";
  return oss.str();
}

template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  const int nOps = Ops_.size();
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out << this->description() << std::endl;
      OSTab tab2(out);
      *out
        <<  "Constituent LinearOpBase objects for M = Op[0]*...*Op[numOps-1]:\n";
      OSTab tab3(out);
      for( int k = 0; k < nOps; ++k ) {
        *out << "Op["<<k<<"] = " << Teuchos::describe(*getOp(k),verbLevel);
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// protected


// Overridden from LinearOpBase


template<class Scalar>
bool DefaultMultipliedLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  bool overallOpSupported = true;
  for( int k = 0; k < static_cast<int>(Ops_.size()); ++k )
    if(!Thyra::opSupported(*getOp(k),M_trans)) overallOpSupported = false;
  return overallOpSupported;
  // ToDo: Cache these?
}

template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::allocateVecs(const Ordinal dim) const {
  const int nOps = Ops_.size();
  if ((T_k_.size() != Teuchos::as<size_t>(nOps+1)) || ((nOps > 0) && (T_k_[0]->domain()->dim() != dim))) {
    // op[0]->range
    // op[0]->domain == op[1]->range
    // ...
    // op[nOps-2]->domain == op[nOps-1]->range
    // op[nOps-1]->domain
    T_k_.resize(0);
    for( int k = 0; k < nOps; ++k ) {
      T_k_.push_back(createMembers(getOp(k)->range(), dim));
    }
    T_k_.push_back(createMembers(getOp(nOps-1)->domain(), dim));
  }
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  using Teuchos::rcpFromPtr;
  using Teuchos::rcpFromRef;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    getClassName()+"::apply(...)", *this, M_trans, X, &*Y
    );
#endif // TEUCHOS_DEBUG  
  const int nOps = Ops_.size();
  const Ordinal m = X.domain()->dim();
  allocateVecs(m);
  if( real_trans(M_trans)==NOTRANS ) {
    //
    // Y = alpha * M * X + beta*Y
    // =>
    // Y = alpha * op(Op[0]) * op(Op[1]) * ... * op(Op[numOps-1]) * X + beta*Y
    //
    RCP<MultiVectorBase<Scalar> > T_kp1, T_k; // Temporary propagated between loops 
    for( int k = nOps-1; k >= 0; --k ) {
      RCP<MultiVectorBase<Scalar> > Y_k;
      RCP<const MultiVectorBase<Scalar> > X_k;
      if(k==0) Y_k = rcpFromPtr(Y);  else Y_k = T_k = T_k_[k];
      if(k==nOps-1) X_k = rcpFromRef(X); else X_k = T_kp1;
      if( k > 0 )
        Thyra::apply(*getOp(k), M_trans, *X_k, Y_k.ptr());
      else
        Thyra::apply(*getOp(k), M_trans, *X_k, Y_k.ptr(), alpha, beta);
      T_kp1 = T_k;
    }
  }
  else {
    //
    // Y = alpha * M' * X + beta*Y
    // =>
    // Y = alpha * Op[numOps-1]' * Op[1]' * ... * Op[0]' * X + beta * Y
    //
    RCP<MultiVectorBase<Scalar> > T_km1, T_k; // Temporary propagated between loops 
    for( int k = 0; k <= nOps-1; ++k ) {
      RCP<MultiVectorBase<Scalar> >         Y_k;
      RCP<const MultiVectorBase<Scalar> >   X_k;
      if(k==nOps-1) Y_k = rcpFromPtr(Y);  else Y_k = T_k = T_k_[k+1];
      if(k==0) X_k = rcpFromRef(X); else X_k = T_km1;
      if( k < nOps-1 )
        Thyra::apply(*getOp(k), M_trans, *X_k, Y_k.ptr());
      else
        Thyra::apply(*getOp(k), M_trans, *X_k, Y_k.ptr(), alpha, beta);
      T_km1 = T_k;
    }
  }
}


// private


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::validateOps()
{
  using Teuchos::toString;
#ifdef TEUCHOS_DEBUG
  try {
    const int nOps = Ops_.size();
    for( int k = 0; k < nOps; ++k ) {
      TEUCHOS_TEST_FOR_EXCEPT( Ops_[k]().get() == NULL );
      if( k < nOps-1 ) {
        THYRA_ASSERT_LINEAR_OP_TIMES_LINEAR_OP_SPACES_NAMES(
          getClassName()+"::initialize(...)"
          ,*Ops_[k],NOTRANS,("Ops["+toString(k)+"]")
          ,*Ops_[k+1],NOTRANS,("Ops["+toString(k+1)+"]")
          );
      }
    }
  }
  catch(...) {
    uninitialize();
    throw;
  }
#endif
}


template<class Scalar>
void DefaultMultipliedLinearOp<Scalar>::setupDefaultObjectLabel()
{
  std::ostringstream label;
  const int nOps = Ops_.size();
  for( int k = 0; k < nOps; ++k ) {
    std::string Op_k_label = Ops_[k]->getObjectLabel();
    if (Op_k_label.length() == 0)
      Op_k_label = "ANYM";
    if (k > 0)
      label << "*";
    label << "("<<Op_k_label<<")";
  }
  this->setObjectLabel(label.str());
  validateOps();
}


}	// end namespace Thyra


//
// Nonmember implementations
//

template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstMultiply(
  const RCP<LinearOpBase<Scalar> > &A,
  const RCP<LinearOpBase<Scalar> > &B,
  const std::string &M_label
  )
{
  using Teuchos::tuple;
  RCP<DefaultMultipliedLinearOp<Scalar> > multOp = 
    defaultMultipliedLinearOp<Scalar>(tuple(A, B)());
  if(M_label.length())
    multOp->setObjectLabel(M_label);
  return multOp;
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::multiply(
  const RCP<const LinearOpBase<Scalar> > &A,
  const RCP<const LinearOpBase<Scalar> > &B,
  const std::string &M_label
  )
{
  using Teuchos::tuple;
  RCP<DefaultMultipliedLinearOp<Scalar> > multOp =
    defaultMultipliedLinearOp<Scalar>(tuple(A, B)());
  if(M_label.length())
    multOp->setObjectLabel(M_label);
  return multOp;
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::multiply(
  const RCP<const LinearOpBase<Scalar> > &A,
  const RCP<const LinearOpBase<Scalar> > &B,
  const RCP<const LinearOpBase<Scalar> > &C,
  const std::string &M_label
  )
{
  using Teuchos::tuple;
  RCP<DefaultMultipliedLinearOp<Scalar> > multOp =
    defaultMultipliedLinearOp<Scalar>(tuple(A, B, C)());
  if(M_label.length())
    multOp->setObjectLabel(M_label);
  return multOp;
}


//
// Explicit instantiation macro
//
// Must be expanded from within the Thyra namespace!
//


#define THYRA_DEFAULT_MULTIPLIED_LINEAR_OP_INSTANT(SCALAR) \
  \
  template class DefaultMultipliedLinearOp<SCALAR >; \
 \
  template RCP<LinearOpBase<SCALAR > > \
  nonconstMultiply( \
    const RCP<LinearOpBase<SCALAR > > &A, \
    const RCP<LinearOpBase<SCALAR > > &B, \
    const std::string &M_label \
    ); \
   \
  template RCP<const LinearOpBase<SCALAR > > \
  multiply( \
    const RCP<const LinearOpBase<SCALAR > > &A, \
    const RCP<const LinearOpBase<SCALAR > > &B, \
    const std::string &M_label \
    ); \
   \
  template RCP<const LinearOpBase<SCALAR > > \
  multiply( \
    const RCP<const LinearOpBase<SCALAR > > &A, \
    const RCP<const LinearOpBase<SCALAR > > &B, \
    const RCP<const LinearOpBase<SCALAR > > &C, \
    const std::string &M_label \
    ); \


#endif	// THYRA_DEFAULT_MULTIPLIED_LINEAR_OP_DEF_HPP
