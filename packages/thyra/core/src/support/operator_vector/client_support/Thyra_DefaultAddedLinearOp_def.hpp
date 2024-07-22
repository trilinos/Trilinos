// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_ADDED_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_ADDED_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultAddedLinearOp_decl.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"


namespace Thyra {


// Inline members only used in this class impl


template<class Scalar>
inline
void DefaultAddedLinearOp<Scalar>::assertInitialized() const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( numOps() > 0 ) );
#endif
}


template<class Scalar>
inline
std::string DefaultAddedLinearOp<Scalar>::getClassName() const
{
  return Teuchos::Describable::description();
}


template<class Scalar>
inline
Ordinal DefaultAddedLinearOp<Scalar>::getRangeDim() const
{
  return (numOps() > 0 ? this->range()->dim() : 0);
}


template<class Scalar>
inline
Ordinal DefaultAddedLinearOp<Scalar>::getDomainDim() const
{
  return (numOps() > 0 ? this->domain()->dim() : 0);
}


// Constructors/initializers/accessors


template<class Scalar>
DefaultAddedLinearOp<Scalar>::DefaultAddedLinearOp()
{}


template<class Scalar>
DefaultAddedLinearOp<Scalar>::DefaultAddedLinearOp(
  const ArrayView<const RCP<LinearOpBase<Scalar> > > &Ops )
{
  initialize(Ops);
}


template<class Scalar>
DefaultAddedLinearOp<Scalar>::DefaultAddedLinearOp(
  const ArrayView<const RCP<const LinearOpBase<Scalar> > > &Ops )
{
  initialize(Ops);
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::initialize(
  const ArrayView<const RCP<LinearOpBase<Scalar> > > &Ops )
{
  const int l_numOps = Ops.size();
  Ops_.resize(l_numOps);
  for( int k = 0; k < l_numOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::initialize(
  const ArrayView<const RCP<const LinearOpBase<Scalar> > > &Ops )
{
  const int l_numOps = Ops.size();
  Ops_.resize(l_numOps);
  for( int k = 0; k < l_numOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::uninitialize()
{
  Ops_.resize(0);
  setupDefaultObjectLabel();
}


// Overridden form AddedLinearOpBase


template<class Scalar>
int DefaultAddedLinearOp<Scalar>::numOps() const
{
  return Ops_.size();
}


template<class Scalar>
bool DefaultAddedLinearOp<Scalar>::opIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].isConst();
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultAddedLinearOp<Scalar>::getNonconstOp(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].getNonconstObj();
}


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultAddedLinearOp<Scalar>::getOp(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k];
}


// Overridden from LinearOpBase


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultAddedLinearOp<Scalar>::range() const
{
  if (numOps()) {
    return getOp(0)->range();
  }
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultAddedLinearOp<Scalar>::domain() const
{
  if (numOps()) {
    return getOp(numOps()-1)->domain();
  }
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultAddedLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be!
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultAddedLinearOp<Scalar>::description() const
{
  std::ostringstream oss;
  oss << getClassName() << "{numOps="<<numOps()
      <<",rangeDim=" << getRangeDim()
      << ",domainDim="<< getDomainDim() <<"}";
  return oss.str();
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  const int l_numOps = Ops_.size();
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
      tab.incrTab();
      for( int k = 0; k < l_numOps; ++k ) {
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
bool DefaultAddedLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  bool isOpSupported = true;
  for( int k = 0; k < static_cast<int>(Ops_.size()); ++k )
    if(!Thyra::opSupported(*getOp(k),M_trans)) isOpSupported = false;
  return isOpSupported;
  // ToDo: Cache these?
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    getClassName()+"::apply(...)", *this, M_trans, X, &*Y
    );
#endif // TEUCHOS_DEBUG  
  //
  // Y = alpha * op(M) * X + beta*Y
  //
  // =>
  //
  // Y = beta*Y + sum(alpha*op(Op[j])*X),j=0...numOps-1)
  //
  const int l_numOps = Ops_.size();
  for( int j = 0; j < l_numOps; ++j )
    Thyra::apply(*getOp(j), M_trans, X, Y, alpha, j==0 ? beta : ST::one());
}


// private


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::validateOps()
{
  using Teuchos::toString;
#ifdef TEUCHOS_DEBUG
  try {
    const int l_numOps = Ops_.size();
    for( int k = 0; k < l_numOps; ++k ) {
      TEUCHOS_TEST_FOR_EXCEPT( Ops_[k]().get() == NULL );
      if( k > 0 ) {
        THYRA_ASSERT_LINEAR_OP_PLUS_LINEAR_OP_SPACES_NAMES(
          getClassName()+"::initialize(...)"
          ,*Ops_[0], NOTRANS, ("Ops[0]")
          ,*Ops_[k], NOTRANS, ("Ops["+toString(k)+"]")
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
void DefaultAddedLinearOp<Scalar>::setupDefaultObjectLabel()
{
  std::ostringstream label;
  const int l_numOps = Ops_.size();
  for( int k = 0; k < l_numOps; ++k ) {
    std::string Op_k_label = Ops_[k]->getObjectLabel();
    if (Op_k_label.length() == 0)
      Op_k_label = "ANYM";
    if (k > 0)
      label << "+";
    label << "("<<Op_k_label<<")";
  }
  this->setObjectLabel(label.str());
  validateOps();
}


}	// end namespace Thyra


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstAdd(
  const RCP<LinearOpBase<Scalar> > &A,
  const RCP<LinearOpBase<Scalar> > &B,
  const std::string &label
  )
{
  using Teuchos::tuple;
  RCP<LinearOpBase<Scalar> > alo =
    defaultAddedLinearOp<Scalar>(
      tuple<RCP<LinearOpBase<Scalar> > >(A, B)()
      );
  if (label.length())
    alo->setObjectLabel(label);
  return alo;
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::add(
  const RCP<const LinearOpBase<Scalar> > &A,
  const RCP<const LinearOpBase<Scalar> > &B,
  const std::string &label
  )
{
  using Teuchos::tuple;
  RCP<LinearOpBase<Scalar> > alo =
    defaultAddedLinearOp<Scalar>(
      tuple<RCP<const LinearOpBase<Scalar> > >(A, B)()
      );
  if (label.length())
    alo->setObjectLabel(label);
  return alo;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstSubtract(
  const RCP<LinearOpBase<Scalar> > &A,
  const RCP<LinearOpBase<Scalar> > &B,
  const std::string &label
  )
{
  typedef ScalarTraits<Scalar> ST;
  using Teuchos::tuple;
  RCP<LinearOpBase<Scalar> > alo =
    defaultAddedLinearOp<Scalar>(
      tuple<RCP<LinearOpBase<Scalar> > >(
        A, nonconstScale<Scalar>(-ST::one(),B) )()
      );
  if (label.length())
    alo->setObjectLabel(label);
  return alo;
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::subtract(
  const RCP<const LinearOpBase<Scalar> > &A,
  const RCP<const LinearOpBase<Scalar> > &B,
  const std::string &label
  )
{
  typedef ScalarTraits<Scalar> ST;
  using Teuchos::tuple;
  RCP<LinearOpBase<Scalar> > alo =
    defaultAddedLinearOp<Scalar>(
      tuple<RCP<const LinearOpBase<Scalar> > >(
          A, scale<Scalar>(-ST::one(),B)
        )()
      );
  if (label.length())
    alo->setObjectLabel(label);
  return alo;
}


//
// Explicit instantiation macro
//

#define THYRA_DEFAULT_ADDED_LINEAR_OP_INSTANT(SCALAR) \
  template class DefaultAddedLinearOp<SCALAR >; \
  \
  template RCP<LinearOpBase<SCALAR > > \
  nonconstAdd( \
    const RCP<LinearOpBase<SCALAR > > &A, \
    const RCP<LinearOpBase<SCALAR > > &B, \
    const std::string &label \
    ); \
   \
  template RCP<const LinearOpBase<SCALAR > > \
  add( \
    const RCP<const LinearOpBase<SCALAR > > &A, \
    const RCP<const LinearOpBase<SCALAR > > &B, \
    const std::string &label \
    ); \
   \
  template RCP<LinearOpBase<SCALAR > > \
  nonconstSubtract( \
    const RCP<LinearOpBase<SCALAR > > &A, \
    const RCP<LinearOpBase<SCALAR > > &B, \
    const std::string &label \
    ); \
   \
  template RCP<const LinearOpBase<SCALAR > > \
  subtract( \
    const RCP<const LinearOpBase<SCALAR > > &A, \
    const RCP<const LinearOpBase<SCALAR > > &B, \
    const std::string &label \
    );


#endif	// THYRA_DEFAULT_ADDED_LINEAR_OP_DEF_HPP
