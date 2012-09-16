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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_ADDED_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_ADDED_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultAddedLinearOp_decl.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"


namespace Thyra {


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
  assertInitialized();
  return getOp(0)->range();
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultAddedLinearOp<Scalar>::domain() const
{
  assertInitialized();
  return getOp(numOps()-1)->domain();
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
  assertInitialized();
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{numOps = "<<numOps()<<"}";
  return oss.str();
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  assertInitialized();
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
      *out
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim="<< this->domain()->dim() << "}\n";
      OSTab tab2(out);
      *out
        <<  "numOps="<< l_numOps << std::endl
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
    "DefaultAddedLinearOp<Scalar>::apply(...)", *this, M_trans, X, &*Y
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
  typedef std::string s;
  using Teuchos::toString;
#ifdef TEUCHOS_DEBUG
  try {
    const int l_numOps = Ops_.size();
    for( int k = 0; k < l_numOps; ++k ) {
      TEUCHOS_TEST_FOR_EXCEPT( Ops_[k]().get() == NULL );
      if( k > 0 ) {
        THYRA_ASSERT_LINEAR_OP_PLUS_LINEAR_OP_SPACES_NAMES(
          "DefaultMultipliedLinearOp<Scalar>::initialize(...)"
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
