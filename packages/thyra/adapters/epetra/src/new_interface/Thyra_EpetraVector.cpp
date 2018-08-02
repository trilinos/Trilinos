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

#include "Thyra_EpetraVectorSpace.hpp"
#include "Thyra_EpetraVector.hpp"
#include "Thyra_EpetraMultiVector.hpp"

#include "Epetra_Comm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

namespace Thyra {


// Constructors/initializers/accessors


EpetraVector::EpetraVector()
{}


void EpetraVector::initialize(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<Epetra_Vector> &epetraVector
  )
{
  initializeImpl(epetraVectorSpace, epetraVector);
}


void EpetraVector::constInitialize(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<const Epetra_Vector> &epetraVector
  )
{
  initializeImpl(epetraVectorSpace, epetraVector);
}


RCP<Epetra_Vector>
EpetraVector::getEpetraVector()
{
  return epetraVector_.getNonconstObj();
}


RCP<const Epetra_Vector>
EpetraVector::getConstEpetraVector() const
{
  return epetraVector_;
}


// Overridden from VectorDefaultBase
RCP<const VectorSpaceBase<double>>
EpetraVector::domain() const
{
  if (domainSpace_.is_null()) {
    domainSpace_ = epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(1,0,epetraVector_.getConstObj()->Map().Comm()) ));
  }
  return domainSpace_;
}


// Overridden from SpmdMultiVectorBase


RCP<const SpmdVectorSpaceBase<double>>
EpetraVector::spmdSpaceImpl() const
{
  return epetraVectorSpace_;
}


// Overridden from SpmdVectorBase


void EpetraVector::getNonconstLocalVectorDataImpl(
  const Ptr<ArrayRCP<double>> &localValues )
{
  Epetra_MultiVector& emv = *epetraVector_.getNonconstObj();
  double* data = emv[0];
  *localValues = Teuchos::arcp(data,0,epetraVector_.getConstObj()->MyLength(),false);

  // The arcp is not owning the data. This could potentially lead to dangling pointers,
  // in case the epetraVector is deleted *before* the arcp. To prevent this, we set
  // the epetraVector rcp as extra data for the arcp, so that epetraVector will be
  // destroyed *after* the arcp. Notice that there is no cyclic dependency, since
  // the arcp is not stored in the epetra vector.
  Teuchos::set_extra_data(epetraVector_,"epetra vector",localValues,Teuchos::POST_DESTROY);
}


void EpetraVector::getLocalVectorDataImpl(
  const Ptr<ArrayRCP<const double>> &localValues ) const
{
  const Epetra_MultiVector& emv = *epetraVector_.getConstObj();
  const double* data = emv[0];
  *localValues = Teuchos::arcp(data,0,epetraVector_->MyLength(),false);

  // The arcp is not owning the data. This could potentially lead to dangling pointers,
  // in case the epetraVector is deleted *before* the arcp. To prevent this, we set
  // the epetraVector rcp as extra data for the arcp, so that epetraVector will be
  // destroyed *after* the arcp. Notice that there is no cyclic dependency, since
  // the arcp is not stored in the epetra vector.
  Teuchos::set_extra_data(epetraVector_,"epetra vector",localValues,Teuchos::POST_DESTROY);
}


// Overridden from VectorBase


void EpetraVector::randomizeImpl(
  double l,
  double u
  )
{
  // Epetra may randomize with different seed for each proc, so need a global
  // reduction to get locally-replicated random vector same on each proc.
  // Additionally, Epetra only offers randomization in [-1,1]. Therefore,
  // set v = (l+u)/2 + (u-l)/2 * rand_vec
  Epetra_Vector ev(epetraVector_.getConstObj()->Map());
  if (!ev.DistributedGlobal()) {
    if (ev.Comm().MyPID() == 0) {
      ev.Random();
    } else {
      ev.PutScalar(0.0);
    }
    ev.Reduce();
  } else {
    ev.Random();
  }
  epetraVector_.getNonconstObj()->PutScalar(0.5*(l+u));
  epetraVector_.getNonconstObj()->Update(0.5*(u-l),ev,1.0);
}


void EpetraVector::absImpl(
  const VectorBase<double>& x
  )
{
  auto ex = this->getConstEpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (Teuchos::nonnull(ex)) {
    epetraVector_.getNonconstObj()->Abs(*ex);
  } else {
    VectorDefaultBase<double>::absImpl(x);
  }
}


void EpetraVector::reciprocalImpl(
  const VectorBase<double>& x
  )
{
  auto ex = this->getConstEpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (Teuchos::nonnull(ex)) {
    epetraVector_.getNonconstObj()->Reciprocal(*ex);
  } else {
    VectorDefaultBase<double>::reciprocalImpl(x);
  }
}


void EpetraVector::eleWiseScaleImpl(
  const VectorBase<double>& x
  )
{
  auto ex = this->getConstEpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (Teuchos::nonnull(ex)) {
    typedef Teuchos::ScalarTraits<double> ST;
    epetraVector_.getNonconstObj()->Multiply(
      ST::one(), *ex, *epetraVector_.getConstObj(), ST::zero());
  } else {
    VectorDefaultBase<double>::eleWiseScaleImpl(x);
  }
}


typename Teuchos::ScalarTraits<double>::magnitudeType
EpetraVector::norm2WeightedImpl(
  const VectorBase<double>& x
  ) const
{
  auto ex = this->getConstEpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (Teuchos::nonnull(ex)) {
    // Weighted 2-norm function for Epetra vector does something different
    // from what Thyra wants: Thyra expects sqrt(sum(w[i]*v[i]*v[i])),
    // while Epetra does sqrt(sum((v[i]/w[i])^2)/n).
    // To get the result Thyra wants, we do element-wise multiplication,
    // then call norm_2 on the result
    typedef Teuchos::ScalarTraits<double> ST;

    // Let tmp = (*this) * x (element-wise)
    Epetra_Vector tmp(epetraVector_.getConstObj()->Map());
    tmp.Multiply(ST::one(),*ex,*epetraVector_.getConstObj(),ST::zero());

    // Compute the the dot product (*this, tmp)
    ST::magnitudeType result;
    tmp.Dot(*epetraVector_.getConstObj(),&result);
    return ST::squareroot(result);
  } else {
    return VectorDefaultBase<double>::norm2WeightedImpl(x);
  }
}


void EpetraVector::applyOpImpl(
  const RTOpPack::RTOpT<double> &op,
  const ArrayView<const Ptr<const VectorBase<double>>> &vecs,
  const ArrayView<const Ptr<VectorBase<double>>> &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset
  ) const
{
  SpmdVectorDefaultBase<double>::applyOpImpl(op, vecs, targ_vecs, reduct_obj, global_offset);
}


void EpetraVector::
acquireDetachedVectorViewImpl(
  const Range1D& rng,
  RTOpPack::ConstSubVectorView<double>* sub_vec
  ) const
{
  SpmdVectorDefaultBase<double>::acquireDetachedVectorViewImpl(rng, sub_vec);
}


void EpetraVector::
acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng,
  RTOpPack::SubVectorView<double>* sub_vec 
  )
{
  SpmdVectorDefaultBase<double>::acquireNonconstDetachedVectorViewImpl(rng, sub_vec);
}


void EpetraVector::
commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<double>* sub_vec
  )
{
  SpmdVectorDefaultBase<double>::commitNonconstDetachedVectorViewImpl(sub_vec);
}


// Overridden protected functions from MultiVectorBase


void EpetraVector::assignImpl(double alpha)
{
  epetraVector_.getNonconstObj()->PutScalar(alpha);
}


void EpetraVector::
assignMultiVecImpl(const MultiVectorBase<double>& mv)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(mv.domain()->dim(), 1);
#endif
  auto emv = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (Teuchos::nonnull(emv)) {
    epetraVector_.getNonconstObj()->Scale(1.0,*emv);
  } else {
    MultiVectorDefaultBase<double>::assignMultiVecImpl(mv);
  }
}


void EpetraVector::scaleImpl(double alpha)
{
  epetraVector_.getNonconstObj()->Scale(alpha);
}


void EpetraVector::updateImpl(
  double alpha,
  const MultiVectorBase<double>& mv
  )
{
  auto emv = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  typedef Teuchos::ScalarTraits<double> ST;
  if (Teuchos::nonnull(emv)) {
    epetraVector_.getNonconstObj()->Update(alpha, *emv, ST::one());
  } else {
    MultiVectorDefaultBase<double>::updateImpl(alpha, mv);
  }
}


void EpetraVector::linearCombinationImpl(
  const ArrayView<const double>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<double>>>& mv,
  const double& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
#endif

  // Try to cast mv to Epetra objects
  Teuchos::Array<RCP<const Epetra_MultiVector>> emvs(mv.size());
  RCP<const Epetra_MultiVector> emv;
  bool allCastsSuccessful = true;
  {
    auto mvIter = mv.begin();
    auto emvIter = emvs.begin();
    for (; mvIter != mv.end(); ++mvIter, ++emvIter) {
      emv = this->getConstEpetraMultiVector(Teuchos::rcpFromPtr(*mvIter));
      if (Teuchos::nonnull(emv)) {
        *emvIter = emv;
      } else {
        allCastsSuccessful = false;
        break;
      }
      TEUCHOS_ASSERT_EQUALITY(emv->NumVectors(), 1);
    }
  }

  // If casts succeeded, or input arrays are size 0, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  auto len = mv.size();
  if (len == 0) {
    epetraVector_.getNonconstObj()->Scale(beta);
  } else if (len == 1 && allCastsSuccessful) {
    epetraVector_.getNonconstObj()->Update(alpha[0], *emvs[0], beta);
  } else if (len == 2 && allCastsSuccessful) {
    epetraVector_.getNonconstObj()->Update(alpha[0], *emvs[0], alpha[1], *emvs[1], beta);
  } else if (allCastsSuccessful) {
    typedef Teuchos::ScalarTraits<double> ST;
    auto emvIter = emvs.begin();
    auto alphaIter = alpha.begin();

    // Check if any entry of emvs aliases this object's wrapped vector.
    // If so, replace that entry in the array with a copy.
    emv = Teuchos::null;
    for (; emvIter != emvs.end(); ++emvIter) {
      if (emvIter->getRawPtr() == epetraVector_.getConstObj().getRawPtr()) {
        if (emv.is_null()) {
          emv = Teuchos::rcp(new Epetra_MultiVector(*epetraVector_.getConstObj()));
        }
        *emvIter = emv;
      }
    }
    emvIter = emvs.begin();

    // We add two MVs at a time, so only scale if even num MVs,
    // and additionally do the first addition if odd num MVs.
    if ((emvs.size() % 2) == 0) {
      epetraVector_.getNonconstObj()->Scale(beta);
    } else {
      epetraVector_.getNonconstObj()->Update(*alphaIter, *(*emvIter), beta);
      ++emvIter;
      ++alphaIter;
    }
    for (; emvIter != emvs.end(); emvIter+=2, alphaIter+=2) {
      epetraVector_.getNonconstObj()->Update(
        *alphaIter, *(*emvIter), *(alphaIter+1), *(*(emvIter+1)), ST::one());
    }
  } else {
    MultiVectorDefaultBase<double>::linearCombinationImpl(alpha, mv, beta);
  }
}


void EpetraVector::dotsImpl(
    const MultiVectorBase<double>& mv,
    const ArrayView<double>& prods
    ) const
{
  auto emv = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (Teuchos::nonnull(emv)) {
    epetraVector_.getConstObj()->Dot(*emv, prods.getRawPtr());
  } else {
    MultiVectorDefaultBase<double>::dotsImpl(mv, prods);
  }
}


void EpetraVector::norms1Impl(
  const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
  ) const
{
  epetraVector_.getConstObj()->Norm1(norms.getRawPtr());
}


void EpetraVector::norms2Impl(
    const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
    ) const
{
  epetraVector_.getConstObj()->Norm2(norms.getRawPtr());
}


void EpetraVector::normsInfImpl(
  const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
  ) const
{
  epetraVector_.getConstObj()->NormInf(norms.getRawPtr());
}


void EpetraVector::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X,
  const Ptr<MultiVectorBase<double>> &Y,
  const double alpha,
  const double beta
  ) const
{
  // Try to extract Epetra objects from X and Y
  Teuchos::RCP<const Epetra_MultiVector> X_epetra = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(X));
  Teuchos::RCP<Epetra_MultiVector> Y_epetra = this->getEpetraMultiVector(Teuchos::rcpFromPtr(Y));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the default implementation.
  if (Teuchos::nonnull(X_epetra) && Teuchos::nonnull(Y_epetra)) {
    auto real_M_trans = real_trans(M_trans);
    char trans;
    switch (real_M_trans) {
      case NOTRANS:
        trans = 'N';
        break;
      case TRANS:
        trans = 'T';
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! Unexpected real transpose enum value.\n");
        break;
    }

    Y_epetra->Multiply(trans, 'N', alpha, *epetraVector_.getConstObj(), *X_epetra, beta);
  } else {
    VectorDefaultBase<double>::applyImpl(M_trans, X, Y, alpha, beta);
  }

}


// private


template<class EpetraVector_t>
void EpetraVector::initializeImpl(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<EpetraVector_t> &epetraVector
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(Teuchos::nonnull(epetraVectorSpace));
  TEUCHOS_ASSERT(Teuchos::nonnull(epetraVector));
#endif
  epetraVectorSpace_ = Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(epetraVectorSpace);
  TEUCHOS_TEST_FOR_EXCEPTION(epetraVectorSpace_.is_null(), std::runtime_error,
                             "Error! Could not cast input vector space to EpetraVectorSpace.\n");
  TEUCHOS_TEST_FOR_EXCEPTION(!epetraVectorSpace_->getEpetraMap()->SameAs(epetraVector->Map()),
                             std::logic_error,
                             "Error! Input multivector has a map that is different "
                             "from that stored in the input EpetraVectorSpace.\n");

  epetraVector_.initialize(epetraVector);
  this->updateSpmdSpace();
}


RCP<Epetra_MultiVector>
EpetraVector::
getEpetraMultiVector(const RCP<MultiVectorBase<double>>& mv) const
{
  using Teuchos::rcp_dynamic_cast;

  RCP<EpetraMultiVector> emv = rcp_dynamic_cast<EpetraMultiVector>(mv);
  if (Teuchos::nonnull(emv)) {
    return emv->getEpetraMultiVector();
  }

  RCP<EpetraVector> ev = rcp_dynamic_cast<EpetraVector>(mv);
  if (Teuchos::nonnull(ev)) {
    return ev->getEpetraVector();
  }

  return Teuchos::null;
}


RCP<const Epetra_MultiVector>
EpetraVector::
getConstEpetraMultiVector(const RCP<const MultiVectorBase<double>>& mv) const
{
  using Teuchos::rcp_dynamic_cast;

  RCP<const EpetraMultiVector> emv = rcp_dynamic_cast<const EpetraMultiVector>(mv);
  if (Teuchos::nonnull(emv)) {
    return emv->getConstEpetraMultiVector();
  }

  RCP<const EpetraVector> ev = rcp_dynamic_cast<const EpetraVector>(mv);
  if (Teuchos::nonnull(ev)) {
    return ev->getConstEpetraVector();
  }

  return Teuchos::null;
}


RCP<Epetra_Vector>
EpetraVector::
getEpetraVector(const RCP<VectorBase<double>>& v) const
{
  RCP<EpetraVector> ev = Teuchos::rcp_dynamic_cast<EpetraVector>(v);
  if (Teuchos::nonnull(ev)) {
    return ev->getEpetraVector();
  } else {
    return Teuchos::null;
  }
}


RCP<const Epetra_Vector>
EpetraVector::
getConstEpetraVector(const RCP<const VectorBase<double>>& v) const
{
  RCP<const EpetraVector> ev = Teuchos::rcp_dynamic_cast<const EpetraVector>(v);
  if (Teuchos::nonnull(ev)) {
    return ev->getConstEpetraVector();
  } else {
    return Teuchos::null;
  }
}


} // end namespace Thyra
