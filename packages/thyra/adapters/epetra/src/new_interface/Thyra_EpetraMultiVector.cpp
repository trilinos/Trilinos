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

#include "Thyra_EpetraMultiVector.hpp"
#include "Thyra_EpetraVectorSpace.hpp"
#include "Thyra_EpetraVector.hpp"

#include "Teuchos_Assert.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

namespace Thyra {


// Constructors/initializers/accessors


EpetraMultiVector::EpetraMultiVector()
{}


void EpetraMultiVector::initialize(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
  const RCP<Epetra_MultiVector> &epetraMultiVector
  )
{
  initializeImpl(epetraVectorSpace, domainSpace, epetraMultiVector);
}


void EpetraMultiVector::constInitialize(
  const RCP<const VectorSpaceBase<double> > &epetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
  const RCP<const Epetra_MultiVector> &epetraMultiVector
  )
{
  initializeImpl(epetraVectorSpace, domainSpace, epetraMultiVector);
}


RCP<Epetra_MultiVector>
EpetraMultiVector::getEpetraMultiVector()
{
  return epetraMultiVector_.getNonconstObj();
}


RCP<const Epetra_MultiVector>
EpetraMultiVector::getConstEpetraMultiVector() const
{
  return epetraMultiVector_;
}


// Overridden public functions form MultiVectorAdapterBase


RCP< const ScalarProdVectorSpaceBase<double> >
EpetraMultiVector::domainScalarProdVecSpc() const
{
  return domainSpace_;
}


// Overridden protected functions from MultiVectorBase


void
EpetraMultiVector::assignImpl(double alpha)
{
  epetraMultiVector_.getNonconstObj()->PutScalar(alpha);
}


void EpetraMultiVector::
assignMultiVecImpl(const MultiVectorBase<double>& mv)
{
  auto emv = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(emv)) {
    epetraMultiVector_.getNonconstObj()->Scale(1.0,*emv);
  } else {
    MultiVectorDefaultBase<double>::assignMultiVecImpl(mv);
  }
}


void
EpetraMultiVector::scaleImpl(double alpha)
{
  epetraMultiVector_.getNonconstObj()->Scale(alpha);
}


void EpetraMultiVector::updateImpl(
  double alpha,
  const MultiVectorBase<double>& mv
  )
{
  auto emv = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(mv));
 
  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(emv)) {
    typedef Teuchos::ScalarTraits<double> ST;
    epetraMultiVector_.getNonconstObj()->Update(alpha, *emv, ST::one());
  } else {
    MultiVectorDefaultBase<double>::updateImpl(alpha, mv);
  }
}


void EpetraMultiVector::linearCombinationImpl(
  const ArrayView<const double>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<double> > >& mv,
  const double& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
#endif

  // Try to cast mv to an array of this type
  Teuchos::Array<RCP<const Epetra_MultiVector> > emvs(mv.size());
  RCP<const Epetra_MultiVector> emv;
  bool allCastsSuccessful = true;
  {
    auto mvIter = mv.begin();
    auto emvIter = emvs.begin();
    for (; mvIter != mv.end(); ++mvIter, ++emvIter) {
      emv = this->getConstEpetraMultiVector(Teuchos::rcpFromPtr(*mvIter));
      if (nonnull(emv)) {
        *emvIter = emv;
      } else {
        allCastsSuccessful = false;
        break;
      }
    }
  }

  // If casts succeeded, or input arrays are size 0, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  auto len = emvs.size();
  if (len == 0) {
    epetraMultiVector_.getNonconstObj()->Scale(beta);
  } else if (len == 1 && allCastsSuccessful) {
    epetraMultiVector_.getNonconstObj()->Update(alpha[0], *emvs[0], beta);
  } else if (len == 2 && allCastsSuccessful) {
    epetraMultiVector_.getNonconstObj()->Update(alpha[0], *emvs[0], alpha[1], *emvs[1], beta);
  } else if (allCastsSuccessful) {
    typedef Teuchos::ScalarTraits<double> ST;
    auto emvIter = emvs.begin();
    auto alphaIter = alpha.begin();

    // Check if any entry of emvs aliases this object's wrapped vector.
    // If so, replace that entry in the array with a copy.
    emv = Teuchos::null;
    for (; emvIter != emvs.end(); ++emvIter) {
      if (emvIter->getRawPtr() == epetraMultiVector_.getConstObj().getRawPtr()) {
        if (emv.is_null()) {
          emv = Teuchos::rcp(new Epetra_MultiVector(*epetraMultiVector_.getConstObj()));
        }
        *emvIter = emv;
      }
    }
    emvIter = emvs.begin();

    // We add two MVs at a time, so only scale if even num MVs,
    // and additionally do the first addition if odd num MVs.
    if ((emvs.size() % 2) == 0) {
      epetraMultiVector_.getNonconstObj()->Scale(beta);
    } else {
      epetraMultiVector_.getNonconstObj()->Update(*alphaIter, *(*emvIter), beta);
      ++emvIter;
      ++alphaIter;
    }
    for (; emvIter != emvs.end(); emvIter+=2, alphaIter+=2) {
      epetraMultiVector_.getNonconstObj()->Update(
        *alphaIter, *(*emvIter), *(alphaIter+1), *(*(emvIter+1)), ST::one());
    }
  } else {
    MultiVectorDefaultBase<double>::linearCombinationImpl(alpha, mv, beta);
  }
}


void EpetraMultiVector::dotsImpl(
    const MultiVectorBase<double>& mv,
    const ArrayView<double>& prods
    ) const
{
  auto emv = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(emv)) {
    epetraMultiVector_.getConstObj()->Dot(*emv, prods.getRawPtr());
  } else {
    MultiVectorDefaultBase<double>::dotsImpl(mv, prods);
  }
}


void EpetraMultiVector::norms1Impl(
  const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
  ) const
{
  epetraMultiVector_.getConstObj()->Norm1(norms.getRawPtr());
}


void EpetraMultiVector::norms2Impl(
    const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
    ) const
{
  epetraMultiVector_.getConstObj()->Norm2(norms.getRawPtr());
}


void EpetraMultiVector::normsInfImpl(
  const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
  ) const
{
  epetraMultiVector_.getConstObj()->NormInf(norms.getRawPtr());
}


RCP<const VectorBase<double> >
EpetraMultiVector::colImpl(Ordinal j) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, this->domain()->dim());
#endif
  const Epetra_Vector* col_j = (*epetraMultiVector_.getConstObj())(j);
  auto col_j_ptr = Teuchos::rcp(col_j,false); // not owning memory

  // Set the underlying epetraMultiVector_ as extra data for the col_j_ptr,
  // so that epetraMultiVector_ will NOT be deleted before col_j_ptr (which
  // would cause col_j_ptr to be a dangling pointer
  Teuchos::set_extra_data(epetraMultiVector_.getConstObj(),"the multi vector",Teuchos::inOutArg(col_j_ptr),Teuchos::POST_DESTROY,false);

  return constEpetraVector(epetraVectorSpace_, col_j_ptr);
}


RCP<VectorBase<double> >
EpetraMultiVector::nonconstColImpl(Ordinal j)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, this->domain()->dim());
#endif
  Epetra_Vector& col_j = *(*epetraMultiVector_.getNonconstObj())(j);
  auto col_j_ptr = Teuchos::rcpFromRef(col_j);

  // Set the underlying epetraMultiVector_ as extra data for the col_j_ptr,
  // so that epetraMultiVector_ will NOT be deleted before col_j_ptr (which
  // would cause col_j_ptr to be a dangling pointer
  Teuchos::set_extra_data(epetraMultiVector_.getConstObj(),"the multi vector",Teuchos::inOutArg(col_j_ptr),Teuchos::POST_DESTROY,false);

  return epetraVector(epetraVectorSpace_, col_j_ptr);
}


RCP<const MultiVectorBase<double> >
EpetraMultiVector::contigSubViewImpl(
  const Range1D& col_rng_in
  ) const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nEpetraMultiVector::subView(Range1D) const called!\n";
#endif
  const Range1D colRng = this->validateColRange(col_rng_in);

  RCP<const Epetra_MultiVector> epetraView;
  const Epetra_MultiVector& emv = *epetraMultiVector_.getConstObj();
  const Ordinal num_cols = colRng.size();
  if (emv.ConstantStride()) {
    double* start = emv[colRng.lbound()];
    epetraView = Teuchos::rcp(new Epetra_MultiVector(View,emv.Map(),start,emv.Stride(),num_cols));
  } else {
    Teuchos::Array<double*> ptrs(num_cols);
    int k = 0;
    for (Ordinal j=colRng.lbound(); j!=colRng.ubound(); ++j, ++k) {
      ptrs[k] = emv[j];
    }
    epetraView = Teuchos::rcp(new Epetra_MultiVector(View,emv.Map(),ptrs.getRawPtr(),num_cols));
  }

  const RCP<const ScalarProdVectorSpaceBase<double> > viewDomainSpace =
    epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(epetraView->NumVectors(),0,emv.Map().Comm()) ));

  return constEpetraMultiVector(
      epetraVectorSpace_,
      viewDomainSpace,
      epetraView
      );
}


RCP<MultiVectorBase<double> >
EpetraMultiVector::nonconstContigSubViewImpl(
  const Range1D& col_rng_in
  )
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nEpetraMultiVector::subView(Range1D) called!\n";
#endif
  const Range1D colRng = this->validateColRange(col_rng_in);

  RCP<Epetra_MultiVector> epetraView;
  Epetra_MultiVector& emv = *epetraMultiVector_.getNonconstObj();
  const Ordinal num_cols = colRng.size();
  if (emv.ConstantStride()) {
    double* start = emv[colRng.lbound()];
    epetraView = Teuchos::rcp(new Epetra_MultiVector(View,emv.Map(),start,emv.Stride(),num_cols));
  } else {
    Teuchos::Array<double*> ptrs(num_cols);
    int k = 0;
    for (Ordinal j=colRng.lbound(); j!=colRng.ubound(); ++j, ++k) {
      ptrs[k] = emv[j];
    }
    epetraView = Teuchos::rcp(new Epetra_MultiVector(View,emv.Map(),ptrs.getRawPtr(),num_cols));
  }

  const RCP<const ScalarProdVectorSpaceBase<double> > viewDomainSpace =
    epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(epetraView->NumVectors(),0,emv.Map().Comm()) ));

  return epetraMultiVector(
      epetraVectorSpace_,
      viewDomainSpace,
      epetraView
      );
}


RCP<const MultiVectorBase<double> >
EpetraMultiVector::nonContigSubViewImpl(
  const ArrayView<const int>& cols_in
  ) const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nEpetraMultiVector::subView(ArrayView) const called!\n";
#endif
  const Epetra_MultiVector& emv = *epetraMultiVector_.getConstObj();

  const int num_cols = cols_in.size();
  Teuchos::Array<double*> ptrs(num_cols);
  for (int k=0; k<num_cols; ++k) {
    ptrs[k] = emv[cols_in[k]];
  }
  RCP<const Epetra_MultiVector> epetraView = Teuchos::rcp(new Epetra_MultiVector(View,emv.Map(),ptrs.getRawPtr(),num_cols));

  const RCP<const ScalarProdVectorSpaceBase<double> > viewDomainSpace =
    epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(epetraView->NumVectors(),0,emv.Map().Comm()) ));

  return constEpetraMultiVector(
      epetraVectorSpace_,
      viewDomainSpace,
      epetraView
      );
}


RCP<MultiVectorBase<double> >
EpetraMultiVector::nonconstNonContigSubViewImpl(
  const ArrayView<const int>& cols_in
  )
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nEpetraMultiVector::subView(ArrayView) called!\n";
#endif
  Epetra_MultiVector& emv = *epetraMultiVector_.getNonconstObj();

  const int num_cols = cols_in.size();
  Teuchos::Array<double*> ptrs(num_cols);
  for (int k=0; k<num_cols; ++k) {
    ptrs[k] = emv[cols_in[k]];
  }
  RCP<Epetra_MultiVector> epetraView = Teuchos::rcp(new Epetra_MultiVector(View,emv.Map(),ptrs.getRawPtr(),num_cols));

  const RCP<const ScalarProdVectorSpaceBase<double> > viewDomainSpace =
    epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(epetraView->NumVectors(),0,emv.Map().Comm()) ));

  return epetraMultiVector(
      epetraVectorSpace_,
      viewDomainSpace,
      epetraView
      );
}

//////////////////////////////////////// CONTINUE HERE /////////////////////////////////////////
void EpetraMultiVector::
mvMultiReductApplyOpImpl(
  const RTOpPack::RTOpT<double> &primary_op,
  const ArrayView<const Ptr<const MultiVectorBase<double> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<double> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal primary_global_offset
  ) const
{
  MultiVectorAdapterBase<double>::mvMultiReductApplyOpImpl(
    primary_op, multi_vecs, targ_multi_vecs, reduct_objs, primary_global_offset);
}


void EpetraMultiVector::
acquireDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::ConstSubMultiVectorView<double>* sub_mv
  ) const
{
  SpmdMultiVectorDefaultBase<double>::
    acquireDetachedMultiVectorViewImpl(rowRng, colRng, sub_mv);
}


void EpetraMultiVector::
acquireNonconstDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::SubMultiVectorView<double>* sub_mv
  )
{
  SpmdMultiVectorDefaultBase<double>::
    acquireNonconstDetachedMultiVectorViewImpl(rowRng, colRng, sub_mv);
}


void EpetraMultiVector::
commitNonconstDetachedMultiVectorViewImpl(
  RTOpPack::SubMultiVectorView<double>* sub_mv
  )
{
  SpmdMultiVectorDefaultBase<double>::
    commitNonconstDetachedMultiVectorViewImpl(sub_mv);
}


// Overridden protected members from SpmdMultiVectorBase


RCP<const SpmdVectorSpaceBase<double> >
EpetraMultiVector::spmdSpaceImpl() const
{
  return epetraVectorSpace_;
}


void EpetraMultiVector::getNonconstLocalMultiVectorDataImpl(
  const Ptr<ArrayRCP<double> > &localValues, const Ptr<Ordinal> &leadingDim
  )
{
  Epetra_MultiVector& emv = *epetraMultiVector_.getNonconstObj();
  double* data = emv[0];
  *leadingDim = emv.Stride();
  *localValues = Teuchos::arcp(data,0,emv.MyLength()*emv.NumVectors(),false);

  // The arcp is not owning the data. This could potentially lead to dangling pointers,
  // in case the epetraMultiVector is deleted *before* the arcp. To prevent this, we set
  // the epetraMultiVector rcp as extra data for the arcp, so that the Epetra MultiVector will be
  // destroyed *after* the arcp. Notice that there is no cyclic dependency, since
  // the arcp is not stored in the epetra multi vector.
  Teuchos::set_extra_data(epetraMultiVector_.getNonconstObj(),"epetra multi vector",localValues,Teuchos::POST_DESTROY);
}


void EpetraMultiVector::getLocalMultiVectorDataImpl(
  const Ptr<ArrayRCP<const double> > &localValues, const Ptr<Ordinal> &leadingDim
  ) const
{
  const Epetra_MultiVector& emv = *epetraMultiVector_.getConstObj();
  double* data = emv[0];
  *leadingDim = emv.Stride();
  *localValues = Teuchos::arcp(data,0,emv.MyLength()*emv.NumVectors(),false);

  // The arcp is not owning the data. This could potentially lead to dangling pointers,
  // in case the epetraMultiVector is deleted *before* the arcp. To prevent this, we set
  // the epetraMultiVector rcp as extra data for the arcp, so that the Epetra MultiVector will be
  // destroyed *after* the arcp. Notice that there is no cyclic dependency, since
  // the arcp is not stored in the epetra multi vector.
  Teuchos::set_extra_data(epetraMultiVector_.getConstObj(),"epetra multi vector",localValues,Teuchos::POST_DESTROY);
}


void EpetraMultiVector::euclideanApply(
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X,
  const Ptr<MultiVectorBase<double> > &Y,
  const double alpha,
  const double beta
  ) const
{
  // Try to extract Epetra objects from X and Y
  Teuchos::RCP<const Epetra_MultiVector> X_epetra = this->getConstEpetraMultiVector(Teuchos::rcpFromRef(X));
  Teuchos::RCP<Epetra_MultiVector> Y_epetra = this->getEpetraMultiVector(Teuchos::rcpFromPtr(Y));

  // If the cast succeeded, call Epetra directly.
  // Otherwise, fall back to the default implementation.
  if (nonnull(X_epetra) && nonnull(Y_epetra)) {
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

    Y_epetra->Multiply(trans, 'N', alpha, *epetraMultiVector_.getConstObj(), *X_epetra, beta);
  } else {
    SpmdMultiVectorDefaultBase<double>::euclideanApply(M_trans, X, Y, alpha, beta);
  }

}

// private


template<class EpetraMultiVector_t>
void EpetraMultiVector::initializeImpl(
  const RCP<const VectorSpaceBase<double> > &epetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
  const RCP<EpetraMultiVector_t> &epetraMultiVector
  )
{
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(epetraVectorSpace));
  TEUCHOS_ASSERT(nonnull(domainSpace));
  TEUCHOS_ASSERT(nonnull(epetraMultiVector));
#endif
  epetraVectorSpace_ = Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(epetraVectorSpace);
  TEUCHOS_TEST_FOR_EXCEPTION(epetraVectorSpace_.is_null(), std::runtime_error,
                             "Error! Could not cast input vector space to EpetraVectorSpace.\n");
  TEUCHOS_TEST_FOR_EXCEPTION(!epetraVectorSpace_->getEpetraMap()->SameAs(epetraMultiVector->Map()),
                             std::logic_error,
                             "Error! Input multivector has a map that is different "
                             "from that stored in the input EpetraVectorSpace.\n");

  domainSpace_ = domainSpace;
  epetraMultiVector_.initialize(epetraMultiVector);
  this->updateSpmdSpace();
}


RCP<Epetra_MultiVector>
EpetraMultiVector::
getEpetraMultiVector(const RCP<MultiVectorBase<double> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::EpetraMultiVector TEMV;
  typedef Thyra::EpetraVector TEV;

  RCP<TEMV> emv = rcp_dynamic_cast<TEMV>(mv);
  if (nonnull(emv)) {
    return emv->getEpetraMultiVector();
  }

  RCP<TEV> ev = rcp_dynamic_cast<TEV>(mv);
  if (nonnull(ev)) {
    return ev->getEpetraVector();
  }

  return Teuchos::null;
}

RCP<const Epetra_MultiVector>
EpetraMultiVector::
getConstEpetraMultiVector(const RCP<const MultiVectorBase<double> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::EpetraMultiVector TEMV;
  typedef Thyra::EpetraVector TEV;

  RCP<const TEMV> emv = rcp_dynamic_cast<const TEMV>(mv);
  if (nonnull(emv)) {
    return emv->getConstEpetraMultiVector();
  }

  RCP<const TEV> ev = rcp_dynamic_cast<const TEV>(mv);
  if (nonnull(ev)) {
    return ev->getConstEpetraVector();
  }

  return Teuchos::null;
}


} // end namespace Thyra
