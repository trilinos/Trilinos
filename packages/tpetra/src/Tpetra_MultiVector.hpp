// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

// FINISH: some of these arrayview objects should be something else, like Ptr
// FINISH: consider requiring that ArrayView objects are the exact size needed, and no larger.
// FINISH: if both MVs are contiguous, a lot of the level one blas routines below can be turned from multiple calls to one call for best efficiency (eliminate loop over numVecs)

#ifndef TPETRA_MULTIVECTOR_HPP
#define TPETRA_MULTIVECTOR_HPP

#include <Teuchos_TestForException.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include "Tpetra_MultiVectorDecl.hpp"
#include "Tpetra_MultiVectorData.hpp"

namespace Tpetra {

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, Ordinal NumVectors, bool zeroOut) 
    : Teuchos::CompObject(), DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::MultiVector")
  {
    using Teuchos::as;
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const Ordinal myLen = myLength();
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->values_ = Teuchos::arcp<Scalar>(NumVectors*myLen);
    if (zeroOut) {
      std::fill(MVData_->values_.begin(),MVData_->values_.end(),Teuchos::ScalarTraits<Scalar>::zero());
    }
    MVData_->pointers_ = Teuchos::arcp<Teuchos::ArrayRCP<Scalar> >(NumVectors);
    for (Ordinal i = as<Ordinal>(0); i < NumVectors; ++i) {
      MVData_->pointers_[i] = MVData_->values_.persistingView(i*myLen,myLen);
    }
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source) 
    : Teuchos::CompObject(), DistObject<Ordinal,Scalar>(source)
  {
    // copy data from the source MultiVector into this multivector
    // this multivector will be allocated with constant stride, even if the source multivector does not have constant stride
    using Teuchos::as;
    const Ordinal myLen   = myLength();
    const Ordinal numVecs = source.numVectors();
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->values_ = Teuchos::arcp<Scalar>(numVecs*myLen);
    MVData_->pointers_ = Teuchos::arcp<Teuchos::ArrayRCP<Scalar> >(numVecs);
    for (Ordinal i = as<Ordinal>(0); i < numVecs; ++i) {
      MVData_->pointers_[i] = MVData_->values_.persistingView(i*myLen,myLen);
      std::copy( source.MVData_->pointers_[i].begin(), source.MVData_->pointers_[i].end(), MVData_->pointers_[i].begin() );
    }
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &A, Ordinal LDA, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::MultiVector")
  {
    using Teuchos::ArrayView;
    using std::copy;
    using Teuchos::as;
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const Ordinal myLen = this->getMap().getNumMyEntries();
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(LDA < myLen, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): LDA must be large enough to accomodate the local entries.");
    // need LDA*(NumVectors-1)+myLen elements in A
    TEST_FOR_EXCEPTION(A.size() < LDA*(NumVectors-1)+myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(): A,LDA must be large enough to accomodate the local entries.");
#endif
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->values_ = Teuchos::arcp<Scalar>(NumVectors*myLen);
    MVData_->pointers_ = Teuchos::arcp<Teuchos::ArrayRCP<Scalar> >(NumVectors);
    for (Ordinal i = as<Ordinal>(0); i < NumVectors; ++i) {
      MVData_->pointers_[i] = MVData_->values_.persistingView(i*myLen,myLen);
      // copy data from A to my internal data structure
      ArrayView<const Scalar> Aptr = A(i*LDA,myLen);
      copy(Aptr.begin(),Aptr.end(),MVData_->pointers_[i].begin());
    }
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &arrayOfArrays, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(map, map.getPlatform()->createComm(), "Tpetra::MultiVector")
  {
    TEST_FOR_EXCEPT(true);
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source, const Teuchos::ArrayView<const Ordinal> &indices)
    : DistObject<Ordinal,Scalar>(source)
  {
    TEST_FOR_EXCEPT(true);
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source, Ordinal startIndex, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(source)
  {
    TEST_FOR_EXCEPT(true);
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::~MultiVector()
  {
  }


  template <typename Ordinal, typename Scalar> 
  bool MultiVector<Ordinal,Scalar>::constantStride() const
  {
    return MVData_->constantStride_;
  }


  template <typename Ordinal, typename Scalar> 
  Ordinal MultiVector<Ordinal,Scalar>::myLength() const
  {
    return this->getMap().getNumMyEntries();
  }


  template <typename Ordinal, typename Scalar> 
  Ordinal MultiVector<Ordinal,Scalar>::globalLength() const
  {
    return this->getMap().getNumGlobalEntries();
  }


  template <typename Ordinal, typename Scalar> 
  Ordinal MultiVector<Ordinal,Scalar>::stride() const
  {
    return MVData_->stride_;
  }


  template <typename Ordinal, typename Scalar> 
  void MultiVector<Ordinal,Scalar>::print(std::ostream &os) const
  {
    using std::endl;
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm = this->getMap().getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    for(int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        if (myImageID == 0) {
          os << "Number of vectors: " << numVectors() << endl;
          os << "Global length: " << globalLength() << endl;
        }
        os << "Local length: " << myLength() << endl;
        os << "Local stride: " << stride() << endl;
        os << "Constant stride: " << (constantStride() ? "true" : "false") << endl;
      }
      // Do a few global ops to give I/O a chance to complete
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }
  }


  template <typename Ordinal, typename Scalar> 
  void MultiVector<Ordinal,Scalar>::printValues(std::ostream &os) const
  {
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  bool MultiVector<Ordinal,Scalar>::checkSizes(const DistObject<Ordinal,Scalar> &/*sourceObj*/) 
  {
    TEST_FOR_EXCEPT(true);
    return true;
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::copyAndPermute(
      const DistObject<Ordinal,Scalar> &/*sourceObj*/,
      Ordinal /*numSameIDs*/, Ordinal /*numPermuteIDs*/,
      const std::vector<Ordinal> &/*permuteToLIDs*/, const std::vector<Ordinal> &/*permuteFromLIDs*/)
  {
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::packAndPrepare(
      const DistObject<Ordinal,Scalar> &/*sourceObj*/,
      Ordinal /*numExportIDs*/,
      const std::vector<Ordinal> &/*exportLIDs*/, std::vector<Scalar> &/*exports*/,
      Ordinal &/*packetSize*/,
      Distributor<Ordinal> &/*distor*/) 
  {
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::unpackAndCombine(
      Ordinal /*numImportIDs*/,
      const std::vector<Ordinal> &/*importLIDs*/,
      const std::vector<Scalar> &/*imports*/,
      Distributor<Ordinal> &/*distor*/,
      CombineMode /*CM*/) 
  {
    TEST_FOR_EXCEPT(true);
  }


  template<typename Ordinal, typename Scalar>
  Ordinal MultiVector<Ordinal,Scalar>::numVectors() const 
  {
    return Teuchos::as<Ordinal>(MVData_->pointers_.size());
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::dot(
      const MultiVector<Ordinal,Scalar> &A, 
      const Teuchos::ArrayView<Scalar> &dots) const 
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    // compute local dot products of *this and A
    // sum these across all nodes
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors must have the same number of vectors.");
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(dots.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::dots(A,dots): dots.size() must be as large as the number of vectors in *this and A.");
#endif
    Teuchos::Array<Scalar> ldots(numVecs);
    for (Ordinal i=ZERO; i<numVecs; ++i) {
      ldots[i] = blas.DOT(MVData_->pointers_[i].size(),MVData_->pointers_[i].getRawPtr(),ONE,A[i].getRawPtr(),ONE);
    }
    Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,ldots.getRawPtr(),dots.getRawPtr());
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::norm1(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    // compute local components of the norms
    // sum these across all nodes
    const Ordinal numVecs = this->numVectors();
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm1(norms): norms.size() must be as large as the number of vectors in *this.");
#endif
    Teuchos::Array<Mag> lnorms(numVecs);
    for (Ordinal i=ZERO; i<numVecs; ++i) {
      lnorms[i] = blas.ASUM(MVData_->pointers_[i].size(),MVData_->pointers_[i].getRawPtr(),ONE);
    }
    Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::norm2(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    // compute local components of the norms
    // sum these across all nodes
    const Ordinal numVecs = this->numVectors();
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm2(norms): norms.size() must be as large as the number of vectors in *this.");
#endif
    Teuchos::Array<Mag> lnorms(numVecs,ScalarTraits<Mag>::zero());
    for (Ordinal j=ZERO; j<numVecs; ++j) {
      Teuchos::ArrayRCP<const Scalar> cpos = MVData_->pointers_[j].getConst();
      for (; cpos != cpos.end(); ++cpos) {
        lnorms[j] += ScalarTraits<Scalar>::magnitude( 
                       (*cpos) * ScalarTraits<Scalar>::conjugate(*cpos)
                     );
      }
    }
    Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n);
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::normWeighted(
      const MultiVector<Ordinal,Scalar> &weights,
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::as;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal numImages = this->getMap().getComm()->getSize();
    bool OneW = false;
    const Ordinal numVecs = this->numVectors();
    if (weights.numVectors() == ONE) {
      OneW = true;
    }
    else {
      TEST_FOR_EXCEPTION(weights.numVectors() != numVecs, std::runtime_error,
          "Tpetra::MultiVector::normWeighted(): MultiVector of weights must contain either one vector or the same number of vectors as this.");
    }
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
        "Tpetra::MultiVector::normWeighted(): MultiVectors must have compatible Maps.");
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::normWeighted(): norms.size() must be as large as the number of vectors in *this.");
#endif
    // compute local components of the norms
    // sum these across all nodes
    Array<Mag> lnorms(numVecs,ScalarTraits<Mag>::zero());
    for (Ordinal j=ZERO; j<numVecs; ++j) {
      ArrayRCP<const Scalar> wpos = (OneW ? weights[0] : weights[j]);
      ArrayRCP<const Scalar> cpos = MVData_->pointers_[j].getConst();
      for (; cpos != cpos.end(); ++cpos, ++wpos) {
        Scalar tmp = *cpos / *wpos;
        lnorms[j] += SCT::magnitude( tmp * SCT::conjugate(tmp) );
      }
    }
    Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n/as<Mag>(numImages));
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::normInf(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    // compute local components of the norms
    // sum these across all nodes
    const Ordinal numVecs = this->numVectors();
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::normInf(norms): norms.size() must be as large as the number of vectors in *this.");
#endif
    Teuchos::Array<Mag> lnorms(numVecs);
    for (Ordinal i=ZERO; i<numVecs; ++i) {
      // careful!!! IAMAX returns FORTRAN-style (i.e., one-based) index. subtract ind by one
      Ordinal ind = blas.IAMAX(MVData_->pointers_[i].size(),MVData_->pointers_[i].getRawPtr(),ONE) - ONE;
      lnorms[i] = Teuchos::ScalarTraits<Scalar>::magnitude( MVData_->pointers_[i][ind] );
    }
    Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::random() 
  {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal myLen   = this->myLength();
    const Ordinal numVecs = this->numVectors();
    for (Ordinal j=ZERO; j<numVecs; ++j) {
      Teuchos::ArrayRCP<Scalar> cpos = MVData_->pointers_[j];
      for (Ordinal i=ZERO; i<myLen; ++i) {
        *cpos = Teuchos::ScalarTraits<Scalar>::random();
        ++cpos;
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::putScalar(const Scalar &alpha) 
  {
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    const Ordinal numVecs = this->numVectors();
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      std::fill(curpos.begin(),curpos.end(),alpha);
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::scale(const Scalar &alpha) 
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    using Teuchos::OrdinalTraits;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    using Teuchos::ArrayRCP;
    const Ordinal numVecs = this->numVectors();
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      blas.SCAL(curpos.size(),alpha,curpos.getRawPtr(),ONE);
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::scale(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A) 
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have the same number of vectors.");
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      ArrayRCP<Scalar> &Apos = A.MVData_->pointers_[i];
      // copy A to *this
      blas.COPY(curpos.size(),Apos.getRawPtr(),ONE,curpos.getRawPtr(),ONE);
      // then scale *this in-situ
      blas.SCAL(curpos.size(),alpha,curpos.getRawPtr(),ONE);
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::reciprocal(const MultiVector<Ordinal,Scalar> &A) 
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayRCP;
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::reciprocal(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::reciprocal(): MultiVectors must have the same number of vectors.");
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      ArrayRCP<Scalar> &Apos = A.MVData_->pointers_[i];
      for (; curpos != curpos.end(); ++curpos, ++Apos) {
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPTION( ScalarTraits<Scalar>::magnitude(*Apos) <= ScalarTraits<Scalar>::sfmin() ||
            *Apos == ScalarTraits<Scalar>::sfmin(), std::runtime_error,
            "Tpetra::MultiVector::reciprocal(): element of A was zero or too small to invert: " << *Apos );
#endif
        *curpos = ScalarTraits<Scalar>::one()/(*Apos);
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::abs(const MultiVector<Ordinal,Scalar> &A) 
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::abs(): MultiVectors must have the same number of vectors.");
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::abs(): MultiVectors must have compatible Maps.");
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      ArrayRCP<Scalar> &Apos = A.MVData_->pointers_[i];
      for (; curpos != curpos.end(); ++curpos, ++Apos) {
        *curpos = Teuchos::ScalarTraits<Scalar>::magnitude(*Apos);
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta) 
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      ArrayRCP<Scalar> &Apos = A.MVData_->pointers_[i];
      blas.SCAL(curpos.size(),beta,curpos.getRawPtr(),ONE);
      blas.AXPY(curpos.size(),alpha,Apos.getRawPtr(),ONE,curpos.getRawPtr(),ONE);
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta, const MultiVector<Ordinal,Scalar> &B, const Scalar &gamma)
  {
    Teuchos::BLAS<Ordinal,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()) || !this->getMap().isCompatible(B.getMap()),
        std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs || B.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayRCP<Scalar> &curpos = MVData_->pointers_[i];
      ArrayRCP<Scalar> &Apos = A.MVData_->pointers_[i];
      ArrayRCP<Scalar> &Bpos = B.MVData_->pointers_[i];
      blas.SCAL(curpos.size(),gamma,curpos.getRawPtr(),ONE);
      blas.AXPY(curpos.size(),alpha,Apos.getRawPtr(),ONE,curpos.getRawPtr(),ONE);
      blas.AXPY(curpos.size(),beta, Bpos.getRawPtr(),ONE,curpos.getRawPtr(),ONE);
    }
  }


  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayRCP<const Scalar> MultiVector<Ordinal,Scalar>::operator[](Ordinal i) const
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return MVData_->pointers_[i].getConst();
  }


  /*

      // Basic constructor
      MultiVector(const VectorSpace<Ordinal, Scalar>& vectorSpace, const Ordinal NumVectors) :
        vectorSpace_(vectorSpace),
        NumVectors_(NumVectors)
      {
        Init(); 

        array_.resize(getNumVectors());

        for (Ordinal i = OrdinalZero_ ; i < NumVectors ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<Ordinal, Scalar> (vectorSpace));
        }
      }

      // Creates a deep copy of the input set of vectors.
      MultiVector(const VectorSpace<Ordinal, Scalar>& vectorSpace, 
                  std::vector<Tpetra::Vector<Ordinal, Scalar> const *> list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init(); 

        array_.resize(getNumVectors());

        for (Ordinal i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // deep copy of each of the vectors
          array_[i] = Teuchos::rcp(new Vector<Ordinal, Scalar> (*(list[i])));
        }
      }

      // Creates a shallow copy of the input set of vectors.
      MultiVector(const VectorSpace<Ordinal, Scalar>& vectorSpace, 
                  std::vector<Teuchos::RCP<Tpetra::Vector<Ordinal, Scalar> > > list) :
        vectorSpace_(vectorSpace),
        NumVectors_(list.size())
      {
        Init();

        array_.resize(NumVectors_);

        for (Ordinal i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          // copy RCP's from the list to this array.
          array_[i] = list[i];
        }
      }

      // Copy constructor.
      MultiVector(const MultiVector& rhs) :
        vectorSpace_(rhs.vectorSpace()),
        NumVectors_(rhs.getNumVectors())
      {
        Init();

        array_.resize(NumVectors_);

        for (Ordinal i = OrdinalZero_ ; i < NumVectors_ ; ++i)
        {
          array_[i] = Teuchos::rcp(new Vector<Ordinal, Scalar> (*(rhs.GetVector(i))));
        }
      }

      //! Returns the global number of entries.
      Ordinal getNumGlobalEntries() const
      {
        return(vectorSpace_.getNumGlobalEntries());
      }

      //! Returns the number of entries on the calling image.
      Ordinal getNumMyEntries() const
      {
        return(vectorSpace_.getNumMyEntries());
      }

      //! Returns the number of vectors in this multivector.
      Ordinal getNumVectors() const
      {
        return(NumVectors_);
      }

      //! Returns a reference to the vector space of this multivector.
      VectorSpace<Ordinal, Scalar> const& vectorSpace () const
      {
        return(vectorSpace_);
      }

      //! Returns a reference to the i-th element of the j-th vector.
      Scalar& operator() (const Ordinal i, const Ordinal j)
      {
        return((*array_[j])[i]);
      }

      //! Returns a reference to the i-th element of the j-th vector (const version)
      Scalar const& operator() (const Ordinal i, const Ordinal j) const
      {
        return((*array_[j])[i]);
      }

      //! Sets all elements of all vectors to the given value.
      void setAllToScalar(Scalar const value)
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          array_[i]->setAllToScalar(value);
      }

      //! Sets all elements of all vectors to random value.
      void setAllToRandom()
      {
        // FIXME: sets only the real part to random
        for (int i = 0 ; i < NumVectors_ ; ++i)
        {
          //array_[i]->setAllToRandom();
          for (int j = 0 ; j < array_[0]->getNumMyEntries() ; ++j)
          {
            // FIXME (*array_[i])[j] = complex<double>(Teuchos::ScalarTraits<double>::random(), 0.0);
            (*array_[i])[j] = Teuchos::ScalarTraits<Scalar>::random();
          }
        }
      }

      //! Prints the vector to cout. FIXME
      void Print() const
      {
        for (int i = 0 ; i < NumVectors_ ; ++i)
          cout << (*array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector.
      Teuchos::RCP<Tpetra::Vector<Ordinal, Scalar> > GetRCP(const int i)
      {
        return(array_[i]);
      }

      //! Returns a RCP pointer to the i-th vector (const version).
      Teuchos::RCP<Tpetra::Vector<Ordinal, Scalar> const > GetRCP(const int i) const
      {
        return(array_[i]);
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector.
      Tpetra::Vector<Ordinal, Scalar>* GetVector(const int i)
      {
        return(array_[i].get());
      }

      //! Returns a Tpetra::Vector pointer to the i-th vector (const version).
      const Tpetra::Vector<Ordinal, Scalar>* GetVector(const int i) const
      {
        return(array_[i].get());
      }

      void norm1(Scalar* Values) const
      {
        for (Ordinal i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<Scalar>::magnitude(array_[i]->norm1());
        }
      }

      void norm2(typename Teuchos::ScalarTraits<Scalar>::magnitudeType* Values) const
      {
        for (Ordinal i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<Scalar>::magnitude(array_[i]->norm2());
        }
      }

      void normInf(Scalar* Values) const
      {
        for (Ordinal i = OrdinalZero_ ; i < getNumVectors() ; ++i)
        {
          Values[i] = Teuchos::ScalarTraits<Scalar>::magnitude(array_[i]->normInf());
        }
      }

      void dotProduct(const MultiVector<Ordinal, Scalar>& A, Scalar* Values) const
      {
        for (int v = 0 ; v < getNumVectors() ; ++v)
        {
          Values[v] = GetVector(v)->dotProduct(*(A.GetVector(v)));
        }
      }

    private:

      void Init()
      {
        OrdinalZero_ = Teuchos::ScalarTraits<Ordinal>::zero();
        OrdinalOne_  = Teuchos::ScalarTraits<Ordinal>::one();

        ScalarZero_ = Teuchos::ScalarTraits<Scalar>::zero();
        ScalarOne_  = Teuchos::ScalarTraits<Scalar>::one();
      }

      std::vector<Teuchos::RCP<Vector<Ordinal, Scalar> > > array_;
      const VectorSpace<Ordinal, Scalar>& vectorSpace_;
      Ordinal NumVectors_;

      Ordinal OrdinalZero_;
      Ordinal OrdinalOne_;

      Scalar ScalarZero_;
      Scalar ScalarOne_;
*/

} // namespace Tpetra

#endif // TPETRA_MULTIVECTOR_HPP
