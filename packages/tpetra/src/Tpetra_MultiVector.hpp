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

// TODO: with contiguous MVs, some of the level one blas routines below can be turned from multiple calls to one call for best efficiency (eliminate loop over numVecs)

#ifndef TPETRA_MULTIVECTOR_HPP
#define TPETRA_MULTIVECTOR_HPP

#include <Teuchos_TestForException.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_BLAS.hpp>

#include "Tpetra_MultiVectorDecl.hpp"
#include "Tpetra_MultiVectorData.hpp"

namespace Tpetra {

  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, Ordinal NumVectors, bool zeroOut) 
    : DistObject<Ordinal,Scalar>(map, map.getComm(), "Tpetra::MultiVector")
  {
    using Teuchos::as;
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const Ordinal myLen = myLength();
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->ptrs_.resize(NumVectors,Teuchos::null);
    if (myLen > as<Ordinal>(0)) {
      MVData_->values_ = Teuchos::arcp<Scalar>(NumVectors*myLen);
      if (zeroOut) {
        std::fill(MVData_->values_.begin(),MVData_->values_.end(),Teuchos::ScalarTraits<Scalar>::zero());
      }
      for (Ordinal i = as<Ordinal>(0); i < NumVectors; ++i) {
        MVData_->ptrs_[i] = MVData_->values_(i*myLen,myLen);
      }
    }
    MVData_->updateConstPointers();
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const MultiVector<Ordinal,Scalar> &source) 
    : DistObject<Ordinal,Scalar>(source)
  {
    // copy data from the source MultiVector into this multivector
    // this multivector will be allocated with constant stride, even if the source multivector does not have constant stride
    using Teuchos::as;
    const Ordinal myLen   = myLength();
    const Ordinal numVecs = source.numVectors();
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->ptrs_.resize(numVecs,Teuchos::null);
    if (myLen > as<Ordinal>(0)) {
      MVData_->values_ = Teuchos::arcp<Scalar>(numVecs*myLen);
      for (Ordinal i = as<Ordinal>(0); i < numVecs; ++i) {
        MVData_->ptrs_[i] = MVData_->values_(i*myLen,myLen);
        std::copy( source.MVData_->ptrs_[i].begin(), source.MVData_->ptrs_[i].end(), MVData_->ptrs_[i].begin() );
      }
    }
    MVData_->updateConstPointers();
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &A, Ordinal LDA, Ordinal NumVectors)
    : DistObject<Ordinal,Scalar>(map, map.getComm(), "Tpetra::MultiVector")
  {
    using Teuchos::ArrayView;
    using Teuchos::as;
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const Ordinal myLen = myLength();
#ifdef TEUCHOS_DEBUG
    // need LDA*(NumVectors-1)+myLen elements in A
    TEST_FOR_EXCEPTION(A.size() < LDA*(NumVectors-1)+myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(): A,LDA must be large enough to accomodate the local entries.");
#endif
    TEST_FOR_EXCEPTION(LDA < myLen, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): LDA must be large enough to accomodate the local entries.");
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->ptrs_.resize(NumVectors,Teuchos::null);
    if (myLen > as<Ordinal>(0)) {
      MVData_->values_ = Teuchos::arcp<Scalar>(NumVectors*myLen);
      for (Ordinal i = as<Ordinal>(0); i < NumVectors; ++i) {
        MVData_->ptrs_[i] = MVData_->values_(i*myLen,myLen);
        // copy data from A to my internal data structure
        ArrayView<const Scalar> Aptr = A(i*LDA,myLen);
        std::copy(Aptr.begin(),Aptr.end(),MVData_->ptrs_[i].begin());
      }
    }
    MVData_->updateConstPointers();
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > &mvdata) 
    : DistObject<Ordinal,Scalar>(map, map.getComm(), "Tpetra::MultiVector"), MVData_(mvdata)
  {}


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &arrayOfArrays)
    : DistObject<Ordinal,Scalar>(map, map.getComm(), "Tpetra::MultiVector")
  {
    using Teuchos::as;
    const Ordinal myLen = myLength();
    Ordinal NumVectors = arrayOfArrays.size();
    TEST_FOR_EXCEPTION(NumVectors < 1, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(map,arrayOfArrays): arrayOfArrays.size() must be strictly positive.");
    MVData_ = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    MVData_->constantStride_ = true;
    MVData_->stride_ = myLen;
    MVData_->ptrs_.resize(NumVectors,Teuchos::null);
    if (myLen > as<Ordinal>(0)) {
      MVData_->values_ = Teuchos::arcp<Scalar>(NumVectors*myLen);
      for (Ordinal i = as<Ordinal>(0); i < NumVectors; ++i) {
        MVData_->ptrs_[i] = MVData_->values_(i*myLen,myLen);
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPTION(arrayOfArrays[i].size() != myLength(), std::runtime_error,
            "Tpetra::MultiVector::MultiVector(map,arrayOfArrays): arrayOfArrays[" << i << "].size() (==" << arrayOfArrays[i].size() 
            << ") != myLength() (==" << myLength() << ")");
#endif
        std::copy(arrayOfArrays[i].begin(),arrayOfArrays[i].end(),MVData_->ptrs_[i].begin());
      }
    }
    MVData_->updateConstPointers();
  }


  template <typename Ordinal, typename Scalar> 
  MultiVector<Ordinal,Scalar>::~MultiVector()
  {}


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
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
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
    using std::endl;
    using std::setw;
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm = this->getMap().getComm();
    const Map<Ordinal> &map = this->getMap();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        for (int i=0; i<myLength(); ++i) {
          os << setw(4) << map.getGlobalIndex(i) << "\t";
          for (int j=0; j<numVectors(); ++j) {
            os << setw(20) << MVData_->ptrs_[j][i] << " ";
          }
          os << endl;
        }
      }
      // Do a few global ops to give I/O a chance to complete
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }
  }


  template<typename Ordinal, typename Scalar>
  bool MultiVector<Ordinal,Scalar>::checkSizes(const DistObject<Ordinal,Scalar> &sourceObj, Ordinal &packetSize) 
  {
    const MultiVector<Ordinal,Scalar> &A = dynamic_cast<const MultiVector<Ordinal,Scalar>&>(sourceObj);
    // objects maps have already been checked. simply check the number of vectors.
    packetSize = this->numVectors();
    return (A.numVectors() == packetSize);
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::copyAndPermute(
      const DistObject<Ordinal,Scalar> & sourceObj,
      Ordinal numSameIDs,
      const Teuchos::ArrayView<const Ordinal> &permuteToLIDs,
      const Teuchos::ArrayView<const Ordinal> &permuteFromLIDs)
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    using Teuchos::ArrayView;
    const MultiVector<Ordinal,Scalar> &sourceMV = dynamic_cast<const MultiVector<Ordinal,Scalar> &>(sourceObj);
    ArrayView<Scalar>          dstView(Teuchos::null);
    ArrayView<const Scalar>    srcView(Teuchos::null);
    typename ArrayView<const Ordinal>::iterator pTo, pFrom;
#ifdef TEUCHOS_DEBUG
    // any other error will be caught by Teuchos
    TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::copyAndPermute(): permuteToLIDs and permuteFromLIDs must have the same size.");
#endif
    // one vector at a time
    for (Ordinal j = OT::zero(); j < numVectors(); ++j) {
      // the first numImportIDs GIDs are the same between source and target,
      // we can just copy them
      dstView = (*this)[j];
      srcView = sourceMV[j];
      std::copy(srcView.begin(),srcView.begin()+numSameIDs,dstView.begin());
      // next, do permutations
      for (pTo = permuteToLIDs.begin(), pFrom = permuteFromLIDs.begin();
           pTo != permuteToLIDs.end(); ++pTo, ++pFrom)
      {
        dstView[*pTo] = srcView[*pFrom];
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::packAndPrepare(
      const DistObject<Ordinal,Scalar> & sourceObj,
      const Teuchos::ArrayView<const Ordinal> &exportLIDs,
      const Teuchos::ArrayView<Scalar> &exports,
      Distributor<Ordinal> &distor)
  {
    const MultiVector<Ordinal,Scalar> &sourceMV = dynamic_cast<const MultiVector<Ordinal,Scalar> &>(sourceObj);
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    using Teuchos::ArrayView;
    (void)distor;    // we don't use these, but we don't want unused parameter warnings
    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ; 
                   ....
                   all of the data from row exportLIDs.back() }
      this doesn't have the best locality, but is necessary because the data for a Packet
      (all data associated with an LID) is required to be contiguous */
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(exports.size() != numVectors()*exportLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::packAndPrepare(): sizing of exports buffer should be appropriate for the amount of data to be exported.");
#endif
    typename ArrayView<       Scalar>::iterator  expptr;
    typename ArrayView<const Ordinal>::iterator  idptr;
    expptr = exports.begin();
    for (idptr = exportLIDs.begin(); idptr != exportLIDs.end(); ++idptr) {
      for (Ordinal j = OT::zero(); j < numVectors(); ++j) {
        *expptr++ = sourceMV[j][*idptr];
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::unpackAndCombine(
      const Teuchos::ArrayView<const Ordinal> &importLIDs,
      const Teuchos::ArrayView<const Scalar> &imports,
      Distributor<Ordinal> &distor,
      CombineMode CM)
  {
    (void)distor; // we don't use this, but we don't want unused parameter warnings
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    using Teuchos::ArrayView;
    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ; 
                   ....
                   all of the data from row exportLIDs.back() }
      this doesn't have the best locality, but is necessary because the data for a Packet
      (all data associated with an LID) is required to be contiguous */
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(imports.size() != numVectors()*importLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::unpackAndCombine(): sizing of imports buffer should be appropriate for the amount of data to be exported.");
#endif
    typename ArrayView<const  Scalar>::iterator  impptr;
    typename ArrayView<const Ordinal>::iterator  idptr;
    impptr = imports.begin();
    if (CM == INSERT || CM == REPLACE) 
    {
      for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
        for (Ordinal j = OT::zero(); j < numVectors(); ++j) {
          (*this)[j][*idptr] = *impptr++;
        }
      }
    }
    else if (CM == ADD) {
      for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
        for (Ordinal j = OT::zero(); j < numVectors(); ++j) {
          (*this)[j][*idptr] += *impptr++;
        }
      }
    }
    else {
      TEST_FOR_EXCEPTION(CM != ADD && CM != REPLACE && CM != INSERT, std::invalid_argument,
          "Tpetra::MultiVector::unpackAndCombine(): Invalid CombineMode: " << CM);
    }
  }


  template<typename Ordinal, typename Scalar>
  Ordinal MultiVector<Ordinal,Scalar>::numVectors() const 
  {
    return Teuchos::as<Ordinal>(MVData_->ptrs_.size());
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::dot(
      const MultiVector<Ordinal,Scalar> &A, 
      const Teuchos::ArrayView<Scalar> &dots) const 
  {
    Teuchos::BLAS<int,Scalar> blas;
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
      ldots[i] = blas.DOT(MVData_->cPtrs_[i].size(),MVData_->cPtrs_[i].getRawPtr(),ONE,A[i].getRawPtr(),ONE);
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,ldots.getRawPtr(),dots.getRawPtr());
    }
    else {
      std::copy(ldots.begin(),ldots.end(),dots.begin());
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::norm1(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    Teuchos::BLAS<int,Scalar> blas;
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
      lnorms[i] = blas.ASUM(MVData_->cPtrs_[i].size(),MVData_->cPtrs_[i].getRawPtr(),ONE);
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    else {
      std::copy(lnorms.begin(),lnorms.end(),norms.begin());
    }
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
      typename Teuchos::ArrayView<const Scalar>::iterator cpos = MVData_->cPtrs_[j].begin(),
                                                          cend = MVData_->cPtrs_[j].end();
      for (; cpos != cend; ++cpos) {
        lnorms[j] += ScalarTraits<Scalar>::magnitude( 
                       (*cpos) * ScalarTraits<Scalar>::conjugate(*cpos)
                     );
      }
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    else {
      std::copy(lnorms.begin(),lnorms.end(),norms.begin());
    }
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
      typename ArrayView<const Scalar>::iterator wpos = (OneW ? weights[0] : weights[j]).begin(),
                                                 cpos = MVData_->cPtrs_[j].begin(),
                                                 cend = MVData_->cPtrs_[j].end();
      for (; cpos != cend; ++cpos, ++wpos) {
        Scalar tmp = *cpos / *wpos;
        lnorms[j] += SCT::magnitude( tmp * SCT::conjugate(tmp) );
      }
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    else {
      std::copy(lnorms.begin(),lnorms.end(),norms.begin());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n/as<Mag>(numImages));
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::normInf(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    Teuchos::BLAS<int,Scalar> blas;
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
      Ordinal ind = blas.IAMAX(MVData_->cPtrs_[i].size(),MVData_->cPtrs_[i].getRawPtr(),ONE) - ONE;
      lnorms[i] = Teuchos::ScalarTraits<Scalar>::magnitude( MVData_->cPtrs_[i][ind] );
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    else {
      std::copy(lnorms.begin(),lnorms.end(),norms.begin());
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::random() 
  {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal numVecs = this->numVectors();
    for (Ordinal j=ZERO; j<numVecs; ++j) {
      typename Teuchos::ArrayView<Scalar>::iterator cpos = MVData_->ptrs_[j].begin(),
                                                    cend = MVData_->ptrs_[j].end();
      for (; cpos != cend; ++cpos) {
        *cpos = Teuchos::ScalarTraits<Scalar>::random();
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::putScalar(const Scalar &alpha) 
  {
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    const Ordinal numVecs = this->numVectors();
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      ArrayView<Scalar> &curpos = MVData_->ptrs_[i];
      std::fill(curpos.begin(),curpos.end(),alpha);
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::scale(const Scalar &alpha) 
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    using Teuchos::ArrayView;
    const Ordinal numVecs = this->numVectors();
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // do nothing
    }
    else if (alpha == Teuchos::ScalarTraits<Scalar>::zero()) {
      putScalar(alpha);
    }
    else {
      for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
        ArrayView<Scalar> &curpos = MVData_->ptrs_[i];
        blas.SCAL(curpos.size(),alpha,curpos.getRawPtr(),ONE);
      }
    }
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::scale(Teuchos::ArrayView<const Scalar> alphas)
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    using Teuchos::ArrayView;
    const Ordinal numVecs = this->numVectors();
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(alphas.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as the number of vectors in *this.");
#endif
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      if (alphas[i] == Teuchos::ScalarTraits<Scalar>::one()) {
        // do nothing
      }
      else if (alphas[i] == Teuchos::ScalarTraits<Scalar>::zero()) {
        ArrayView<Scalar> &curpos = MVData_->ptrs_[i];
        std::fill(curpos.begin(),curpos.end(),Teuchos::ScalarTraits<Scalar>::zero());
      }
      else {
        ArrayView<Scalar> &curpos = MVData_->ptrs_[i];
        blas.SCAL(curpos.size(),alphas[i],curpos.getRawPtr(),ONE);
      }
    }
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::scale(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A) 
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have the same number of vectors.");
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    if (alpha == Teuchos::ScalarTraits<Scalar>::zero()) {
      putScalar(alpha); // set me = 0.0
    }
    else if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      *this = A;        // set me = A
    }
    else {
      // set me == alpha*A
      for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
        ArrayView<Scalar> &curpos =   MVData_->ptrs_[i],
                          &Apos   = A.MVData_->ptrs_[i];
        // copy A to *this
        blas.COPY(curpos.size(),Apos.getRawPtr(),ONE,curpos.getRawPtr(),ONE);
        // then scale *this in-situ
        blas.SCAL(curpos.size(),alpha,curpos.getRawPtr(),ONE);
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::reciprocal(const MultiVector<Ordinal,Scalar> &A) 
  {
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::reciprocal(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::reciprocal(): MultiVectors must have the same number of vectors.");
    for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
      typename Teuchos::ArrayView<Scalar>::iterator curpos =   MVData_->ptrs_[i].begin(),
                                                      Apos = A.MVData_->ptrs_[i].begin();
      for (; curpos != MVData_->ptrs_[i].end(); ++curpos, ++Apos) {
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
    Teuchos::BLAS<int,Scalar> blas;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::abs(): MultiVectors must have the same number of vectors.");
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::abs(): MultiVectors must have compatible Maps.");
    for (Ordinal j = OrdinalTraits<Ordinal>::zero(); j < numVecs; ++j) {
      typename ArrayView<Scalar>::iterator curpos =   MVData_->ptrs_[j].begin(),
                                             Apos = A.MVData_->ptrs_[j].begin();
      for (; curpos != MVData_->ptrs_[j].end(); ++curpos, ++Apos) {
        *curpos = Teuchos::ScalarTraits<Scalar>::magnitude(*Apos);
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta) 
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename Teuchos::ArrayView<Scalar>::iterator avi;
    typedef typename Teuchos::ArrayView<const Scalar>::iterator avci;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayView;
    if (alpha == ST::zero()) {
      scale(beta);
      return;
    }
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");

    if (beta == ST::zero()) { // this = alpha*A
      scale(alpha,A);
      return;
    }
    else if (beta == ST::one()) { // this = this + alpha*A
      if (alpha == ST::one()) { // this = this + A
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = A.MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = (*curpos) + (*Apos); }
        }
      }
      else { // this = this + alpha*A
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = A.MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = (*curpos) + alpha*(*Apos); }
        }
      }
    }
    else { // this = beta*this + alpha*A
      if (alpha == ST::one()) { // this = beta*this + A
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = A.MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = beta*(*curpos) + (*Apos); }
        }
      }
      else { // this = beta*this + alpha*A
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = A.MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = beta*(*curpos) + alpha*(*Apos); }
        }
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta, const MultiVector<Ordinal,Scalar> &B, const Scalar &gamma)
  {
    typedef typename Teuchos::ArrayView<Scalar>::iterator avi;
    typedef typename Teuchos::ArrayView<const Scalar>::iterator avci;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    if (alpha == ST::zero()) {
      update(beta,B,gamma);
      return;
    }
    else if (beta == ST::zero()) {
      update(alpha,A,gamma);
      return;
    }
    const Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()) || !this->getMap().isCompatible(B.getMap()),
        std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have compatible Maps.");
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs || B.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    // determine if alpha==1 xor beta==1
    // if only one of these is 1.0, make it alpha
    Teuchos::Ptr<const MultiVector<Ordinal,Scalar> > Aptr = Teuchos::ptr(&A), Bptr = Teuchos::ptr(&B);
    Teuchos::Ptr<const Scalar> lalpha = Teuchos::ptr(&alpha),
                               lbeta  = Teuchos::ptr(&beta);
    if (alpha!=ST::one() && beta==ST::one()) {
      // switch them
      Aptr = Teuchos::ptr(&B);
      Bptr = Teuchos::ptr(&A);
      lalpha = Teuchos::ptr(&beta);
      lbeta  = Teuchos::ptr(&alpha);
    }

    if (gamma == ST::zero()) { // this = lalpha*A + lbeta*B
      if (*lalpha == ST::one()) {
        if (*lbeta == ST::one()) { // this = gamma*this + A + B
          for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
            avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*Apos) + (*Bpos); }
          }
        }
        else { // this = A + lbeta*B
          for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
            avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*Apos) + (*lbeta)*(*Bpos); }
          }
        }
      }
      else { // this = lalpha*A + lbeta*B
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
        }
      }
    }
    else if (gamma == ST::one()) { // this = this + lalpha*A + lbeta*B
      if ((*lalpha) == ST::one()) {
        if ((*lbeta) == ST::one()) { // this = this + A + B
          for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
            avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*Apos) + (*Bpos); }
          }
        }
        else { // this = this + A + lbeta*B
          for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
            avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*Apos) + (*lbeta)*(*Bpos); }
          }
        }
      }
      else { // this = this + lalpha*A + lbeta*B
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
        }
      }
    }
    else { // this = gamma*this + lalpha*A + lbeta*B
      if ((*lalpha) == ST::one()) {
        if ((*lbeta) == ST::one()) { // this = gamma*this + A + B
          for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
            avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*Apos) + (*Bpos); }
          }
        }
        else { // this = gamma*this + A + lbeta*B
          for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
            avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*Apos) + (*lbeta)*(*Bpos); }
          }
        }
      }
      else { // this = gamma*this + lalpha*A + lbeta*B
        for (Ordinal i = OrdinalTraits<Ordinal>::zero(); i < numVecs; ++i) {
          avi curpos = MVData_->ptrs_[i].begin(), cend = MVData_->ptrs_[i].end(); avci Apos = Aptr->MVData_->cPtrs_[i].begin(), Bpos = Bptr->MVData_->cPtrs_[i].begin();
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
        }
      }
    }
  }


  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayView<const Scalar> MultiVector<Ordinal,Scalar>::operator[](Ordinal i) const
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return MVData_->cPtrs_[i];
  }

  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayView<Scalar> MultiVector<Ordinal,Scalar>::operator[](Ordinal i)
  {
    // teuchos does the bounds checking here, if TEUCHOS_DEBUG
    return MVData_->ptrs_[i];
  }

  template<typename Ordinal, typename Scalar>
  MultiVector<Ordinal,Scalar>& MultiVector<Ordinal,Scalar>::operator=(const MultiVector<Ordinal,Scalar> &source) {
    // Check for special case of this=Source
    if (this != &source) {
      TEST_FOR_EXCEPTION( !this->getMap().isCompatible(source.getMap()), std::runtime_error,
          "Tpetra::MultiVector::operator=(): MultiVectors must have compatible Maps.");
      if (constantStride() && source.constantStride() && myLength()==stride() && source.myLength()==source.stride()) {
        // can copy in one call
        std::copy( source.MVData_->values_.begin(), source.MVData_->values_.begin() + source.numVectors()*source.stride(),
                   MVData_->values_.begin() );
      }
      else {
        for (Ordinal j=0; j<numVectors(); ++j) {
          std::copy( source.MVData_->ptrs_[j].begin(), source.MVData_->ptrs_[j].end(), 
                     MVData_->ptrs_[j].begin() );
        }
      }
    }
    return(*this);
  }

  /*
  template<typename Ordinal, typename Scalar>
  Teuchos::RCP<MultiVector<Ordinal,Scalar> > MultiVector<Ordinal,Scalar>::subCopy(const Teuchos::Range1D &colRng) const
  {
    // FINISH
    TEST_FOR_EXCEPT(true);
    return Teuchos::null;
  }
  */

  template<typename Ordinal, typename Scalar>
  Teuchos::RCP<MultiVector<Ordinal,Scalar> > MultiVector<Ordinal,Scalar>::subCopy(const Teuchos::ArrayView<const Teuchos_Index> &cols) const
  {
    // TODO: in debug mode, do testing that cols[j] are distinct
    Ordinal numCols = cols.size();
    // allocate new MV
    const bool zeroData = false;
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > mv = rcp( new MultiVector<Ordinal,Scalar>(this->getMap(),numCols,zeroData) );
    // copy data from *this into mv
    for (Ordinal j=0; j<numCols; ++j)
    {
      std::copy( MVData_->ptrs_[cols[j]].begin(), MVData_->ptrs_[cols[j]].end(), mv->MVData_->ptrs_[j].begin() );
    }
    return mv;
  }

  /*
  template<typename Ordinal, typename Scalar>
  Teuchos::RCP<MultiVector<Ordinal,Scalar> > MultiVector<Ordinal,Scalar>::subView(const Teuchos::Range1D &colRng) 
  {
    // the range is contiguous, so if this multivector is contiguous, so will be the resulting view
    // make sure that its ArrayRCP contains only the data pertaining to that multivector
    // FINISH
    TEST_FOR_EXCEPT(true);
    return Teuchos::null;
  }
  */

  template<typename Ordinal, typename Scalar>
  Teuchos::RCP<MultiVector<Ordinal,Scalar> > MultiVector<Ordinal,Scalar>::subView(const Teuchos::ArrayView<const Teuchos_Index> &cols) 
  {
    using Teuchos::as;
    const Ordinal numVecs = cols.size();
    Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > mvdata = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    // FINISH: if cols are contiguous, we keep constantStride. in that case, values_ needs to point at the beginning.
    mvdata->constantStride_ = false;
    mvdata->values_ = MVData_->values_;
    mvdata->stride_ = stride();
    mvdata->ptrs_.resize(numVecs,Teuchos::null);
    for (Ordinal j = as<Ordinal>(0); j < numVecs; ++j) {
      mvdata->ptrs_[j] = MVData_->ptrs_[cols[j]];
    }
    mvdata->updateConstPointers();
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > mv = Teuchos::rcp( new MultiVector<Ordinal,Scalar>(this->getMap(),mvdata) );
    return mv;
  }

  /*
  template<typename Ordinal, typename Scalar>
  Teuchos::RCP<const MultiVector<Ordinal,Scalar> > MultiVector<Ordinal,Scalar>::subViewConst(const Teuchos::Range1D &colRng) const 
  {
    // the range is contiguous, so if this multivector is contiguous, so will be the resulting view
    // make sure that its ArrayRCP contains only the data pertaining to that multivector
    // FINISH
    TEST_FOR_EXCEPT(true);
    return Teuchos::null;
  }
  */

  template<typename Ordinal, typename Scalar>
  Teuchos::RCP<const MultiVector<Ordinal,Scalar> > MultiVector<Ordinal,Scalar>::subViewConst(const Teuchos::ArrayView<const Teuchos_Index> &cols) const
  {
    using Teuchos::as;
    const Ordinal numVecs = cols.size();
    Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > mvdata = Teuchos::rcp( new MultiVectorData<Ordinal,Scalar>() );
    mvdata->constantStride_ = false;
    mvdata->stride_ = stride();
    mvdata->values_ = MVData_->values_;
    mvdata->ptrs_.resize(numVecs,Teuchos::null);
    for (Ordinal j = as<Ordinal>(0); j < numVecs; ++j) {
      mvdata->ptrs_[j] = MVData_->ptrs_[cols[j]];
    }
    mvdata->updateConstPointers();
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > mv = Teuchos::rcp( new MultiVector<Ordinal,Scalar>(this->getMap(),mvdata) );
    return mv;
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::extractCopy(Teuchos::ArrayView<Scalar> A, Ordinal &MyLDA) const 
  {
    TEST_FOR_EXCEPTION(constantStride() == false, std::runtime_error,
      "MultiVector::extractCopy(A,LDA): only supported for constant stride multivectors.");
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(A.size() != stride()*numVectors(), std::runtime_error,
      "MultiVector::extractCopy(A,LDA): A must be large enough to hold contents of MultiVector.");
#endif
    MyLDA = stride();
    std::copy(MVData_->values_.begin(), MVData_->values_.begin()+stride()*numVectors(),
              A.begin());
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::extractCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > arrayOfArrays) const
  {
#ifdef TEUCHOS_DEBUG 
    TEST_FOR_EXCEPTION(arrayOfArrays.size() != numVectors(), std::runtime_error,
        "Tpetra::MultiVector::extractCopy(arrayOfArrays): arrayOfArrays.size() (==" << arrayOfArrays.size() << ") must match"
        "numVectors() (==" << numVectors() << ")");
    for (typename Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> >::iterator it=arrayOfArrays.begin(); 
         it != arrayOfArrays.end(); ++it) {
      TEST_FOR_EXCEPTION(it->size() != myLength(), std::runtime_error, 
          "Tpetra::MultiVector::extractCopy(arrayOfArrays): ArrayView's size (==" 
          << it->size() << ") must match local MultiVector length (==" << myLength() << ")");
    }
#endif
    // copy each multivector into the user provided 2-D array
    for (Ordinal i=0; i<numVectors(); ++i) {
      std::copy(MVData_->ptrs_[i].begin(), MVData_->ptrs_[i].end(), arrayOfArrays[i].begin());
    }
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::extractView(Teuchos::ArrayView<Scalar> &A, Ordinal &MyLDA)
  {
    TEST_FOR_EXCEPTION(constantStride() == false, std::runtime_error,
      "MultiVector::extractView(A,LDA): only supported for constant stride multivectors.");
    A = MVData_->values_;
    MyLDA = MVData_->stride_;
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::extractConstView(Teuchos::ArrayView<const Scalar> &A, Ordinal &MyLDA) const
  {
    TEST_FOR_EXCEPTION(constantStride() == false, std::runtime_error,
      "MultiVector::extractConstView(A,LDA): only supported for constant stride multivectors.");
    A = MVData_->values_.getConst();
    MyLDA = MVData_->stride_;
  }

  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > MultiVector<Ordinal,Scalar>::extractView()
  {
    return MVData_->ptrs_().getConst();
  }

  template<typename Ordinal, typename Scalar>
  Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > MultiVector<Ordinal,Scalar>::extractConstView() const
  {
    return MVData_->cPtrs_().getConst();
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta) 
  {
    // This routine performs a variety of matrix-matrix multiply operations, interpreting
    // the MultiVector (this-aka C , A and B) as 2D matrices.  Variations are due to
    // the fact that A, B and C can be local replicated or global distributed
    // MultiVectors and that we may or may not operate with the transpose of 
    // A and B.  Possible cases are:
    using Teuchos::NO_TRANS;      // enums
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::null;
    using Teuchos::ScalarTraits;  // traits
    using Teuchos::OrdinalTraits;
    using Teuchos::as;
    using Teuchos::RCP;           // data structures
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::rcp;           // initializers for data structures
    using Teuchos::arcp;

    //                                       Num
    //      OPERATIONS                        cases  Notes
    //  1) C(local) = A^X(local) * B^X(local)  4    (X=Trans or Not, No comm needed) 
    //  2) C(local) = A^T(distr) * B  (distr)  1    (2D dot product, replicate C)
    //  3) C(distr) = A  (distr) * B^X(local)  2    (2D vector update, no comm needed)
    //
    // The following operations are not meaningful for 1D distributions:
    //
    // u1) C(local) = A^T(distr) * B^T(distr)  1
    // u2) C(local) = A  (distr) * B^X(distr)  2
    // u3) C(distr) = A^X(local) * B^X(local)  4
    // u4) C(distr) = A^X(local) * B^X(distr)  4
    // u5) C(distr) = A^T(distr) * B^X(local)  2
    // u6) C(local) = A^X(distr) * B^X(local)  4
    // u7) C(distr) = A^X(distr) * B^X(local)  4
    // u8) C(local) = A^X(local) * B^X(distr)  4
    //
    // Total of 32 case (2^5).

    std::string errPrefix("Tpetra::MultiVector::multiply(transOpA,transOpB,A,B): ");

    TEST_FOR_EXCEPTION( ScalarTraits<Scalar>::isComplex && (transA == TRANS || transB == TRANS), std::invalid_argument,
        errPrefix << "non-conjugate transpose not supported for complex types.");
    transA = (transA == NO_TRANS ? NO_TRANS : CONJ_TRANS);
    transB = (transB == NO_TRANS ? NO_TRANS : CONJ_TRANS);

    // Compute effective dimensions, w.r.t. transpose operations on 
    Ordinal A_nrows = (transA==CONJ_TRANS) ? A.numVectors() : A.myLength();
    Ordinal A_ncols = (transA==CONJ_TRANS) ? A.myLength() : A.numVectors();
    Ordinal B_nrows = (transB==CONJ_TRANS) ? B.numVectors() : B.myLength();
    Ordinal B_ncols = (transB==CONJ_TRANS) ? B.myLength() : B.numVectors();

    Scalar beta_local = beta; // local copy of beta; might be reassigned below

    TEST_FOR_EXCEPTION( myLength() != A_nrows || numVectors() != B_ncols || A_ncols != B_nrows, std::runtime_error,
        errPrefix << "dimension of *this, op(A) and op(B) must be consistent.");

    bool A_is_local = !A.isDistributed();
    bool B_is_local = !B.isDistributed();
    bool C_is_local = !this->isDistributed();
    bool Case1 = ( C_is_local &&  A_is_local &&  B_is_local);                                           // Case 1: C(local) = A^X(local) * B^X(local)
    bool Case2 = ( C_is_local && !A_is_local && !B_is_local && transA==CONJ_TRANS && transB==NO_TRANS); // Case 2: C(local) = A^T(distr) * B  (distr)
    bool Case3 = (!C_is_local && !A_is_local &&  B_is_local && transA==NO_TRANS  );                     // Case 3: C(distr) = A  (distr) * B^X(local)

    // Test that we are considering a meaningful cases
    TEST_FOR_EXCEPTION( !Case1 && !Case2 && !Case3, std::runtime_error,
        errPrefix << "multiplication of op(A) and op(B) into *this is not a supported use case.");

    if (beta != ScalarTraits<Scalar>::zero() && Case2) // 
    {
      // if Case2, then C is local and contributions must be summed across all nodes
      // however, if beta != 0, then accumulate beta*C into the sum
      // when summing across all nodes, we only want to accumulate this once, so 
      // set beta == 0 on all nodes except node 0
      int MyPID = this->getMap().getComm()->getRank();
      if (MyPID!=0) beta_local = ScalarTraits<Scalar>::zero();
    }

    // Check if A, B, C have constant stride, if not then make temp copy (strided)
    RCP<const MultiVector<Ordinal,Scalar> > Atmp, Btmp; 
    RCP<MultiVector<Ordinal,Scalar> > Ctmp;
    if (constantStride() == false) Ctmp = rcp(new MultiVector<Ordinal,Scalar>(*this));
    else Ctmp = rcp(this,false);

    if (A.constantStride() == false) Atmp = rcp(new MultiVector<Ordinal,Scalar>(A));
    else Atmp = rcp(&A,false);

    if (B.constantStride() == false) Btmp = rcp(new MultiVector<Ordinal,Scalar>(B));
    else Btmp = rcp(&B,false);

#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(!Ctmp->constantStride() || !Btmp->constantStride() || !Atmp->constantStride(), std::logic_error,
        errPrefix << "failed making temporary strided copies of input multivectors.");
#endif

    Ordinal m = this->myLength();
    Ordinal n = this->numVectors();
    Ordinal k = A_ncols;
    Ordinal lda, ldb, ldc;
    ArrayView<const Scalar> Ap(Teuchos::null), Bp(Teuchos::null);
    ArrayView<Scalar> Cp(Teuchos::null);
    Atmp->extractConstView(Ap,lda);
    Btmp->extractConstView(Bp,ldb);
    Ctmp->extractView(Cp,ldc);

    Teuchos::BLAS<int,Scalar> blas;
    // do the arithmetic now
    blas.GEMM(transA,transB,m,n,k,alpha,Ap.getRawPtr(),lda,Bp.getRawPtr(),ldb,beta_local,Cp.getRawPtr(),ldc);

    // Dispose of (possibly) extra copies of A, B
    Atmp = null;
    Btmp = null;

    // If *this was not strided, copy the data from the strided version and then delete it
    if (constantStride() == false) {
      Array<ArrayView<Scalar> > aoa(MVData_->ptrs_.size(),null);
      for (Ordinal i=0; i<as<Ordinal>(aoa.size()); ++i) {
        aoa[i] = MVData_->ptrs_[i]();
      }
      Ctmp->extractCopy(aoa());
    }
    Ctmp = null;

    // If Case 2 then sum up C and distribute it to all processors.
    if (Case2) 
    {
      RCP<const Teuchos::Comm<Ordinal> > comm = this->getMap().getComm();
      // Global reduction on each entry of a Replicated Local MultiVector
      // Comm requires that local and global buffers be congruous and distinct
      // Therefore, we must allocate storage for the local values
      // Furthermore, if the storage in C (our destination for the global results)
      //   is not packed, we must allocate storage for the result as well.
      ArrayRCP<Scalar> source = arcp<Scalar>(m*n), target;
      bool packed = constantStride() && (stride() == m);
      if (packed) {
        // copy local info into source buffer
        // target buffer will be multivector storage
        std::copy(MVData_->values_.begin(),MVData_->values_.begin()+m*n,
                  source.begin());
        // local storage is packed. can use it for target buffer.
        target = MVData_->values_;
      }
      else {
        // copy local info into source buffer
        ArrayRCP<Scalar> sptr = source;
        for (Ordinal j=OrdinalTraits<Ordinal>::zero(); j<n; ++j) 
        {
          // copy j-th local MV data into source buffer
          std::copy(MVData_->ptrs_[j].begin(),MVData_->ptrs_[j].begin()+m,
                    sptr.begin());
          // advance ptr into source buffer
          sptr += m;
        }
        // must allocate packed storage for target buffer
        target = arcp<Scalar>(m*n);
      }
      // reduce 
      Teuchos::reduceAll<Ordinal,Scalar>(*comm,Teuchos::REDUCE_SUM,m*n,source.getRawPtr(),target.getRawPtr());
      if (!packed) {
        // copy target buffer into multivector storage buffer
        ArrayRCP<Scalar> tptr = target;
        for (Ordinal j=OrdinalTraits<Ordinal>::zero(); j<n; ++j) 
        {
          std::copy(tptr.begin(),tptr.begin()+m,
                    MVData_->ptrs_[j].begin()
                   );
          tptr += m;
        }
      }
      // clear allocated buffers
      source = null;
      target = null;
    } // Case2 reduction
  }


  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::reduce() 
  {
    using Teuchos::as;
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
    // this should be called only for "local" MultiVectors (!isDistributed())
    TEST_FOR_EXCEPTION(this->isDistributed() == true, std::runtime_error,
        "Tpetra::MultiVector::reduce() should only be called for non-distributed MultiVectors.");
    // sum the data across all multivectors
    // need to have separate (contiguous) buffers for send and receive
    // if we're contiguous, we'll just set the receive buffer as our data, the send as a copy
    // if we're non-contig, we'll use separate buffers for both. 
    Ordinal bufsize = myLength() * numVectors();
    Teuchos::ArrayRCP<Scalar> sndbuf, rcvbuf;
    if (constantStride() == true) {
      sndbuf = Teuchos::arcpClone<Scalar>(MVData_->values_());
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(sndbuf.size() != MVData_->values_.size(), std::logic_error,
          "Tpetra::MultiVector::reduce(): Error in Tpetra. Please contact Tpetra developers.");
#endif
      rcvbuf = MVData_->values_;
    }
    else {
      sndbuf = Teuchos::arcp<Scalar>(bufsize);
      rcvbuf = Teuchos::arcp<Scalar>(bufsize);
      Teuchos::ArrayRCP<Scalar> sndptr = sndbuf;
      for (Ordinal j=OT::zero(); j<numVectors(); ++j) {
        Teuchos::ArrayView<const Scalar> jvec = (*this)[j];
        std::copy(jvec.begin(),jvec.end(),sndptr);
        sndptr += myLength();
      }
    }
    Teuchos::reduceAll(*(this->getMap().getComm()), Teuchos::REDUCE_SUM, bufsize, sndbuf().getConst().getRawPtr(), rcvbuf.getRawPtr() );
    if (constantStride() == false) {
      Teuchos::ArrayRCP<const Scalar> rcvptr = rcvbuf;
      for (Ordinal j=OT::zero(); j<numVectors(); ++j) {
        Teuchos::ArrayView<Scalar> jvec = (*this)[j];
        std::copy(rcvptr.begin(),rcvptr.begin()+myLength(),jvec.begin());
        rcvptr += myLength();
      }
    }
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::replaceMap(const Map<Ordinal> &map)
  {
    TEST_FOR_EXCEPT(true);
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::replaceMyValue(Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(MyRow < this->getMap().getMinLocalIndex() || MyRow > this->getMap().getMaxLocalIndex(), std::runtime_error,
        "Tpetra::MultiVector::replaceMyValue(): row index is invalid.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::replaceMyValue(): vector index is invalid.");
#endif
    MVData_->ptrs_[VectorIndex][MyRow] = ScalarValue;
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::sumIntoMyValue(Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(MyRow < this->getMap().getMinLocalIndex() || MyRow > this->getMap().getMaxLocalIndex(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoMyValue(): row index is invalid.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoMyValue(): vector index is invalid.");
#endif
    MVData_->ptrs_[VectorIndex][MyRow] += ScalarValue;
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::replaceGlobalValue(Ordinal GlobalRow, Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(!this->getMap().isMyGlobalIndex(GlobalRow), std::runtime_error,
        "Tpetra::MultiVector::replaceGlobalValue(): row index is not present on this processor.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::replaceGlobalValue(): vector index is invalid.");
#endif
    Ordinal MyRow = this->getMap().getLocalIndex(GlobalRow);
    MVData_->ptrs_[VectorIndex][MyRow] = ScalarValue;
  }

  template<typename Ordinal, typename Scalar>
  void MultiVector<Ordinal,Scalar>::sumIntoGlobalValue(Ordinal GlobalRow, Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(!this->getMap().isMyGlobalIndex(GlobalRow), std::runtime_error,
        "Tpetra::MultiVector::sumIntoGlobalValue(): row index is not present on this processor.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoGlobalValue(): vector index is invalid.");
#endif
    Ordinal MyRow = this->getMap().getLocalIndex(GlobalRow);
    MVData_->ptrs_[VectorIndex][MyRow] += ScalarValue;
  }

} // namespace Tpetra

#endif // TPETRA_MULTIVECTOR_HPP
