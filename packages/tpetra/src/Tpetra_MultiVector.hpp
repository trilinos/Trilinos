// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_MULTIVECTOR_HPP
#define TPETRA_MULTIVECTOR_HPP

#include <Kokkos_DefaultArithmetic.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Tpetra_MultiVectorDecl.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos_Ordinal NumVectors, bool zeroOut) 
  : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map), lclMV_(node) {
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const LocalOrdinal myLen = getMyLength();
    if (myLen > 0) {
      Teuchos::ArrayRCP<Scalar> data = node.template allocBuffer<Scalar>(myLen*NumVectors);
      lclMV_.initializeValues(myLen,NumVectors,data,myLen);
      if (zeroOut) {
        DMVA::Init(lclMV_,0.0);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) 
  : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(source), lclMV_(source.lclMV_.getNode()) {
    // copy data from the source MultiVector into this multivector
    Node &node = source.lclMV_.getNode();
    const LocalOrdinal myLen = getMyLength();
    const Kokkos::size_type numVecs = source.getNumVectors();
    if (myLen > 0) {
      // allocate data
      Teuchos::ArrayRCP<Scalar> data = node.template allocBuffer<Scalar>(myLen*numVecs);
      lclMV_.initializeValues(myLen,numVecs,data,myLen);
      // copy data
      {
        Teuchos::ArrayRCP<Scalar> dstdata = data;
        Teuchos::ArrayRCP<const Scalar> srcdata = source.lclMV_.getValues();
        for (Kokkos::size_type j = 0; j < numVecs; ++j) {
          Teuchos::ArrayRCP<const Scalar> srcj = source.getSubArrayRCP(srcdata,j);
          node.template copyBuffers<Scalar>(myLen,srcj,dstdata);
          dstdata += myLen;
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(
                        Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, 
                        const Teuchos::ArrayView<const Scalar> &A, Teuchos_Ordinal LDA, 
                        Teuchos_Ordinal NumVectors)
  : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map), lclMV_(node) {
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const Teuchos_Ordinal myLen = getMyLength();
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(LDA < myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(): LDA must be large enough to accomodate the local entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(A.size() < LDA*(NumVectors-1)+myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(A,LDA): A does not contain enough data to specify the entries in this.");
#endif
    if (myLen > 0) {
      Teuchos::ArrayRCP<Scalar> mydata = node.template allocBuffer<Scalar>(myLen*NumVectors);
      lclMV_.initializeValues(myLen,NumVectors,mydata,myLen);
      Teuchos::ArrayRCP<Scalar> myview = node.template viewBufferNonConst<Scalar>(true,myLen*NumVectors,mydata);
      typename Teuchos::ArrayView<const Scalar>::iterator srcit = A.begin();
      for (Teuchos_Ordinal j = 0; j < NumVectors; ++j) {
        std::copy(srcit,srcit+myLen,myview);
        srcit += LDA;
        myview += myLen;
      }
      mydata = Teuchos::null;
      myview = Teuchos::null;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, 
                                                                   const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, 
                                                                   Teuchos_Ordinal NumVectors)
    : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map), lclMV_(node)
  {
    TEST_FOR_EXCEPTION(NumVectors < 1 || NumVectors != ArrayOfPtrs.size(), std::runtime_error,
        "Tpetra::MultiVector::MultiVector(ArrayOfPtrs): ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
    const Teuchos_Ordinal myLen = getMyLength();
    if (myLen > 0) {
      Teuchos::ArrayRCP<Scalar> mydata = node.template allocBuffer<Scalar>(myLen*NumVectors);
      lclMV_.initializeValues(myLen,NumVectors,mydata,myLen);
      Teuchos::ArrayRCP<Scalar> myview = node.template viewBufferNonConst<Scalar>(true,myLen*NumVectors,mydata);
      for (Teuchos_Ordinal j = 0; j < NumVectors; ++j) {
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION(ArrayOfPtrs[j].size() != getMyLength(), std::runtime_error,
          "Tpetra::MultiVector::MultiVector(ArrayOfPtrs): ArrayOfPtrs[" << j << "].size() (== " << ArrayOfPtrs[j].size() << 
          ") is not equal to getMyLength() (== " << getMyLength());
#endif
        typename Teuchos::ArrayView<const Scalar>::iterator src = ArrayOfPtrs[j].begin();
        std::copy(src,src+myLen,myview);
        myview += myLen;
      }
      myview = Teuchos::null;
      mydata = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map,
              Teuchos::ArrayRCP<Scalar> data, Teuchos_Ordinal LDA, Teuchos::ArrayView<const Teuchos_Ordinal> WhichVectors)
  : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map), lclMV_(node), whichVectors_(WhichVectors) {
    const Teuchos_Ordinal myLen = getMyLength();
    Teuchos_Ordinal maxVector = *std::max_element(WhichVectors.begin(), WhichVectors.end());
    TEST_FOR_EXCEPTION(LDA < myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(): LDA must be large enough to accomodate the local entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(data.size() < LDA * maxVector + myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(data,LDA,WhichVectors): data does not contain enough data to specify the entries in this.");
#endif
    if (WhichVectors.size() == 1) {
      // shift data so that desired vector is vector 0
      maxVector = 0;
      data += LDA*WhichVectors[0];
      // kill whichVectors_; we are constant stride
      whichVectors_.clear();
    }
    if (myLen > 0) {
      lclMV_.initializeValues(myLen,maxVector+1,data,LDA);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map,
                Teuchos::ArrayRCP<Scalar> data, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors)
  : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map), lclMV_(node) {
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const LocalOrdinal myLen = getMyLength();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(data.size() < LDA*(NumVectors-1)+myLen, std::runtime_error,
        "Tpetra::MultiVector::MultiVector(): data does not contain enough data to specify the entries in this.");
#endif
    lclMV_.initializeValues(myLen,NumVectors,data,LDA);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~MultiVector() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isConstantStride() const {
    return whichVectors_.empty();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getMyLength() const {
    return this->getMap().getNumMyEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalLength() const {
    return this->getMap().getNumGlobalEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos_Ordinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getStride() const {
    if (isConstantStride()) {
      return lclMV_.getStride();
    }
    return 0;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj, Teuchos_Ordinal &packetSize) 
  {
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(sourceObj);
    // objects maps have already been checked. simply check the number of vectors.
    packetSize = this->getNumVectors();
    return (A.getNumVectors() == packetSize);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::copyAndPermute(
                          const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & sourceObj,
                          Teuchos_Ordinal numSameIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs) {
    TEST_FOR_EXCEPT(!isConstantStride());
    using Teuchos::ArrayView;
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceMV = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &>(sourceObj);
    typename ArrayView<const LocalOrdinal>::iterator pTo, pFrom;
    // any other error will be caught by Teuchos
    TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::copyAndPermute(): permuteToLIDs and permuteFromLIDs must have the same size.");
    // one vector at a time
    Node &node = lclMV_.getNode();
    const int numCols   = getNumVectors(),
              srcStride = sourceMV.lclMV_.getStride(),
              dstStride =          lclMV_.getStride();
    Teuchos::ArrayRCP<Scalar> dstptr;
    Teuchos::ArrayRCP<const Scalar> srcptr;
    // Get a host view of the local Kokkos::MultiVector data:
    //   Kokkos::MultiVector::getValues(j) returns a parallel buffer pointing to 
    //   the j-th column; viewBuffer/Const() will return a host pointer
    // TODO: determine whether this viewBuffer is write-only or not; for now, safe option is not
    dstptr =  node.template viewBufferNonConst<Scalar>(false,dstStride*numCols,lclMV_.getValuesNonConst());
    srcptr =  node.template viewBuffer<Scalar>( srcStride*numCols,sourceMV.lclMV_.getValues());
    for (Teuchos_Ordinal j = 0; j < numCols; ++j) {
      // The first numImportIDs GIDs are the same between source and target,
      // We can just copy them
      std::copy(srcptr,srcptr+numSameIDs,dstptr);
      // next, do permutations
      for (pTo = permuteToLIDs.begin(), pFrom = permuteFromLIDs.begin();
           pTo != permuteToLIDs.end(); ++pTo, ++pFrom)
      {
        dstptr[*pTo] = srcptr[*pFrom];
      }
      dstptr += dstStride;
      srcptr += srcStride;
    }
    dstptr = Teuchos::null;
    srcptr = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::packAndPrepare(
          const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & sourceObj,
          const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
          const Teuchos::ArrayView<Scalar> &exports,
          Distributor &distor) {
    TEST_FOR_EXCEPT(!isConstantStride());
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceMV = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &>(sourceObj);
    using Teuchos::ArrayView;
    (void)distor;    // we don't use these, but we don't want unused parameter warnings
    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ; 
                   ....
                   all of the data from row exportLIDs.back() }
      this doesn't have the best locality, but is necessary because the data for a Packet
      (all data associated with an LID) is required to be contiguous */
    TEST_FOR_EXCEPTION(exports.size() != getNumVectors()*exportLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::packAndPrepare(): sizing of exports buffer should be appropriate for the amount of data to be exported.");
    const KMV &srcData = sourceMV.lclMV_;
    const int myStride = srcData.getStride(),
              numCols   = getNumVectors();
    typename ArrayView<const LocalOrdinal>::iterator idptr;
    typename ArrayView<Scalar>::iterator expptr;
    expptr = exports.begin();
    Node &node = srcData.getNode();
    Teuchos::ArrayRCP<const Scalar> myptr;
    myptr = node.template viewBuffer<Scalar>(myStride*numCols,srcData.getValues());
    for (idptr = exportLIDs.begin(); idptr != exportLIDs.end(); ++idptr) {
      for (Teuchos_Ordinal j = 0; j < numCols; ++j) {
        *expptr++ = myptr[j*myStride + (*idptr)];
      }
    }
    myptr = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::unpackAndCombine(
                  const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                  const Teuchos::ArrayView<const Scalar> &imports,
                  Distributor &distor,
                  CombineMode CM) {
    TEST_FOR_EXCEPT(!isConstantStride());
    (void)distor; // we don't use this, but we don't want unused parameter warnings
    using Teuchos::ArrayView;
    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ; 
                   ....
                   all of the data from row exportLIDs.back() }
      this doesn't have the best locality, but is necessary because the data for a Packet
      (all data associated with an LID) is required to be contiguous */
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(imports.size() != getNumVectors()*importLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::unpackAndCombine(): sizing of imports buffer should be appropriate for the amount of data to be exported.");
#endif
    const int myStride = lclMV_.getStride(),
              numCols   = getNumVectors();
    Node &node = lclMV_.getNode();
    Teuchos::ArrayRCP<Scalar> myptr;
    // TODO: determine whether this viewBuffer is write-only or not; for now, safe option is not
    myptr = node.template viewBufferNonConst<Scalar>(false,myStride*numCols,lclMV_.getValuesNonConst());
    typename ArrayView<const       Scalar>::iterator impptr;
    typename ArrayView<const LocalOrdinal>::iterator  idptr;
    impptr = imports.begin();
    if (CM == INSERT || CM == REPLACE) 
    {
      for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
        for (Teuchos_Ordinal j = 0; j < getNumVectors(); ++j) {
          myptr[j*myStride + *idptr] = *impptr++;
        }
      }
    }
    else if (CM == ADD) {
      for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
        for (Teuchos_Ordinal j = 0; j < getNumVectors(); ++j) {
          myptr[j*myStride + *idptr] += *impptr++;
        }
      }
    }
    else {
      TEST_FOR_EXCEPTION(CM != ADD && CM != REPLACE && CM != INSERT, std::invalid_argument,
          "Tpetra::MultiVector::unpackAndCombine(): Invalid CombineMode: " << CM);
    }
    myptr = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  inline Teuchos_Ordinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumVectors() const {
    if (isConstantStride()) {
      return lclMV_.getNumCols();
    }
    else {
      return whichVectors_.size();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot(
      const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, 
      const Teuchos::ArrayView<Scalar> &dots) const 
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    const Teuchos_Ordinal numVecs = getNumVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != A.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.getNumVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors must have the same number of vectors.");
    TEST_FOR_EXCEPTION(dots.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::dots(A,dots): dots.size() must be as large as the number of vectors in *this and A.");
    DMVA::Dot(lclMV_,A.lclMV_,dots);
    if (this->isDistributed()) {
      Teuchos::Array<Scalar> ldots(dots);
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,ldots.getRawPtr(),dots.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2(
                  const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Teuchos_Ordinal numVecs = this->getNumVectors();
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm2(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride()) {
      DMVA::Norm2Squared(lclMV_,norms);
    }
    else {
      KMV v(lclMV_.getNode());
      Teuchos::ArrayRCP<Scalar> vi;
      for (Teuchos_Ordinal i=0; i < numVecs; ++i) {
        vi = Teuchos::arcp_const_cast<Scalar>( lclMV_.getValues(whichVectors_[i]) );
        v.initializeValues(lclMV_.getNumRows(), 1, vi, lclMV_.getStride());
        norms[i] = DMVA::Norm2Squared(v);
      }
    }
    if (this->isDistributed()) {
      Teuchos::Array<Mag> lnorms(norms);
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      (*n) = ScalarTraits<Mag>::squareroot(*n);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted(
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights,
          const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
    TEST_FOR_EXCEPT(!isConstantStride());
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;
    const Mag OneOverN = ScalarTraits<Mag>::one() / Teuchos::as<Mag>(getGlobalLength());
    bool OneW = false;
    const Teuchos_Ordinal numVecs = this->getNumVectors();
    if (weights.getNumVectors() == 1) {
      OneW = true;
    }
    else {
      TEST_FOR_EXCEPTION(weights.getNumVectors() != numVecs, std::runtime_error,
          "Tpetra::MultiVector::normWeighted(): MultiVector of weights must contain either one vector or the same number of vectors as this.");
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
        "Tpetra::MultiVector::normWeighted(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "weights.getMap(): " << std::endl << weights.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != weights.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::normWeighted(): MultiVectors do not have the same local length.");
#endif
    // 
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::normWeighted(): norms.size() must be as large as the number of vectors in *this.");
    DMVA::WeightedNorm(lclMV_,weights.lclMV_,norms);
    if (this->isDistributed()) {
      Teuchos::Array<Mag> lnorms(norms);
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n * OneOverN);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1(
                  const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Teuchos_Ordinal numVecs = this->getNumVectors();
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm1(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride()) {
      DMVA::Norm1(lclMV_,norms);
    }
    else {
      TEST_FOR_EXCEPT(true);
      KMV v(lclMV_.getNode());
      Teuchos::ArrayRCP<Scalar> vj;
      for (Teuchos_Ordinal j=0; j < numVecs; ++j) {
        vj = Teuchos::arcp_const_cast<Scalar>( lclMV_.getValues(whichVectors_[j]) );
        v.initializeValues(lclMV_.getNumRows(), 1, vj, lclMV_.getStride());
        norms[j] = DMVA::Norm1((const KMV&)v);
      }
    }
    if (this->isDistributed()) {
      Teuchos::Array<Mag> lnorms(norms);
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf(
        const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Teuchos_Ordinal numVecs = this->getNumVectors();
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::normInf(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride()) {
      DMVA::NormInf(lclMV_,norms);
    }
    else {
      KMV v(lclMV_.getNode());
      Teuchos::ArrayRCP<Scalar> vj;
      for (Teuchos_Ordinal j=0; j < numVecs; ++j) {
        vj = Teuchos::arcp_const_cast<Scalar>( lclMV_.getValues(whichVectors_[j]) );
        v.initializeValues(lclMV_.getNumRows(), 1, vj, lclMV_.getStride());
        norms[j] = DMVA::NormInf((const KMV&)v);
      }
    }
    if (this->isDistributed()) {
      Teuchos::Array<Mag> lnorms(norms);
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::randomize() {
    if (isConstantStride()) {
      DMVA::Random(lclMV_);
    }
    else {
      TEST_FOR_EXCEPT(true);
      const Teuchos_Ordinal numVecs = this->getNumVectors();
      KMV v(lclMV_.getNode());
      Teuchos::ArrayRCP<Scalar> vj;
      for (Teuchos_Ordinal j=0; j < numVecs; ++j) {
        vj = lclMV_.getValuesNonConst(whichVectors_[j]);
        v.initializeValues(lclMV_.getNumRows(), 1, vj, lclMV_.getStride());
        DMVA::Random(v);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::putScalar(const Scalar &alpha) {
    const Teuchos_Ordinal numVecs = getNumVectors();
    if (isConstantStride()) {
      DMVA::Init(lclMV_,alpha);
    }
    else {
      KMV v(lclMV_.getNode());
      Teuchos::ArrayRCP<Scalar> vj;
      for (Teuchos_Ordinal j=0; j < numVecs; ++j) {
        vj = lclMV_.getValuesNonConst(whichVectors_[j]);
        v.initializeValues(lclMV_.getNumRows(), 1, vj, lclMV_.getStride());
        DMVA::Init(v,alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(const Scalar &alpha) {
    // NOTE: can't substitute putScalar(0.0) for scale(0.0), because 
    //       the former will overwrite NaNs present in the MultiVector, while the 
    //       semantics of this call require multiplying them by 0, which IEEE requires to be NaN
    const Teuchos_Ordinal numVecs = getNumVectors();
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // do nothing
    }
    else if (isConstantStride()) {
      DMVA::Scale(lclMV_,alpha);
    }
    else {
      TEST_FOR_EXCEPT(true);
      KMV v(lclMV_.getNode());
      Teuchos::ArrayRCP<Scalar> vj;
      for (Teuchos_Ordinal j=0; j < numVecs; ++j) {
        vj = Teuchos::arcp_const_cast<Scalar>( lclMV_.getValues(whichVectors_[j]) );
        v.initializeValues(lclMV_.getNumRows(), 1, vj, lclMV_.getStride());
        DMVA::Scale(v,alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(Teuchos::ArrayView<const Scalar> alphas)
  {
    const Teuchos_Ordinal numVecs = this->getNumVectors();
    TEST_FOR_EXCEPTION(alphas.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as the number of vectors in *this.");
    KMV vec(lclMV_.getNode());
    for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
      if (alphas[j] == Teuchos::ScalarTraits<Scalar>::one()) {
        // do nothing: NaN * 1.0 == NaN, Number*1.0 == Number
      }
      else {
        vec.initializeValues(getMyLength(),1,lclMV_.getValuesNonConst(j),getMyLength());
        DMVA::Scale(vec,alphas[j]);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) {
    TEST_FOR_EXCEPT(!isConstantStride());
    const Teuchos_Ordinal numVecs = this->getNumVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != A.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.getNumVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have the same number of vectors.");
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // set me = A
      DMVA::Assign(lclMV_,(const KMV&)A.lclMV_);
    }
    else {
      // set me == alpha*A
      DMVA::Scale(lclMV_,alpha,(const KMV&)A.lclMV_);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) {
    TEST_FOR_EXCEPT(!isConstantStride());
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != A.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::reciprocal(): MultiVectors must have the same number of vectors.");
    try {
      DMVA::Divide(lclMV_,(const KMV&)A.lclMV_);
    }
    catch (std::runtime_error &e) {
      TEST_FOR_EXCEPTION(true,std::runtime_error,
          "Tpetra::MultiVector::reciprocal(A): caught exception from Kokkos:" << std::endl
          << e.what() << std::endl);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) {
    TEST_FOR_EXCEPT(!isConstantStride());
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != A.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have the same number of vectors.");
    DMVA::Abs(lclMV_,(const KMV&)A.lclMV_);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, 
                      const Scalar &beta) {
    TEST_FOR_EXCEPT(!isConstantStride());
    // this = beta*this + alpha*A
    // must support case where &this == &A
    // can't short circuit on alpha==0.0 or beta==0.0, because 0.0*NaN != 0.0
    using Teuchos::ArrayView;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != A.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    DMVA::GESUM(lclMV_,alpha,(const KMV&)A.lclMV_,beta);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, 
                      const Scalar &beta,  const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, 
                      const Scalar &gamma) {
    TEST_FOR_EXCEPT(!isConstantStride());
    // this = alpha*A + beta*B + gamma*this
    // must support case where &this == &A or &this == &B
    // can't short circuit on alpha==0.0 or beta==0.0 or gamma==0.0, because 0.0*NaN != 0.0
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()) || !this->getMap().isCompatible(B.getMap()),
        std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl
        << "B.getMap(): " << std::endl << B.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( getMyLength() != A.getMyLength() || getMyLength() != B.getMyLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.getNumVectors() != this->getNumVectors() || B.getNumVectors() != this->getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    DMVA::GESUM(lclMV_,alpha,(const KMV&)A.lclMV_,beta,(const KMV&)B.lclMV_,gamma);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getData(Teuchos_Ordinal j) const
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    Node &node = lclMV_.getNode();
    return node.template viewBuffer<Scalar>(getMyLength(), lclMV_.getValues(j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDataNonConst(Teuchos_Ordinal j)
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    Node &node = lclMV_.getNode();
    return node.template viewBufferNonConst<Scalar>(false, getMyLength(), lclMV_.getValuesNonConst(j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator=(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) 
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    TEST_FOR_EXCEPT(!source.isConstantStride());
    // Check for special case of this=Source, in which case we do nothing
    if (this != &source) {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( !this->getMap().isCompatible(source.getMap()), std::runtime_error,
          "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
          << "this->getMap(): " << std::endl << this->getMap() 
          << "source.getMap(): " << std::endl << source.getMap() << std::endl);
#else
      TEST_FOR_EXCEPTION( getMyLength() != source.getMyLength(), std::runtime_error,
          "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
      TEST_FOR_EXCEPTION(source.getNumVectors() != getNumVectors(), std::runtime_error,
          "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
      Node &node = lclMV_.getNode();
      if (isConstantStride() && source.isConstantStride() && getMyLength()==getStride() && source.getMyLength()==source.getStride()) {
        // we're both packed, we can copy in one call
        node.template copyBuffers<Scalar>(getMyLength()*getNumVectors(), source.lclMV_.getValues(), lclMV_.getValuesNonConst() );
      }
      else {
        for (Kokkos::size_type j=0; j < lclMV_.getNumCols(); ++j) {
          node.template copyBuffers<Scalar>(getMyLength(), source.lclMV_.getValues(j), lclMV_.getValuesNonConst(j) );
        }
      }
    }
    return(*this);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subCopy(const Teuchos::ArrayView<const Teuchos_Ordinal> &cols) const
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    TEST_FOR_EXCEPTION(cols.size() < 1, std::runtime_error,
        "Tpetra::MultiVector::subCopy(cols): cols must contain at least one column.");
    Teuchos_Ordinal numCopyVecs = cols.size();
    const bool zeroData = false;
    Node &node = lclMV_.getNode();
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv; 
    mv = Teuchos::rcp( new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(node,this->getMap(),numCopyVecs,zeroData) );
    // copy data from *this into mv
    for (Teuchos_Ordinal j=0; j<numCopyVecs; ++j) {
      node.template copyBuffers<Scalar>( getMyLength(), lclMV_.getValues(cols[j]), mv->lclMV_.getValuesNonConst(j) );
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subCopy(const Teuchos::Range1D &colRng) const {
    TEST_FOR_EXCEPT(!isConstantStride());
    TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subCopy(Range1D): range must include at least one vector.");
    Teuchos_Ordinal numCopyVecs = colRng.size();
    const bool zeroData = false;
    Node &node = lclMV_.getNode();
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv; 
    mv = Teuchos::rcp( new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(node,this->getMap(),numCopyVecs,zeroData) );
    // copy data from *this into mv
    for (Teuchos_Ordinal js=colRng.lbound(), jd=0; jd<numCopyVecs; ++jd, ++js) {
      node.template copyBuffers<Scalar>( getMyLength(), lclMV_.getValues(js), mv->lclMV_.getValuesNonConst(jd) );
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subView(Teuchos::ArrayView<const Teuchos_Ordinal> cols) const {
    TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subView(ArrayView): range must include at least one vector.");
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because 
    //   a) the constructor doesn't modify the data, and 
    //   b) we are encapsulating in a const MV before returning
    Teuchos::ArrayRCP<const Scalar> cbuf = lclMV_.getValues();
    Teuchos::ArrayRCP<Scalar>      ncbuf = Teuchos::arcp_const_cast<Scalar>(cbuf);
    if (isConstantStride()) {
      return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),
                                                                                  ncbuf,lclMV_.getStride(),
                                                                                  cols) );
    }
    // else, lookup current whichVectors_ using cols
    Teuchos::Array<Teuchos_Ordinal> newcols(cols.size());
    for (Teuchos_Ordinal j=0; j < cols.size(); ++j) {
      newcols[j] = whichVectors_[cols[j]];
    }
    return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),
                                                                                ncbuf,lclMV_.getStride(),
                                                                                newcols()) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subView(const Teuchos::Range1D &colRng) const {
    TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
    Teuchos_Ordinal numViewVecs = colRng.size();
    using Teuchos::ArrayRCP;
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because 
    //   a) the constructor doesn't modify the data, and 
    //   b) we are encapsulating in a const MV before returning
    ArrayRCP<const Scalar> cbuf = lclMV_.getValues();
    ArrayRCP<Scalar>      ncbuf = Teuchos::arcp_const_cast<Scalar>(cbuf);
    // resulting MultiVector is constant stride only if *this is 
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      ArrayRCP<Scalar> subdata = ncbuf.persistingView( getStride() * colRng.lbound(),
                                                       getStride() * (numViewVecs-1) + getMyLength() );
      return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),
                                                                                  subdata,lclMV_.getStride(),numViewVecs) );
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Teuchos::Array<Teuchos_Ordinal> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),  
                                                                                ncbuf,lclMV_.getStride(),whchvecs) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subViewNonConst(Teuchos::ArrayView<const Teuchos_Ordinal> cols) {
    TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subViewNonConst(ArrayView): range must include at least one vector.");
    if (isConstantStride()) {
      return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),
                                                                                  lclMV_.getValuesNonConst(),lclMV_.getStride(),
                                                                                  cols) );
    }
    // else, lookup current whichVectors_ using cols
    Teuchos::Array<Teuchos_Ordinal> newcols(cols.size());
    for (Teuchos_Ordinal j=0; j < cols.size(); ++j) {
      newcols[j] = whichVectors_[cols[j]];
    }
    return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),
                                                                                lclMV_.getValuesNonConst(),lclMV_.getStride(),
                                                                                newcols()) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subViewNonConst(const Teuchos::Range1D &colRng) {
    TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subViewNonConst(Range1D): range must include at least one vector.");
    Teuchos_Ordinal numViewVecs = colRng.size();
    using Teuchos::ArrayRCP;
    // resulting MultiVector is constant stride only if *this is 
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      ArrayRCP<Scalar> data = lclMV_.getValuesNonConst();
      ArrayRCP<Scalar> subdata = data.persistingView( getStride() * colRng.lbound(),
                                                      getStride() * (numViewVecs-1) + getMyLength() );
      return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),
                                                                                  subdata,getStride(),numViewVecs) );
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Teuchos::Array<Teuchos_Ordinal> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    ArrayRCP<Scalar> data = lclMV_.getValuesNonConst();
    return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),  
                                                                                data,getStride(),whchvecs) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getVector(Teuchos_Ordinal j) const {
    TEST_FOR_EXCEPT(!isConstantStride());
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(j < 0 || j >= this->getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::operator()(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif
    // this is const, so lclMV_ is const, so we get const buff
    // it is safe to cast away the const because we will wrap it in a const Vector below
    Teuchos::ArrayRCP<Scalar> ncbuff;
    if (getMyLength() > 0) {
      Teuchos::ArrayRCP<const Scalar> cbuff = lclMV_.getValues(j);
      ncbuff = Teuchos::arcp_const_cast<Scalar>(cbuff);
    }
    return Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),ncbuff));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getVectorNonConst(Teuchos_Ordinal j) {
    TEST_FOR_EXCEPT(!isConstantStride());
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(j < 0 || j >= this->getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::operator()(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif
    Teuchos::ArrayRCP<Scalar> ncbuff;
    if (getMyLength() > 0) {
      ncbuff = lclMV_.getValuesNonConst(j);
    }
    return Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lclMV_.getNode(),this->getMap(),ncbuff));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy(Teuchos::ArrayView<Scalar> A, Teuchos_Ordinal LDA) const
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(LDA < getMyLength(), std::runtime_error,
      "Tpetra::MultiVector::get1dCopy(A,LDA): specified stride is not large enough for local vector length.");
    TEST_FOR_EXCEPTION(A.size() < LDA*(getNumVectors()-1)+getMyLength(), std::runtime_error,
      "Tpetra::MultiVector::get1dCopy(A,LDA): specified stride/storage is not large enough for the number of vectors.");
    Node &node = lclMV_.getNode();
    const int myStride = lclMV_.getStride(),
              numCols = getNumVectors(),
              myLen   = getMyLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mydata = lclMV_.getValues();
      ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,mydata);
      typename Teuchos::ArrayView<Scalar>::iterator Aptr = A.begin();
      for (Teuchos_Ordinal j=0; j<numCols; j++) {
        std::copy(myview,myview+myLen,Aptr);
        Aptr += LDA;
        myview += myStride;
      }
      myview = Teuchos::null;
      mydata = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const
  {
    TEST_FOR_EXCEPTION(ArrayOfPtrs.size() != getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::get2dCopy(ArrayOfPtrs): Array of pointers must contain as many pointers as the MultiVector has rows.");
    Node &node = lclMV_.getNode();
    if (isConstantStride()) {
      const int myStride = lclMV_.getStride(),
                 numCols = getNumVectors(),
                 myLen   = getMyLength();
      if (myLen > 0) {
        Teuchos::ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,lclMV_.getValues());
        for (Teuchos_Ordinal j=0; j<numCols; ++j) {
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(ArrayOfPtrs[j].size() != getMyLength(), std::runtime_error,
              "Tpetra::MultiVector::get2dCopy(ArrayOfPtrs): The ArrayView provided in ArrayOfPtrs[" << j << "] was not large enough to contain the local entries.");
#endif
          std::copy(myview,myview+myLen,ArrayOfPtrs[j].begin());
          myview += myStride;
        }
        myview = Teuchos::null;
      }
    }
    else {
      const int myStride = lclMV_.getStride(),
                 numCols = lclMV_.getNumCols(),
                  myCols = getNumVectors(),
                  myLen  = lclMV_.getNumRows();
      if (myLen > 0) {
        Teuchos::ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,lclMV_.getValues());
        for (Teuchos_Ordinal j=0; j<myCols; ++j) {
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(ArrayOfPtrs[j].size() != getMyLength(), std::runtime_error,
              "Tpetra::MultiVector::get2dCopy(ArrayOfPtrs): The ArrayView provided in ArrayOfPtrs[" << j << "] was not large enough to contain the local entries.");
#endif
          Teuchos::ArrayRCP<const Scalar> viewj = myview + whichVectors_[j]*myStride;
          std::copy(viewj,viewj+myLen,ArrayOfPtrs[j].begin());
        }
        myview = Teuchos::null;
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar> MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dView() const
  {
    TEST_FOR_EXCEPTION(!isConstantStride(), std::runtime_error,
      "Tpetra::MultiVector::get1dView() requires that this MultiVector have constant stride.");
    Node & node = lclMV_.getNode();
    return node.template viewBuffer<Scalar>( getStride()*(getNumVectors()-1)+getMyLength(), lclMV_.getValues() );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar> MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dViewNonConst()
  {
    TEST_FOR_EXCEPTION(!isConstantStride(), std::runtime_error,
      "Tpetra::MultiVector::get1dViewNonConst(): requires that this MultiVector have constant stride.");
    Node & node = lclMV_.getNode();
    return node.template viewBufferNonConst<Scalar>( false, getStride()*(getNumVectors()-1)+getMyLength(), lclMV_.getValuesNonConst() );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get2dViewNonConst()
  {
    Node & node = lclMV_.getNode();
    Teuchos::ArrayRCP< Teuchos::ArrayRCP<Scalar> > views = Teuchos::arcp<Teuchos::ArrayRCP<Scalar> >(getNumVectors());
    if (isConstantStride()) {
      const int myStride = lclMV_.getStride(),
                 numCols = getNumVectors(),
                 myLen   = getMyLength();
      if (myLen > 0) {
        Teuchos::ArrayRCP<Scalar> myview = node.template viewBufferNonConst<Scalar>(false,myStride*(numCols-1)+myLen,lclMV_.getValuesNonConst());
        for (Teuchos_Ordinal j=0; j<numCols; ++j) {
          views[j] = myview.persistingView(0,myLen);
          myview += myStride;
        }
      }
    }
    else {
      const int myStride = lclMV_.getStride(),
                 numCols = lclMV_.getNumCols(),
                  myCols = getNumVectors(),
                  myLen  = lclMV_.getNumRows();
      if (myLen > 0) {
        Teuchos::ArrayRCP<Scalar> myview = node.template viewBufferNonConst<Scalar>(false,myStride*(numCols-1)+myLen,lclMV_.getValuesNonConst());
        for (Teuchos_Ordinal j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > 
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get2dView() const
  {
    Node & node = lclMV_.getNode();
    Teuchos::ArrayRCP< Teuchos::ArrayRCP<const Scalar> > views = Teuchos::arcp<Teuchos::ArrayRCP<const Scalar> >(getNumVectors());
    if (isConstantStride()) {
      const int myStride = lclMV_.getStride(),
                 numCols = getNumVectors(),
                 myLen   = getMyLength();
      if (myLen > 0) {
        Teuchos::ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,lclMV_.getValues());
        for (Teuchos_Ordinal j=0; j<numCols; ++j) {
          views[j] = myview.persistingView(0,myLen);
          myview += myStride;
        }
      }
    }
    else {
      const int myStride = lclMV_.getStride(),
                 numCols = lclMV_.getNumCols(),
                  myCols = getNumVectors(),
                  myLen  = lclMV_.getNumRows();
      if (myLen > 0) {
        Teuchos::ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,lclMV_.getValues());
        for (Teuchos_Ordinal j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta) 
  {
    using Teuchos::NO_TRANS;      // enums
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::null;
    using Teuchos::ScalarTraits;  // traits
    using Teuchos::as;
    using Teuchos::RCP;           // data structures
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::rcp;           // initializers for data structures

    // This routine performs a variety of matrix-matrix multiply operations, interpreting
    // the MultiVector (this-aka C , A and B) as 2D matrices.  Variations are due to
    // the fact that A, B and C can be local replicated or global distributed
    // MultiVectors and that we may or may not operate with the transpose of 
    // A and B.  Possible cases are:
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
    Teuchos_Ordinal A_nrows = (transA==CONJ_TRANS) ? A.getNumVectors() : A.getMyLength();
    Teuchos_Ordinal A_ncols = (transA==CONJ_TRANS) ? A.getMyLength() : A.getNumVectors();
    Teuchos_Ordinal B_nrows = (transB==CONJ_TRANS) ? B.getNumVectors() : B.getMyLength();
    Teuchos_Ordinal B_ncols = (transB==CONJ_TRANS) ? B.getMyLength() : B.getNumVectors();

    Scalar beta_local = beta; // local copy of beta; might be reassigned below

    TEST_FOR_EXCEPTION( getMyLength() != A_nrows || getNumVectors() != B_ncols || A_ncols != B_nrows, std::runtime_error,
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

    if (beta != ScalarTraits<Scalar>::zero() && Case2)
    {
      // if Case2, then C is local and contributions must be summed across all nodes
      // however, if beta != 0, then accumulate beta*C into the sum
      // when summing across all nodes, we only want to accumulate this once, so 
      // set beta == 0 on all nodes except node 0
      int MyPID = this->getMap().getComm()->getRank();
      if (MyPID!=0) beta_local = ScalarTraits<Scalar>::zero();
    }

    // Check if A, B, C have constant stride, if not then make temp copy (strided)
    RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Atmp, Btmp; 
    RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >       Ctmp;
    if (isConstantStride() == false) Ctmp = rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(*this));
    else Ctmp = rcp(this,false);

    if (A.isConstantStride() == false) Atmp = rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A));
    else Atmp = rcp(&A,false);

    if (B.isConstantStride() == false) Btmp = rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(B));
    else Btmp = rcp(&B,false);

#ifdef HAVE_TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(!Ctmp->isConstantStride() || !Btmp->isConstantStride() || !Atmp->isConstantStride(), std::logic_error,
        errPrefix << "failed making temporary strided copies of input multivectors.");
#endif

    KMV &C_mv = Ctmp->lclMV_;
    {
      // get local multivectors
      const KMV &A_mv = Atmp->lclMV_;
      const KMV &B_mv = Btmp->lclMV_;
      // do the multiply (GEMM)
      DMVA::GEMM(C_mv,transA,transB,alpha,A_mv,B_mv,beta_local);
    }

    // Dispose of (possibly) extra copies of A, B
    Atmp = null;
    Btmp = null;

    Node &node = lclMV_.getNode();
    // If *this was not strided, copy the data from the strided version and then delete it
    if (isConstantStride() == false) {
      // *this is not strided, we must put data from Ctmp into *this
      TEST_FOR_EXCEPT(&C_mv != &lclMV_);
      const Kokkos::size_type numVecs = lclMV_.getNumCols();
      for (Kokkos::size_type j=0; j < numVecs; ++j) {
        node.template copyBuffers<Scalar>(getMyLength(),C_mv.getValues(j),lclMV_.getValuesNonConst(j));
      }
    }

    // If Case 2 then sum up *this and distribute it to all processors.
    if (Case2) {
      this->reduce();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reduce() {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    // this should be called only for "local" MultiVectors (!isDistributed())
    TEST_FOR_EXCEPTION(this->isDistributed() == true, std::runtime_error,
        "Tpetra::MultiVector::reduce() should only be called for non-distributed MultiVectors.");
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap().getComm();
    if (comm->getSize() == 1) return;
    Node &node = lclMV_.getNode();
    // sum the data across all multivectors
    // need to have separate packed buffers for send and receive
    // if we're packed, we'll just set the receive buffer as our data, the send as a copy
    // if we're not packed, we'll use allocated buffers for both. 
    ArrayView<Scalar> target;
    const int myStride = lclMV_.getStride(),
              numCols  = lclMV_.getNumCols(),
              myLen    = lclMV_.getNumRows();
    Teuchos::Array<Scalar> sourceBuffer(numCols*myLen), tmparr(0);
    bool packed = isConstantStride() && (myStride == myLen);
    ArrayRCP<Scalar> bufView = node.template viewBufferNonConst<Scalar>(
                                          false,myStride*(numCols-1)+myLen,
                                          lclMV_.getValuesNonConst() );
    if (packed) {
      // copy data from view to sourceBuffer, reduce into view below
      target = bufView(0,myLen*numCols);
      std::copy(target.begin(),target.end(),sourceBuffer.begin());
    }
    else {
      // copy data into sourceBuffer, reduce from sourceBuffer into tmparr below, copy back to view after that
      tmparr.resize(myLen*numCols);
      Scalar *sptr = sourceBuffer.getRawPtr();
      ArrayRCP<const Scalar> vptr = bufView;
      for (Teuchos_Ordinal j=0; j<numCols; ++j) 
      {
        std::copy(vptr,vptr+myLen,sptr);
        sptr += myLen;
        vptr += myStride;
      }
      target = tmparr();
    }
    // reduce 
    Teuchos::reduceAll<int,Scalar>(*comm,Teuchos::REDUCE_SUM,numCols*myLen,sourceBuffer.getRawPtr(),target.getRawPtr());
    if (!packed) {
      // copy tmparr back into view
      const Scalar *sptr = tmparr.getRawPtr();
      ArrayRCP<Scalar> vptr = bufView;
      for (Teuchos_Ordinal j=0; j<numCols; ++j) 
      {
        std::copy(sptr,sptr+myLen,vptr);
        sptr += myLen;
        vptr += myStride;
      }
    }
    bufView = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceMyValue(LocalOrdinal MyRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(MyRow < this->getMap().getMinLocalIndex() || MyRow > this->getMap().getMaxLocalIndex(), std::runtime_error,
        "Tpetra::MultiVector::replaceMyValue(): row index is invalid.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::replaceMyValue(): vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoMyValue(LocalOrdinal MyRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(MyRow < this->getMap().getMinLocalIndex() || MyRow > this->getMap().getMaxLocalIndex(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoMyValue(): row index is invalid.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoMyValue(): vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal GlobalRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap().getLocalIndex(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        "Tpetra::MultiVector::replaceGlobalValue(): row index is not present on this processor.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::replaceGlobalValue(): vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal GlobalRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap().getLocalIndex(GlobalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoGlobalValue(): row index is not present on this processor.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= getNumVectors(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoGlobalValue(): vector index is invalid.");
#endif
    getDataNonConst(VectorIndex)[MyRow] += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue(const Teuchos::ArrayView<Scalar> &means) const
  {
    TEST_FOR_EXCEPT(!isConstantStride());
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    const Teuchos_Ordinal numVecs = this->getNumVectors();
    TEST_FOR_EXCEPTION(means.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::meanValue(): means.size() must be as large as the number of vectors in *this.");
    // compute local components of the means
    // sum these across all nodes
    DMVA::Sum(lclMV_,means);
    if (this->isDistributed()) {
      Teuchos::Array<Scalar> lmeans(means);
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lmeans.getRawPtr(),means.getRawPtr());
    }
    const Scalar OneOverN = Teuchos::ScalarTraits<Scalar>::one() / Teuchos::as<Scalar>(getGlobalLength());
    for (typename Teuchos::ArrayView<Scalar>::iterator i = means.begin(); i != means.begin()+numVecs; ++i) {
      (*i) = (*i)*OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class T>
  Teuchos::ArrayRCP<T> MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getSubArrayRCP(Teuchos::ArrayRCP<T> arr, Teuchos_Ordinal j) const {
    Kokkos::size_type stride = lclMV_.getStride(),
                       myLen = getMyLength();
    Teuchos::ArrayRCP<T> ret;
    if (isConstantStride()) {
      ret = arr.persistingView(j*stride,myLen);
    }
    else {
      ret = arr.persistingView(whichVectors_[j]*stride,myLen);
    }
    return ret;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
  {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<getGlobalLength()
        << ",getNumVectors="<<getNumVectors()
        << ",isConstantStride()="<<isConstantStride()
        << "}";
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap().getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    int width = 1;
    for (int dec=10; dec<getGlobalLength(); dec *= 10) {
      ++width;
    }
    Teuchos::OSTab tab(out);
    if (vl != VERB_NONE) {
      // VERB_LOW and higher prints description()
      if (myImageID == 0) out << this->description() << std::endl; 
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // VERB_MEDIUM and higher prints getMyLength()
            out << "node " << setw(width) << myImageID << ": local length=" << getMyLength();
            if (vl != VERB_MEDIUM) {
              // VERB_HIGH and higher prints isConstantStride() and getStride()
              if (isConstantStride()) out << ", constant stride=" << getStride() << endl;
              else out << ", non-constant stride" << endl;
              if (vl == VERB_EXTREME && getMyLength() > 0) {
                Node &node = lclMV_.getNode();
                if (isConstantStride()) {
                  Teuchos::ArrayRCP<const Scalar> myview = node.template viewBuffer<Scalar>(
                                          getMyLength()+getStride()*(getNumVectors()-1), 
                                          lclMV_.getValues() );
                  // VERB_EXTREME prints values
                  for (Teuchos_Ordinal i=0; i<getMyLength(); ++i) {
                    out << setw(width) << this->getMap().getGlobalIndex(i) << ": ";
                    for (Teuchos_Ordinal j=0; j<getNumVectors(); ++j) {
                      out << myview[j*getStride()] << "  ";
                    }
                    ++myview;
                    out << endl;
                  }
                  myview = Teuchos::null;
                }
                else {
                  Kokkos::size_type stride = lclMV_.getStride(),
                                    rows   = lclMV_.getNumRows(),
                                    cols   = lclMV_.getNumCols();
                  Teuchos::ArrayRCP<const Scalar> myview = 
                    node.template viewBuffer<Scalar>( rows + stride * (cols - 1), lclMV_.getValues() );
                  // VERB_EXTREME prints values
                  for (Teuchos_Ordinal i=0; i<getMyLength(); ++i) {
                    out << setw(width) << this->getMap().getGlobalIndex(i) << ": ";
                    for (Teuchos_Ordinal j=0; j<getNumVectors(); ++j) {
                      out << myview[whichVectors_[j]*stride + i] << "  ";
                    }
                    out << endl;
                  }
                  myview = Teuchos::null;
                }
              }
            }
            else {
              out << endl;
            }
          }
        }
      }
    }
  }

} // namespace Tpetra

#endif // TPETRA_MULTIVECTOR_HPP
