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
#include "Tpetra_MultiVectorData.hpp"
#include "Tpetra_Vector.hpp"

// TODO: add real constructors for MultiVectorData so we don't have to setup its data and call setupPointers()

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos_Ordinal NumVectors, bool zeroOut) 
    : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
  {
    TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
        "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const LocalOrdinal myLen = myLength();
    lclMV_ = Teuchos::rcp( new KMV(node) );
    if (myLen > 0) {
      typedef Node::buffer<Scalar>::buffer_t data = node.template allocBuffer<Scalar>(myLen*NumVectors);
      lclMV_->initializeValues(myLen,NumVectors,data,myLen);
      if (zeroOut) {
        DMVA::Init(*lclMV_,0.0);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) 
    : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(source)
  {
    // copy data from the source MultiVector into this multivector
    Node &node = source.lclMV_->getNode();
    const LocalOrdinal myLen = myLength();
    const Teuchos_Ordinal numVecs = source.numVectors();
    {
      // hijack standard constructor; allocates contiguous strided data
      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> victim(node,source.getMap(),numVecs,false);
      lclMV_ = victim.lclMV_;
    }
    if (myLen > 0) {
      // copy data
      Teuchos::RCP<KMV> sourceLclMV = source.lclMV_;
      typename Node::buffer<Scalar>::buffer_t dest = lclMV_->getValues(0);
      int offset = 0;
      for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
        typename Node::buffer<Scalar>::buffer_t src = sourceLclMv->getValues(j);
        node.copyBuffers(myLen,src,0,dest,offset);
        offset += mylen;
      }
    }
  }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos::DataAccess CopyView, Scalar **ArrayOfPtrs, Teuchos_Ordinal NumVectors, bool OwnsMem)
  //REFACTOR//   : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPT(true);
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos::DataAccess CopyView, Scalar *A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors, bool OwnsMem)
  //REFACTOR//   : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
  //REFACTOR// {
  //REFACTOR//   if (CopyView == Teuchos::Copy) {
  //REFACTOR//     MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> victim(map, Teuchos::arrayView(A,LDA*(NumVectors-1)+myLength()),NumVectors);
  //REFACTOR//     MVData_ = victim.MVData_;
  //REFACTOR//   }
  //REFACTOR//   else {
  //REFACTOR//     TEST_FOR_EXCEPT(true);
  //REFACTOR//   }
  //REFACTOR// }


//REFACTOR//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//REFACTOR//   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayRCP<Scalar> &A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors)
//REFACTOR//     : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
//REFACTOR//   {
//REFACTOR//     TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
//REFACTOR//     const Teuchos_Ordinal myLen = myLength();
//REFACTOR//     TEST_FOR_EXCEPTION(LDA < myLen, std::runtime_error,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(): LDA must be large enough to accomodate the local entries.");
//REFACTOR// #ifdef HAVE_TPETRA_DEBUG
//REFACTOR//     TEST_FOR_EXCEPTION(A.size() < LDA*(NumVectors-1)+myLen, std::runtime_error,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(A,LDA): A.size() is not large enough for the specified NumVectors and LDA.");
//REFACTOR// #endif
//REFACTOR//     MVData_ = Teuchos::rcp( new MultiVectorData<Scalar>() );
//REFACTOR//     MVData_->constantStride_ = true;
//REFACTOR//     MVData_->contigValues_ = A;
//REFACTOR//     MVData_->stride_ = LDA;
//REFACTOR//     MVData_->setupPointers(myLen,NumVectors);
//REFACTOR//   }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayRCP<Scalar> > &ArrayOfPtrs, Teuchos_Ordinal NumVectors)
  //REFACTOR//   : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPT(true);
  //REFACTOR// }


//REFACTOR//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//REFACTOR//   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors)
//REFACTOR//     : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
//REFACTOR//   {
//REFACTOR//     TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
//REFACTOR//     const Teuchos_Ordinal myLen = myLength();
//REFACTOR//     using Teuchos::ArrayRCP;
//REFACTOR//     TEST_FOR_EXCEPTION(LDA < myLen, std::runtime_error,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(): LDA must be large enough to accomodate the local entries.");
//REFACTOR// #ifdef HAVE_TPETRA_DEBUG
//REFACTOR//     TEST_FOR_EXCEPTION(A.size() < LDA*(NumVectors-1)+myLen, std::runtime_error,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(A,LDA): A.size() is not large enough for the specified NumVectors and LDA.");
//REFACTOR// #endif
//REFACTOR//     {
//REFACTOR//       MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> victim(map,NumVectors,false);
//REFACTOR//       MVData_ = victim.MVData_;
//REFACTOR//     }
//REFACTOR//     if (myLen > 0) {
//REFACTOR//       const_pointer srcit = A.begin();
//REFACTOR//       for (Teuchos_Ordinal j = 0; j < NumVectors; ++j) {
//REFACTOR//         std::copy(srcit,srcit+myLen,(*this)[j]);
//REFACTOR//         srcit += LDA;
//REFACTOR//       }
//REFACTOR//     }
//REFACTOR//   }


//REFACTOR//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//REFACTOR//   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, 
//REFACTOR//                                                               const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, 
//REFACTOR//                                                               Teuchos_Ordinal NumVectors)
//REFACTOR//     : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm())
//REFACTOR//   {
//REFACTOR//     TEST_FOR_EXCEPTION(NumVectors < 1 || NumVectors != ArrayOfPtrs.size(), std::runtime_error,
//REFACTOR//         "Tpetra::MultiVector::MultiVector(map,ArrayOfPtrs): ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
//REFACTOR//     const Teuchos_Ordinal myLen = myLength();
//REFACTOR//     {
//REFACTOR//       MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> victim(map,NumVectors,false);
//REFACTOR//       MVData_ = victim.MVData_;
//REFACTOR//     }
//REFACTOR//     if (myLen > 0) {
//REFACTOR//       for (Teuchos_Ordinal j = 0; j < NumVectors; ++j) {
//REFACTOR// #ifdef HAVE_TPETRA_DEBUG
//REFACTOR//         TEST_FOR_EXCEPTION(ArrayOfPtrs[j].size() != myLength(), std::runtime_error,
//REFACTOR//           "Tpetra::MultiVector::MultiVector(map,ArrayOfPtrs): ArrayOfPtrs[" << j << "].size() (== " << ArrayOfPtrs[j].size() << 
//REFACTOR//           ") is not equal to myLength() (== " << myLength());
//REFACTOR// #endif
//REFACTOR//         typename Teuchos::ArrayView<const Scalar>::iterator src = ArrayOfPtrs[j].begin();
//REFACTOR//         std::copy(src,src+myLen,(*this)[j]);
//REFACTOR//       }
//REFACTOR//     }
//REFACTOR//   }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MultiVector(const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::RCP<Kokkos::MultiVector<Scalar,LocalOrdinal,Node> > &lclmv) 
    : DistObject<Scalar,LocalOrdinal,GlobalOrdinal>(map, map.getComm()), lclMV_(lclmv)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~MultiVector()
  {
    // free allocated memory 
    if (myLength() > 0) {
      typename Node::buffer<Scalar>::buffer_t buf = lclMV_.getValues(0);
      node.template freeBuffer<Scalar>(buf);
    }
    // delete Kokkos::MultiVector computational object
    lclMV_ = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::constantStride() const
  {
    return true;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::myLength() const
  {
    return this->getMap().getNumMyEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::globalLength() const
  {
    return this->getMap().getNumGlobalEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos_Ordinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::stride() const
  {
    return lclMV_->getStride();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj, Teuchos_Ordinal &packetSize) 
  {
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(sourceObj);
    // objects maps have already been checked. simply check the number of vectors.
    packetSize = this->numVectors();
    return (A.numVectors() == packetSize);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::copyAndPermute(
      const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & sourceObj,
      Teuchos_Ordinal numSameIDs,
      const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
      const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::ArrayView;
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceMV = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &>(sourceObj);
    typename ArrayView<const LocalOrdinal>::iterator pTo, pFrom;
    // any other error will be caught by Teuchos
    TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::copyAndPermute(): permuteToLIDs and permuteFromLIDs must have the same size.");
    // one vector at a time
    Node &node = lclMV_->getNode();
    const int numCols   = numVectors(),
              srcStride = sourceMV.lclMV_->getStride(),
              dstStride =          lclMV_->getStride();
    Scalar *dstptr;
    const Scalar *srcptr;
    // Get a host view of the local Kokkos::MultiVector data:
    //   Kokkos::MultiVector::getValues(j) returns a parallel buffer pointing to 
    //   the j-th column; viewBuffer/Const() will return a host pointer
    // TODO: determine whether this viewBuffer is write-only or not; for now, safe option is not
    dstptr =  node.viewBuffer(false,dstStride*numCols,lclMV_->getValues(0),0);
    srcptr =  node.viewBufferConst( srcStride*numCols,soruceMV.lclMV_->getValues(0),0);
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
    node.template releaseView<Scalar>(dstptr - numCols*dstStride);
    node.template releaseView<const Scalar>(srcptr - numCols*srcStride);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::packAndPrepare(
      const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & sourceObj,
      const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
      const Teuchos::ArrayView<Scalar> &exports,
      Distributor &distor)
  {
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceMV = dynamic_cast<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &>(sourceObj);
    using Teuchos::ArrayView;
    (void)distor;    // we don't use these, but we don't want unused parameter warnings
    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ; 
                   ....
                   all of the data from row exportLIDs.back() }
      this doesn't have the best locality, but is necessary because the data for a Packet
      (all data associated with an LID) is required to be contiguous */
    TEST_FOR_EXCEPTION(exports.size() != numVectors()*exportLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::packAndPrepare(): sizing of exports buffer should be appropriate for the amount of data to be exported.");
    const int myStride = lclMV_->getStride(),
              numCols   = numVectors();
    typename ArrayView<const LocalOrdinal>::iterator idptr;
    typename ArrayView<Scalar>::iterator expptr;
    expptr = exports.begin();
    Node &node = lclMV_->getNode();
    const Scalar * myptr;
    myptr = node.viewBufferConst(myStride*numCols,lclMV_->getValues(0),0);
    for (idptr = exportLIDs.begin(); idptr != exportLIDs.end(); ++idptr) {
      for (Teuchos_Ordinal j = 0; j < numCols; ++j) {
        *expptr++ = myptr[j*myStride + (*idptr)];
      }
    }
    node.template releaseView<const Scalar>(myptr);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::unpackAndCombine(
      const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
      const Teuchos::ArrayView<const Scalar> &imports,
      Distributor &distor,
      CombineMode CM)
  {
    (void)distor; // we don't use this, but we don't want unused parameter warnings
    using Teuchos::ArrayView;
    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ; 
                   ....
                   all of the data from row exportLIDs.back() }
      this doesn't have the best locality, but is necessary because the data for a Packet
      (all data associated with an LID) is required to be contiguous */
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(imports.size() != numVectors()*importLIDs.size(), std::runtime_error,
        "Tpetra::MultiVector::unpackAndCombine(): sizing of imports buffer should be appropriate for the amount of data to be exported.");
#endif
    const int myStride = lclMV_->getStride(),
              numCols   = numVectors();
    Node &node = lclMV_->getNode();
    const Scalar * myptr;
    // TODO: determine whether this viewBuffer is write-only or not; for now, safe option is not
    myptr = node.viewBuffer(false,myStride*numCols,lclMV_->getValues(0),0);
    typename ArrayView<const       Scalar>::iterator impptr;
    typename ArrayView<const LocalOrdinal>::iterator  idptr;
    impptr = imports.begin();
    if (CM == INSERT || CM == REPLACE) 
    {
      for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
        for (Teuchos_Ordinal j = 0; j < numVectors(); ++j) {
          (myptr+j*myStride)[*idptr] = *impptr++;
        }
      }
    }
    else if (CM == ADD) {
      for (idptr = importLIDs.begin(); idptr != importLIDs.end(); ++idptr) {
        for (Teuchos_Ordinal j = 0; j < numVectors(); ++j) {
          (myptr+j*myStride)[*idptr] += *impptr++;
        }
      }
    }
    else {
      TEST_FOR_EXCEPTION(CM != ADD && CM != REPLACE && CM != INSERT, std::invalid_argument,
          "Tpetra::MultiVector::unpackAndCombine(): Invalid CombineMode: " << CM);
    }
    node.template releaseView<Scalar>(myptr);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  inline Teuchos_Ordinal MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::numVectors() const 
  {
    return lclMV_->getNumCols();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot(
      const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, 
      const Teuchos::ArrayView<Scalar> &dots) const 
  {
    const Teuchos_Ordinal numVecs = numVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != A.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors must have the same number of vectors.");
    TEST_FOR_EXCEPTION(dots.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::dots(A,dots): dots.size() must be as large as the number of vectors in *this and A.");
    if (!this->isDistributed()) {
      MVDA::Dot(*lclMV_,*A.lclMV_,dots);
    }
    else {
      Teuchos::Array<Scalar> ldots(numVecs);
      MVDA::Dot(*lclMV_,*A.lclMV_,ldots());
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,ldots.getRawPtr(),dots.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Teuchos_Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm2(norms): norms.size() must be as large as the number of vectors in *this.");
    if (!this->isDistributed()) {
      MVDA::Norm2Squared(*lclMV_,norms);
    }
    else {
      Teuchos::Array<Mag> lnorms(numVecs);
      MVDA::Norm2Squared(*lclMV_,lnorms());
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      (*n) = ScalarTraits<Mag>::squareroot(*n);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted(
      const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights,
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::ScalarTraits;
    using Teuchos::ArrayView;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;
    const Mag OneOverN = ScalarTraits<Mag>::one() / Teuchos::as<Mag>(globalLength());
    bool OneW = false;
    const Teuchos_Ordinal numVecs = this->numVectors();
    if (weights.numVectors() == 1) {
      OneW = true;
    }
    else {
      TEST_FOR_EXCEPTION(weights.numVectors() != numVecs, std::runtime_error,
          "Tpetra::MultiVector::normWeighted(): MultiVector of weights must contain either one vector or the same number of vectors as this.");
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(weights.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "weights.getMap(): " << std::endl << weights.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != weights.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    // 
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::normWeighted(): norms.size() must be as large as the number of vectors in *this.");
    if (!this->isDistributed()) {
      MVDA::WeightedNorm(*weights.lclMV_,*lclMV_,norms);
    }
    else {
      Teuchos::Array<Mag> lnorms(numVecs);
      MVDA::WeightedNorm(*weights.lclMV_,*lclMV_,lnorms());
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n * OneOverN);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Teuchos_Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm1(norms): norms.size() must be as large as the number of vectors in *this.");
    if (!this->isDistributed()) {
      MVDA::Norm1(*lclMV_,*A.lclMV_,norms);
    }
    else {
      Teuchos::Array<Mag> lnorms(numVecs);
      MVDA::Norm1(*lclMV_,*A.lclMV_,lnorms());
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf(
      const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    const Teuchos_Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(norms.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::normInf(norms): norms.size() must be as large as the number of vectors in *this.");
    if (!this->isDistributed()) {
      MVDA::Norm1(*lclMV_,*A.lclMV_,norms);
    }
    else {
      Teuchos::Array<Mag> lnorms(numVecs);
      MVDA::Norm1(*lclMV_,*A.lclMV_,lnorms());
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_MAX,numVecs,lnorms.getRawPtr(),norms.getRawPtr());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::random() 
  {
    MVDA::Random(*lclMV_);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::putScalar(const Scalar &alpha) 
  {
    MVDA::Init(*lclMV_,alpha);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(const Scalar &alpha) 
  {
    // NOTE: can't substitute putScalar(0.0) for scale(0.0), because 
    //       the former will overwrite NaNs present in the MultiVector, while the 
    //       semantics of this call require multiplying them by 0, which IEEE requires to be NaN
    const Teuchos_Ordinal numVecs = this->numVectors();
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // do nothing
    }
    else {
      MVDA::Scale(*lclMV_,alpha);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(Teuchos::ArrayView<const Scalar> alphas)
  {
    const Teuchos_Ordinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(alphas.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as the number of vectors in *this.");
    KMV vec(lclMV_->getNode());
    for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
      if (alphas[j] == Teuchos::ScalarTraits<Scalar>::one()) {
        // do nothing: NaN * 1.0 == NaN, Number*1.0 == Number
      }
      else {
        vec.initializeValues(myLength(),1,lclMV_->getValues(j),myLength());
        MVDA::Scale(vec,alphas[j]);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) 
  {
    const Teuchos_Ordinal numVecs = this->numVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != A.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have the same number of vectors.");
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // set me = A
      MVDA::Assign(*lclMV_,(const KMV&)(*A.lclMV_));
    }
    else {
      // set me == alpha*A
      MVDA::Scale(*lclMV_,alpha,(const KMV&)(*A.lclMV_));
    }
  }

  //REFACTOR// FINISH HERE

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) 
  {
    typedef Teuchos::ScalarTraits<Scalar>                     ST;
    typedef Teuchos::ScalarTraits<typename ST::magnitudeType> MT;
    const Teuchos_Ordinal numVecs = this->numVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != A.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::reciprocal(): MultiVectors must have the same number of vectors.");
    int ierr = 0;
    for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
      pointer myptr = (*this)[j],
              myend = (*this)[j]+myLength();
      const_pointer Aptr = A[j];
      for (; myptr != myend; ++myptr, ++Aptr) {
        if (ST::magnitude(*Aptr) <= ST::sfmin()) {
          if (*Aptr == ST::zero()) {
            ierr = 1;
            (*myptr) = ST::nan();
          }
          else {
            ierr = 2;
            if (ST::magnitude(*Aptr) < MT::zero()) { // negative
              (*myptr) = -ST::rmax();
            }
            else {
              (*myptr) =  ST::rmax();
            }
          }
        }
        else {
          (*myptr) = ST::one()/(*Aptr);
        }
      }
    }
    if (ierr != 0) {
      TEST_FOR_EXCEPTION(ierr == 1, std::runtime_error, 
          "Tpetra::MultiVector::reciprocal(): element of A was zero.");
      TEST_FOR_EXCEPTION(ierr == 2, std::runtime_error, 
          "Tpetra::MultiVector::reciprocal(): element of A was too small to invert.");
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) 
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Teuchos_Ordinal numVecs = this->numVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != A.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::scale(): MultiVectors must have the same number of vectors.");
    for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
      pointer myptr = (*this)[j],
              myend = (*this)[j]+myLength();
      const_pointer Aptr = A[j];
      for (; myptr != myend;) {
        (*myptr++) = ST::magnitude(*Aptr++);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta) 
  {
    // this = beta*this + alpha*A
    // do atomically, in case &this == &A, and to increase memory efficiency
    // can't short circuit on alpha==0.0 or beta==0.0, because 0.0*NaN != 0.0
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::ArrayView;
    const Teuchos_Ordinal numVecs = this->numVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != A.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    if (beta == ST::one()) {
      if (alpha == ST::one()) { // this = this + A
        for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
          pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = A[j];
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = (*curpos) + (*Apos); }
        }
      }
      else {                    // this = this + alpha*A
        for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
          pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = A[j];
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = (*curpos) + alpha*(*Apos); }
        }
      }
    }
    else {
      if (alpha == ST::one()) { // this = beta*this + A
        for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
          pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = A[j];
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = beta*(*curpos) + (*Apos); }
        }
      }
      else {                    // this = beta*this + alpha*A
        for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
          pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = A[j];
          for (; curpos != cend; ++curpos, ++Apos) { *curpos = beta*(*curpos) + alpha*(*Apos); }
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma)
  {
    // this = alpha*A + beta*B + gamma*this
    // do atomically, in case &this == &A or &this == &B, and to increase memory efficiency
    // can't short circuit on alpha==0.0 or beta==0.0 or gamma==0.0, because 0.0*NaN != 0.0
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Teuchos_Ordinal numVecs = this->numVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( !this->getMap().isCompatible(A.getMap()) || !this->getMap().isCompatible(B.getMap()),
        std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap() 
        << "A.getMap(): " << std::endl << A.getMap() << std::endl
        << "B.getMap(): " << std::endl << B.getMap() << std::endl);
#else
    TEST_FOR_EXCEPTION( myLength() != A.myLength() || myLength() != B.myLength(), std::runtime_error,
        "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
#endif
    TEST_FOR_EXCEPTION(A.numVectors() != numVecs || B.numVectors() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
    // if only one of alpha or beta is 1.0, make it alpha. then below, alpha != 1.0 implies beta != 1.0
    Teuchos::Ptr<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Aptr = Teuchos::ptr(&A), Bptr = Teuchos::ptr(&B);
    Teuchos::Ptr<const Scalar> lalpha = Teuchos::ptr(&alpha),
                               lbeta  = Teuchos::ptr(&beta);
    if (alpha!=ST::one() && beta==ST::one()) {
      // switch them
      Aptr = Teuchos::ptr(&B);
      Bptr = Teuchos::ptr(&A);
      lalpha = Teuchos::ptr(&beta);
      lbeta  = Teuchos::ptr(&alpha);
    }
    if (gamma == ST::one()) {
      if ((*lalpha) == ST::one()) {
        if ((*lbeta) == ST::one()) {    // this = this + A + B
          for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
            pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = (*Aptr)[j], Bpos = (*Bptr)[j];
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*Apos) + (*Bpos); }
          }
        }
        else {                          // this = this + A + beta*B
          for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
            pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = (*Aptr)[j], Bpos = (*Bptr)[j];
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*Apos) + (*lbeta)*(*Bpos); }
          }
        }
      }
      else {                            // this = this + alpha*A + beta*B
        for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
          pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = (*Aptr)[j], Bpos = (*Bptr)[j];
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = (*curpos) + (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
        }
      }
    }
    else {
      if ((*lalpha) == ST::one()) {
        if ((*lbeta) == ST::one()) {    // this = gamma*this + A + B
          for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
            pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = (*Aptr)[j], Bpos = (*Bptr)[j];
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*Apos) + (*Bpos); }
          }
        }
        else {                          // this = gamma*this + A + beta*B
          for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
            pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = (*Aptr)[j], Bpos = (*Bptr)[j];
            for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*Apos) + (*lbeta)*(*Bpos); }
          }
        }
      }
      else {                            // this = gamma*this + alpha*A + beta*B
        for (Teuchos_Ordinal j = 0; j < numVecs; ++j) {
          pointer curpos = (*this)[j], cend = (*this)[j]+myLength(); const_pointer Apos = (*Aptr)[j], Bpos = (*Bptr)[j];
          for (; curpos != cend; ++curpos, ++Apos, ++Bpos) { *curpos = gamma*(*curpos) + (*lalpha)*(*Apos) + (*lbeta)*(*Bpos); }
        }
      }
    }
  }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::const_pointer 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator[](Teuchos_Ordinal j) const
  //REFACTOR// {
  //REFACTOR//   return MVData_->ptrs_[j];
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::pointer 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator[](Teuchos_Ordinal i)
  //REFACTOR// {
  //REFACTOR//   return MVData_->ptrs_[i];
  //REFACTOR// }


//REFACTOR//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//REFACTOR//   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator=(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) 
//REFACTOR//   {
//REFACTOR//     // Check for special case of this=Source, in which case we do nothing
//REFACTOR//     if (this != &source) {
//REFACTOR// #ifdef HAVE_TPETRA_DEBUG
//REFACTOR//       TEST_FOR_EXCEPTION( !this->getMap().isCompatible(source.getMap()), std::runtime_error,
//REFACTOR//           "Tpetra::MultiVector::dots(): MultiVectors do not have compatible Maps:" << std::endl
//REFACTOR//           << "this->getMap(): " << std::endl << this->getMap() 
//REFACTOR//           << "source.getMap(): " << std::endl << source.getMap() << std::endl);
//REFACTOR// #else
//REFACTOR//       TEST_FOR_EXCEPTION( myLength() != source.myLength(), std::runtime_error,
//REFACTOR//           "Tpetra::MultiVector::dots(): MultiVectors do not have the same local length.");
//REFACTOR// #endif
//REFACTOR//       TEST_FOR_EXCEPTION(source.numVectors() != numVectors(), std::runtime_error,
//REFACTOR//           "Tpetra::MultiVector::update(): MultiVectors must have the same number of vectors.");
//REFACTOR//       if (constantStride() && source.constantStride() && myLength()==stride() && source.myLength()==source.stride()) {
//REFACTOR//         // we're both packed, we can copy in one call
//REFACTOR//         std::copy( source.MVData_->contigValues_.begin(), source.MVData_->contigValues_.end(), MVData_->contigValues_.begin() );
//REFACTOR//       }
//REFACTOR//       else {
//REFACTOR//         for (Teuchos_Ordinal j=0; j<numVectors(); ++j) {
//REFACTOR//           std::copy( source[j], source[j]+myLength(), (*this)[j] );
//REFACTOR//         }
//REFACTOR//       }
//REFACTOR//     }
//REFACTOR//     return(*this);
//REFACTOR//   }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subCopy(const Teuchos::ArrayView<const Teuchos_Index> &cols) const
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(cols.size() < 1, std::runtime_error,
  //REFACTOR//       "Tpetra::MultiVector::subCopy(cols): cols must contain at least one column.");
  //REFACTOR//   Teuchos_Ordinal numCopyVecs = cols.size();
  //REFACTOR//   const bool zeroData = false;
  //REFACTOR//   Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv = Teuchos::rcp( new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),numCopyVecs,zeroData) );
  //REFACTOR//   // copy data from *this into mv
  //REFACTOR//   for (Teuchos_Ordinal j=0; j<numCopyVecs; ++j)
  //REFACTOR//   {
  //REFACTOR//     std::copy( (*this)[cols[j]], (*this)[cols[j]]+myLength(), (*mv)[j] );
  //REFACTOR//   }
  //REFACTOR//   return mv;
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subView(const Teuchos::ArrayView<const Teuchos_Index> &cols) 
  //REFACTOR// {
  //REFACTOR//   const Teuchos_Ordinal numViewVecs = cols.size();
  //REFACTOR//   TEST_FOR_EXCEPTION(cols.size() < 1, std::runtime_error,
  //REFACTOR//       "Tpetra::MultiVector::subView(cols): cols must contain at least one column.");
  //REFACTOR//   Teuchos::RCP<MultiVectorData<Scalar> > mvdata = Teuchos::rcp( new MultiVectorData<Scalar>() );
  //REFACTOR//   mvdata->constantStride_ = false;
  //REFACTOR//   mvdata->stride_ = 0;
  //REFACTOR//   mvdata->contigValues_ = Teuchos::null;
  //REFACTOR//   mvdata->nonContigValues_.resize(numViewVecs);
  //REFACTOR//   if (myLength() > 0) {
  //REFACTOR//     if (constantStride()) {
  //REFACTOR//       for (Teuchos_Ordinal j = 0; j < numViewVecs; ++j) {
  //REFACTOR//         mvdata->nonContigValues_[j] = MVData_->contigValues_.persistingView( stride()*cols[j], myLength() );
  //REFACTOR//       }
  //REFACTOR//     }
  //REFACTOR//     else {
  //REFACTOR//       for (Teuchos_Ordinal j = 0; j < numViewVecs; ++j) {
  //REFACTOR//         mvdata->nonContigValues_[j] = MVData_->nonContigValues_[cols[j]];
  //REFACTOR//       }
  //REFACTOR//     }
  //REFACTOR//   }
  //REFACTOR//   mvdata->setupPointers(myLength(),numViewVecs);
  //REFACTOR//   Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv = Teuchos::rcp( new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),mvdata) );
  //REFACTOR//   return mv;
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subViewConst(const Teuchos::ArrayView<const Teuchos_Index> &cols) const
  //REFACTOR// {
  //REFACTOR//   const Teuchos_Ordinal numVecs = cols.size();
  //REFACTOR//   TEST_FOR_EXCEPTION(cols.size() < 1, std::runtime_error,
  //REFACTOR//       "Tpetra::MultiVector::subViewConst(cols): cols must contain at least one column.");
  //REFACTOR//   Teuchos::RCP<MultiVectorData<Scalar> > mvdata = Teuchos::rcp( new MultiVectorData<Scalar>() );
  //REFACTOR//   mvdata->constantStride_ = false;
  //REFACTOR//   mvdata->stride_ = stride();
  //REFACTOR//   mvdata->contigValues_ = MVData_->contigValues_;
  //REFACTOR//   mvdata->ptrs_.resize(numVecs);
  //REFACTOR//   for (GlobalOrdinal j = 0; j < numVecs; ++j) {
  //REFACTOR//     mvdata->ptrs_[j] = MVData_->ptrs_[cols[j]];
  //REFACTOR//   }
  //REFACTOR//   Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv = Teuchos::rcp( new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),mvdata) );
  //REFACTOR//   return mv;
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subCopy(const Teuchos::Range1D &colRng) const
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
  //REFACTOR//       "Tpetra::MultiVector::subCopy(Range1D): range must include at least one vector.");
  //REFACTOR//   Teuchos_Ordinal numCopyVecs = colRng.size();
  //REFACTOR//   const bool zeroData = false;
  //REFACTOR//   Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mv = Teuchos::rcp( new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),numCopyVecs,zeroData) );
  //REFACTOR//   // copy data from *this into mv
  //REFACTOR//   for (Teuchos_Ordinal js=colRng.lbound(), jd=0; jd<numCopyVecs; ++jd, ++js)
  //REFACTOR//   {
  //REFACTOR//     std::copy( (*this)[js], (*this)[js]+myLength(), (*mv)[jd] );
  //REFACTOR//   }
  //REFACTOR//   return mv;
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subView(const Teuchos::Range1D &colRng)
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
  //REFACTOR//       "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
  //REFACTOR//   Teuchos_Ordinal numViewVecs = colRng.size();
  //REFACTOR//   using Teuchos::ArrayRCP;
  //REFACTOR//   // resulting MultiVector is constant stride only if *this is 
  //REFACTOR//   if (constantStride()) {
  //REFACTOR//     // view goes from first entry of first vector to last entry of last vector
  //REFACTOR//     ArrayRCP<Scalar> data = MVData_->contigValues_.persistingView(stride()*colRng.lbound(),
  //REFACTOR//                                                                   stride()*(numViewVecs-1) + myLength());
  //REFACTOR//     return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),data,stride(),numViewVecs));
  //REFACTOR//   }
  //REFACTOR//   else if (numViewVecs == 1) {
  //REFACTOR//     // not constant stride, so nonContigValues_ is filled
  //REFACTOR//     ArrayRCP<Scalar> data = MVData_->nonContigValues_[colRng.lbound()].persistingView(0,myLength());
  //REFACTOR//     return Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal>(this->getMap(),data));
  //REFACTOR//   }
  //REFACTOR//   else {
  //REFACTOR//     Teuchos::Array<ArrayRCP<Scalar> > rcps(numViewVecs);
  //REFACTOR//     for (Teuchos_Ordinal js=colRng.lbound(), jd=0; jd<numViewVecs; ++jd, ++js) {
  //REFACTOR//       rcps[jd] = MVData_->nonContigValues_[js];
  //REFACTOR//     }
  //REFACTOR//     return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),rcps(),numViewVecs));
  //REFACTOR//   }
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::subViewConst(const Teuchos::Range1D &colRng) const
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
  //REFACTOR//       "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
  //REFACTOR//   Teuchos_Ordinal numViewVecs = colRng.size();
  //REFACTOR//   using Teuchos::ArrayRCP;
  //REFACTOR//   // resulting MultiVector is constant stride only if *this is 
  //REFACTOR//   if (constantStride()) {
  //REFACTOR//     // view goes from first entry of first vector to last entry of last vector
  //REFACTOR//     ArrayRCP<Scalar> data = MVData_->contigValues_.persistingView(stride()*colRng.lbound(),
  //REFACTOR//                                                                   stride()*(numViewVecs-1) + myLength());
  //REFACTOR//     return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),data,stride(),numViewVecs));
  //REFACTOR//   }
  //REFACTOR//   else if (numViewVecs == 1) {
  //REFACTOR//     // not constant stride, so nonContigValues_ is filled
  //REFACTOR//     ArrayRCP<Scalar> data = MVData_->nonContigValues_[colRng.lbound()].persistingView(0,myLength());
  //REFACTOR//     return Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal>(this->getMap(),data));
  //REFACTOR//   }
  //REFACTOR//   else {
  //REFACTOR//     Teuchos::Array<ArrayRCP<Scalar> > rcps(numViewVecs);
  //REFACTOR//     for (Teuchos_Ordinal js=colRng.lbound(), jd=0; jd<numViewVecs; ++jd, ++js) {
  //REFACTOR//       rcps[jd] = MVData_->nonContigValues_[js];
  //REFACTOR//     }
  //REFACTOR//     return Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(this->getMap(),rcps(),numViewVecs));
  //REFACTOR//   }
  //REFACTOR// }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal> > MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator()(Teuchos_Ordinal j)
  {
    Teuchos::ArrayRCP<Scalar> data;
    if (constantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      data = MVData_->contigValues_.persistingView(stride()*j,myLength());
    }
    else {
      data = MVData_->nonContigValues_[j].persistingView(0,myLength());
    }
    return Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal>(this->getMap(),data));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal> > MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::operator() (Teuchos_Ordinal j) const
  {
    Teuchos::ArrayRCP<Scalar> data;
    if (constantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      data = MVData_->contigValues_.persistingView(stride()*j,myLength());
    }
    else {
      data = MVData_->nonContigValues_[j].persistingView(0,myLength());
    }
    return Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal>(this->getMap(),data));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractCopy1D(typename Teuchos::ArrayView<Scalar> A, Teuchos_Ordinal LDA) const
  {
    TEST_FOR_EXCEPTION(LDA < myLength(), std::runtime_error,
      "Tpetra::MultiVector::extractCopy1D(A,LDA): specified stride is not large enough for local vector length.");
    TEST_FOR_EXCEPTION(A.size() < LDA*(numVectors()-1)+myLength(), std::runtime_error,
      "Tpetra::MultiVector::extractCopy1D(A,LDA): specified stride is not large enough for local vector length.");
    pointer Aptr = A.begin();
    for (Teuchos_Ordinal j=0; j<numVectors(); j++) {
      std::copy((*this)[j], (*this)[j]+myLength(), Aptr);
      Aptr += LDA;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractCopy1D(Scalar *A, Teuchos_Ordinal LDA) const
  {
    for (Teuchos_Ordinal j=0; j<numVectors(); j++) {
      std::copy((*this)[j], (*this)[j]+myLength(), A);
      A += LDA;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractCopy2D(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const
  {
    TEST_FOR_EXCEPTION(ArrayOfPtrs.size() != numVectors(), std::runtime_error,
        "Tpetra::MultiVector::extractCopy2D(ArrayOfPtrs): Array of pointers must contain as many pointers as the MultiVector has rows.");
    for (Teuchos_Ordinal j=0; j<numVectors(); ++j) {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(ArrayOfPtrs[j].size() < myLength(), std::runtime_error,
        "Tpetra::MultiVector::extractCopy2D(ArrayOfPtrs): The ArrayView provided in ArrayOfPtrs[" << j << "] was not large enough to contain the local entries.");
#endif
      std::copy((*this)[j], (*this)[j]+myLength(), ArrayOfPtrs[j].begin());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractCopy2D(Scalar * const * ArrayOfPtrs) const
  {
    for (Teuchos_Ordinal j=0; j<numVectors(); ++j) {
      std::copy((*this)[j], (*this)[j]+myLength(), ArrayOfPtrs[j]);
    }
  }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractView1D(Scalar * &A, Teuchos_Ordinal &LDA)
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(!constantStride(), std::runtime_error,
  //REFACTOR//     "Tpetra::MultiVector::extractView1D(): cannot call for MultiVectors with non-constant stride.");
  //REFACTOR//   A = MVData_->contigValues_.getRawPtr();
  //REFACTOR//   LDA = stride();
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractView1D(Teuchos::ArrayView<Scalar> &A, Teuchos_Ordinal &LDA)
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(!constantStride(), std::runtime_error,
  //REFACTOR//     "Tpetra::MultiVector::extractView1D(): cannot call for MultiVectors with non-constant stride.");
  //REFACTOR//   A = MVData_->contigValues_();
  //REFACTOR//   LDA = stride();
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractConstView1D(const Scalar * &A, Teuchos_Ordinal &LDA) const
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(!constantStride(), std::runtime_error,
  //REFACTOR//     "Tpetra::MultiVector::extractConstView1D(): cannot call for MultiVectors with non-constant stride.");
  //REFACTOR//   A = MVData_->contigValues_.getRawPtr();
  //REFACTOR//   LDA = stride();
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractConstView1D(Teuchos::ArrayView<const Scalar> &A, Teuchos_Ordinal &LDA) const
  //REFACTOR// {
  //REFACTOR//   TEST_FOR_EXCEPTION(!constantStride(), std::runtime_error,
  //REFACTOR//     "Tpetra::MultiVector::extractConstView1D(): cannot call for MultiVectors with non-constant stride.");
  //REFACTOR//   A = MVData_->contigValues_();
  //REFACTOR//   LDA = stride();
  //REFACTOR// }


  //REFACTOR// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //REFACTOR// typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::double_pointer
  //REFACTOR// MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractView2D()
  //REFACTOR// {
  //REFACTOR//   return MVData_->ptrs_().begin();
  //REFACTOR// }


//REFACTOR//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//REFACTOR//   typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::const_double_pointer
//REFACTOR//   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::extractConstView2D() const
//REFACTOR//   {
//REFACTOR// #ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
//REFACTOR//       // MVData_->ptrs_                 Array<ArrayView<Scalar>::iterator>
//REFACTOR//       // MVData_->ptrs_()           ArrayView<ArrayView<Scalar>::iterator>
//REFACTOR//       // MVData_->ptrs_().begin()   ArrayView<ArrayView<Scalar>::iterator>::iterator
//REFACTOR//       // in debug mode, this is      ArrayRCP<ArrayRCP<Scalar> >
//REFACTOR//       return Teuchos::arcp_reinterpret_cast<Teuchos::ArrayRCP<const Scalar> >(MVData_->ptrs_().begin());
//REFACTOR// #else
//REFACTOR//       return MVData_->ptrs_().begin();
//REFACTOR// #endif
//REFACTOR//   }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta) 
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
    using Teuchos::as;
    using Teuchos::RCP;           // data structures
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::rcp;           // initializers for data structures

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
    Teuchos_Ordinal A_nrows = (transA==CONJ_TRANS) ? A.numVectors() : A.myLength();
    Teuchos_Ordinal A_ncols = (transA==CONJ_TRANS) ? A.myLength() : A.numVectors();
    Teuchos_Ordinal B_nrows = (transB==CONJ_TRANS) ? B.numVectors() : B.myLength();
    Teuchos_Ordinal B_ncols = (transB==CONJ_TRANS) ? B.myLength() : B.numVectors();

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
    RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Ctmp;
    if (constantStride() == false) Ctmp = rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(*this));
    else Ctmp = rcp(this,false);

    if (A.constantStride() == false) Atmp = rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A));
    else Atmp = rcp(&A,false);

    if (B.constantStride() == false) Btmp = rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(B));
    else Btmp = rcp(&B,false);

#ifdef HAVE_TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(!Ctmp->constantStride() || !Btmp->constantStride() || !Atmp->constantStride(), std::logic_error,
        errPrefix << "failed making temporary strided copies of input multivectors.");
#endif

    Teuchos_Ordinal m = this->myLength();
    Teuchos_Ordinal n = this->numVectors();
    Teuchos_Ordinal k = A_ncols;
    Teuchos_Ordinal lda, ldb, ldc;
    ArrayView<const Scalar> Ap(Teuchos::null), Bp(Teuchos::null);
    ArrayView<Scalar> Cp(Teuchos::null);
    Atmp->extractConstView1D(Ap,lda);
    Btmp->extractConstView1D(Bp,ldb);
    Ctmp->extractView1D(Cp,ldc);

    Teuchos::BLAS<int,Scalar> blas;
    // do the arithmetic now
    blas.GEMM(transA,transB,m,n,k,alpha,Ap.getRawPtr(),lda,Bp.getRawPtr(),ldb,beta_local,Cp.getRawPtr(),ldc);

    // Dispose of (possibly) extra copies of A, B
    Atmp = null;
    Btmp = null;

    // If *this was not strided, copy the data from the strided version and then delete it
    if (constantStride() == false) {
      // *this is not strided, so we must use extract data into our 2D structure
      Array<ArrayView<Scalar> > aoa(MVData_->ptrs_.size(),null);
      for (Teuchos_Ordinal j=0; j<Teuchos::as<Teuchos_Ordinal>(aoa.size()); ++j) {
        aoa[j] = MVData_->nonContigValues_[j]();
      }
      Ctmp->extractCopy2D(aoa());
    }
    Ctmp = null;

    // If Case 2 then sum up *this and distribute it to all processors.
    if (Case2) {
      this->reduce();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reduce() 
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;
    // this should be called only for "local" MultiVectors (!isDistributed())
    TEST_FOR_EXCEPTION(this->isDistributed() == true, std::runtime_error,
        "Tpetra::MultiVector::reduce() should only be called for non-distributed MultiVectors.");
    // sum the data across all multivectors
    // need to have separate packed buffers for send and receive
    // if we're packed, we'll just set the receive buffer as our data, the send as a copy
    // if we're not packed, we'll use allocated buffers for both. 
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap().getComm();
    if (comm->getSize() == 1) return;
    ArrayRCP<Scalar> source = arcp<Scalar>(myLength()*numVectors()), target;
    bool packed = constantStride() && (stride() == myLength());
    if (packed) {
      // copy local info into source buffer
      // target buffer will be multivector storage
      std::copy(MVData_->contigValues_.begin(),MVData_->contigValues_.end(), source.begin());
      // local storage is packed. can use it for target buffer.
      target = MVData_->contigValues_;
    }
    else {
      // copy local data into source buffer
      ArrayRCP<Scalar> sptr = source;
      for (Teuchos_Ordinal j=0; j<numVectors(); ++j) 
      {
        // copy j-th local MV data into source buffer
        std::copy((*this)[j],(*this)[j]+myLength(),sptr);
        // advance ptr into source buffer
        sptr += myLength();
      }
      // must allocate packed storage for target buffer
      target = arcp<Scalar>(myLength()*numVectors());
    }
    // reduce 
    Teuchos::reduceAll<int,Scalar>(*comm,Teuchos::REDUCE_SUM,myLength()*numVectors(),source.getRawPtr(),target.getRawPtr());
    if (!packed) {
      // copy target buffer into multivector storage buffer
      ArrayRCP<Scalar> tptr = target;
      for (Teuchos_Ordinal j=0; j<numVectors(); ++j) 
      {
        std::copy(tptr,tptr+myLength(),(*this)[j]);
        tptr += myLength();
      }
    }
    // clear allocated buffers
    source = Teuchos::null;
    target = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceMap(const Map<LocalOrdinal,GlobalOrdinal> &map)
  {
    TEST_FOR_EXCEPT(true);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceMyValue(LocalOrdinal MyRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(MyRow < this->getMap().getMinLocalIndex() || MyRow > this->getMap().getMaxLocalIndex(), std::runtime_error,
        "Tpetra::MultiVector::replaceMyValue(): row index is invalid.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::replaceMyValue(): vector index is invalid.");
#endif
    MVData_->ptrs_[VectorIndex][MyRow] = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoMyValue(LocalOrdinal MyRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(MyRow < this->getMap().getMinLocalIndex() || MyRow > this->getMap().getMaxLocalIndex(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoMyValue(): row index is invalid.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoMyValue(): vector index is invalid.");
#endif
    MVData_->ptrs_[VectorIndex][MyRow] += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal GlobalRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap().getLocalIndex(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        "Tpetra::MultiVector::replaceGlobalValue(): row index is not present on this processor.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::replaceGlobalValue(): vector index is invalid.");
#endif
    (*this)[VectorIndex][MyRow] = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal GlobalRow, Teuchos_Ordinal VectorIndex, const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap().getLocalIndex(GlobalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoGlobalValue(): row index is not present on this processor.");
    TEST_FOR_EXCEPTION(VectorIndex < 0 || VectorIndex >= numVectors(), std::runtime_error,
        "Tpetra::MultiVector::sumIntoGlobalValue(): vector index is invalid.");
#endif
    (*this)[VectorIndex][MyRow] += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue(const Teuchos::ArrayView<Scalar> &means) const
  {
    using Teuchos::ArrayView;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    const int numImages = this->getMap().getComm()->getSize();
    const GlobalOrdinal numVecs = this->numVectors();
    TEST_FOR_EXCEPTION(means.size() != numVecs, std::runtime_error,
        "Tpetra::MultiVector::meanValue(): means.size() must be as large as the number of vectors in *this.");
    // compute local components of the means
    // sum these across all nodes
    Teuchos::Array<Scalar> lmeans(numVecs,SCT::zero());
    for (Teuchos_Ordinal j=0; j<numVecs; ++j) {
      const_pointer cpos = (*this)[j],
                    cend = (*this)[j]+myLength();
      for (; cpos != cend; ++cpos) {
        lmeans[j] += *cpos;
      }
    }
    if (this->isDistributed()) {
      // only combine if we are a distributed MV
      Teuchos::reduceAll(*this->getMap().getComm(),Teuchos::REDUCE_SUM,numVecs,lmeans.getRawPtr(),means.getRawPtr());
    }
    else {
      std::copy(lmeans.begin(),lmeans.end(),means.begin());
    }
    const Scalar OneOverN = Teuchos::as<Scalar>(globalLength());
    for (pointer n = means.begin(); n != means.begin()+numVecs; ++n) {
      (*n) = (*n)*OneOverN;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
  {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{length="<<globalLength()
        << ",numVectors="<<numVectors()
        << "}";
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
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
    for (int dec=10; dec<globalLength(); dec *= 10) {
      ++width;
    }
    Teuchos::OSTab tab(out);
    if (vl != VERB_NONE) {
      // VERB_LOW and higher prints description()
      if (myImageID == 0) out << this->description() << std::endl; 
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // VERB_MEDIUM and higher prints myLength()
            out << "node " << setw(width) << myImageID << ": local length=" << myLength();
            if (vl != VERB_MEDIUM) {
              // VERB_HIGH and higher prints constantStride() and stride()
              if (constantStride()) out << ", constant stride=" << stride() << endl;
              else out << ", non-constant stride" << endl;
              if (vl == VERB_EXTREME) {
                // VERB_EXTREME prints values
                for (Teuchos_Ordinal i=0; i<myLength(); ++i) {
                  out << setw(width) << this->getMap().getGlobalIndex(i) << ": ";
                  for (Teuchos_Ordinal j=0; j<numVectors(); ++j) {
                    out << MVData_->ptrs_[j][i] << "  ";
                  }
                  out << endl;
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
