// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_MULTIVECTOR_DEF_HPP
#define TPETRA_MULTIVECTOR_DEF_HPP

#include <Kokkos_DefaultArithmetic.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Vector.hpp>

#if TPETRA_USE_KOKKOS_DISTOBJECT
#  include <Tpetra_Details_MultiVectorDistObjectKernels.hpp>
#endif

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_MultiVector_decl.hpp"
#endif

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  vectorIndexOutOfRange(size_t VectorIndex) const {
    return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector () :
    DO (Teuchos::rcp (new Map<LocalOrdinal, GlobalOrdinal, Node> ())),
    lclMV_ (this->getMap ()->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               size_t NumVectors,
               bool zeroOut) : /* default is true */
    DO (map),
    lclMV_ (map->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    TEUCHOS_TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
      "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
    if (myLen > 0) {
      RCP<Node> node = map->getNode();
      // On host-type Kokkos Nodes, allocBuffer() just calls the
      // one-argument version of arcp to allocate memory.  This should
      // not fill the memory by default, otherwise we would lose the
      // first-touch allocation optimization.
      ArrayRCP<Scalar> data = node->template allocBuffer<Scalar>(myLen*NumVectors);
      MVT::initializeValues(lclMV_,myLen,NumVectors,data,myLen);
      if (zeroOut) {
        // MVT uses the Kokkos Node's parallel_for in this case, for
        // first-touch allocation (across rows).
        MVT::Init(lclMV_, Teuchos::ScalarTraits<Scalar>::zero());
      }
    }
    else {
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source) :
    DO (source),
    lclMV_ (MVT::getNode (source.lclMV_)),
    hasViewSemantics_ (source.hasViewSemantics_)
  {
    if (source.hasViewSemantics_) {
      // Shallow copy.  The const_cast is ugly but does no harm:
      // ArrayRCP is a pointer, so there's no reason why a method that
      // returns it can't be marked const.  (If you're not changing
      // the pointer itself, e.g., by reallocating, the pointer is
      // const.  It doesn't matter whether or not you can change the
      // values to which the pointer points.)
      MVT::initializeValues (lclMV_, source.lclMV_.getNumRows (),
                             source.lclMV_.getNumCols (),
                             const_cast<KMV&> (source.lclMV_).getValuesNonConst (),
                             source.lclMV_.getStride (),
                             source.lclMV_.getOrigNumRows (),
                             source.lclMV_.getOrigNumCols ());
      // Don't forget to copy whichVectors_ (for non-constant-stride view).
      if (source.whichVectors_.size () > 0) {
        whichVectors_ = source.whichVectors_; // Array::op= does a deep copy
      }
    }
    else {
      using Teuchos::ArrayRCP;
      using Teuchos::RCP;

      // copy data from the source MultiVector into this multivector
      RCP<Node> node = MVT::getNode(source.lclMV_);
      const LocalOrdinal myLen = source.getLocalLength();
      const size_t numVecs = source.getNumVectors();

      // On host-type Kokkos Nodes, allocBuffer() just calls the
      // one-argument version of arcp to allocate memory.  This should
      // not fill the memory by default, otherwise we would lose the
      // first-touch allocation optimization.
      ArrayRCP<Scalar> data = (myLen > 0) ?
        node->template allocBuffer<Scalar> (myLen * numVecs) :
        Teuchos::null;
      const size_t stride = (myLen > 0) ? myLen : size_t (0);

      // This just sets the dimensions, pointer, and stride of lclMV_.
      MVT::initializeValues (lclMV_, myLen, numVecs, data, stride);
      // This actually copies the data.  It uses the Node's
      // parallel_for to copy, which should ensure first-touch
      // allocation on systems that support it.
      if (source.isConstantStride ()) {
        MVT::Assign (lclMV_, source.lclMV_);
      }
      else {
        MVT::Assign (lclMV_, source.lclMV_, source.whichVectors_);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source,
               const Teuchos::DataAccess copyOrView) :
    DO (source),
    lclMV_ (MVT::getNode (source.lclMV_)),
    hasViewSemantics_ (copyOrView == Teuchos::View)
  {
    if (copyOrView == Teuchos::View) {
      // Shallow copy.  The const_cast is ugly but does no harm:
      // ArrayRCP is a pointer, so there's no reason why a method that
      // returns it can't be marked const.  (If you're not changing
      // the pointer itself, e.g., by reallocating, the pointer is
      // const.  It doesn't matter whether or not you can change the
      // values to which the pointer points.)
      MVT::initializeValues (lclMV_, source.lclMV_.getNumRows (),
                             source.lclMV_.getNumCols (),
                             const_cast<KMV&> (source.lclMV_).getValuesNonConst (),
                             source.lclMV_.getStride (),
                             source.lclMV_.getOrigNumRows (),
                             source.lclMV_.getOrigNumCols ());
      // Don't forget to copy whichVectors_ (for non-constant-stride view).
      if (source.whichVectors_.size () > 0) {
        whichVectors_ = source.whichVectors_; // Array::op= does a deep copy
      }
    }
    else if (copyOrView == Teuchos::Copy) {
      using Teuchos::ArrayRCP;
      using Teuchos::RCP;

      // copy data from the source MultiVector into this multivector
      RCP<Node> node = MVT::getNode(source.lclMV_);
      const LocalOrdinal myLen = source.getLocalLength();
      const size_t numVecs = source.getNumVectors();

      ArrayRCP<Scalar> data = (myLen > 0) ?
        node->template allocBuffer<Scalar> (myLen * numVecs) :
        Teuchos::null;
      const size_t stride = (myLen > 0) ? myLen : size_t (0);

      // This just sets the dimensions, pointer, and stride of lclMV_.
      MVT::initializeValues (lclMV_, myLen, numVecs, data, stride);
      // This actually copies the data.  It uses the Node's
      // parallel_for to copy, which should ensure first-touch
      // allocation on systems that support it.
      if (source.isConstantStride ()) {
        MVT::Assign (lclMV_, source.lclMV_);
      }
      else {
        MVT::Assign (lclMV_, source.lclMV_, source.whichVectors_);
      }
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::MultiVector copy constructor: "
        "The second argument 'copyOrView' has an invalid value " << copyOrView
        << ".  Valid values include Teuchos::Copy = " << Teuchos::Copy <<
        " and Teuchos::View = " << Teuchos::View << ".");
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const KokkosClassic::MultiVector<Scalar,Node>& localMultiVector,
               EPrivateComputeViewConstructor /* dummy */) :
    DO (map),
    lclMV_ (localMultiVector),
    hasViewSemantics_ (false) // for now, not for long though
  {
    const size_t localNumElts = map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      localMultiVector.getNumRows () < localNumElts,
      std::invalid_argument,
      "Tpetra::MultiVector constructor: The Map contains " << localNumElts
      << " on process " << map->getComm ()->getRank () << ", but the given "
      "KokkosClassic::MultiVector only has " << localMultiVector.getNumRows ()
      << " rows.");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayRCP<Scalar>& view,
               size_t LDA,
               size_t NumVectors,
               EPrivateHostViewConstructor /* dummy */) :
    DO (map),
    lclMV_ (map->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    const char tfecfFuncName[] = "MultiVector(view,LDA,NumVector)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive, but you specified NumVectors = "
      << NumVectors << ".");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (view.size ()) < LDA*(NumVectors - 1) + myLen,
      std::invalid_argument,
      ": view does not contain enough data to specify the entries in this.");
    MVT::initializeValues (lclMV_, myLen, NumVectors, view, myLen);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
               Teuchos::ArrayRCP<Scalar> data,
               size_t LDA,
               size_t NumVectors,
               EPrivateComputeViewConstructor /* dummy */) :
    DO (map),
    lclMV_ (map->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    const char tfecfFuncName[] = "MultiVector(data,LDA,NumVector)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (data.size ()) < LDA*(NumVectors - 1) + myLen,
      std::runtime_error,
      ": data does not contain enough data to specify the entries in this.");
#endif
    MVT::initializeValues(lclMV_,myLen,NumVectors,data,LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
               Teuchos::ArrayRCP<Scalar> data,
               size_t LDA,
               Teuchos::ArrayView<const size_t> WhichVectors,
               EPrivateComputeViewConstructor /* dummy */) :
    DO (map),
    lclMV_ (map->getNode ()),
    whichVectors_ (WhichVectors),
    hasViewSemantics_ (false) // for now, not for long though
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    const char tfecfFuncName[] = "MultiVector(data,LDA,WhichVectors)";
    const size_t myLen = getLocalLength();
    size_t maxVector = *std::max_element(WhichVectors.begin(), WhichVectors.end());
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::runtime_error,
      ": LDA must be large enough to accomodate the local entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (data.size ()) < LDA * maxVector + myLen, std::runtime_error,
      ": data does not contain enough data to specify the entries in this.");
#endif
    if (WhichVectors.size () == 1) {
      // shift data so that desired vector is vector 0
      maxVector = 0;
      data += LDA * WhichVectors[0];
      // kill whichVectors_; we are constant stride
      whichVectors_.clear ();
    }
    if (myLen > 0) {
      MVT::initializeValues(lclMV_,myLen,maxVector+1,data,LDA);
    }
    else {
      MVT::initializeValues(lclMV_,0,WhichVectors.size(),Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayRCP<Scalar>& data,
               const size_t LDA,
               const size_t numVecs) :
    DO (map),
    lclMV_ (map->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {
    const char tfecfFuncName[] = "MultiVector(map,data,LDA,numVecs)";
    const size_t numRows = this->getLocalLength ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < numRows, std::runtime_error,
      ": LDA = " << LDA << " < numRows = " << numRows << ".");
    MVT::initializeValues (lclMV_, numRows, numVecs, data, LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
               const KokkosClassic::MultiVector<Scalar,Node>& localMultiVector,
               Teuchos::ArrayView<const size_t> WhichVectors,
               EPrivateComputeViewConstructor /* dummy */) :
    DO (map),
    lclMV_ (localMultiVector),
    whichVectors_ (WhichVectors),
    hasViewSemantics_ (false) // for now, not for long though
  {
    const size_t localNumElts = map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      localMultiVector.getNumRows () < localNumElts,
      std::invalid_argument,
      "Tpetra::MultiVector constructor: The Map contains " << localNumElts
      << " on process " << map->getComm ()->getRank () << ", but the given "
      "KokkosClassic::MultiVector only has " << localMultiVector.getNumRows ()
      << " rows.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      WhichVectors.size () == 0,
      std::invalid_argument,
      "Tpetra::MultiVector constructor: WhichVectors.size() == 0.  "
      "The noncontiguous column view constructor only makes sense for a "
      "nonzero number of columns to view.");

    // Optimization for the case of a single column: just make a
    // contiguous ("constant stride") view of the one column.
    if (WhichVectors.size () == 1) {
      const size_t offsetRow = 0;
      const size_t offsetCol = WhichVectors[0];
      lclMV_ = lclMV_.offsetViewNonConst (localNumElts, 1, offsetRow, offsetCol);
      whichVectors_.clear (); // This being empty defines "constant stride."
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const Scalar>& A,
               size_t LDA,
               size_t NumVectors) :
    DO (map),
    lclMV_ (map->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    const char tfecfFuncName[] = "MultiVector(A,LDA)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::runtime_error,
      ": LDA must be large enough to accomodate the local entries.");
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (A.size ()) < LDA*(NumVectors - 1) + myLen,
      std::runtime_error,
      ": A does not contain enough data to specify the entries in this.");
#endif
    if (myLen > 0) {
      RCP<Node> node = MVT::getNode(lclMV_);
      ArrayRCP<Scalar> mydata = node->template allocBuffer<Scalar>(myLen*NumVectors);
      MVT::initializeValues(lclMV_,myLen,NumVectors,mydata,myLen);
      // FIXME (mfh 13 Sep 2012) It would be better to use the Kokkos
      // Node's copyToBuffer method to push data directly to the
      // device pointer, rather than using an intermediate host
      // buffer.  Also, we should have an optimization for the
      // contiguous storage case (constant stride, and the stride
      // equals the local number of rows).
      ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(KokkosClassic::WriteOnly,myLen*NumVectors,mydata);
      typename ArrayView<const Scalar>::iterator srcit = A.begin();
      for (size_t j = 0; j < NumVectors; ++j) {
        std::copy(srcit,srcit+myLen,myview);
        srcit += LDA;
        myview += myLen;
      }
      mydata = Teuchos::null;
      myview = Teuchos::null;
    }
    else {
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const ArrayView<const Scalar> >& ArrayOfPtrs,
               size_t NumVectors) :
    DO (map),
    lclMV_ (map->getNode ()),
    hasViewSemantics_ (false) // for now, not for long though
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    const char tfecfFuncName[] = "MultiVector(ArrayOfPtrs)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      NumVectors < 1 || NumVectors != static_cast<size_t>(ArrayOfPtrs.size()),
      std::runtime_error,
      ": ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
    const size_t myLen = getLocalLength();
    if (myLen > 0) {
      RCP<Node> node = MVT::getNode(lclMV_);
      ArrayRCP<Scalar> mydata = node->template allocBuffer<Scalar>(myLen*NumVectors);
      MVT::initializeValues(lclMV_,myLen,NumVectors,mydata,myLen);
      ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(KokkosClassic::WriteOnly,myLen*NumVectors,mydata);
      for (size_t j = 0; j < NumVectors; ++j) {
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          static_cast<size_t>(ArrayOfPtrs[j].size()) != getLocalLength(),
          std::runtime_error,
          ": ArrayOfPtrs[" << j << "].size() (== " << ArrayOfPtrs[j].size()
          << ") is not equal to getLocalLength() (== " << getLocalLength());
#endif
        typename ArrayView<const Scalar>::iterator src = ArrayOfPtrs[j].begin();
        std::copy(src,src+myLen,myview);
        myview += myLen;
      }
      myview = Teuchos::null;
      mydata = Teuchos::null;
    }
    else {
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~MultiVector() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isConstantStride() const {
    return whichVectors_.empty();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalLength() const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getNodeNumElements ();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalLength() const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getGlobalNumElements ();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getStride() const {
    if (isConstantStride()) {
      return MVT::getStride(lclMV_);
    }
    return 0;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  checkSizes (const SrcDistObject& sourceObj)
  {
    // Check whether the source object is a MultiVector.  If not, then
    // we can't even compare sizes, so it's definitely not OK to
    // Import or Export from it.
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
    const this_type* src = dynamic_cast<const this_type*> (&sourceObj);
    if (src == NULL) {
      return false;
    } else {
      // The target of the Import or Export calls checkSizes() in
      // DistObject::doTransfer().  By that point, we've already
      // constructed an Import or Export object using the two
      // multivectors' Maps, which means that (hopefully) we've
      // already checked other attributes of the multivectors.  Thus,
      // all we need to do here is check the number of columns.
      return src->getNumVectors () == this->getNumVectors ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  constantNumberOfPackets () const {
    return this->getNumVectors ();
  }

#if TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  copyAndPermute (
    const SrcDistObject& sourceObj,
    size_t numSameIDs,
    const Kokkos::View<const LocalOrdinal*, device_type> &permuteToLIDs,
    const Kokkos::View<const LocalOrdinal*, device_type> &permuteFromLIDs)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
    const char tfecfFuncName[] = "copyAndPermute";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      ": permuteToLIDs and permuteFromLIDs must have the same size."
      << std::endl << "permuteToLIDs.size() = " << permuteToLIDs.size ()
      << " != permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    const size_t numCols = this->getNumVectors ();
    const size_t stride = MVT::getStride (lclMV_);

    // Copy rows [0, numSameIDs-1] of the local multivectors.
    //
    // For GPU Nodes: All of this happens using device pointers; this
    // does not require host views of either source or destination.
    if (isConstantStride ()) {
      const size_t numSrcCols = MVT::getNumCols (sourceMV.lclMV_);
      const size_t numDestCols = MVT::getNumCols (lclMV_);

      const KMV src = sourceMV.lclMV_.offsetView (numSameIDs, numSrcCols, 0, 0);
      KMV dest = lclMV_.offsetViewNonConst (numSameIDs, numDestCols, 0, 0);
      if (sourceMV.isConstantStride ()) {
        MVT::Assign (dest, src);
      }
      else {
        MVT::Assign (dest, src, sourceMV.whichVectors_ ());
      }
    }
    else {
      // Copy the columns one at a time, since MVT doesn't have an
      // Assign method for noncontiguous access to the destination
      // multivector.
      for (size_t j = 0; j < numCols; ++j) {
        const size_t destCol = isConstantStride () ? j : whichVectors_[j];
        const size_t srcCol = sourceMV.isConstantStride () ?
          j : sourceMV.whichVectors_[j];
        KMV dest_j = lclMV_.offsetViewNonConst (numSameIDs, 1, 0, destCol);
        const KMV src_j = sourceMV.lclMV_.offsetView (numSameIDs, 1, 0, srcCol);
        MVT::Assign (dest_j, src_j); // Copy src_j into dest_j
      }
    }

    // For the remaining GIDs, execute the permutations.  This may
    // involve noncontiguous access of both source and destination
    // vectors, depending on the LID lists.
    //
    // FIXME (mfh 20 June 2012) For an Export with duplicate GIDs on
    // the same process, this merges their values by replacement of
    // the last encountered GID, not by the specified merge rule
    // (such as ADD).

    // If there are no permutations, we are done
    if (permuteFromLIDs.size() == 0 || permuteToLIDs.size() == 0)
      return;

    if (this->isConstantStride ()) {
      Details::PermuteArrayMultiColumnConstantStride<Scalar,LocalOrdinal,device_type> perm;
      perm.permuteToLIDs = permuteToLIDs;
      perm.permuteFromLIDs = permuteFromLIDs;
      perm.src = sourceMV.getKokkosView();
      perm.dest = getKokkosViewNonConst();
      perm.src_stride = MVT::getStride (sourceMV.lclMV_);
      perm.dest_stride = stride;
      perm.numCols = numCols;

      perm.permute();
    }
    else {
      Details::PermuteArrayMultiColumnVariableStride<Scalar,LocalOrdinal,device_type> perm;
      perm.permuteToLIDs = permuteToLIDs;
      perm.permuteFromLIDs = permuteFromLIDs;
      perm.src = sourceMV.getKokkosView();
      perm.dest = getKokkosViewNonConst();
      perm.src_whichVectors =
        getKokkosViewDeepCopy<device_type>(sourceMV.whichVectors_ ());
      perm.dest_whichVectors =
        getKokkosViewDeepCopy<device_type>(whichVectors_ ());
      perm.src_stride = MVT::getStride (sourceMV.lclMV_);
      perm.dest_stride = stride;
      perm.numCols = numCols;

      perm.permute();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  packAndPrepare (
    const SrcDistObject& sourceObj,
    const Kokkos::View<const LocalOrdinal*, device_type> &exportLIDs,
    Kokkos::View<Scalar*, device_type> &exports,
    const Kokkos::View<size_t*, device_type> &numExportPacketsPerLID,
    size_t& constantNumPackets,
    Distributor & /* distor */ )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Array<size_t>::size_type size_type;

    // If we have no exports, there is nothing to do
    if (exportLIDs.size() == 0)
      return;

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    // We don't need numExportPacketsPerLID; forestall "unused
    // variable" compile warnings.
    (void) numExportPacketsPerLID;

    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    const size_t numCols = sourceMV.getNumVectors ();
    const size_t stride = MVT::getStride (sourceMV.lclMV_);
    // This spares us from needing to fill numExportPacketsPerLID.
    // Setting constantNumPackets to a nonzero value signals that
    // all packets have the same number of entries.
    constantNumPackets = numCols;

    const size_t numExportLIDs = exportLIDs.size ();
    const size_t newExportsSize = numCols * numExportLIDs;
    if (exports.size () != newExportsSize) {
      Kokkos::Compat::realloc(exports, newExportsSize);
    }

    if (numCols == 1) { // special case for one column only
      // MultiVector always represents a single column with constant
      // stride, but it doesn't hurt to implement both cases anyway.
      if (sourceMV.isConstantStride ()) {
        Details::PackArraySingleColumnConstantStride<Scalar,LocalOrdinal,device_type> pack;
        pack.exportLIDs = exportLIDs;
        pack.exports = exports;
        pack.src = sourceMV.getKokkosView();
        pack.pack();
      }
      else {
        Details::PackArraySingleColumnOffset<Scalar,LocalOrdinal,device_type> pack;
        pack.exportLIDs = exportLIDs;
        pack.exports = exports;
        pack.src = sourceMV.getKokkosView();
        pack.offset = sourceMV.whichVectors_[0] * stride;
        pack.pack();
      }
    }
    else { // the source MultiVector has multiple columns
      if (sourceMV.isConstantStride ()) {
        Details::PackArrayMultiColumnConstantStride<Scalar,LocalOrdinal,device_type> pack;
        pack.exportLIDs = exportLIDs;
        pack.exports = exports;
        pack.src = sourceMV.getKokkosView();
        pack.stride = stride;
        pack.numCols = numCols;
        pack.pack ();
      }
      else {
        Details::PackArrayMultiColumnVariableStride<Scalar,LocalOrdinal,device_type> pack;
        pack.exportLIDs = exportLIDs;
        pack.exports = exports;
        pack.src = sourceMV.getKokkosView();
        pack.srcWhichVectors =
          getKokkosViewDeepCopy<device_type>(sourceMV.whichVectors_ ());
        pack.stride = stride;
        pack.numCols = numCols;
        pack.pack ();
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  unpackAndCombine (
    const Kokkos::View<const LocalOrdinal*, device_type> &importLIDs,
    const Kokkos::View<const Scalar*, device_type> &imports,
    const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
    size_t constantNumPackets,
    Distributor & /* distor */,
    CombineMode CM)
  {
    using Teuchos::ArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
    const char tfecfFuncName[] = "unpackAndCombine";

    // If we have no imports, there is nothing to do
    if (importLIDs.size() == 0)
      return;

    // We don't need numPacketsPerLID; forestall "unused variable"
    // compile warnings.
    (void) numPacketsPerLID;

    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    const size_t numVecs = getNumVectors ();

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (imports.size()) != getNumVectors()*importLIDs.size(),
      std::runtime_error,
      ": 'imports' buffer size must be consistent with the amount of data to "
      "be sent.  " << std::endl << "imports.size() = " << imports.size()
      << " != getNumVectors()*importLIDs.size() = " << getNumVectors() << "*"
      << importLIDs.size() << " = " << getNumVectors() * importLIDs.size()
      << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      constantNumPackets == static_cast<size_t> (0), std::runtime_error,
      ": constantNumPackets input argument must be nonzero.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numVecs) != static_cast<size_t> (constantNumPackets),
      std::runtime_error, ": constantNumPackets must equal numVecs.");
#endif // HAVE_TPETRA_DEBUG

    if (numVecs > 0 && importLIDs.size () > 0) {
      const size_t myStride = MVT::getStride (lclMV_);

      // NOTE (mfh 10 Mar 2012) If you want to implement custom
      // combine modes, start editing here.  Also, if you trust
      // inlining, it would be nice to condense this code by using a
      // binary function object f:
      //
      // ncview_[...] = f (ncview_[...], *impptr++);
      if (CM == INSERT || CM == REPLACE) {
        if (isConstantStride()) {
          Details::UnpackArrayMultiColumnConstantStride<Scalar,LocalOrdinal,Details::InsertOp,device_type> unpack;
          unpack.importLIDs = importLIDs;
          unpack.imports = imports;
          unpack.dest = getKokkosViewNonConst();
          unpack.stride = myStride;
          unpack.numCols = numVecs;
          unpack.op = Details::InsertOp();

          unpack.unpack();
        }
        else {
          Details::UnpackArrayMultiColumnVariableStride<Scalar,LocalOrdinal,Details::InsertOp,device_type> unpack;
          unpack.importLIDs = importLIDs;
          unpack.imports = imports;
          unpack.dest = getKokkosViewNonConst();
          unpack.whichVectors =
            getKokkosViewDeepCopy<device_type>(whichVectors_ ());
          unpack.stride = myStride;
          unpack.numCols = numVecs;
          unpack.op = Details::InsertOp();

          unpack.unpack();
        }
      }
      else if (CM == ADD) {
        if (isConstantStride()) {
          Details::UnpackArrayMultiColumnConstantStride<Scalar,LocalOrdinal,Details::AddOp,device_type> unpack;
          unpack.importLIDs = importLIDs;
          unpack.imports = imports;
          unpack.dest = getKokkosViewNonConst();
          unpack.stride = myStride;
          unpack.numCols = numVecs;
          unpack.op = Details::AddOp();

          unpack.unpack();
        }
        else {
          Details::UnpackArrayMultiColumnVariableStride<Scalar,LocalOrdinal,Details::AddOp,device_type> unpack;
          unpack.importLIDs = importLIDs;
          unpack.imports = imports;
          unpack.dest = getKokkosViewNonConst();
          unpack.whichVectors =
            getKokkosViewDeepCopy<device_type>(whichVectors_ ());
          unpack.stride = myStride;
          unpack.numCols = numVecs;
          unpack.op = Details::AddOp();

          unpack.unpack();
        }
      }
      else if (CM == ABSMAX) {
        if (isConstantStride()) {
          Details::UnpackArrayMultiColumnConstantStride<Scalar,LocalOrdinal,Details::AbsMaxOp,device_type> unpack;
          unpack.importLIDs = importLIDs;
          unpack.imports = imports;
          unpack.dest = getKokkosViewNonConst();
          unpack.stride = myStride;
          unpack.numCols = numVecs;
          unpack.op = Details::AbsMaxOp();

          unpack.unpack();
        }
        else {
          Details::UnpackArrayMultiColumnVariableStride<Scalar,LocalOrdinal,Details::AbsMaxOp,device_type> unpack;
          unpack.importLIDs = importLIDs;
          unpack.imports = imports;
          unpack.dest = getKokkosViewNonConst();
          unpack.whichVectors =
            getKokkosViewDeepCopy<device_type>(whichVectors_ ());
          unpack.stride = myStride;
          unpack.numCols = numVecs;
          unpack.op = Details::AbsMaxOp();

          unpack.unpack();
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          CM != ADD && CM != REPLACE && CM != INSERT && CM != ABSMAX,
          std::invalid_argument, ": Invalid CombineMode: " << CM << ".  Valid "
          "CombineMode values are ADD, REPLACE, INSERT, and ABSMAX.");
      }
    }
  }

#else

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  copyAndPermute (const SrcDistObject& sourceObj,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                  const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
    const char tfecfFuncName[] = "copyAndPermute";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      ": permuteToLIDs and permuteFromLIDs must have the same size."
      << std::endl << "permuteToLIDs.size() = " << permuteToLIDs.size ()
      << " != permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    const size_t numCols = this->getNumVectors ();

    // Copy rows [0, numSameIDs-1] of the local multivectors.
    //
    // For GPU Nodes: All of this happens using device pointers; this
    // does not require host views of either source or destination.
    //
    // mfh 10 Jan 2013: In the current implementation, device data has
    // already been copied to the host views.  If we do a
    // device-device copy here, it won't synch the host views
    // correctly.  Thus, we replace this step with a host copy, if
    // Node is a GPU Node.  That ensures correctness for now, though
    // it's really something we need to fix if we plan to run
    // efficiently on GPUs.
    if (Node::isHostNode) {
      if (isConstantStride ()) {
        const size_t numSrcCols = MVT::getNumCols (sourceMV.lclMV_);
        const size_t numDestCols = MVT::getNumCols (lclMV_);

        const KMV src = sourceMV.lclMV_.offsetView (numSameIDs, numSrcCols, 0, 0);
        KMV dest = lclMV_.offsetViewNonConst (numSameIDs, numDestCols, 0, 0);
        if (sourceMV.isConstantStride ()) {
          MVT::Assign (dest, src);
        }
        else {
          MVT::Assign (dest, src, sourceMV.whichVectors_ ());
        }
      }
      else {
        // Copy the columns one at a time, since MVT doesn't have an
        // Assign method for noncontiguous access to the destination
        // multivector.
        for (size_t j = 0; j < numCols; ++j) {
          const size_t destCol = isConstantStride () ? j : whichVectors_[j];
          const size_t srcCol = sourceMV.isConstantStride () ?
            j : sourceMV.whichVectors_[j];
          KMV dest_j = lclMV_.offsetViewNonConst (numSameIDs, 1, 0, destCol);
          const KMV src_j = sourceMV.lclMV_.offsetView (numSameIDs, 1, 0, srcCol);
          MVT::Assign (dest_j, src_j); // Copy src_j into dest_j
        }
      }
    }

    // Copy the vectors one at a time.  This works whether or not the
    // multivectors have constant stride.
    for (size_t j = 0; j < numCols; ++j) {
      // Get host views of the current source and destination column.
      //
      // FIXME (mfh 10 Mar 2012) Copying should really be done on the
      // device.  There's no need to bring everything back to the
      // host, copy there, and then copy back.
      ArrayView<const Scalar> srcView =
        sourceMV.getSubArrayRCP (sourceMV.cview_, j) ();
      ArrayView<Scalar> dstView = getSubArrayRCP (ncview_, j) ();

      if (! Node::isHostNode) {
        // The first numSameIDs IDs are the same between source and
        // target, so we can just copy the data.  (This favors faster
        // contiguous access whenever we can do it.)
        std::copy (srcView.begin (), srcView.begin () + numSameIDs,
                   dstView.begin ());
      }

      // For the remaining GIDs, execute the permutations.  This may
      // involve noncontiguous access of both source and destination
      // vectors, depending on the LID lists.
      //
      // FIXME (mfh 09 Apr 2012) Permutation should be done on the
      // device, not on the host.
      //
      // FIXME (mfh 20 June 2012) For an Export with duplicate GIDs on
      // the same process, this merges their values by replacement of
      // the last encountered GID, not by the specified merge rule
      // (such as ADD).
      const size_type numPermuteLIDs =
        std::min (permuteToLIDs.size (), permuteFromLIDs.size ());
      for (size_type k = 0; k < numPermuteLIDs; ++k) {
        dstView[permuteToLIDs[k]] = srcView[permuteFromLIDs[k]];
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  packAndPrepare (const SrcDistObject& sourceObj,
                  const ArrayView<const LocalOrdinal> &exportLIDs,
                  Array<Scalar> &exports,
                  const ArrayView<size_t> &numExportPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor & /* distor */ )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Array<size_t>::size_type size_type;

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    // We don't need numExportPacketsPerLID; forestall "unused
    // variable" compile warnings.
    (void) numExportPacketsPerLID;

    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    const size_t numCols = sourceMV.getNumVectors ();
    const size_t stride = MVT::getStride (sourceMV.lclMV_);
    // This spares us from needing to fill numExportPacketsPerLID.
    // Setting constantNumPackets to a nonzero value signals that
    // all packets have the same number of entries.
    constantNumPackets = numCols;

    const size_type numExportLIDs = exportLIDs.size ();
    const size_type newExportsSize = numCols * numExportLIDs;
    if (exports.size () != newExportsSize) {
      exports.resize (newExportsSize);
    }

    ArrayView<const Scalar> srcView = sourceMV.cview_ ();
    if (numCols == 1) { // special case for one column only
      // MultiVector always represents a single column with constant
      // stride, but it doesn't hurt to implement both cases anyway.
      if (sourceMV.isConstantStride ()) {
        for (size_type k = 0; k < numExportLIDs; ++k) {
          exports[k] = srcView[exportLIDs[k]];
        }
      }
      else {
        const size_t offset = sourceMV.whichVectors_[0] * stride;
        for (size_type k = 0; k < numExportLIDs; ++k) {
          exports[k] = srcView[exportLIDs[k] + offset];
        }
      }
    }
    else { // the source MultiVector has multiple columns
      typename Array<Scalar>::iterator expptr = exports.begin ();
      if (sourceMV.isConstantStride ()) {
        for (size_type k = 0; k < numExportLIDs; ++k) {
          const size_t localRow = static_cast<size_t> (exportLIDs[k]);
          for (size_t j = 0; j < numCols; ++j) {
            *expptr++ = srcView[localRow + j*stride];
          }
        }
      }
      else {
        ArrayView<const size_t> srcWhichVectors = sourceMV.whichVectors_ ();
        for (size_type k = 0; k < numExportLIDs; ++k) {
          const size_t localRow = static_cast<size_t> (exportLIDs[k]);
          for (size_t j = 0; j < numCols; ++j) {
            *expptr++ = srcView[localRow + srcWhichVectors[j]*stride];
          }
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  unpackAndCombine (const ArrayView<const LocalOrdinal> &importLIDs,
                    const ArrayView<const Scalar> &imports,
                    const ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor & /* distor */,
                    CombineMode CM)
  {
    using Teuchos::ArrayView;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
    const char tfecfFuncName[] = "unpackAndCombine";

    // We don't need numPacketsPerLID; forestall "unused variable"
    // compile warnings.
    (void) numPacketsPerLID;

    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    const size_t numVecs = getNumVectors ();

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (imports.size()) != getNumVectors()*importLIDs.size(),
      std::runtime_error,
      ": 'imports' buffer size must be consistent with the amount of data to "
      "be sent.  " << std::endl << "imports.size() = " << imports.size()
      << " != getNumVectors()*importLIDs.size() = " << getNumVectors() << "*"
      << importLIDs.size() << " = " << getNumVectors() * importLIDs.size()
      << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      constantNumPackets == static_cast<size_t> (0), std::runtime_error,
      ": constantNumPackets input argument must be nonzero.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numVecs) != static_cast<size_t> (constantNumPackets),
      std::runtime_error, ": constantNumPackets must equal numVecs.");
#endif // HAVE_TPETRA_DEBUG

    if (numVecs > 0 && importLIDs.size () > 0) {
      typename ArrayView<const Scalar>::iterator impptr = imports.begin ();
      ArrayView<Scalar> destView = ncview_ ();
      const size_type numImportLIDs = importLIDs.size ();
      const size_t myStride = MVT::getStride (lclMV_);

      // NOTE (mfh 10 Mar 2012) If you want to implement custom
      // combine modes, start editing here.  Also, if you trust
      // inlining, it would be nice to condense this code by using a
      // binary function object f:
      //
      // ncview_[...] = f (ncview_[...], *impptr++);
      if (CM == INSERT || CM == REPLACE) {
        if (isConstantStride()) {
          for (size_type k = 0; k < numImportLIDs; ++k) {
            const size_t localRow = static_cast<size_t> (importLIDs[k]);
            for (size_t j = 0; j < numVecs; ++j) {
              destView[localRow + myStride*j] = *impptr++;
            }
          }
        }
        else {
          for (size_type k = 0; k < numImportLIDs; ++k) {
            const size_t localRow = static_cast<size_t> (importLIDs[k]);
            for (size_t j = 0; j < numVecs; ++j) {
              destView[localRow + myStride*whichVectors_[j]] = *impptr++;
            }
          }
        }
      }
      else if (CM == ADD) {
        if (isConstantStride()) {
          for (size_type k = 0; k < numImportLIDs; ++k) {
            const size_t localRow = static_cast<size_t> (importLIDs[k]);
            for (size_t j = 0; j < numVecs; ++j) {
              destView[localRow + myStride*j] += *impptr++;
            }
          }
        }
        else {
          for (size_type k = 0; k < numImportLIDs; ++k) {
            const size_t localRow = static_cast<size_t> (importLIDs[k]);
            for (size_t j = 0; j < numVecs; ++j) {
              destView[localRow + myStride*whichVectors_[j]] += *impptr++;
            }
          }
        }
      }
      else if (CM == ABSMAX) {
        if (isConstantStride()) {
          for (size_type k = 0; k < numImportLIDs; ++k) {
            const size_t localRow = static_cast<size_t> (importLIDs[k]);
            for (size_t j = 0; j < numVecs; ++j) {
              Scalar &curval       = destView[localRow + myStride*j];
              const Scalar &newval = *impptr++;
              curval = std::max( SCT::magnitude(curval), SCT::magnitude(newval) );
            }
          }
        }
        else {
          for (size_type k = 0; k < numImportLIDs; ++k) {
            const size_t localRow = static_cast<size_t> (importLIDs[k]);
            for (size_t j = 0; j < numVecs; ++j) {
              Scalar &curval       = destView[localRow + myStride*whichVectors_[j]];
              const Scalar &newval = *impptr++;
              curval = std::max( SCT::magnitude(curval), SCT::magnitude(newval) );
            }
          }
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          CM != ADD && CM != REPLACE && CM != INSERT && CM != ABSMAX,
          std::invalid_argument, ": Invalid CombineMode: " << CM << ".  Valid "
          "CombineMode values are ADD, REPLACE, INSERT, and ABSMAX.");
      }
    }
  }

#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  inline size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNumVectors() const {
    if (isConstantStride()) {
      return MVT::getNumCols(lclMV_);
    }
    else {
      return whichVectors_.size();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
       const Teuchos::ArrayView<Scalar> &dots) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;

    const char tfecfFuncName[] = "dot()";
    const size_t myLen   = getLocalLength();
    const size_t numVecs = getNumVectors();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != numVecs, std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (dots.size ()) != numVecs, std::runtime_error,
      ": dots.size() must be as large as the number of vectors in *this and A.");
    if (isConstantStride() && A.isConstantStride()) {
      MVT::Dot(lclMV_,A.lclMV_,dots);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr  = arcp_const_cast<Scalar>(MVT::getValues(lclMV_));
      ArrayRCP<Scalar> avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj  =   getSubArrayRCP( vptr,j);
        ArrayRCP<Scalar> avj = A.getSubArrayRCP(avptr,j);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        dots[j] = MVT::Dot((const KMV&)v,(const KMV &)a);
      }
    }
    if (this->isDistributed ()) {
      Array<Scalar> ldots(dots);
      Teuchos::reduceAll (* (this->getMap ()->getComm ()), Teuchos::REDUCE_SUM,
                          static_cast<int> (numVecs), ldots.getRawPtr (),
                          dots.getRawPtr ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  norm2 (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    const size_t numVecs = this->getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      numNorms < numVecs, std::runtime_error, "Tpetra::MultiVector::norm2: "
      "'norms' must have at least as many entries as the number of vectors in "
      "*this.  norms.size() = " << numNorms << " < this->getNumVectors() = "
      << numVecs << ".");

    if (isConstantStride ()) {
      MVT::Norm2Squared (lclMV_, norms);
    }
    else {
      KMV v (MVT::getNode (lclMV_));
      ArrayRCP<Scalar> vi;
      for (size_t j = 0; j < numVecs; ++j) {
        vi = arcp_const_cast<Scalar> (MVT::getValues (lclMV_, whichVectors_[j]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vi,
                               MVT::getStride (lclMV_));
        norms[j] = MVT::Norm2Squared (v);
      }
    }
    if (this->isDistributed ()) {
      RCP<const Comm<int> > comm = this->getMap ().is_null () ? null :
        this->getMap ()->getComm ();
      // The calling process only participates in the collective if
      // both the Map and its Comm on that process are nonnull.
      if (! comm.is_null ()) {
        Array<MT> lnorms (norms);
        reduceAll<int, MT> (*comm, REDUCE_SUM, static_cast<int> (numVecs),
                            lnorms.getRawPtr (), norms.getRawPtr ());
      }
    }

    for (typename Teuchos::ArrayView<MT>::iterator n = norms.begin ();
         n != norms.begin () + numVecs; ++n) {
      (*n) = STM::squareroot (*n);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& weights,
                const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;

    const char tfecfFuncName[] = "normWeighted: ";

    const MT OneOverN =
      STM::one () / static_cast<MT> (this->getGlobalLength ());
    const size_t numVecs = this->getNumVectors ();
    const size_t myLen = this->getLocalLength ();
    const bool OneW = (weights.getNumVectors () == 1);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! OneW && weights.getNumVectors () != numVecs, std::runtime_error,
      "The input MultiVector of weights must either have only one vector "
      "(column), or it must have the same number of vectors as this.  "
      "weights.getNumVectors() = " << weights.getNumVectors () << ", and "
      "this->getNumVectors() = " << numVecs << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (norms.size ()) != numVecs, std::invalid_argument,
      "norms.size() must be as large as the number of vectors in *this.  "
      "norms.size() = " << norms.size () << " != this->getNumVectors() = "
      << numVecs << ".");

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! this->getMap ()->isCompatible (* (weights.getMap ())),
      std::runtime_error, "MultiVectors do not have compatible Maps.");
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      myLen != weights.getLocalLength (), std::runtime_error,
      "This MultiVector and the input argument 'weights' do not have the same "
      "local length.  this->getLocalLength() = " << myLen << " != weights."
      "getLocalLength() = " << weights.getLocalLength () << ".");
#endif // HAVE_TPETRA_DEBUG

    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const Scalar> X_j = this->getData (j);
      ArrayRCP<const Scalar> W_j =
        OneW ? weights.getData (0) : weights.getData (j);
      MT norm_j = STM::zero ();
      for (size_t i = 0; i < myLen; ++i) {
        const Scalar tmp = X_j[i] / W_j[i];
        const MT tmp_real = STS::real (tmp);
        const MT tmp_imag = STS::imag (tmp);
        norm_j += tmp_real * tmp_real + tmp_imag * tmp_imag;
      }
      norms[j] = norm_j;
    }

    if (this->isDistributed ()) {
      Array<MT> lnorms (norms);
      reduceAll<int, MT> (* (this->getMap ()->getComm ()), REDUCE_SUM,
                          static_cast<int> (numVecs), lnorms.getRawPtr (),
                          norms.getRawPtr ());
    }
    for (size_t k = 0; k < numVecs; ++k) {
      norms[k] = STM::squareroot (norms[k] * OneOverN);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  norm1 (const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    const size_t numVecs = this->getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      numNorms < numVecs, std::runtime_error, "Tpetra::MultiVector::norm1: "
      "'norms' must have at least as many entries as the number of vectors in "
      "*this.  norms.size() = " << numNorms << " < this->getNumVectors() = "
      << numVecs << ".");

    if (isConstantStride ()) {
      MVT::Norm1 (lclMV_, norms);
    }
    else {
      KMV v (MVT::getNode (lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j = 0; j < numVecs; ++j) {
        vj = arcp_const_cast<Scalar> (MVT::getValues (lclMV_, whichVectors_[j]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj,
                               MVT::getStride (lclMV_));
        norms[j] = MVT::Norm1 (v);
      }
    }
    if (this->isDistributed ()) {
      RCP<const Comm<int> > comm = this->getMap ().is_null () ? null :
        this->getMap ()->getComm ();
      // The calling process only participates in the collective if
      // both the Map and its Comm on that process are nonnull.
      if (! comm.is_null ()) {
        Array<MT> lnorms (norms);
        reduceAll<int, MT> (*comm, REDUCE_SUM, static_cast<int> (numVecs),
                            lnorms.getRawPtr (), norms.getRawPtr ());
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  normInf (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_MAX;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(norms.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::normInf(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride ()) {
      MVT::NormInf(lclMV_,norms);
    }
    else {
      KMV v (MVT::getNode (lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j = 0; j < numVecs; ++j) {
        vj = arcp_const_cast<Scalar> (MVT::getValues (lclMV_, whichVectors_[j]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj,
                               MVT::getStride (lclMV_));
        norms[j] = MVT::NormInf (v);
      }
    }
    if (this->isDistributed ()) {
      RCP<const Comm<int> > comm = this->getMap ().is_null () ? null :
        this->getMap ()->getComm ();
      // The calling process only participates in the collective if
      // both the Map and its Comm on that process are nonnull.
      if (! comm.is_null ()) {
        Array<MT> lnorms (norms);
        reduceAll<int, MT> (*comm, REDUCE_MAX, static_cast<int> (numVecs),
                            lnorms.getRawPtr (), norms.getRawPtr ());
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  meanValue (const Teuchos::ArrayView<Scalar> &means) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Teuchos::ScalarTraits<Scalar> SCT;

    const size_t numVecs = getNumVectors();
    const size_t myLen   = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(means.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::meanValue(): means.size() must be as large as the number of vectors in *this.");
    // compute local components of the means
    // sum these across all nodes
    if (isConstantStride()) {
      MVT::Sum(lclMV_,means);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = arcp_const_cast<Scalar>(MVT::getValues(lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        means[j] = MVT::Sum((const KMV &)v);
      }
    }
    if (this->isDistributed ()) {
      Array<Scalar> lmeans(means);
      // only combine if we are a distributed MV
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM,
                 static_cast<int> (numVecs), lmeans.getRawPtr (),
                 means.getRawPtr ());
    }
    // mfh 12 Apr 2012: Don't take out the cast from the ordinal type
    // to the magnitude type, since operator/ (std::complex<T>, int)
    // isn't necessarily defined.
    const Scalar OneOverN =
      SCT::one() / static_cast<typename SCT::magnitudeType> (getGlobalLength ());
    for (typename ArrayView<Scalar>::iterator i = means.begin();
         i != means.begin() + numVecs;
         ++i)
    {
      (*i) = (*i) * OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  randomize()
  {
    if (isConstantStride ()) {
      MVT::Random (lclMV_);
    }
    else {
      const size_t numVecs = this->getNumVectors ();
      KMV v (MVT::getNode (lclMV_));
      Teuchos::ArrayRCP<Scalar> vj;
      for (size_t j = 0; j < numVecs; ++j) {
        vj = MVT::getValuesNonConst (lclMV_, whichVectors_[j]);
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj, MVT::getStride (lclMV_));
        MVT::Random (v);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  putScalar (const Scalar &alpha)
  {
    const size_t numVecs = getNumVectors();
    if (isConstantStride()) {
      MVT::Init(lclMV_,alpha);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      Teuchos::ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = MVT::getValuesNonConst(lclMV_,whichVectors_[j]);
        MVT::initializeValues(v,MVT::getNumRows(lclMV_), 1, vj, MVT::getStride(lclMV_));
        MVT::Init(v,alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  replaceMap (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& newMap)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::Comm;
    using Teuchos::RCP;

    // mfh 28 Mar 2013: This method doesn't forget whichVectors_, so
    // it might work if the MV is a column view of another MV.
    // However, things might go wrong when restoring the original
    // Map, so we don't allow this case for now.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->isConstantStride (), std::logic_error,
      "Tpetra::MultiVector::replaceMap: This method does not currently work "
      "if the MultiVector is a column view of another MultiVector (that is, if "
      "isConstantStride() == false).");

    // Case 1: current Map and new Map are both nonnull on this process.
    // Case 2: current Map is nonnull, new Map is null.
    // Case 3: current Map is null, new Map is nonnull.
    // Case 4: both Maps are null: forbidden.
    //
    // Case 1 means that we don't have to do anything on this process,
    // other than assign the new Map.  (We always have to do that.)
    // It's an error for the user to supply a Map that requires
    // resizing in this case.
    //
    // Case 2 means that the calling process is in the current Map's
    // communicator, but will be excluded from the new Map's
    // communicator.  We don't have to do anything on the calling
    // process; just leave whatever data it may have alone.
    //
    // Case 3 means that the calling process is excluded from the
    // current Map's communicator, but will be included in the new
    // Map's communicator.  This means we need to (re)allocate the
    // local (KokkosClassic::)MultiVector if it does not have the right
    // number of rows.  If the new number of rows is nonzero, we'll
    // fill the newly allocated local data with zeros, as befits a
    // projection operation.
    //
    // The typical use case for Case 3 is that the MultiVector was
    // first created with the Map with more processes, then that Map
    // was replaced with a Map with fewer processes, and finally the
    // original Map was restored on this call to replaceMap.

#ifdef HAVE_TEUCHOS_DEBUG
    // mfh 28 Mar 2013: We can't check for compatibility across the
    // whole communicator, unless we know that the current and new
    // Maps are nonnull on _all_ participating processes.
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   origNumProcs == newNumProcs && ! this->getMap ()->isCompatible (*map),
    //   std::invalid_argument, "Tpetra::MultiVector::project: "
    //   "If the input Map's communicator is compatible (has the same number of "
    //   "processes as) the current Map's communicator, then the two Maps must be "
    //   "compatible.  The replaceMap() method is not for data redistribution; "
    //   "use Import or Export for that purpose.");

    // TODO (mfh 28 Mar 2013) Add compatibility checks for projections
    // of the Map, in case the process counts don't match.
#endif // HAVE_TEUCHOS_DEBUG

    if (this->getMap ().is_null ()) { // current Map is null
      TEUCHOS_TEST_FOR_EXCEPTION(
        newMap.is_null (), std::invalid_argument,
        "Tpetra::MultiVector::replaceMap: both current and new Maps are null.  "
        "This probably means that the input Map is incorrect.");
      // Case 3: current Map is null, new Map is nonnull.
      const size_t newNumRows = newMap->getNodeNumElements ();
      const size_t origNumRows = lclMV_.getNumRows ();
      const size_t numCols = getNumVectors ();

      if (origNumRows != newNumRows) {
        RCP<Node> node = newMap->getNode ();
        ArrayRCP<Scalar> data = newNumRows == 0 ? Teuchos::null :
          node->template allocBuffer<Scalar> (newNumRows * numCols);
        const size_t stride = newNumRows;
        MVT::initializeValues (lclMV_, newNumRows, numCols, data, stride);
        if (newNumRows > 0) {
          MVT::Init (lclMV_, Teuchos::ScalarTraits<Scalar>::zero ());
        }
      }
    }
    else if (newMap.is_null ()) { // Case 2: current Map is nonnull, new Map is null
      // I am an excluded process.  Reinitialize my data so that I
      // have 0 rows.  Keep the number of columns as before.
      const size_t numVecs = getNumVectors ();
      MVT::initializeValues (lclMV_, 0, numVecs, Teuchos::null, 0);
    }
    this->map_ = newMap;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  scale (const Scalar &alpha)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    // NOTE: can't substitute putScalar(0.0) for scale(0.0), because
    //       the former will overwrite NaNs present in the MultiVector, while the
    //       semantics of this call require multiplying them by 0, which IEEE requires to be NaN
    const size_t numVecs = getNumVectors();
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // do nothing
    }
    else if (isConstantStride()) {
      MVT::Scale(lclMV_,alpha);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = arcp_const_cast<Scalar> (MVT::getValues(lclMV_, whichVectors_[j]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj, MVT::getStride (lclMV_));
        MVT::Scale (v, alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  scale (Teuchos::ArrayView<const Scalar> alphas)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(alphas.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as the number of vectors in *this.");

    // Unfortunately, MVT doesn't have a Scale method that works on a
    // whole MultiVector.  Thus, we have to scale a column at a time.
    if (this->isConstantStride ()) {
      for (size_t j = 0; j < numVecs; ++j) {
        const Scalar alpha_j = alphas[j];
        const size_t curCol = j;

        // If alphas[j] is zero, don't scale; just fill with zeros.
        // This follows the BLAS convention.  (zero*NaN is NaN.)
        //
        // Similarly, if alphas[j] is one, don't scale; just do nothing.
        // This follows the BLAS convention.  (one*NaN is NaN.)
        if (alpha_j == Teuchos::ScalarTraits<Scalar>::zero ()) {
          KMV MV_j = lclMV_.offsetViewNonConst (lclMV_.getNumRows (), 1, 0, curCol);
          MVT::Init (MV_j, alpha_j);
        }
        else if (alpha_j != Teuchos::ScalarTraits<Scalar>::one ()) {
          KMV MV_j = lclMV_.offsetViewNonConst (lclMV_.getNumRows (), 1, 0, curCol);
          MVT::Scale (MV_j, alpha_j);
        }
      }
    } else { // not constant stride
      for (size_t j = 0; j < numVecs; ++j) {
        const Scalar alpha_j = alphas[j];
        const size_t curCol = whichVectors_[j];

        // If alphas[j] is zero, don't scale; just fill with zeros.
        // This follows the BLAS convention.  (zero*NaN is NaN.)
        //
        // Similarly, if alphas[j] is one, don't scale; just do nothing.
        // This follows the BLAS convention.  (one*NaN is NaN.)
        if (alpha_j == Teuchos::ScalarTraits<Scalar>::zero ()) {
          KMV MV_j = lclMV_.offsetViewNonConst (lclMV_.getNumRows (), 1, 0, curCol);
          MVT::Init (MV_j, alpha_j);
        }
        else if (alpha_j != Teuchos::ScalarTraits<Scalar>::one ()) {
          KMV MV_j = lclMV_.offsetViewNonConst (lclMV_.getNumRows (), 1, 0, curCol);
          MVT::Scale (MV_j, alpha_j);
        }
      }
    }

    // KMV vec(MVT::getNode(lclMV_));
    // const size_t myLen = MVT::getNumRows(lclMV_);
    // if (myLen == 0) return;
    // ArrayRCP<Scalar> mybuf = MVT::getValuesNonConst(lclMV_);
    // for (size_t j = 0; j < numVecs; ++j) {
    //   if (alphas[j] == Teuchos::ScalarTraits<Scalar>::one()) {
    //     // do nothing: NaN * 1.0 == NaN, Number*1.0 == Number
    //   }
    //   else {
    //     ArrayRCP<Scalar> mybufj = getSubArrayRCP(mybuf,j);
    //     MVT::initializeValues(vec,myLen,1,mybufj,myLen);
    //     MVT::Scale(vec,alphas[j]);
    //   }
    // }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  scale (const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const char tfecfFuncName[] = "scale(alpha,A)";

    const size_t numVecs = getNumVectors(),
                 myLen   = getLocalLength();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! this->getMap ()->isCompatible (*A.getMap ()),
      std::runtime_error, ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *(this->getMap ())
        << "A.getMap(): " << std::endl << *(A.getMap ()) << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != numVecs, std::runtime_error,
      ": MultiVectors must have the same number of vectors.");
    if (isConstantStride() && A.isConstantStride()) {
      // set me == alpha*A
      MVT::Scale(lclMV_,alpha,(const KMV&)A.lclMV_);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr  = MVT::getValuesNonConst (lclMV_);
      ArrayRCP<Scalar> avptr = arcp_const_cast<Scalar> (MVT::getValues (A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar>  vj =   getSubArrayRCP (vptr,  j);
        ArrayRCP<Scalar> avj = A.getSubArrayRCP (avptr, j);
        MVT::initializeValues (a,myLen, 1, avj, myLen);
        MVT::initializeValues (v,myLen, 1,  vj, myLen);
        MVT::Scale (v, alpha, (const KMV &)a);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  reciprocal (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const char tfecfFuncName[] = "reciprocal()";

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    const size_t numVecs = getNumVectors();
    const size_t myLen = getLocalLength();
    try {
      if (isConstantStride() && A.isConstantStride()) {
        MVT::Recip(lclMV_,(const KMV&)A.lclMV_);
      }
      else {
        KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
        ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                        avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
        for (size_t j=0; j < numVecs; ++j) {
          ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                          avj = A.getSubArrayRCP(avptr,j);
          MVT::initializeValues(a,myLen, 1, avj, myLen);
          MVT::initializeValues(v,myLen, 1,  vj, myLen);
          MVT::Recip(v,(const KMV &)a);
        }
      }
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true,std::runtime_error,
          ": caught exception from Kokkos:" << std::endl
          << e.what() << std::endl);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const char tfecfFuncName[] = "abs()";
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    const size_t myLen = getLocalLength();
    const size_t numVecs = getNumVectors();
    if (isConstantStride() && A.isConstantStride()) {
      MVT::Abs(lclMV_,(const KMV&)A.lclMV_);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                      avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                        avj = A.getSubArrayRCP(avptr,j);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        MVT::Abs(v,(const KMV &)a);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
                      const Scalar &beta)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const char tfecfFuncName[] = "update()";
    // this = beta*this + alpha*A
    // must support case where &this == &A
    // can't short circuit on alpha==0.0 or beta==0.0, because 0.0*NaN != 0.0
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()), std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    const size_t myLen = getLocalLength();
    const size_t numVecs = getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    if (isConstantStride() && A.isConstantStride()) {
      MVT::GESUM(lclMV_,alpha,(const KMV&)A.lclMV_,beta);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                      avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                        avj = A.getSubArrayRCP(avptr,j);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        MVT::GESUM(v,alpha,(const KMV &)a,beta);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
                      const Scalar &beta,  const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                      const Scalar &gamma)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const char tfecfFuncName[] = "update()";
    // this = alpha*A + beta*B + gamma*this
    // must support case where &this == &A or &this == &B
    // can't short circuit on alpha==0.0 or beta==0.0 or gamma==0.0, because 0.0*NaN != 0.0
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*A.getMap()) || !this->getMap()->isCompatible(*B.getMap()),
        std::runtime_error,
        ": MultiVectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "A.getMap(): " << std::endl << *A.getMap() << std::endl
        << "B.getMap(): " << std::endl << *B.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength() || getLocalLength() != B.getLocalLength(), std::runtime_error,
        ": MultiVectors do not have the same local length.");
#endif
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors() || B.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");
    const size_t myLen = getLocalLength();
    const size_t numVecs = getNumVectors();
    if (isConstantStride() && A.isConstantStride() && B.isConstantStride()) {
      MVT::GESUM(lclMV_,alpha,(const KMV&)A.lclMV_,beta,(const KMV&)B.lclMV_,gamma);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), a(MVT::getNode(lclMV_)), b(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = MVT::getValuesNonConst(lclMV_),
                      avptr = arcp_const_cast<Scalar>(MVT::getValues(A.lclMV_)),
                      bvptr = arcp_const_cast<Scalar>(MVT::getValues(B.lclMV_));
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj =   getSubArrayRCP( vptr,j),
                        avj = A.getSubArrayRCP(avptr,j),
                        bvj = B.getSubArrayRCP(bvptr,j);
        MVT::initializeValues(b,myLen, 1, bvj, myLen);
        MVT::initializeValues(a,myLen, 1, avj, myLen);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        MVT::GESUM(v,alpha,(const KMV&)a,beta,(const KMV&)b,gamma);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getData (size_t j) const
  {
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    return node->template viewBuffer<Scalar> (getLocalLength (), getSubArrayRCP (MVT::getValues(lclMV_), j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getDataNonConst(size_t j)
  {
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    return node->template viewBufferNonConst<Scalar> (KokkosClassic::ReadWrite, getLocalLength (), getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  operator= (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source)
  {
    if (source.hasViewSemantics_) {
      this->map_ = source.map_; // "op=" for DistObject
      lclMV_ = source.lclMV_; // Shallow copy

      // Don't forget to copy whichVectors_ (for non-constant-stride view).
      if (source.whichVectors_.size () > 0) {
        whichVectors_ = source.whichVectors_; // Array::op= does a deep copy
      }
      // View semantics are "contagious" for the "classic" version of
      // MultiVector.  That is, the right-hand side of the assignment
      // transmits view semantics to the left-hand side.  (The Kokkos
      // refactor version of MultiVector always has view semantics.)
      hasViewSemantics_ = true;
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      using Teuchos::outArg;
      using Teuchos::REDUCE_MIN;
      using Teuchos::reduceAll;
      using std::endl;
#endif // HAVE_TPETRA_DEBUG

      const char tfecfFuncName[] = "operator=";
      // Check for special case of this=Source, in which case we do nothing.
      if (this != &source) {
        // Whether the input and *this are compatible on the calling process.
        const int locallyCompat =
          (this->getLocalLength () == source.getLocalLength ()) ? 1 : 0;
#ifdef HAVE_TPETRA_DEBUG
        int globallyCompat = 1; // innocent until proven guilty
        reduceAll<int, int> (* (this->getMap ()->getComm ()), REDUCE_MIN,
                             locallyCompat, outArg (globallyCompat));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          globallyCompat == 0, std::invalid_argument,
          ": MultiVectors do not have the same local length on all processes.");
#else
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          locallyCompat == 0, std::invalid_argument,
          ": MultiVectors do not have the same local length on the calling "
          "process " << this->getMap ()->getComm ()->getSize () << ".  *this "
          "has " << this->getLocalLength () << " local rows, and source has "
          << source.getLocalLength () << " rows.");
#endif // HAVE_TPETRA_DEBUG

        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          source.getNumVectors() != getNumVectors(), std::invalid_argument,
          ": MultiVectors must have the same number of vectors.");

        Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
        const size_t numVecs = getNumVectors();
        if (isConstantStride () && source.isConstantStride () &&
            getLocalLength () == getStride () &&
            source.getLocalLength ()== source.getStride ()) {
          // Both multivectors' data are stored contiguously, so we can
          // copy in one call.
          node->template copyBuffers<Scalar> (getLocalLength () * numVecs,
                                              MVT::getValues (source.lclMV_),
                                              MVT::getValuesNonConst (lclMV_));
        }
        else {
          // We have to copy the columns one at a time.
          for (size_t j=0; j < numVecs; ++j) {
            node->template copyBuffers<Scalar> (getLocalLength (),
                                                source.getSubArrayRCP (MVT::getValues (source.lclMV_), j),
                                                getSubArrayRCP (MVT::getValuesNonConst(lclMV_), j));
          }
        }
      }
    }
    return *this;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subCopy (const Teuchos::ArrayView<const size_t>& cols) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Teuchos::ArrayView<const size_t>::size_type size_type;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    const size_type numCopyVecs = static_cast<size_type> (cols.size ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numCopyVecs < 1, std::runtime_error, "Tpetra::MultiVector::subCopy"
      "(const Teuchos::ArrayView<const size_t>&): The input array cols must "
      "contain at least one entry, but cols.size() = " << cols.size ()
      << " == 0.");

    // Check whether the index set in cols is contiguous.  If it is,
    // use the more efficient Range1D version of subCopy.
    bool contiguous = true;
    for (size_type j = 1; j < numCopyVecs; ++j) {
      if (cols[j] != cols[j-1] + static_cast<size_t> (1)) {
        contiguous = false;
        break;
      }
    }

    const bool zeroData = false;
    if (contiguous) {
      if (numCopyVecs == 0) {
        // The output MV has no columns, so there is nothing to copy.
        return rcp (new MV (this->getMap (), numCopyVecs, zeroData));
      } else {
        // Use the more efficient contiguous-index-range version.
        return this->subCopy (Teuchos::Range1D (cols[0], cols[numCopyVecs-1]));
      }
    }
    else { // Copy data one column at a time.
      RCP<MV> mv = rcp (new MV (this->getMap (), numCopyVecs, zeroData));
      RCP<Node> node = this->getMap ()->getNode ();
      for (size_type j = 0; j < numCopyVecs; ++j) {
        node->template copyBuffers<Scalar> (getLocalLength (),
                                            getSubArrayRCP (MVT::getValues (lclMV_), cols[j]),
                                            MVT::getValuesNonConst (mv->lclMV_, j));
      }
      return mv;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subCopy (const Teuchos::Range1D& colRng) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(
      colRng.size () == 0, std::runtime_error, "Tpetra::MultiVector::"
      "subCopy(Range1D): Range must include at least one vector.");

    const size_t numCopyVecs = colRng.size ();
    const bool zeroData = false;
    RCP<Node> node = MVT::getNode (lclMV_);
    RCP<MV> mv;
    // mv is allocated with constant stride
    mv = rcp (new MV (this->getMap (), numCopyVecs, zeroData));
    // copy data from *this into mv
    for (size_t js = colRng.lbound (), jd = 0; jd < numCopyVecs; ++jd, ++js) {
      node->template copyBuffers<Scalar> (getLocalLength (),
                                          getSubArrayRCP (MVT::getValues (lclMV_), js),
                                          MVT::getValuesNonConst (mv->lclMV_, jd));
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
              size_t offset) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    const size_t newNumRows = subMap->getNodeNumElements();
    const bool tooManyElts = newNumRows + offset > lclMV_.getOrigNumRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > MVT::getNumRows (lclMV_),
        std::runtime_error,
        "Tpetra::MultiVector::offsetView: Invalid input Map.  Input Map owns "
        << subMap->getNodeNumElements() << " elements on process " << myRank
        << ".  offset = " << offset << ".  Yet, the MultiVector contains only "
        << lclMV_.getOrigNumRows () << " on this process.");
    }
    RCP<const MV> constViewMV;
    if (isConstantStride()) {
      KokkosClassic::MultiVector<Scalar, Node> newLocalMV =
        lclMV_.offsetView (newNumRows, lclMV_.getNumCols (), offset, 0);
      constViewMV = rcp (new MV (subMap, newLocalMV, COMPUTE_VIEW_CONSTRUCTOR));
    }
    else {
      // Compute the max column being viewed.  This tells us where to stop the view.
      const size_t maxCol = *std::max_element (whichVectors_.begin(), whichVectors_.end());
      KokkosClassic::MultiVector<Scalar, Node> newLocalMV =
        lclMV_.offsetView (newNumRows, maxCol+1, offset, 0);
      constViewMV = rcp (new MV (subMap, newLocalMV, whichVectors_, COMPUTE_VIEW_CONSTRUCTOR));
    }
    return constViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                      size_t offset)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    const size_t newNumRows = subMap->getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > lclMV_.getOrigNumRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > MVT::getNumRows (lclMV_),
        std::runtime_error,
        "Tpetra::MultiVector::offsetViewNonconst: Invalid input Map.  Input Map owns "
        << subMap->getNodeNumElements() << " elements on process " << myRank
        << ".  offset = " << offset << ".  Yet, the MultiVector contains only "
        << lclMV_.getOrigNumRows () << " on this process.");
    }
    RCP<MV> subViewMV;
    if (isConstantStride()) {
      KokkosClassic::MultiVector<Scalar, Node> newLocalMV =
        lclMV_.offsetViewNonConst (newNumRows, lclMV_.getNumCols (), offset, 0);
      subViewMV = rcp (new MV (subMap, newLocalMV, COMPUTE_VIEW_CONSTRUCTOR));
    }
    else {
      // Compute the max column being viewed.  This tells us where to stop the view.
      const size_t maxCol = *std::max_element (whichVectors_.begin(), whichVectors_.end());
      KokkosClassic::MultiVector<Scalar, Node> newLocalMV =
        lclMV_.offsetViewNonConst (newNumRows, maxCol+1, offset, 0);
      subViewMV = rcp (new MV (subMap, newLocalMV, whichVectors_, COMPUTE_VIEW_CONSTRUCTOR));
    }
    return subViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subView (const ArrayView<const size_t>& cols) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Teuchos::ArrayView<size_t>::size_type size_type;

    const size_type numViewCols = static_cast<size_type> (cols.size ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numViewCols < 1, std::runtime_error, "Tpetra::MultiVector::subView"
      "(const Teuchos::ArrayView<const size_t>&): The input array 'cols' "
      "must contain at least one entry, but cols.size() = " << cols.size ()
      << " == 0.");

    // Check whether the index set in cols is contiguous.  If it is,
    // use the more efficient Range1D version of subView.
    bool contiguous = true;
    for (size_type j = 1; j < numViewCols; ++j) {
      if (cols[j] != cols[j-1] + static_cast<size_t> (1)) {
        contiguous = false;
        break;
      }
    }
    if (contiguous) {
      if (numViewCols == 0) {
        // The output MV has no columns, so there is nothing to view.
        return rcp (new MV (this->getMap (), numViewCols));
      } else {
        // Use the more efficient contiguous-index-range version.
        return this->subView (Teuchos::Range1D (cols[0], cols[numViewCols-1]));
      }
    }

    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning
    const size_t myStride = MVT::getStride (lclMV_);
    const size_t myLen = MVT::getNumRows (lclMV_);
    // use the smallest view possible of the buffer: from the first
    // element of the minInd vector to the last element of the maxInd
    // vector.  this minimizes overlap between views, and keeps view
    // of the minimum amount necessary, in order to allow node to
    // achieve maximum efficiency.  adjust the indices appropriately;
    // shift so that smallest index is 0
    ArrayRCP<const Scalar> cbuf = MVT::getValues(lclMV_);
    ArrayRCP<Scalar>      ncbuf = arcp_const_cast<Scalar> (cbuf);
    Array<size_t> newCols (numViewCols);
    size_t minInd = Teuchos::OrdinalTraits<size_t>::max();
    size_t maxInd = Teuchos::OrdinalTraits<size_t>::zero();
    if (isConstantStride()) {
      for (size_type j = 0; j < numViewCols; ++j) {
        newCols[j] = cols[j];
        if (newCols[j] < minInd) minInd = newCols[j];
        if (maxInd < newCols[j]) maxInd = newCols[j];
      }
    }
    else {
      for (size_type j = 0; j < numViewCols; ++j) {
        newCols[j] = whichVectors_[cols[j]];
        if (newCols[j] < minInd) minInd = newCols[j];
        if (maxInd < newCols[j]) maxInd = newCols[j];
      }
    }
    ArrayRCP<Scalar> minbuf =
      ncbuf.persistingView (minInd * myStride,
                            myStride * (maxInd - minInd) + myLen);
    for (size_type j = 0; j < numViewCols; ++j) {
      newCols[j] -= minInd;
    }
    RCP<const MV> constViewMV;
    return rcp_const_cast<const MV> (rcp (new MV (this->getMap (), minbuf,
                                                  myStride, newCols (),
                                                  COMPUTE_VIEW_CONSTRUCTOR)));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subView (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
    size_t numViewVecs = colRng.size();
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning
    ArrayRCP<const Scalar> cbuf = MVT::getValues(lclMV_);
    ArrayRCP<Scalar>      ncbuf = arcp_const_cast<Scalar>(cbuf);
    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      ArrayRCP<Scalar> subdata = ncbuf.persistingView( MVT::getStride(lclMV_) * colRng.lbound(),
                                                       MVT::getStride(lclMV_) * (numViewVecs-1) + getLocalLength() );
      return rcp (new MV (this->getMap (), subdata, MVT::getStride (lclMV_), numViewVecs, COMPUTE_VIEW_CONSTRUCTOR));
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Array<size_t> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    return rcp (new MV (this->getMap (), ncbuf, MVT::getStride (lclMV_), whchvecs, COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subViewNonConst (const ArrayView<const size_t> &cols)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    const size_t numViewCols = static_cast<size_t> (cols.size ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numViewCols < 1, std::runtime_error, "Tpetra::MultiVector::"
      "subViewNonConst(const Teuchos::ArrayView<const size_t>&): The input "
      "array 'cols' must contain at least one entry, but cols.size() = "
      << cols.size () << " == 0.");

    // Check whether the index set in cols is contiguous.  If it is,
    // use the more efficient Range1D version of subViewNonConst.
    bool contiguous = true;
    for (size_t j = 1; j < numViewCols; ++j) {
      if (cols[j] != cols[j-1] + static_cast<size_t> (1)) {
        contiguous = false;
        break;
      }
    }
    if (contiguous) {
      if (numViewCols == 0) {
        // The output MV has no columns, so there is nothing to view.
        return rcp (new MV (this->getMap (), numViewCols));
      } else {
        // Use the more efficient contiguous-index-range version.
        return this->subViewNonConst (Range1D (cols[0], cols[numViewCols-1]));
      }
    }

    if (isConstantStride ()) {
      return rcp (new MV (this->getMap (), MVT::getValuesNonConst (lclMV_),
                          MVT::getStride (lclMV_), cols,
                          COMPUTE_VIEW_CONSTRUCTOR));
    }
    // else, lookup current whichVectors_ using cols
    Array<size_t> newcols (numViewCols);
    for (size_t j = 0; j < numViewCols; ++j) {
      newcols[j] = whichVectors_[cols[j]];
    }
    return rcp (new MV (this->getMap (), MVT::getValuesNonConst (lclMV_),
                        MVT::getStride (lclMV_), newcols (),
                        COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  subViewNonConst (const Teuchos::Range1D &colRng)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subViewNonConst(Range1D): range must include at least one vector.");
    size_t numViewVecs = colRng.size();
    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride ()) {
      // view goes from first entry of first vector to last entry of last vector
      const size_t stride = MVT::getStride(lclMV_);
      ArrayRCP<Scalar> data = MVT::getValuesNonConst(lclMV_);
      ArrayRCP<Scalar> subdata = data.persistingView( stride * colRng.lbound(),
                                                      stride * (numViewVecs-1) + getLocalLength() );
      return rcp (new MV (this->getMap (), subdata, stride, numViewVecs, COMPUTE_VIEW_CONSTRUCTOR));
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Array<size_t> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    const size_t stride = MVT::getStride(lclMV_);
    ArrayRCP<Scalar> data = MVT::getValuesNonConst(lclMV_);
    return rcp (new MV (this->getMap (), data, stride, whchvecs, COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getVector (size_t j) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVector(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif
    // this is const, so lclMV_ is const, so we get const buff
    // it is safe to cast away the const because we will wrap it in a const Vector below
    ArrayRCP<Scalar> ncbuff;
    if (getLocalLength() > 0) {
      ArrayRCP<const Scalar> cbuff = getSubArrayRCP(MVT::getValues(lclMV_),j);
      ncbuff = arcp_const_cast<Scalar>(cbuff);
    }
    return rcp_const_cast<const V> (rcp (new V (this->getMap (), ncbuff, COMPUTE_VIEW_CONSTRUCTOR)));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getVectorNonConst(size_t j)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVectorNonConst(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif
    ArrayRCP<Scalar> ncbuff;
    if (getLocalLength() > 0) {
      ncbuff = getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j);
    }
    return rcp (new V (this->getMap (), ncbuff, COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  get1dCopy (Teuchos::ArrayView<Scalar> A, size_t LDA) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const char tfecfFuncName[] = "get1dCopy(A,LDA)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < getLocalLength(), std::runtime_error,
      ": specified stride is not large enough for local vector length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (A.size ()) < LDA * (getNumVectors() - 1) + getLocalLength (),
      std::runtime_error, ": specified stride/storage is not large enough for "
      "the number of vectors.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t myStride = MVT::getStride(lclMV_),
                  numCols = getNumVectors(),
                  myLen   = getLocalLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mydata = MVT::getValues(lclMV_);
      ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,mydata);
      typename ArrayView<Scalar>::iterator Aptr = A.begin();
      for (size_t j=0; j<numCols; j++) {
        ArrayRCP<const Scalar> myviewj = getSubArrayRCP(myview,j);
        std::copy(myviewj,myviewj+myLen,Aptr);
        Aptr += LDA;
      }
      myview = Teuchos::null;
      mydata = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  get2dCopy (Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const char tfecfFuncName[] = "get2dCopy(ArrayOfPtrs)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (ArrayOfPtrs.size ()) != getNumVectors(), std::runtime_error,
      ": Array of pointers must contain as many pointers as the MultiVector has rows.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t numCols = getNumVectors(),
                   myLen = getLocalLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mybuff = MVT::getValues(lclMV_);
      ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(mybuff.size(), mybuff);
      for (size_t j=0; j<numCols; ++j) {
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          static_cast<size_t> (ArrayOfPtrs[j].size ()) != getLocalLength (),
          std::runtime_error, ": The ArrayView provided in ArrayOfPtrs[" << j
          << "] was not large enough to contain the local entries.");
#endif
        ArrayRCP<const Scalar> myviewj = getSubArrayRCP(myview,j);
        std::copy(myviewj,myviewj+myLen,ArrayOfPtrs[j].begin());
      }
      myview = Teuchos::null;
      mybuff = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dView () const
  {
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      return Teuchos::null;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
        "get1dView: This MultiVector does not have constant stride, so it is "
        "not possible to view its data as a single array.  You may check "
        "whether a MultiVector has constant stride by calling "
        "isConstantStride().");
      Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
      return node->template viewBuffer<Scalar> (getStride () * (getNumVectors () - 1) +
                                                getLocalLength (),
                                                MVT::getValues (lclMV_));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dViewNonConst ()
  {
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      return Teuchos::null;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
        "get1dViewNonConst: This MultiVector does not have constant stride, so "
        "it is not possible to view its data as a single array.  You may check "
        "whether a MultiVector has constant stride by calling "
        "isConstantStride().");
      Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
      return node->template viewBufferNonConst<Scalar> (KokkosClassic::ReadWrite,
                                                        getStride () * (getNumVectors () - 1) +
                                                        getLocalLength (),
                                                        MVT::getValuesNonConst (lclMV_));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get2dViewNonConst()
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    RCP<Node> node = MVT::getNode(lclMV_);
    ArrayRCP<ArrayRCP<Scalar> > views = arcp<ArrayRCP<Scalar> > (getNumVectors ());
    if (isConstantStride ()) {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = getNumVectors(),
                    myLen   = getLocalLength();
      if (myLen > 0) {
        ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(KokkosClassic::ReadWrite,myStride*(numCols-1)+myLen,MVT::getValuesNonConst(lclMV_));
        for (size_t j=0; j<numCols; ++j) {
          views[j] = myview.persistingView(0,myLen);
          myview += myStride;
        }
      }
    }
    else {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = MVT::getNumCols(lclMV_),
                     myCols = getNumVectors(),
                     myLen  = MVT::getNumRows(lclMV_);
      if (myLen > 0) {
        ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(KokkosClassic::ReadWrite,myStride*(numCols-1)+myLen,MVT::getValuesNonConst(lclMV_));
        for (size_t j=0; j<myCols; ++j) {
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
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    RCP<Node> node = MVT::getNode(lclMV_);
    ArrayRCP< ArrayRCP<const Scalar> > views = arcp<ArrayRCP<const Scalar> >(getNumVectors());
    if (isConstantStride()) {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = getNumVectors(),
                    myLen   = getLocalLength();
      if (myLen > 0) {
        ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,MVT::getValues(lclMV_));
        for (size_t j=0; j<numCols; ++j) {
          views[j] = myview.persistingView(0,myLen);
          myview += myStride;
        }
      }
    }
    else {
      const size_t myStride = MVT::getStride(lclMV_),
                    numCols = MVT::getNumCols(lclMV_),
                     myCols = getNumVectors(),
                     myLen  = MVT::getNumRows(lclMV_);
      if (myLen > 0) {
        ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,MVT::getValues(lclMV_));
        for (size_t j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  multiply (Teuchos::ETransp transA,
            Teuchos::ETransp transB,
            const Scalar &alpha,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
            const Scalar &beta)
  {
    using Teuchos::NO_TRANS;      // enums
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::null;
    using Teuchos::ScalarTraits;  // traits
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

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

    const char errPrefix[] = "Tpetra::MultiVector::multiply(transOpA,transOpB,alpha,A,B,beta): ";

    TEUCHOS_TEST_FOR_EXCEPTION( ScalarTraits<Scalar>::isComplex && (transA == TRANS || transB == TRANS), std::invalid_argument,
        errPrefix << "non-conjugate transpose not supported for complex types.");
    transA = (transA == NO_TRANS ? NO_TRANS : CONJ_TRANS);
    transB = (transB == NO_TRANS ? NO_TRANS : CONJ_TRANS);

    // Compute effective dimensions, w.r.t. transpose operations on
    size_t A_nrows = (transA==CONJ_TRANS) ? A.getNumVectors() : A.getLocalLength();
    size_t A_ncols = (transA==CONJ_TRANS) ? A.getLocalLength() : A.getNumVectors();
    size_t B_nrows = (transB==CONJ_TRANS) ? B.getNumVectors() : B.getLocalLength();
    size_t B_ncols = (transB==CONJ_TRANS) ? B.getLocalLength() : B.getNumVectors();

    Scalar beta_local = beta; // local copy of beta; might be reassigned below

    TEUCHOS_TEST_FOR_EXCEPTION(
      getLocalLength() != A_nrows || getNumVectors() != B_ncols || A_ncols != B_nrows,
      std::runtime_error,
      errPrefix << "dimension of *this, op(A) and op(B) must be consistent.  "
      << std::endl << "The local part of *this is "
      << getLocalLength() << " x " << getNumVectors()
      << ", A is " << A_nrows << " x " << A_ncols
      << ", and B is " << B_nrows << " x " << B_ncols << ".");

    bool A_is_local = ! A.isDistributed ();
    bool B_is_local = ! B.isDistributed ();
    bool C_is_local = ! this->isDistributed ();
    // Case 1: C(local) = A^X(local) * B^X(local)
    bool Case1 = (C_is_local &&  A_is_local && B_is_local);
    // Case 2: C(local) = A^T(distr) * B  (distr)
    bool Case2 = (C_is_local && ! A_is_local && ! B_is_local &&
                  transA == CONJ_TRANS && transB == NO_TRANS);
    // Case 3: C(distr) = A  (distr) * B^X(local)
    bool Case3 = (! C_is_local && ! A_is_local && B_is_local &&
                  transA == NO_TRANS);

    // Test that we are considering a meaningful cases
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! Case1 && ! Case2 && ! Case3, std::runtime_error,
      errPrefix << "multiplication of op(A) and op(B) into *this is not a "
      "supported use case.");

    if (beta != ScalarTraits<Scalar>::zero() && Case2) {
      // If Case2, then C is local and contributions must be summed
      // across all processes.  However, if beta != 0, then accumulate
      // beta*C into the sum.  When summing across all processes, we
      // only want to accumulate this once, so set beta == 0 on all
      // processes except Process 0.
      const int MyPID = this->getMap ()->getComm ()->getRank ();
      if (MyPID!=0) beta_local = ScalarTraits<Scalar>::zero();
    }

    // Check if A, B, C have constant stride, if not then make temp copy (strided)
    RCP<const MV> Atmp, Btmp;
    RCP<MV>       Ctmp;
    if (isConstantStride() == false) Ctmp = rcp (new MV (*this));
    else Ctmp = rcp(this,false);

    if (A.isConstantStride() == false) Atmp = rcp (new MV (A));
    else Atmp = rcp(&A,false);

    if (B.isConstantStride() == false) Btmp = rcp (new MV (B));
    else Btmp = rcp(&B,false);

#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(!Ctmp->isConstantStride() || !Btmp->isConstantStride() || !Atmp->isConstantStride(), std::logic_error,
        errPrefix << "failed making temporary strided copies of input multivectors.");
#endif

    KMV &C_mv = Ctmp->lclMV_;
    {
      // get local multivectors
      const KMV &A_mv = Atmp->lclMV_;
      const KMV &B_mv = Btmp->lclMV_;
      // do the multiply (GEMM)
      MVT::GEMM(C_mv,transA,transB,alpha,A_mv,B_mv,beta_local);
    }

    // Dispose of (possibly) extra copies of A, B
    Atmp = null;
    Btmp = null;

    RCP<Node> node = MVT::getNode(lclMV_);
    // If *this was not strided, copy the data from the strided version and then delete it
    if (! isConstantStride ()) {
      // *this is not strided, we must put data from Ctmp into *this
      TEUCHOS_TEST_FOR_EXCEPT(&C_mv != &lclMV_);
      const size_t numVecs = MVT::getNumCols(lclMV_);
      for (size_t j=0; j < numVecs; ++j) {
        node->template copyBuffers<Scalar>(getLocalLength(),MVT::getValues(C_mv,j),MVT::getValuesNonConst(lclMV_,whichVectors_[j]));
      }
    }

    // If Case 2 then sum up *this and distribute it to all processors.
    if (Case2) {
      this->reduce();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  elementWiseMultiply (Scalar scalarAB,
                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                       const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
                       Scalar scalarThis)
  {
    using Teuchos::arcp_const_cast;
    const char tfecfFuncName[] = "elementWiseMultiply()";

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength() != A.getLocalLength() ||
      getLocalLength() != B.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif // HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      B.getNumVectors() != this->getNumVectors(), std::runtime_error,
      ": MultiVectors 'this' and B must have the same number of vectors.");
    try {
      MVT::ElemMult (lclMV_, scalarThis, scalarAB, (const KMV&) A.lclMV_,
                     (const KMV&) B.lclMV_);
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true,std::runtime_error,
          ": caught exception from Kokkos:" << std::endl
          << e.what() << std::endl);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reduce()
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;

    // This method should be called only for locally replicated
    // MultiVectors (! isDistributed ()).
    TEUCHOS_TEST_FOR_EXCEPTION(this->isDistributed (), std::runtime_error,
      "Tpetra::MultiVector::reduce() should only be called with locally "
      "replicated or otherwise not distributed MultiVector objects.");
    RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
    if (comm->getSize() == 1) return;
    RCP<Node> node = MVT::getNode(lclMV_);
    // sum the data across all multivectors
    // need to have separate packed buffers for send and receive
    // if we're packed, we'll just set the receive buffer as our data, the send as a copy
    // if we're not packed, we'll use allocated buffers for both.
    ArrayView<Scalar> target;
    const size_t myStride = MVT::getStride(lclMV_),
                 numCols  = MVT::getNumCols(lclMV_),
                 myLen    = MVT::getNumRows(lclMV_);
    Array<Scalar> sourceBuffer(numCols*myLen), tmparr(0);
    bool packed = isConstantStride() && (myStride == myLen);
    ArrayRCP<Scalar> bufView = node->template viewBufferNonConst<Scalar>(
                                          KokkosClassic::ReadWrite,myStride*(numCols-1)+myLen,
                                          MVT::getValuesNonConst(lclMV_) );
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
      for (size_t j=0; j<numCols; ++j) {
        std::copy(vptr,vptr+myLen,sptr);
        sptr += myLen;
        vptr += myStride;
      }
      target = tmparr();
    }
    // reduce
    reduceAll<int,Scalar> (*comm, REDUCE_SUM, static_cast<int> (numCols*myLen),
                           sourceBuffer.getRawPtr (), target.getRawPtr ());
    if (! packed) {
      // copy tmparr back into view
      const Scalar *sptr = tmparr.getRawPtr();
      ArrayRCP<Scalar> vptr = bufView;
      for (size_t j=0; j<numCols; ++j) {
        std::copy (sptr, sptr+myLen, vptr);
        sptr += myLen;
        vptr += myStride;
      }
    }
    bufView = Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  replaceLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    const LocalOrdinal minLocalIndex = this->getMap()->getMinLocalIndex();
    const LocalOrdinal maxLocalIndex = this->getMap()->getMaxLocalIndex();
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow < minLocalIndex || MyRow > maxLocalIndex,
      std::runtime_error,
      "Tpetra::MultiVector::replaceLocalValue: row index " << MyRow
      << " is invalid.  The range of valid row indices on this process "
      << this->getMap()->getComm()->getRank() << " is [" << minLocalIndex
      << ", " << maxLocalIndex << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::replaceLocalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    const size_t colInd = isConstantStride () ?
      VectorIndex : whichVectors_[VectorIndex];
    lclMV_.replaceLocalValue (static_cast<size_t> (MyRow), colInd, ScalarValue);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  sumIntoLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    const LocalOrdinal minLocalIndex = this->getMap()->getMinLocalIndex();
    const LocalOrdinal maxLocalIndex = this->getMap()->getMaxLocalIndex();
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow < minLocalIndex || MyRow > maxLocalIndex,
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoLocalValue: row index " << MyRow
      << " is invalid.  The range of valid row indices on this process "
      << this->getMap()->getComm()->getRank() << " is [" << minLocalIndex
      << ", " << maxLocalIndex << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoLocalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    const size_t colInd = isConstantStride () ?
      VectorIndex : whichVectors_[VectorIndex];
    lclMV_.sumIntoLocalValue (static_cast<size_t> (MyRow), colInd, ScalarValue);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  replaceGlobalValue (GlobalOrdinal GlobalRow,
                      size_t VectorIndex,
                      const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap()->getLocalElement(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
      std::runtime_error,
      "Tpetra::MultiVector::replaceGlobalValue: global row index " << GlobalRow
      << "is not present on this process " << this->getMap()->getComm()->getRank()
      << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::replaceGlobalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    replaceLocalValue (MyRow, VectorIndex, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  sumIntoGlobalValue (GlobalOrdinal GlobalRow,
                      size_t VectorIndex,
                      const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap()->getLocalElement(GlobalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: global row index " << GlobalRow
      << "is not present on this process " << this->getMap()->getComm()->getRank()
      << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    sumIntoLocalValue (MyRow, VectorIndex, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class T>
  Teuchos::ArrayRCP<T>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  getSubArrayRCP (Teuchos::ArrayRCP<T> arr,
                  size_t j) const
  {
    const size_t stride = MVT::getStride (lclMV_);
    const size_t myLen = getLocalLength ();
    Teuchos::ArrayRCP<T> ret;
    if (isConstantStride()) {
      ret = arr.persistingView (j*stride, myLen);
    }
    else {
      ret = arr.persistingView (whichVectors_[j]*stride, myLen);
    }
    return ret;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  KokkosClassic::MultiVector<Scalar,Node>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMV() const {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TEUCHOS_DEPRECATED
  KokkosClassic::MultiVector<Scalar,Node>&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMVNonConst() {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
  {
    using std::endl;
    std::ostringstream oss;
    oss << Teuchos::typeName (*this) << " {" << endl
        << "  label: \"" << this->getObjectLabel () << "\"" << endl
        << "  numRows: " << getGlobalLength () << endl
        << "  numCols: " << getNumVectors () << endl
        << "  isConstantStride: " << isConstantStride () << endl;
    if (isConstantStride ()) {
      oss << "  columnStride: " << getStride () << endl;
    }
    oss << "}" << endl;
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using std::endl;
    using std::setw;

    // Set default verbosity if applicable.
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();

    if (vl != VERB_NONE) {
      // Don't set the tab level unless we're printing something.
      Teuchos::OSTab tab (out);

      if (myImageID == 0) { // >= VERB_LOW prints description()
        out << this->description() << endl;
      }
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // At verbosity > VERB_LOW, each process prints something.
            out << "Process " << myImageID << ":" << endl;

            Teuchos::OSTab procTab (out);
            // >= VERB_MEDIUM: print the local vector length.
            out << "local length=" << getLocalLength();
            if (vl != VERB_MEDIUM) {
              // >= VERB_HIGH: print isConstantStride() and getStride()
              if (isConstantStride()) {
                out << ", constant stride=" << getStride() << endl;
              }
              else {
                out << ", not constant stride" << endl;
              }
              if (vl == VERB_EXTREME) {
                // VERB_EXTREME: print all the values in the multivector.
                out << "Values:" << endl;
                ArrayRCP<ArrayRCP<const Scalar> > X = this->get2dView();
                for (size_t i = 0; i < getLocalLength(); ++i) {
                  for (size_t j = 0; j < getNumVectors(); ++j) {
                    out << X[j][i];
                    if (j < getNumVectors() - 1) {
                      out << " ";
                    }
                  } // for each column
                  out << endl;
                } // for each row
              } // if vl == VERB_EXTREME
            } // if (vl != VERB_MEDIUM)
            else { // vl == VERB_LOW
              out << endl;
            }
          } // if vl != VERB_LOW
        } // if it is my process' turn to print
        comm->barrier();
      } // for each process in the communicator
    } // if vl != VERB_NONE
  }

#if TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  createViews() const
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  createViewsNonConst (KokkosClassic::ReadWriteOption rwo)
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::releaseViews () const
  {
    // Do nothing in Kokkos::View implementation
  }

#else

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  createViews() const
  {
    Teuchos::RCP<Node> node = this->getMap ()->getNode ();
    if (cview_.is_null () && getLocalLength () > 0) {
      Teuchos::ArrayRCP<const Scalar> buff = MVT::getValues (lclMV_);
      cview_ = node->template viewBuffer<Scalar> (buff.size (), buff);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  createViewsNonConst (KokkosClassic::ReadWriteOption rwo)
  {
    Teuchos::RCP<Node> node = this->getMap ()->getNode ();
    if (ncview_.is_null () && getLocalLength () > 0) {
      Teuchos::ArrayRCP<Scalar> buff = MVT::getValuesNonConst (lclMV_);
      ncview_ = node->template viewBufferNonConst<Scalar> (rwo, buff.size (), buff);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::releaseViews () const
  {
    const int constViewCount = cview_.total_count ();
    const int nonconstViewCount = ncview_.total_count ();

    // This prevents unused variable compiler warnings, in case Tpetra
    // efficiency warnings aren't enabled.
    (void) constViewCount;
    (void) nonconstViewCount;

    // Release the views.
    cview_ = Teuchos::null;
    ncview_ = Teuchos::null;
  }

#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
  {
    replaceMap (newMap);
  }

  template <class ST, class LO, class GO, class NT>
  MultiVector<ST, LO, GO, NT>
  createCopy (const MultiVector<ST, LO, GO, NT>& src) {
    // This always creates a deep copy.
    MultiVector<ST, LO, GO, NT> dst (src, Teuchos::Copy);
    // Returning by value will invoke the copy constructor, so we need
    // to set the result to have view semantics, so that the copy
    // constructor only does a shallow copy.
    dst.setCopyOrView (Teuchos::View);
    return dst;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src)
  {
    const char prefix[] = "Tpetra::MultiVector::assign (called by deep_copy): ";
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getGlobalLength () != src.getGlobalLength (), std::invalid_argument,
      prefix << "this->getGlobalLength() = " << this->getGlobalLength ()
      << " != src.getGlobalLength() = " << src.getGlobalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getNumVectors () != src.getNumVectors(), std::invalid_argument,
      prefix << "this->getNumVectors() = " << this->getNumVectors () << " != "
      "src.getNumVectors() = " << src.getNumVectors () << ".");
    // FIXME (mfh 11 Sep 2014) This may cause deadlock on error.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getMap ().is_null () || src.getMap ().is_null (), std::runtime_error,
      prefix << "You may not call this method if either the source or target "
      "(this) MultiVector's Map is null.  That means that the calling process "
      "is not participating in the operation.")
    // FIXME (mfh 11 Sep 2014) This may cause deadlock on error.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getLocalLength () != src.getLocalLength (), std::invalid_argument,
      prefix << "this->getLocalLength() = " << this->getLocalLength () << " != "
      "src.getLocalLength() = " << src.getLocalLength () << ".");

    if (this != &src) {
      if (this->isConstantStride ()) {
        if (src.isConstantStride ()) {
          MVT::Assign (this->lclMV_, src.lclMV_);
        }
        else { // src is not constant stride
          MVT::Assign (lclMV_, src.lclMV_, src.whichVectors_);
        }
      }
      else { // we have to copy the columns one at a time
        Teuchos::RCP<Node> node = MVT::getNode (this->lclMV_);
        const size_t lclNumRows = this->getLocalLength ();
        const size_t numVecs = this->getNumVectors ();
        for (size_t j = 0; j < numVecs; ++j) {
          node->template copyBuffers<Scalar> (lclNumRows,
                                              src.getSubArrayRCP (MVT::getValues (src.lclMV_), j),
                                              this->getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j));
        }
      }
    }
  }
} // namespace Tpetra

// Include KokkosRefactor partial specialization if enabled
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Tpetra_KokkosRefactor_MultiVector_def.hpp"
#endif

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class MultiVector< SCALAR , LO , GO , NODE >; \
  template MultiVector< SCALAR , LO , GO , NODE > createCopy( const MultiVector< SCALAR , LO , GO , NODE >& src); \



#endif // TPETRA_MULTIVECTOR_DEF_HPP
