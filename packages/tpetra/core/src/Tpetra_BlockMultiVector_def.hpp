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
// ************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKMULTIVECTOR_DEF_HPP
#define TPETRA_BLOCKMULTIVECTOR_DEF_HPP

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_BlockView.hpp"
#include "Teuchos_OrdinalTraits.hpp"

namespace { // anonymous

  /// \brief Get a raw pointer to the (host) data in a
  ///   Tpetra::MultiVector.
  /// \tparam MultiVectorType A specialization of Tpetra::MultiVector.
  ///
  /// \warning This is an implementation detail of Tpetra.  It is not
  ///   part of the public interface and may change or go away at any
  ///   time.
  ///
  /// \note To Tpetra developers: This struct implements the function
  ///   getRawHostPtrFromMultiVector(); see below.  Call that function
  ///   instead.
  template<class MultiVectorType>
  struct RawHostPtrFromMultiVector {
    typedef typename MultiVectorType::impl_scalar_type impl_scalar_type;

    static impl_scalar_type* getRawPtr (MultiVectorType& X) {
      // NOTE (mfh 09 Jun 2016) This does NOT sync to host, or mark
      // host as modified.  This is on purpose, because we don't want
      // the BlockMultiVector sync'd to host unnecessarily.
      // Otherwise, all the MultiVector and BlockMultiVector kernels
      // would run on host instead of device.  See Github Issue #428.
      auto X_view_host = X.template getLocalView<typename MultiVectorType::dual_view_type::t_host::device_type> ();
      impl_scalar_type* X_raw = X_view_host.data ();
      return X_raw;
    }
  };

  /// \brief Get a raw pointer to the (host) data in a
  ///   Tpetra::MultiVector.
  ///
  /// \warning This is an implementation detail of Tpetra.  It is not
  ///   part of the public interface and may change or go away at any
  ///   time.
  ///
  /// \note To Tpetra developers: This function exists to smooth over
  ///   differences between the "classic" and current ("Kokkos
  ///   refactor," circa 2014/5) versions of Tpetra::MultiVector.  It
  ///   also makes the Tpetra::BlockMultiVector
  ///   implementation below a bit easier to read.
  template<class S, class LO, class GO, class N>
  typename Tpetra::MultiVector<S, LO, GO, N>::impl_scalar_type*
  getRawHostPtrFromMultiVector (Tpetra::MultiVector<S, LO, GO, N>& X) {
    typedef Tpetra::MultiVector<S, LO, GO, N> MV;
    return RawHostPtrFromMultiVector<MV>::getRawPtr (X);
  }

} // namespace (anonymous)

namespace Tpetra {

template<class Scalar, class LO, class GO, class Node>
typename BlockMultiVector<Scalar, LO, GO, Node>::mv_type
BlockMultiVector<Scalar, LO, GO, Node>::
getMultiVectorView () const
{
  return mv_;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const BlockMultiVector<Scalar, LO, GO, Node> >
BlockMultiVector<Scalar, LO, GO, Node>::
getBlockMultiVectorFromSrcDistObject (const Tpetra::SrcDistObject& src)
{
  typedef BlockMultiVector<Scalar, LO, GO, Node> BMV;
  const BMV* src_bmv = dynamic_cast<const BMV*> (&src);
  TEUCHOS_TEST_FOR_EXCEPTION(
    src_bmv == nullptr, std::invalid_argument, "Tpetra::"
    "BlockMultiVector: The source object of an Import or Export to a "
    "BlockMultiVector, must also be a BlockMultiVector.");
  return Teuchos::rcp (src_bmv, false); // nonowning RCP
}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const BlockMultiVector<Scalar, LO, GO, Node>& in,
                  const Teuchos::DataAccess copyOrView) :
  dist_object_type (in),
  meshMap_ (in.meshMap_),
  pointMap_ (in.pointMap_),
  mv_ (in.mv_, copyOrView),
  mvData_ (getRawHostPtrFromMultiVector (mv_)),
  blockSize_ (in.blockSize_)
{}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const map_type& meshMap,
                  const LO blockSize,
                  const LO numVecs) :
  dist_object_type (Teuchos::rcp (new map_type (meshMap))), // shallow copy
  meshMap_ (meshMap),
  pointMap_ (makePointMap (meshMap, blockSize)),
  mv_ (Teuchos::rcpFromRef (pointMap_), numVecs), // nonowning RCP is OK, since pointMap_ won't go away
  mvData_ (getRawHostPtrFromMultiVector (mv_)),
  blockSize_ (blockSize)
{}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const map_type& meshMap,
                  const map_type& pointMap,
                  const LO blockSize,
                  const LO numVecs) :
  dist_object_type (Teuchos::rcp (new map_type (meshMap))), // shallow copy
  meshMap_ (meshMap),
  pointMap_ (pointMap),
  mv_ (Teuchos::rcpFromRef (pointMap_), numVecs),
  mvData_ (getRawHostPtrFromMultiVector (mv_)),
  blockSize_ (blockSize)
{}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const mv_type& X_mv,
                  const map_type& meshMap,
                  const LO blockSize) :
  dist_object_type (Teuchos::rcp (new map_type (meshMap))), // shallow copy
  meshMap_ (meshMap),
  mvData_ (nullptr), // just for now
  blockSize_ (blockSize)
{
  using Teuchos::RCP;

  if (X_mv.getCopyOrView () == Teuchos::View) {
    // The input MultiVector has view semantics, so assignment just
    // does a shallow copy.
    mv_ = X_mv;
  }
  else if (X_mv.getCopyOrView () == Teuchos::Copy) {
    // The input MultiVector has copy semantics.  We can't change
    // that, but we can make a view of the input MultiVector and
    // change the view to have view semantics.  (That sounds silly;
    // shouldn't views always have view semantics? but remember that
    // "view semantics" just governs the default behavior of the copy
    // constructor and assignment operator.)
    RCP<const mv_type> X_view_const;
    const size_t numCols = X_mv.getNumVectors ();
    if (numCols == 0) {
      Teuchos::Array<size_t> cols (0); // view will have zero columns
      X_view_const = X_mv.subView (cols ());
    } else { // Range1D is an inclusive range
      X_view_const = X_mv.subView (Teuchos::Range1D (0, numCols-1));
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      X_view_const.is_null (), std::logic_error, "Tpetra::"
      "BlockMultiVector constructor: X_mv.subView(...) returned null.  This "
      "should never happen.  Please report this bug to the Tpetra developers.");

    // It's perfectly OK to cast away const here.  Those view methods
    // should be marked const anyway, because views can't change the
    // allocation (just the entries themselves).
    RCP<mv_type> X_view = Teuchos::rcp_const_cast<mv_type> (X_view_const);
    TEUCHOS_TEST_FOR_EXCEPTION(
      X_view->getCopyOrView () != Teuchos::View, std::logic_error, "Tpetra::"
      "BlockMultiVector constructor: We just set a MultiVector "
      "to have view semantics, but it claims that it doesn't have view "
      "semantics.  This should never happen.  "
      "Please report this bug to the Tpetra developers.");
    mv_ = *X_view; // MultiVector::operator= does a shallow copy here
  }

  // At this point, mv_ has been assigned, so we can ignore X_mv.
  Teuchos::RCP<const map_type> pointMap = mv_.getMap ();
  if (! pointMap.is_null ()) {
    pointMap_ = *pointMap; // Map::operator= also does a shallow copy
  }
  mvData_ = getRawHostPtrFromMultiVector (mv_);
}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                  const map_type& newMeshMap,
                  const map_type& newPointMap,
                  const size_t offset) :
  dist_object_type (Teuchos::rcp (new map_type (newMeshMap))), // shallow copy
  meshMap_ (newMeshMap),
  pointMap_ (newPointMap),
  mv_ (X.mv_, newPointMap, offset * X.getBlockSize ()), // MV "offset view" constructor
  mvData_ (getRawHostPtrFromMultiVector (mv_)),
  blockSize_ (X.getBlockSize ())
{}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                  const map_type& newMeshMap,
                  const size_t offset) :
  dist_object_type (Teuchos::rcp (new map_type (newMeshMap))), // shallow copy
  meshMap_ (newMeshMap),
  pointMap_ (makePointMap (newMeshMap, X.getBlockSize ())),
  mv_ (X.mv_, pointMap_, offset * X.getBlockSize ()), // MV "offset view" constructor
  mvData_ (getRawHostPtrFromMultiVector (mv_)),
  blockSize_ (X.getBlockSize ())
{}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector () :
  dist_object_type (Teuchos::null),
  mvData_ (nullptr),
  blockSize_ (0)
{}

template<class Scalar, class LO, class GO, class Node>
typename BlockMultiVector<Scalar, LO, GO, Node>::map_type
BlockMultiVector<Scalar, LO, GO, Node>::
makePointMap (const map_type& meshMap, const LO blockSize)
{
  typedef Tpetra::global_size_t GST;
  typedef typename Teuchos::ArrayView<const GO>::size_type size_type;

  const GST gblNumMeshMapInds =
    static_cast<GST> (meshMap.getGlobalNumElements ());
  const size_t lclNumMeshMapIndices =
    static_cast<size_t> (meshMap.getNodeNumElements ());
  const GST gblNumPointMapInds =
    gblNumMeshMapInds * static_cast<GST> (blockSize);
  const size_t lclNumPointMapInds =
    lclNumMeshMapIndices * static_cast<size_t> (blockSize);
  const GO indexBase = meshMap.getIndexBase ();

  if (meshMap.isContiguous ()) {
    return map_type (gblNumPointMapInds, lclNumPointMapInds, indexBase,
                     meshMap.getComm ());
  }
  else {
    // "Hilbert's Hotel" trick: multiply each process' GIDs by
    // blockSize, and fill in.  That ensures correctness even if the
    // mesh Map is overlapping.
    Teuchos::ArrayView<const GO> lclMeshGblInds = meshMap.getNodeElementList ();
    const size_type lclNumMeshGblInds = lclMeshGblInds.size ();
    Teuchos::Array<GO> lclPointGblInds (lclNumPointMapInds);
    for (size_type g = 0; g < lclNumMeshGblInds; ++g) {
      const GO meshGid = lclMeshGblInds[g];
      const GO pointGidStart = indexBase +
        (meshGid - indexBase) * static_cast<GO> (blockSize);
      const size_type offset = g * static_cast<size_type> (blockSize);
      for (LO k = 0; k < blockSize; ++k) {
        const GO pointGid = pointGidStart + static_cast<GO> (k);
        lclPointGblInds[offset + static_cast<size_type> (k)] = pointGid;
      }
    }
    return map_type (gblNumPointMapInds, lclPointGblInds (), indexBase,
                     meshMap.getComm ());
  }
}


template<class Scalar, class LO, class GO, class Node>
void
BlockMultiVector<Scalar, LO, GO, Node>::
replaceLocalValuesImpl (const LO localRowIndex,
                        const LO colIndex,
                        const Scalar vals[]) const
{
  auto X_dst = getLocalBlock (localRowIndex, colIndex);
  typename const_little_vec_type::HostMirror::const_type X_src (reinterpret_cast<const impl_scalar_type*> (vals),
                                                                getBlockSize ());
  Kokkos::deep_copy (X_dst, X_src);
}


template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
replaceLocalValues (const LO localRowIndex,
                    const LO colIndex,
                    const Scalar vals[]) const
{
  if (! meshMap_.isNodeLocalElement (localRowIndex)) {
    return false;
  } else {
    replaceLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
replaceGlobalValues (const GO globalRowIndex,
                     const LO colIndex,
                     const Scalar vals[]) const
{
  const LO localRowIndex = meshMap_.getLocalElement (globalRowIndex);
  if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid ()) {
    return false;
  } else {
    replaceLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
void
BlockMultiVector<Scalar, LO, GO, Node>::
sumIntoLocalValuesImpl (const LO localRowIndex,
                        const LO colIndex,
                        const Scalar vals[]) const
{
  auto X_dst = getLocalBlock (localRowIndex, colIndex);
  typename const_little_vec_type::HostMirror::const_type X_src (reinterpret_cast<const impl_scalar_type*> (vals),
                                                                getBlockSize ());
  AXPY (static_cast<impl_scalar_type> (STS::one ()), X_src, X_dst);
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
sumIntoLocalValues (const LO localRowIndex,
                    const LO colIndex,
                    const Scalar vals[]) const
{
  if (! meshMap_.isNodeLocalElement (localRowIndex)) {
    return false;
  } else {
    sumIntoLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
sumIntoGlobalValues (const GO globalRowIndex,
                     const LO colIndex,
                     const Scalar vals[]) const
{
  const LO localRowIndex = meshMap_.getLocalElement (globalRowIndex);
  if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid ()) {
    return false;
  } else {
    sumIntoLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
getLocalRowView (const LO localRowIndex, const LO colIndex, Scalar*& vals) const
{
  if (! meshMap_.isNodeLocalElement (localRowIndex)) {
    return false;
  } else {
    auto X_ij = getLocalBlock (localRowIndex, colIndex);
    vals = reinterpret_cast<Scalar*> (X_ij.data ());
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
getGlobalRowView (const GO globalRowIndex, const LO colIndex, Scalar*& vals) const
{
  const LO localRowIndex = meshMap_.getLocalElement (globalRowIndex);
  if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid ()) {
    return false;
  } else {
    auto X_ij = getLocalBlock (localRowIndex, colIndex);
    vals = reinterpret_cast<Scalar*> (X_ij.data ());
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
typename BlockMultiVector<Scalar, LO, GO, Node>::little_vec_type::HostMirror
BlockMultiVector<Scalar, LO, GO, Node>::
getLocalBlock (const LO localRowIndex,
               const LO colIndex) const
{
  // NOTE (mfh 07 Jul 2016) It should be correct to add the
  // commented-out test below.  However, I've conservatively commented
  // it out, since users might not realize that they need to have
  // things sync'd correctly.

// #ifdef HAVE_TPETRA_DEBUG
//   TEUCHOS_TEST_FOR_EXCEPTION
//     (mv_.need_sync_host (), std::runtime_error,
//      "Tpetra::BlockMultiVector::getLocalBlock: This method "
//      "accesses host data, but the object is not in sync on host." );
// #endif // HAVE_TPETRA_DEBUG

  if (! isValidLocalMeshIndex (localRowIndex)) {
    return typename little_vec_type::HostMirror ();
  } else {
    const size_t blockSize = getBlockSize ();
    const size_t offset = colIndex * this->getStrideY () +
      localRowIndex * blockSize;
    impl_scalar_type* blockRaw = this->getRawPtr () + offset;
    return typename little_vec_type::HostMirror (blockRaw, blockSize);
  }
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const typename BlockMultiVector<Scalar, LO, GO, Node>::mv_type>
BlockMultiVector<Scalar, LO, GO, Node>::
getMultiVectorFromSrcDistObject (const Tpetra::SrcDistObject& src)
{
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  // The source object of an Import or Export must be a
  // BlockMultiVector or MultiVector (a Vector is a MultiVector).  Try
  // them in that order; one must succeed.  Note that the target of
  // the Import or Export calls checkSizes in DistObject's doTransfer.
  typedef BlockMultiVector<Scalar, LO, GO, Node> this_type;
  const this_type* srcBlkVec = dynamic_cast<const this_type*> (&src);
  if (srcBlkVec == nullptr) {
    const mv_type* srcMultiVec = dynamic_cast<const mv_type*> (&src);
    if (srcMultiVec == nullptr) {
      // FIXME (mfh 05 May 2014) Tpetra::MultiVector's operator=
      // currently does a shallow copy by default.  This is why we
      // return by RCP, rather than by value.
      return rcp (new mv_type ());
    } else { // src is a MultiVector
      return rcp (srcMultiVec, false); // nonowning RCP
    }
  } else { // src is a BlockMultiVector
    return rcpFromRef (srcBlkVec->mv_); // nonowning RCP
  }
}

template<class Scalar, class LO, class GO, class Node>
bool BlockMultiVector<Scalar, LO, GO, Node>::
checkSizes (const Tpetra::SrcDistObject& src)
{
  return ! getMultiVectorFromSrcDistObject (src).is_null ();
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
copyAndPermute
(const SrcDistObject& src,
 const size_t numSameIDs,
 const Kokkos::DualView<const local_ordinal_type*,
 buffer_device_type>& permuteToLIDs,
 const Kokkos::DualView<const local_ordinal_type*,
 buffer_device_type>& permuteFromLIDs)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error,
     "Tpetra::BlockMultiVector::copyAndPermute: Do NOT use this "
     "instead, create a point importer using makePointMap function.");
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
packAndPrepare
(const SrcDistObject& src,
 const Kokkos::DualView<const local_ordinal_type*,
 buffer_device_type>& exportLIDs,
 Kokkos::DualView<packet_type*,
 buffer_device_type>& exports,
 Kokkos::DualView<size_t*,
 buffer_device_type> numPacketsPerLID,
 size_t& constantNumPackets,
 Distributor& distor)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error,
     "Tpetra::BlockMultiVector::copyAndPermute: Do NOT use this; "
     "instead, create a point importer using makePointMap function.");
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
unpackAndCombine
(const Kokkos::DualView<const local_ordinal_type*,
 buffer_device_type>& importLIDs,
 Kokkos::DualView<packet_type*,
 buffer_device_type> imports,
 Kokkos::DualView<size_t*,
 buffer_device_type> numPacketsPerLID,
 const size_t constantNumPackets,
 Distributor& distor,
 const CombineMode combineMode)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error,
     "Tpetra::BlockMultiVector::copyAndPermute: Do NOT use this; "
     "instead, create a point importer using makePointMap function.");
}

template<class Scalar, class LO, class GO, class Node>
bool BlockMultiVector<Scalar, LO, GO, Node>::
isValidLocalMeshIndex (const LO meshLocalIndex) const
{
  return meshLocalIndex != Teuchos::OrdinalTraits<LO>::invalid () &&
    meshMap_.isNodeLocalElement (meshLocalIndex);
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
putScalar (const Scalar& val)
{
  mv_.putScalar (val);
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
scale (const Scalar& val)
{
  mv_.scale (val);
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
update (const Scalar& alpha,
        const BlockMultiVector<Scalar, LO, GO, Node>& X,
        const Scalar& beta)
{
  mv_.update (alpha, X.mv_, beta);
}

namespace Impl {
// Y := alpha * D * X
template <typename Scalar, typename ViewY, typename ViewD, typename ViewX>
struct BlockWiseMultiply {
  typedef typename ViewD::size_type Size;

private:
  typedef typename ViewD::device_type Device;
  typedef typename ViewD::non_const_value_type ImplScalar;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> Unmanaged;

  template <typename View>
  using UnmanagedView = Kokkos::View<typename View::data_type, typename View::array_layout,
                                     typename View::device_type, Unmanaged>;
  typedef UnmanagedView<ViewY> UnMViewY;
  typedef UnmanagedView<ViewD> UnMViewD;
  typedef UnmanagedView<ViewX> UnMViewX;

  const Size block_size_;
  Scalar alpha_;
  UnMViewY Y_;
  UnMViewD D_;
  UnMViewX X_;

public:
  BlockWiseMultiply (const Size block_size, const Scalar& alpha,
                     const ViewY& Y, const ViewD& D, const ViewX& X)
    : block_size_(block_size), alpha_(alpha), Y_(Y), D_(D), X_(X)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const Size k) const {
    const auto zero = Kokkos::Details::ArithTraits<Scalar>::zero();
    auto D_curBlk = Kokkos::subview(D_, k, Kokkos::ALL (), Kokkos::ALL ());
    const auto num_vecs = X_.extent(1);
    for (Size i = 0; i < num_vecs; ++i) {
      Kokkos::pair<Size, Size> kslice(k*block_size_, (k+1)*block_size_);
      auto X_curBlk = Kokkos::subview(X_, kslice, i);
      auto Y_curBlk = Kokkos::subview(Y_, kslice, i);
      // Y_curBlk := alpha * D_curBlk * X_curBlk.
      // Recall that GEMV does an update (+=) of the last argument.
      Tpetra::FILL(Y_curBlk, zero);
      Tpetra::GEMV(alpha_, D_curBlk, X_curBlk, Y_curBlk);
    }
  }
};

template <typename Scalar, typename ViewY, typename ViewD, typename ViewX>
inline BlockWiseMultiply<Scalar, ViewY, ViewD, ViewX>
createBlockWiseMultiply (const int block_size, const Scalar& alpha,
                         const ViewY& Y, const ViewD& D, const ViewX& X) {
  return BlockWiseMultiply<Scalar, ViewY, ViewD, ViewX>(block_size, alpha, Y, D, X);
}

template <typename ViewY,
          typename Scalar,
          typename ViewD,
          typename ViewZ,
          typename LO = typename ViewY::size_type>
class BlockJacobiUpdate {
private:
  typedef typename ViewD::device_type Device;
  typedef typename ViewD::non_const_value_type ImplScalar;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> Unmanaged;

  template <typename ViewType>
  using UnmanagedView = Kokkos::View<typename ViewType::data_type,
                                     typename ViewType::array_layout,
                                     typename ViewType::device_type,
                                     Unmanaged>;
  typedef UnmanagedView<ViewY> UnMViewY;
  typedef UnmanagedView<ViewD> UnMViewD;
  typedef UnmanagedView<ViewZ> UnMViewZ;

  const LO blockSize_;
  UnMViewY Y_;
  const Scalar alpha_;
  UnMViewD D_;
  UnMViewZ Z_;
  const Scalar beta_;

public:
  BlockJacobiUpdate (const ViewY& Y,
                     const Scalar& alpha,
                     const ViewD& D,
                     const ViewZ& Z,
                     const Scalar& beta) :
    blockSize_ (D.extent (1)),
    // numVecs_ (static_cast<int> (ViewY::rank) == 1 ? static_cast<size_type> (1) : static_cast<size_type> (Y_.extent (1))),
    Y_ (Y),
    alpha_ (alpha),
    D_ (D),
    Z_ (Z),
    beta_ (beta)
  {
    static_assert (static_cast<int> (ViewY::rank) == 1,
                   "Y must have rank 1.");
    static_assert (static_cast<int> (ViewD::rank) == 3, "D must have rank 3.");
    static_assert (static_cast<int> (ViewZ::rank) == 1,
                   "Z must have rank 1.");
    // static_assert (static_cast<int> (ViewZ::rank) ==
    //                static_cast<int> (ViewY::rank),
    //                "Z must have the same rank as Y.");
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const LO& k) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::pair<LO, LO> range_type;
    typedef Kokkos::Details::ArithTraits<Scalar> KAT;

    // We only have to implement the alpha != 0 case.

    auto D_curBlk = subview (D_, k, ALL (), ALL ());
    const range_type kslice (k*blockSize_, (k+1)*blockSize_);

    // Z.update (STS::one (), X, -STS::one ()); // assume this is done

    auto Z_curBlk = subview (Z_, kslice);
    auto Y_curBlk = subview (Y_, kslice);
    // Y_curBlk := beta * Y_curBlk + alpha * D_curBlk * Z_curBlk
    if (beta_ == KAT::zero ()) {
      Tpetra::FILL (Y_curBlk, KAT::zero ());
    }
    else if (beta_ != KAT::one ()) {
      Tpetra::SCAL (beta_, Y_curBlk);
    }
    Tpetra::GEMV (alpha_, D_curBlk, Z_curBlk, Y_curBlk);
  }
};

template<class ViewY,
         class Scalar,
         class ViewD,
         class ViewZ,
         class LO = typename ViewD::size_type>
void
blockJacobiUpdate (const ViewY& Y,
                   const Scalar& alpha,
                   const ViewD& D,
                   const ViewZ& Z,
                   const Scalar& beta)
{
  static_assert (Kokkos::Impl::is_view<ViewY>::value, "Y must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<ViewD>::value, "D must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<ViewZ>::value, "Z must be a Kokkos::View.");
  static_assert (static_cast<int> (ViewY::rank) == static_cast<int> (ViewZ::rank),
                 "Y and Z must have the same rank.");
  static_assert (static_cast<int> (ViewD::rank) == 3, "D must have rank 3.");

  const auto lclNumMeshRows = D.extent (0);

#ifdef HAVE_TPETRA_DEBUG
  // D.extent(0) is the (local) number of mesh rows.
  // D.extent(1) is the block size.  Thus, their product should be
  // the local number of point rows, that is, the number of rows in Y.
  const auto blkSize = D.extent (1);
  const auto lclNumPtRows = lclNumMeshRows * blkSize;
  TEUCHOS_TEST_FOR_EXCEPTION
    (Y.extent (0) != lclNumPtRows, std::invalid_argument,
     "blockJacobiUpdate: Y.extent(0) = " << Y.extent (0) << " != "
     "D.extent(0)*D.extent(1) = " << lclNumMeshRows << " * " << blkSize
     << " = " << lclNumPtRows << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (Y.extent (0) != Z.extent (0), std::invalid_argument,
     "blockJacobiUpdate: Y.extent(0) = " << Y.extent (0) << " != "
     "Z.extent(0) = " << Z.extent (0) << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (Y.extent (1) != Z.extent (1), std::invalid_argument,
     "blockJacobiUpdate: Y.extent(1) = " << Y.extent (1) << " != "
     "Z.extent(1) = " << Z.extent (1) << ".");
#endif // HAVE_TPETRA_DEBUG

  BlockJacobiUpdate<ViewY, Scalar, ViewD, ViewZ, LO> functor (Y, alpha, D, Z, beta);
  typedef Kokkos::RangePolicy<typename ViewY::execution_space, LO> range_type;
  // lclNumMeshRows must fit in LO, else the Map would not be correct.
  range_type range (0, static_cast<LO> (lclNumMeshRows));
  Kokkos::parallel_for (range, functor);
}

} // namespace Impl

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
blockWiseMultiply (const Scalar& alpha,
                   const Kokkos::View<const impl_scalar_type***,
                     device_type, Kokkos::MemoryUnmanaged>& D,
                   const BlockMultiVector<Scalar, LO, GO, Node>& X)
{
  using Kokkos::ALL;
  typedef typename device_type::execution_space execution_space;
  typedef typename device_type::memory_space memory_space;
  const LO lclNumMeshRows = meshMap_.getNodeNumElements ();

  if (alpha == STS::zero ()) {
    this->putScalar (STS::zero ());
  }
  else { // alpha != 0
    const LO blockSize = this->getBlockSize ();
    const impl_scalar_type alphaImpl = static_cast<impl_scalar_type> (alpha);
    auto X_lcl = X.mv_.template getLocalView<memory_space> ();
    auto Y_lcl = this->mv_.template getLocalView<memory_space> ();
    auto bwm = Impl::createBlockWiseMultiply (blockSize, alphaImpl, Y_lcl, D, X_lcl);

    // Use an explicit RangePolicy with the desired execution space.
    // Otherwise, if you just use a number, it runs on the default
    // execution space.  For example, if the default execution space
    // is Cuda but the current execution space is Serial, using just a
    // number would incorrectly run with Cuda.
    Kokkos::RangePolicy<execution_space, LO> range (0, lclNumMeshRows);
    Kokkos::parallel_for (range, bwm);
  }
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
blockJacobiUpdate (const Scalar& alpha,
                   const Kokkos::View<const impl_scalar_type***,
                     device_type, Kokkos::MemoryUnmanaged>& D,
                   const BlockMultiVector<Scalar, LO, GO, Node>& X,
                   BlockMultiVector<Scalar, LO, GO, Node>& Z,
                   const Scalar& beta)
{
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef typename device_type::memory_space memory_space;
  typedef impl_scalar_type IST;

  const IST alphaImpl = static_cast<IST> (alpha);
  const IST betaImpl = static_cast<IST> (beta);
  const LO numVecs = mv_.getNumVectors ();

  auto X_lcl = X.mv_.template getLocalView<memory_space> ();
  auto Y_lcl = this->mv_.template getLocalView<memory_space> ();
  auto Z_lcl = Z.mv_.template getLocalView<memory_space> ();

  if (alpha == STS::zero ()) { // Y := beta * Y
    this->scale (beta);
  }
  else { // alpha != 0
    Z.update (STS::one (), X, -STS::one ());
    for (LO j = 0; j < numVecs; ++j) {
      auto X_lcl_j = subview (X_lcl, ALL (), j);
      auto Y_lcl_j = subview (Y_lcl, ALL (), j);
      auto Z_lcl_j = subview (Z_lcl, ALL (), j);
      Impl::blockJacobiUpdate (Y_lcl_j, alphaImpl, D, Z_lcl_j, betaImpl);
    }
  }
}

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_BLOCKMULTIVECTOR_INSTANT(S,LO,GO,NODE) \
  template class BlockMultiVector< S, LO, GO, NODE >;

#endif // TPETRA_BLOCKMULTIVECTOR_DEF_HPP
