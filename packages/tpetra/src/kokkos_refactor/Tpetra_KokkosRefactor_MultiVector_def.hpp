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

#ifndef TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP
#define TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP

#include <Kokkos_DefaultArithmetic.hpp>
#include <Kokkos_NodeTrace.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_as.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Vector.hpp>

#include <Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_KokkosRefactor_MultiVector_decl.hpp"
#endif

#include <KokkosCompat_View.hpp>
#include <Kokkos_MV.hpp>

namespace Tpetra {


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  vectorIndexOutOfRange(size_t VectorIndex) const {
    return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               size_t NumVectors,
               bool zeroOut) : /* default is true */
    DO (map),
    lclMV_ (map->getNode ())
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
      //ArrayRCP<Scalar> data = node->template allocBuffer<Scalar>(myLen*NumVectors);

      // Allocate a DualView from new Kokkos, wrap its device data into an ArrayRCP
      view_ = view_type("MV::dual_view",myLen,NumVectors);
      ArrayRCP<Scalar> data = Kokkos::Compat::persistingView(view_.d_view);
      MVT::initializeValues(lclMV_,myLen,NumVectors,data,myLen);

      // First touch was done by view allocation
      /*if (zeroOut) {
        // MVT uses the Kokkos Node's parallel_for in this case, for
        // first-touch allocation (across rows).
        MVT::Init(lclMV_, Teuchos::ScalarTraits<Scalar>::zero());
      }*/
    }
    else {
      view_ = view_type("MV::dual_view",0,NumVectors);
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& source) :
    DO (source),
    lclMV_ (MVT::getNode (source.lclMV_)),
    view_ (source.view_),
    whichVectors_(source.whichVectors_)
  {
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
/*    ArrayRCP<Scalar> data = (myLen > 0) ?
      node->template allocBuffer<Scalar> (myLen * numVecs) :
      Teuchos::null;*/
    ArrayRCP<Scalar> data = (myLen > 0) ?
      Kokkos::Compat::persistingView(view_.d_view) :
      Teuchos::null;
    const size_t stride = (myLen > 0) ? myLen : size_t (0);
    // This just sets the dimensions, pointer, and stride of lclMV_.
    MVT::initializeValues (lclMV_, myLen, numVecs, data, stride);

    // Refactor: don't copy we use view semantics
    //
    // mfh 27 Nov 2013: Hm, it will break a lot of things if
    // MultiVector uses view semantics.  I prefer view semantics
    // myself, but it definitely won't preserve backwards
    // compatibility.

    // This actually copies the data.  It uses the Node's
    // parallel_for to copy, which should ensure first-touch
    // allocation on systems that support it.
    /*if (source.isConstantStride ()) {
      MVT::Assign (lclMV_, source.lclMV_);
    }
    else {
      MVT::Assign (lclMV_, source.lclMV_, source.whichVectors_);
    }*/
  }



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const view_type view) :
    DO (map),
    lclMV_ (map->getNode ()),
    view_(view)
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    // getting stride of view: if second dimension is 0 stride might be 0, so take view_dimension instead
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = view_.dimension_1() > 1?stride[1]:view_.dimension_0();
    const size_t NumVectors = view_.dimension_1();

    const char tfecfFuncName[] = "MultiVector(view,LDA,NumVector)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive, but you specified NumVectors = "
      << NumVectors << ".");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");

    // This is a shallow copy into the KokkosClassic::MultiVector.
    MVT::initializeValues(lclMV_,myLen,NumVectors,Kokkos::Compat::persistingView(view_.d_view),LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const view_type view,
               const Teuchos::ArrayView<const size_t>& whichVectors) :
    DO (map),
    lclMV_ (map->getNode ()),
    view_(view),
    whichVectors_ (whichVectors.begin (), whichVectors.end ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    // getting stride of view: if second dimension is 0 stride might be 0, so take view_dimension instead
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = view_.dimension_1() > 1?stride[1]:view_.dimension_0();
    size_t NumVectors = view_.dimension_1();

    const char tfecfFuncName[] = "MultiVector(view,LDA,NumVector)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(NumVectors < 1, std::invalid_argument,
      ": NumVectors must be strictly positive, but you specified NumVectors = "
      << NumVectors << ".");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");

    if (whichVectors.size () == 1) {
      // shift data so that desired vector is vector 0
      // kill whichVectors_; we are constant stride
      view_type tmpview = Kokkos::subview<view_type>(view_,Kokkos::ALL(),whichVectors[0]);
      view_ = tmpview;
      NumVectors = 1;
      whichVectors_.clear ();
    }
    // This is a shallow copy into the KokkosClassic::MultiVector.
    // KokkosClassic::MultiVector doesn't know about whichVectors_.
    MVT::initializeValues (lclMV_, myLen, NumVectors, Kokkos::Compat::persistingView (view_.d_view), LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const Scalar>& data,
               const size_t LDA,
               const size_t numVecs) :
    DO (map),
    lclMV_ (map->getNode ())
  {
    // Deep copy constructor, constant stride (NO whichVectors_).
    // There is no need for a deep copy constructor with nonconstant stride.

    const char tfecfFuncName[] = "MultiVector(map,data,LDA,numVecs)";
    const size_t numRows = this->getLocalLength ();
    view_ = view_type("MV::view_",numRows,numVecs);
    for(size_t i = 0; i<numRows;i++)
      for(size_t j = 0; j<numVecs;j++)
        view_.h_view(i,j) = data[j*LDA+i];
    view_.template modify<typename view_type::host_mirror_device_type>();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < numRows, std::runtime_error,
      ": LDA = " << LDA << " < numRows = " << numRows << ".");
    MVT::initializeValues (lclMV_, numRows, numVecs, Kokkos::Compat::persistingView(view_.d_view), LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const ArrayView<const Scalar> >& ArrayOfPtrs,
               const size_t NumVectors) :
    DO (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;

    const char tfecfFuncName[] = "MultiVector(ArrayOfPtrs)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      NumVectors < 1 || NumVectors != as<size_t>(ArrayOfPtrs.size()),
      std::runtime_error,
      ": ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
    const size_t myLen = getLocalLength();
    RCP<Node> node = MVT::getNode(lclMV_);
    view_ = view_type("MV::view_",myLen,NumVectors);

    // TODO: write a functor and use parallel_for

    for(size_t i = 0; i<myLen;i++)
      for(size_t j = 0; j<NumVectors;j++)
        view_.h_view(i,j) = ArrayOfPtrs[j][i];
      view_.template modify<typename view_type::t_host::device_type>();
    MVT::initializeValues(lclMV_,myLen,NumVectors,Kokkos::Compat::persistingView(view_.d_view),myLen);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::~MultiVector() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::isConstantStride() const {
    return whichVectors_.empty();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getLocalLength() const {
    return this->getMap()->getNodeNumElements();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getGlobalLength() const {
    return this->getMap()->getGlobalNumElements();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getStride() const {
    if (isConstantStride()) {
      return MVT::getStride(lclMV_);
    }
    return 0;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  constantNumberOfPackets () const {
    return this->getNumVectors ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  copyAndPermuteNew (
    const SrcDistObject& sourceObj,
    size_t numSameIDs,
    const Kokkos::View<const LocalOrdinal*, device_type> &permuteToLIDs,
    const Kokkos::View<const LocalOrdinal*, device_type> &permuteFromLIDs)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
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

    // TODO (mfh 15 Sep 2013) When we replace
    // KokkosClassic::MultiVector with a Kokkos::View, there are two
    // ways to copy the data:
    //
    // 1. Get a (sub)view of each column and call deep_copy on that.
    // 2. Write a custom kernel to copy the data.
    //
    // The first is easier, but the second might be more performant in
    // case we decide to use layouts other than LayoutLeft.  It might
    // even make sense to hide whichVectors_ in an entirely new layout
    // for Kokkos Views.

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
      KokkosRefactor::Details::permute_array_multi_column(
        getKokkosView(),
        sourceMV.getKokkosView(),
        permuteToLIDs,
        permuteFromLIDs,
        numCols);
    }
    else {
      KokkosRefactor::Details::permute_array_multi_column_variable_stride(
        getKokkosView(),
        sourceMV.getKokkosView(),
        permuteToLIDs,
        permuteFromLIDs,
        getKokkosViewDeepCopy<device_type> (whichVectors_ ()),
        getKokkosViewDeepCopy<device_type> (sourceMV.whichVectors_ ()),
        numCols);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  packAndPrepareNew (
    const SrcDistObject& sourceObj,
    const Kokkos::View<const LocalOrdinal*, device_type> &exportLIDs,
    Kokkos::View<Scalar*, device_type> &exports,
    const Kokkos::View<size_t*, device_type> &numExportPacketsPerLID,
    size_t& constantNumPackets,
    Distributor & /* distor */ )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    typedef Array<size_t>::size_type size_type;

    // If we have no exports, there is nothing to do
    if (exportLIDs.size () == 0) {
      return;
    }

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

    // FIXME (mfh 15 Sep 2013) Would it make sense to rethink the
    // packing scheme in the above comment?  The data going to a
    // particular process must be contiguous, of course, but those
    // data could include entries from multiple LIDs.  DistObject just
    // needs to know how to index into that data.  Kokkos is good at
    // decoupling storage intent from data layout choice.

    const size_t numCols = sourceMV.getNumVectors ();

    // This spares us from needing to fill numExportPacketsPerLID.
    // Setting constantNumPackets to a nonzero value signals that
    // all packets have the same number of entries.
    constantNumPackets = numCols;

    const size_t numExportLIDs = exportLIDs.size ();
    const size_t newExportsSize = numCols * numExportLIDs;
    if (exports.size () != newExportsSize) {
      Kokkos::Compat::realloc (exports, newExportsSize);
    }

    if (numCols == 1) { // special case for one column only
      // MultiVector always represents a single column with constant
      // stride, but it doesn't hurt to implement both cases anyway.
      //
      // ETP:  I'm not sure I agree with the above statement.  Can't a single-
      // column multivector be a subview of another multi-vector, in which case
      // sourceMV.whichVectors_[0] != 0 ?  I think we have to handle that case
      // separately.
      if (sourceMV.isConstantStride ()) {
        KokkosRefactor::Details::pack_array_single_column(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          0);
      }
      else {
        KokkosRefactor::Details::pack_array_single_column(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          sourceMV.whichVectors_[0]);
      }
    }
    else { // the source MultiVector has multiple columns
      if (sourceMV.isConstantStride ()) {
        KokkosRefactor::Details::pack_array_multi_column(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          numCols);
      }
      else {
        KokkosRefactor::Details::pack_array_multi_column_variable_stride(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          getKokkosViewDeepCopy<device_type> (sourceMV.whichVectors_ ()),
          numCols);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  unpackAndCombineNew (
    const Kokkos::View<const LocalOrdinal*, device_type> &importLIDs,
    const Kokkos::View<const Scalar*, device_type> &imports,
    const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
    size_t constantNumPackets,
    Distributor & /* distor */,
    CombineMode CM)
  {
    using Teuchos::ArrayView;
    using Teuchos::as;
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
      as<size_t> (imports.size()) != getNumVectors()*importLIDs.size(),
      std::runtime_error,
      ": 'imports' buffer size must be consistent with the amount of data to "
      "be sent.  " << std::endl << "imports.size() = " << imports.size()
      << " != getNumVectors()*importLIDs.size() = " << getNumVectors() << "*"
      << importLIDs.size() << " = " << getNumVectors() * importLIDs.size()
      << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (constantNumPackets) == as<size_t> (0), std::runtime_error,
      ": constantNumPackets input argument must be nonzero.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (numVecs) != as<size_t> (constantNumPackets),
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
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::InsertOp<Scalar>(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<device_type>(whichVectors_ ()),
            KokkosRefactor::Details::InsertOp<Scalar>(),
            numVecs);
        }
      }
      else if (CM == ADD) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::AddOp<Scalar>(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<device_type>(whichVectors_ ()),
            KokkosRefactor::Details::AddOp<Scalar>(),
            numVecs);
        }
      }
      else if (CM == ABSMAX) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::AbsMaxOp<Scalar>(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<device_type>(whichVectors_ ()),
            KokkosRefactor::Details::AbsMaxOp<Scalar>(),
            numVecs);
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  inline size_t MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getNumVectors() const {
    if (isConstantStride()) {
      return MVT::getNumCols(lclMV_);
    }
    else {
      return whichVectors_.size();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A,
       const Teuchos::ArrayView<Scalar> &dots) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::arcp_const_cast;

    const char tfecfFuncName[] = "dot()";
    //const size_t myLen   = getLocalLength();
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
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(dots.size()) != numVecs, std::runtime_error,
        ": dots.size() must be as large as the number of vectors in *this and A.");
    /*
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
    }*/
    if (isConstantStride() && A.isConstantStride() ) {
      Kokkos::MV_Dot(&dots[0],view_.d_view,A.view_.d_view,getLocalLength());
    } else {
      for (size_t k = 0; k < numVecs; ++k) {
        Kokkos::View<Scalar*,DeviceType> vector_k,vector_Ak;
        if (!isConstantStride() )
          vector_k = Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), whichVectors_[k]);
        else
          vector_k = Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), k);
        if (!A.isConstantStride() )
          vector_Ak = Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (A.view_.d_view, Kokkos::ALL (), A.whichVectors_[k]);
        else
          vector_Ak = Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (A.view_.d_view, Kokkos::ALL (), k);
        dots[k] = Kokkos::V_Dot (vector_k,vector_Ak);
      }
    }
    if (this->isDistributed()) {
      Array<Scalar> ldots(dots);
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,Teuchos::as<int>(numVecs),ldots.getRawPtr(),dots.getRawPtr());
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  norm2 (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    using Teuchos::ScalarTraits;
    typedef typename ScalarTraits<Scalar>::magnitudeType MT;
    typedef ScalarTraits<MT> STM;

    // FIXME (mfh 15 Sep 2013) When migrating to use Kokkos::View, we
    // should have the post-reduce kernel do the square root(s) on the
    // device.  If we make norm2() fill in a Kokkos::View instead of a
    // Teuchos::ArrayView, we could have norm2() leave the results of
    // the norm(s) on the device, instead of bringing them back to the
    // host.  I don't think it makes sense to template norm2() on the
    // output View type; users should use a View type compatible with
    // the MultiVector's native device.

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(
      as<size_t> (norms.size ()) != numVecs,
      std::runtime_error,
      "Tpetra::MultiVector::norm2(norms): norms.size() must be as large as the "
      "number of vectors in *this.  norms.size() = " << norms.size () << ", but "
      "*this.getNumVectors() = " << numVecs << ".");
    /*if (isConstantStride ()) {
      MVT::Norm2Squared (lclMV_,norms);
    }
    else {
      KMV v (MVT::getNode (lclMV_));
      ArrayRCP<Scalar> vi;
      for (size_t i=0; i < numVecs; ++i) {
        vi = arcp_const_cast<Scalar> (MVT::getValues (lclMV_, whichVectors_[i]));
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vi, MVT::getStride (lclMV_));
        norms[i] = MVT::Norm2Squared (v);
      }
    }*/

    if (isConstantStride()) {
      Kokkos::MV_Dot(&norms[0],view_.d_view,view_.d_view,getLocalLength());
    } else {
      for (size_t k = 0; k < numVecs; ++k) {
        Kokkos::View<Scalar*,DeviceType> vector_k = Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), whichVectors_[k]);
        norms[k] = Kokkos::V_Dot (vector_k,vector_k);
      }
    }

    if (this->isDistributed ()) {
      Array<MT> lnorms (norms);
      // FIXME (mfh 25 Apr 2012) Somebody explain to me why we're
      // calling Teuchos::reduceAll when MultiVector has a perfectly
      // good reduce() function.
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs),
                 lnorms.getRawPtr (), norms.getRawPtr ());
    }
    for (typename ArrayView<MT>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      (*n) = STM::squareroot (*n);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& weights,
                const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    using Teuchos::ScalarTraits;
    typedef ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;

    const char tfecfFuncName[] = "normWeighted()";

    const Mag OneOverN = ScalarTraits<Mag>::one() / getGlobalLength();
    bool OneW = false;
    const size_t numVecs = this->getNumVectors();
    if (weights.getNumVectors() == 1) {
      OneW = true;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(weights.getNumVectors() != numVecs, std::runtime_error,
        ": MultiVector of weights must contain either one vector or the same number of vectors as this.");
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*weights.getMap()), std::runtime_error,
      ": MultiVectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << *this->getMap()
      << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != weights.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif
    const size_t myLen = getLocalLength();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(norms.size()) != numVecs, std::runtime_error,
      ": norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride() && weights.isConstantStride()) {
      MVT::WeightedNorm(lclMV_,weights.lclMV_,norms);
    }
    else {
      KMV v(MVT::getNode(lclMV_)), w(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vptr = arcp_const_cast<Scalar> (MVT::getValues (lclMV_));
      ArrayRCP<Scalar> wptr = arcp_const_cast<Scalar> (MVT::getValues (weights.lclMV_));
      ArrayRCP<Scalar> wj   = wptr.persistingView (0,myLen);
      MVT::initializeValues (w,myLen, 1, wj, myLen);
      for (size_t j=0; j < numVecs; ++j) {
        ArrayRCP<Scalar> vj = getSubArrayRCP (vptr, j);
        MVT::initializeValues(v,myLen, 1,  vj, myLen);
        if (! OneW) {
          wj = weights.getSubArrayRCP(wptr,j);
          MVT::initializeValues (w, myLen, 1, wj, myLen);
        }
        norms[j] = MVT::WeightedNorm ((const KMV&)v, (const KMV &)w);
      }
    }
    if (this->isDistributed ()) {
      Array<Mag> lnorms(norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int>(numVecs),
                 lnorms.getRawPtr (), norms.getRawPtr ());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n * OneOverN);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  norm1 (const ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(norms.size()) != numVecs, std::runtime_error,
        "Tpetra::MultiVector::norm1(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride()) {
      MVT::Norm1(lclMV_,norms);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = arcp_const_cast<Scalar> (MVT::getValues(lclMV_,whichVectors_[j]) );
        MVT::initializeValues (v, MVT::getNumRows(lclMV_), 1, vj, MVT::getStride (lclMV_));
        norms[j] = MVT::Norm1 ((const KMV&)v);
      }
    }
    if (this->isDistributed ()) {
      Array<Mag> lnorms (norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs), lnorms.getRawPtr (), norms.getRawPtr ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  normInf (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_MAX;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(norms.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::normInf(norms): norms.size() must be as large as the number of vectors in *this.");
    if (isConstantStride ()) {
      MVT::NormInf(lclMV_,norms);
    }
    else {
      KMV v(MVT::getNode(lclMV_));
      ArrayRCP<Scalar> vj;
      for (size_t j=0; j < numVecs; ++j) {
        vj = Teuchos::arcp_const_cast<Scalar>( MVT::getValues(lclMV_,whichVectors_[j]) );
        MVT::initializeValues(v,MVT::getNumRows(lclMV_), 1, vj, MVT::getStride(lclMV_));
        norms[j] = MVT::NormInf((const KMV&)v);
      }
    }
    if (this->isDistributed()) {
      Array<Mag> lnorms(norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_MAX, as<int> (numVecs), lnorms.getRawPtr (), norms.getRawPtr ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  meanValue (const Teuchos::ArrayView<Scalar> &means) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Teuchos::ScalarTraits<Scalar> SCT;

    const size_t numVecs = getNumVectors();
    const size_t myLen   = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(means.size()) != numVecs, std::runtime_error,
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
    if (this->isDistributed()) {
      Array<Scalar> lmeans(means);
      // only combine if we are a distributed MV
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs),
                 lmeans.getRawPtr (), means.getRawPtr ());
    }
    // mfh 12 Apr 2012: Don't take out the cast from the ordinal type
    // to the magnitude type, since operator/ (std::complex<T>, int)
    // isn't necessarily defined.
    const Scalar OneOverN =
      SCT::one() / as<typename SCT::magnitudeType> (getGlobalLength ());
    for (typename ArrayView<Scalar>::iterator i = means.begin();
         i != means.begin() + numVecs;
         ++i)
    {
      (*i) = (*i) * OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  scale (Teuchos::ArrayView<const Scalar> alphas)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(alphas.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as the number of vectors in *this.");
    KMV vec(MVT::getNode(lclMV_));
    const size_t myLen = MVT::getNumRows(lclMV_);
    if (myLen == 0) return;
    ArrayRCP<Scalar> mybuf = MVT::getValuesNonConst(lclMV_);
    for (size_t j = 0; j < numVecs; ++j) {
      if (alphas[j] == Teuchos::ScalarTraits<Scalar>::one()) {
        // do nothing: NaN * 1.0 == NaN, Number*1.0 == Number
      }
      else {
        ArrayRCP<Scalar> mybufj = getSubArrayRCP(mybuf,j);
        MVT::initializeValues(vec,myLen,1,mybufj,myLen);
        MVT::Scale(vec,alphas[j]);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  scale (const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  reciprocal (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A) {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A,
                      const Scalar &beta)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

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
    const size_t numVecs = getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
        ": MultiVectors must have the same number of vectors.");


    if (isConstantStride() && A.isConstantStride ()) {
      Kokkos::MV_Add(view_.d_view,alpha,A.view_.d_view,beta,view_.d_view,getLocalLength());
    } else {
      for (size_t k = 0; k < numVecs; ++k) {
        Kokkos::V_Add (Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), whichVectors_[k]),
                        alpha,
                        Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (A.view_.d_view, Kokkos::ALL (), A.whichVectors_[k]),
                        beta,
                        Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (view_.d_view, Kokkos::ALL (), whichVectors_[k]));
      }
    }
    /*using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

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
    }*/
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::update(
                      const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A,
                      const Scalar &beta,  const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &B,
                      const Scalar &gamma)
  {
    const char tfecfFuncName[] = "update()";
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      A.getNumVectors() != this->getNumVectors() || B.getNumVectors() != this->getNumVectors(),
      std::runtime_error,
      "Tpetra::MultiVector::update(three MVs): MultiVectors must have the same number of vectors.");

    if (isConstantStride() && A.isConstantStride () && B.isConstantStride ()) {
      if(gamma == Teuchos::ScalarTraits<Scalar>::zero ()) {
        Kokkos::MV_Add(view_.d_view,alpha,A.view_.d_view,beta,B.view_.d_view);
      } else {
        Kokkos::MV_Add(view_.d_view,alpha,A.view_.d_view,gamma,view_.d_view);
        Kokkos::MV_Add(view_.d_view,beta,B.view_.d_view,1.0,view_.d_view);
      }
    } else { // some input (or *this) is not constant stride
      const size_t numVecs = getNumVectors ();
      for (size_t k = 0; k < numVecs; ++k) {
        if (gamma == Teuchos::ScalarTraits<Scalar>::zero ()) {
          // TODO: make sure it only uses LocalLength for add.
          Kokkos::V_Add (Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), whichVectors_[k]),
                          alpha,
                          Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (A.view_.d_view, Kokkos::ALL (), A.whichVectors_[k]),
                          beta,
                          Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (B.view_.d_view, Kokkos::ALL (), B.whichVectors_[k]));
        } else {
          Kokkos::V_Add (Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), whichVectors_[k]),
                          alpha,
                          Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (A.view_.d_view, Kokkos::ALL (), A.whichVectors_[k]),
                          gamma,
                          Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (view_.d_view, Kokkos::ALL (), whichVectors_[k]));
          Kokkos::V_Add (Kokkos::subview<Kokkos::View<Scalar*,DeviceType> > (view_.d_view, Kokkos::ALL (), whichVectors_[k]),
                          beta,
                          Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (B.view_.d_view, Kokkos::ALL (), B.whichVectors_[k]),
                          Teuchos::ScalarTraits<Scalar>::one (),
                          Kokkos::subview<Kokkos::View<Scalar*,DeviceType> >  (view_.d_view, Kokkos::ALL (), whichVectors_[k]));
        }
      }
    }
    /*
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

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
    }*/
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getData (size_t j) const
  {
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::getData()")
    return node->template viewBuffer<Scalar> (getLocalLength (), getSubArrayRCP (MVT::getValues(lclMV_), j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getDataNonConst(size_t j)
  {
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::getDataNonConst()")
    return node->template viewBufferNonConst<Scalar> (KokkosClassic::ReadWrite, getLocalLength (), getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  operator= (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &source)
  {
    //const char tfecfFuncName[] = "operator=";
    DO::operator=(source);
    RCP<Node> node = MVT::getNode(source.lclMV_);
    const LocalOrdinal myLen = source.getLocalLength();
    const size_t numVecs = source.getNumVectors();

    // On host-type Kokkos Nodes, allocBuffer() just calls the
    // one-argument version of arcp to allocate memory.  This should
    // not fill the memory by default, otherwise we would lose the
    // first-touch allocation optimization.
/*    ArrayRCP<Scalar> data = (myLen > 0) ?
      node->template allocBuffer<Scalar> (myLen * numVecs) :
      Teuchos::null;*/

    view_ = source.view_; // OK View Semantics from Kokkos
    whichVectors_ = source.whichVectors_; // Probably not ok (probably constitutes deep copy)

    ArrayRCP<Scalar> data = (myLen > 0) ?
      Kokkos::Compat::persistingView(view_.d_view) :
      Teuchos::null;
    const size_t stride = (myLen > 0) ? myLen : size_t (0);
    // This just sets the dimensions, pointer, and stride of lclMV_.
    MVT::initializeValues (lclMV_, myLen, numVecs, data, stride);

    // Check for special case of this=Source, in which case we do nothing.
    /*if (this != &source) {
      // Whether the input and *this are compatible on the calling process.
      const int locallyCompat =
        (this->getLocalLength () == source.getLocalLength ()) ? 1 : 0;
#ifdef HAVE_TPETRA_DEBUG
      int globallyCompat = 1; // innocent until proven guilty
      Teuchos::reduceAll<int, int> (* (this->getMap ()->getComm ()), Teuchos::REDUCE_MIN,
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
        KOKKOS_NODE_TRACE("MultiVector::operator=()")
        node->template copyBuffers<Scalar> (getLocalLength () * numVecs,
                                            MVT::getValues (source.lclMV_),
                                            MVT::getValuesNonConst (lclMV_));
      }
      else {
        // We have to copy the columns one at a time.
        for (size_t j=0; j < numVecs; ++j) {
          KOKKOS_NODE_TRACE("MultiVector::operator=()")
          node->template copyBuffers<Scalar> (getLocalLength (),
                                              source.getSubArrayRCP (MVT::getValues (source.lclMV_), j),
                                              getSubArrayRCP (MVT::getValuesNonConst(lclMV_), j));
        }
      }
    }*/
    return *this;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subCopy (const Teuchos::ArrayView<const size_t> &cols) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(cols.size() < 1, std::runtime_error,
        "Tpetra::MultiVector::subCopy(cols): cols must contain at least one column.");
    size_t numCopyVecs = cols.size();
    const bool zeroData = false;
    RCP<Node> node = MVT::getNode(lclMV_);
    RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > mv;
    // mv is allocated with constant stride
    mv = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > (this->getMap (), numCopyVecs, zeroData));
    // copy data from *this into mv
    for (size_t j=0; j<numCopyVecs; ++j) {
      KOKKOS_NODE_TRACE("MultiVector::subCopy()")
      node->template copyBuffers<Scalar> (getLocalLength (),
                                          getSubArrayRCP (MVT::getValues (lclMV_), cols[j]),
                                          MVT::getValuesNonConst (mv->lclMV_, j));
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subCopy (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
        "Tpetra::MultiVector::subCopy(Range1D): range must include at least one vector.");
    size_t numCopyVecs = colRng.size();
    const bool zeroData = false;
    RCP<Node> node = MVT::getNode(lclMV_);
    RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > mv;
    // mv is allocated with constant stride
    mv = rcp (new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > (this->getMap (), numCopyVecs, zeroData));
    // copy data from *this into mv
    for (size_t js=colRng.lbound(), jd=0; jd<numCopyVecs; ++jd, ++js) {
      KOKKOS_NODE_TRACE("MultiVector::subCopy()")
      node->template copyBuffers<Scalar> (getLocalLength (),
                                          getSubArrayRCP (MVT::getValues (lclMV_), js),
                                          MVT::getValuesNonConst (mv->lclMV_, jd));
    }
    return mv;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
              size_t offset) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

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
    RCP<const MV> subViewMV;
    if (isConstantStride()) {
      subViewMV = rcp (new MV (subMap,
          Kokkos::subview<view_type> (view_,std::make_pair(offset,offset+newNumRows),Kokkos::ALL())));
    }
    else {
      subViewMV = rcp (new MV (subMap,
          Kokkos::subview<view_type> (view_,std::make_pair(offset,offset+newNumRows),Kokkos::ALL()),
          whichVectors_()));
    }
    return subViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                      size_t offset)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    const size_t newNumRows = subMap->getNodeNumElements();
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
      subViewMV = rcp (new MV (subMap,
          Kokkos::subview<view_type> (view_,std::make_pair(offset,offset+newNumRows),Kokkos::ALL())));
    }
    else {
      subViewMV = rcp (new MV (subMap,
          Kokkos::subview<view_type> (view_,std::make_pair(offset,offset+newNumRows),Kokkos::ALL()),
          whichVectors_()));
    }
    return subViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subView (const ArrayView<const size_t> &cols) const
  {
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    TEUCHOS_TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subViewNonConst(ArrayView): range must include at least one vector.");
    if (isConstantStride()) {
      return rcp (new MV (this->getMap (), view_, cols));
    }
    // else, lookup current whichVectors_ using cols
    Array<size_t> newcols(cols.size());
    for (size_t j = 0; j < as<size_t> (cols.size ()); ++j) {
      newcols[j] = whichVectors_[cols[j]];
    }
    return rcp (new MV (this->getMap (), view_, newcols ()));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subView (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
    //size_t numViewVecs = colRng.size();
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning

    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      return rcp (new MV (this->getMap (), Kokkos::subview<view_type>(view_,Kokkos::ALL(),
                          std::pair<unsigned int,unsigned int>(colRng.lbound(),colRng.ubound()+1))));
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Array<size_t> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    return rcp (new MV (this->getMap (), view_, whchvecs));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subViewNonConst (const ArrayView<const size_t> &cols)
  {
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    TEUCHOS_TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subViewNonConst(ArrayView): range must include at least one vector.");
    if (isConstantStride()) {
      return rcp (new MV (this->getMap (), view_, cols));
    }
    // else, lookup current whichVectors_ using cols
    Array<size_t> newcols(cols.size());
    for (size_t j = 0; j < as<size_t> (cols.size ()); ++j) {
      newcols[j] = whichVectors_[cols[j]];
    }
    return rcp (new MV (this->getMap (), view_, newcols ()));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subViewNonConst (const Teuchos::Range1D &colRng)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    TEUCHOS_TEST_FOR_EXCEPTION(colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(Range1D): range must include at least one vector.");
    //size_t numViewVecs = colRng.size();
    // this is const, so the lclMV_ is const, so that we can only get const buffers
    // we will cast away the const; this is okay, because
    //   a) the constructor doesn't modify the data, and
    //   b) we are encapsulating in a const MV before returning

    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride()) {
      // view goes from first entry of first vector to last entry of last vector
      return rcp (new MV (this->getMap (), Kokkos::subview<view_type>(view_,Kokkos::ALL(),
                          std::pair<unsigned int,unsigned int>(colRng.lbound(),colRng.ubound()+1))));
    }
    // otherwise, use a subset of this whichVectors_ to construct new multivector
    Array<size_t> whchvecs( whichVectors_.begin()+colRng.lbound(), whichVectors_.begin()+colRng.ubound()+1 );
    return rcp (new MV (this->getMap (), view_, whchvecs));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getVector (size_t j) const
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVector(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif

    unsigned int jj;
    if (isConstantStride ()) {
      jj = j;
    } else {
      jj = whichVectors_[j];
    }
    return rcp_const_cast<const V> (rcp (new V (this->getMap (), Kokkos::subview<view_type> (view_,Kokkos::ALL(),jj))));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getVectorNonConst(size_t j)
  {
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVector(j): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif

    unsigned int jj;
    if (isConstantStride ()) {
      jj = j;
    } else {
      jj = whichVectors_[j];
    }
    return rcp_const_cast<V> (rcp (new V (this->getMap (), Kokkos::subview<view_type> (view_,Kokkos::ALL(),jj))));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get1dCopy (Teuchos::ArrayView<Scalar> A, size_t LDA) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const char tfecfFuncName[] = "get1dCopy(A,LDA)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < getLocalLength(), std::runtime_error,
      ": specified stride is not large enough for local vector length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(A.size()) < LDA*(getNumVectors()-1)+getLocalLength(), std::runtime_error,
      ": specified stride/storage is not large enough for the number of vectors.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t myStride = MVT::getStride(lclMV_),
                  numCols = getNumVectors(),
                  myLen   = getLocalLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mydata = MVT::getValues(lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::getVectorNonConst()")
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get2dCopy (Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;

    const char tfecfFuncName[] = "get2dCopy(ArrayOfPtrs)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<size_t>(ArrayOfPtrs.size()) != getNumVectors(), std::runtime_error,
        ": Array of pointers must contain as many pointers as the MultiVector has rows.");
    RCP<Node> node = MVT::getNode(lclMV_);
    const size_t numCols = getNumVectors(),
                   myLen = getLocalLength();
    if (myLen > 0) {
      ArrayRCP<const Scalar> mybuff = MVT::getValues(lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::get2dCopy()")
      ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(mybuff.size(), mybuff);
      for (size_t j=0; j<numCols; ++j) {
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(ArrayOfPtrs[j].size()) != getLocalLength(), std::runtime_error,
          ": The ArrayView provided in ArrayOfPtrs[" << j << "] was not large enough to contain the local entries.");
#endif
        ArrayRCP<const Scalar> myviewj = getSubArrayRCP(myview,j);
        std::copy(myviewj,myviewj+myLen,ArrayOfPtrs[j].begin());
      }
      myview = Teuchos::null;
      mybuff = Teuchos::null;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::get1dView () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!isConstantStride(), std::runtime_error,
      "Tpetra::MultiVector::get1dView() requires that this MultiVector have constant stride.");
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::get1dView()")
    return node->template viewBuffer<Scalar>( getStride()*(getNumVectors()-1)+getLocalLength(), MVT::getValues(lclMV_) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::get1dViewNonConst ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!isConstantStride(), std::runtime_error,
      "Tpetra::MultiVector::get1dViewNonConst(): requires that this MultiVector have constant stride.");
    Teuchos::RCP<Node> node = MVT::getNode(lclMV_);
    KOKKOS_NODE_TRACE("MultiVector::get1dViewNonConst()")
    return node->template viewBufferNonConst<Scalar>( KokkosClassic::ReadWrite, getStride()*(getNumVectors()-1)+getLocalLength(), MVT::getValuesNonConst(lclMV_) );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::get2dViewNonConst()
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
        KOKKOS_NODE_TRACE("MultiVector::get2dViewNonConst()")
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
        KOKKOS_NODE_TRACE("MultiVector::get2dViewNonConst()")
        ArrayRCP<Scalar> myview = node->template viewBufferNonConst<Scalar>(KokkosClassic::ReadWrite,myStride*(numCols-1)+myLen,MVT::getValuesNonConst(lclMV_));
        for (size_t j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::get2dView() const
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
        KOKKOS_NODE_TRACE("MultiVector::get2dView()")
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
        KOKKOS_NODE_TRACE("MultiVector::get2dView()")
        ArrayRCP<const Scalar> myview = node->template viewBuffer<Scalar>(myStride*(numCols-1)+myLen,MVT::getValues(lclMV_));
        for (size_t j=0; j<myCols; ++j) {
          views[j] = myview.persistingView(whichVectors_[j]*myStride,myLen);
        }
      }
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  multiply (Teuchos::ETransp transA,
            Teuchos::ETransp transB,
            const Scalar &alpha,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& B,
            const Scalar &beta)
  {
    using Teuchos::NO_TRANS;      // enums
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::null;
    using Teuchos::ScalarTraits;  // traits
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

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

    bool A_is_local = !A.isDistributed();
    bool B_is_local = !B.isDistributed();
    bool C_is_local = !this->isDistributed();
    bool Case1 = ( C_is_local &&  A_is_local &&  B_is_local);                                           // Case 1: C(local) = A^X(local) * B^X(local)
    bool Case2 = ( C_is_local && !A_is_local && !B_is_local && transA==CONJ_TRANS && transB==NO_TRANS); // Case 2: C(local) = A^T(distr) * B  (distr)
    bool Case3 = (!C_is_local && !A_is_local &&  B_is_local && transA==NO_TRANS  );                     // Case 3: C(distr) = A  (distr) * B^X(local)

    // Test that we are considering a meaningful cases
    TEUCHOS_TEST_FOR_EXCEPTION( !Case1 && !Case2 && !Case3, std::runtime_error,
        errPrefix << "multiplication of op(A) and op(B) into *this is not a supported use case.");

    if (beta != ScalarTraits<Scalar>::zero() && Case2)
    {
      // if Case2, then C is local and contributions must be summed across all nodes
      // however, if beta != 0, then accumulate beta*C into the sum
      // when summing across all nodes, we only want to accumulate this once, so
      // set beta == 0 on all nodes except node 0
      int MyPID = this->getMap()->getComm()->getRank();
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
        KOKKOS_NODE_TRACE("MultiVector::multiply()")
        node->template copyBuffers<Scalar>(getLocalLength(),MVT::getValues(C_mv,j),MVT::getValuesNonConst(lclMV_,whichVectors_[j]));
      }
    }

    // If Case 2 then sum up *this and distribute it to all processors.
    if (Case2) {
      this->reduce();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  elementWiseMultiply (Scalar scalarAB,
                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
                       const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& B,
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::reduce()
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
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
    KOKKOS_NODE_TRACE("MultiVector::reduce()")
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
    reduceAll<int,Scalar> (*comm, REDUCE_SUM, as<int> (numCols*myLen),
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
    lclMV_.replaceLocalValue (Teuchos::as<size_t> (MyRow), colInd, ScalarValue);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
    lclMV_.sumIntoLocalValue (Teuchos::as<size_t> (MyRow), colInd, ScalarValue);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class T>
  Teuchos::ArrayRCP<T>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  const KokkosClassic::MultiVector<Scalar,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getLocalMV() const {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::view_type
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getLocalView() const {
    return view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  KokkosClassic::MultiVector<Scalar,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getLocalMVNonConst() {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::description() const
  {
    using std::endl;
    std::ostringstream oss;
    oss << Teuchos::typeName (*this) << " {"
        << "label: \"" << this->getObjectLabel () << "\""
        << ", numRows: " << getGlobalLength ()
        << ", numCols: " << getNumVectors ()
        << ", isConstantStride: " << isConstantStride ();
    if (isConstantStride ()) {
      oss << ", columnStride: " << getStride ();
    }
    oss << "}";
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
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
      Teuchos::OSTab tab0 (out);

      if (myImageID == 0) { // >= VERB_LOW prints description()
        out << "Tpetra::MultiVector:" << endl;
        Teuchos::OSTab tab1 (out);
        out << "Template parameters:" << endl;
        {
          Teuchos::OSTab tab2 (out);
          out << "Scalar: " << Teuchos::TypeNameTraits<Scalar>::name () << endl
              << "LocalOrdinal: " << Teuchos::TypeNameTraits<LocalOrdinal>::name () << endl
              << "GlobalOrdinal: " << Teuchos::TypeNameTraits<GlobalOrdinal>::name () << endl
              << "Node: " << Teuchos::TypeNameTraits<Node>::name () << endl;
        }
        out << "label: \"" << this->getObjectLabel () << "\"" << endl
            << "numRows: " << getGlobalLength () << endl
            << "numCols: " << getNumVectors () << endl
            << "isConstantStride: " << isConstantStride () << endl;
        if (isConstantStride ()) {
          out << "columnStride: " << getStride () << endl;
        }
      }
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // At verbosity > VERB_LOW, each process prints something.
            out << "Process " << myImageID << ":" << endl;
            Teuchos::OSTab tab2 (out);

            // >= VERB_MEDIUM: print the local vector length.
            out << "local length: " << getLocalLength();
            if (vl != VERB_MEDIUM) {
              // >= VERB_HIGH: print isConstantStride() and getStride()
              if (isConstantStride()) {
                out << "constant stride: " << getStride() << endl;
              }
              else {
                out << "not constant stride" << endl;
              }
              if (vl == VERB_EXTREME) {
                // VERB_EXTREME: print all the values in the multivector.
                out << "values: " << endl;
                ArrayRCP<ArrayRCP<const Scalar> > X = this->get2dView();
                out << "[";
                for (size_t i = 0; i < getLocalLength(); ++i) {
                  for (size_t j = 0; j < getNumVectors(); ++j) {
                    out << X[j][i];
                    if (j + 1 < getNumVectors()) {
                      out << ", ";
                    }
                  } // for each column
                  if (i + 1 < getLocalLength ()) {
                    out << "; ";
                  } else {
                    out << endl;
                  }
                } // for each row
                out << "]" << endl;
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViews() const
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViewsNonConst (KokkosClassic::ReadWriteOption rwo)
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::releaseViews () const
  {
    // Do nothing in Kokkos::View implementation
  }

#else // NOT TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViews() const
  {
    Teuchos::RCP<Node> node = this->getMap ()->getNode ();
    if (cview_.is_null () && getLocalLength () > 0) {
      Teuchos::ArrayRCP<const Scalar> buff = MVT::getValues (lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::createViews()")
      cview_ = node->template viewBuffer<Scalar> (buff.size (), buff);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViewsNonConst (KokkosClassic::ReadWriteOption rwo)
  {
    Teuchos::RCP<Node> node = this->getMap ()->getNode ();
    if (ncview_.is_null () && getLocalLength () > 0) {
      Teuchos::ArrayRCP<Scalar> buff = MVT::getValuesNonConst (lclMV_);
      KOKKOS_NODE_TRACE("MultiVector::createViewsNonConst()")
      ncview_ = node->template viewBufferNonConst<Scalar> (rwo, buff.size (), buff);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::releaseViews () const
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

#endif // TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
  {
    replaceMap (newMap);
  }

  /// \brief Nonmember MultiVector constructor with view semantics using user-allocated data.
  /// \relatesalso MultiVector
  /// \relatesalso Vector
  ///
  /// \warning This function is not supported for all Kokkos Node
  ///   types.  Specifically, it is not typically supported for
  ///   GPU accelerator-based nodes like KokkosClassic::ThrustGPUNode.
  ///
  /// \param map [in] The Map describing the distribution of rows of
  ///   the multivector.
  /// \param view [in/out] A pointer to column-major dense matrix
  ///   data.  This will be the multivector's data on the calling
  ///   process.  The multivector will use the pointer directly,
  ///   without copying.
  /// \param LDA [in] The leading dimension (a.k.a. "stride") of the
  ///   column-major input data.
  /// \param numVectors [in] The number of columns in the input data.
  ///   This will be the number of vectors in the returned
  ///   multivector.
  ///
  /// \node To Kokkos and Tpetra developers: If you add a new Kokkos
  ///   Node type that is a host Node type (where memory lives in user
  ///   space, not in a different space as on a GPU), you will need to
  ///   add a specialization of Tpetra::details::ViewAccepter for your
  ///   new Node type.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  createMultiVectorFromView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >& map,
                             const Teuchos::ArrayRCP<Scalar>& view,
                             size_t LDA,
                             size_t numVectors)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Tpetra::createMultiVectorFromView: Not implemented for KokkosDeviceWrapperNode");

    /*
    using Teuchos::rcp;
    typedef Tpetra::details::ViewAccepter<Node> VAN;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    // This uses a protected MultiVector constructor, but this
    // nonmember function was declared a friend of MultiVector.
    //
    // The ViewAccepter expression will fail to compile for
    // unsupported Kokkos Node types.
    return rcp (new MV (map, VAN::template acceptView<Scalar> (view),
                        LDA, numVectors, HOST_VIEW_CONSTRUCTOR));*/
  }

  template<class DstType, class SrcType, class DeviceType,bool DstConstStride,bool SrcConstStride>
  struct DeepCopySelectedVectors {
    typedef DeviceType device_type;
    DstType dst;
    SrcType src;
    Kokkos::View<int*,DeviceType> whichVectorSrc;
    Kokkos::View<int*,DeviceType> whichVectorDst;
    int n;
    DeepCopySelectedVectors(DstType dst_, SrcType src_,
                            Kokkos::View<int*,DeviceType> whichVectorDst_,
                            Kokkos::View<int*,DeviceType> whichVectorSrc_):
                              dst(dst_),src(src_),whichVectorDst(whichVectorDst_),whichVectorSrc(whichVectorSrc_),n(whichVectorSrc_.dimension_0()) {};
    void KOKKOS_INLINE_FUNCTION operator()(int i) const {
      if(DstConstStride ) {
        if(SrcConstStride) {
          for(int j = 0; j<n ; j++)
            dst(i,j) = src(i,j);
        } else {
          for(int j = 0; j<n ; j++)
            dst(i,j) = src(i,whichVectorSrc(j));
        }
      } else {
        if(SrcConstStride) {
          for(int j = 0; j<n ; j++)
            dst(i,whichVectorDst(j)) = src(i,j);
        } else {
          for(int j = 0; j<n ; j++)
            dst(i,whichVectorDst(j)) = src(i,whichVectorSrc(j));
        }
      }
    }
  };

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
    createCopy( const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& src) {
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    MV cpy(src.getMap(),src.getNumVectors());
    if(src.isConstantStride())
      Kokkos::deep_copy(cpy.getLocalView(),src.getLocalView());
    else {
      if(src.getLocalView().modified_device>=src.getLocalView().modified_host) {
        Kokkos::View<int*,DeviceType> whichVectors("MultiVector::createCopy::WhichVectors",src.whichVectors_.size());
        for(int i = 0; i < src.whichVectors_.size(); i++)
          whichVectors(i)=src.whichVectors_[i];
        Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MV::view_type::t_dev,typename MV::view_type::t_dev,DeviceType,true,false>
                                                  (cpy.getLocalView().template view<DeviceType>(),
                                                   src.getLocalView().template view<DeviceType>(),
                                                   whichVectors,whichVectors));
      } else {
        Kokkos::View<int*,typename DeviceType::host_mirror_device_type> whichVectors("MultiVector::createCopy::WhichVectors",src.whichVectors_.size());
        for(int i = 0; i < src.whichVectors_.size(); i++)
          whichVectors(i)=src.whichVectors_[i];
        Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MV::view_type::t_host,typename MV::view_type::t_host,typename DeviceType::host_mirror_device_type,true,false>
                                                  (cpy.getLocalView().template view<typename DeviceType::host_mirror_device_type>(),
                                                   src.getLocalView().template view<typename DeviceType::host_mirror_device_type>(),
                                                   whichVectors,whichVectors));
      }
    }
    return cpy;
  }

  template <class DS, class DL, class DG, class DD, class SS, class SL, class SG, class SD>
  void deep_copy( MultiVector<DS,DL,DG,Kokkos::Compat::KokkosDeviceWrapperNode<DD> >& dst,
                  const MultiVector<SS,SL,SG,Kokkos::Compat::KokkosDeviceWrapperNode<SD> >& src) {
    typedef MultiVector<DS,DL,DG,Kokkos::Compat::KokkosDeviceWrapperNode<DD> > MVD;
    typedef const MultiVector<SS,SL,SG,Kokkos::Compat::KokkosDeviceWrapperNode<SD> > MVS;
    if(src.isConstantStride() && dst.isConstantStride()) {
      Kokkos::deep_copy(dst.getLocalView(),src.getLocalView());
    }
    else {
      if(dst.isConstantStride()) {
        if(src.getLocalView().modified_device>=src.getLocalView().modified_host) {
          Kokkos::View<int*,DD> whichVectors("MultiVector::createCopy::WhichVectors",src.whichVectors_.size());
          for(int i = 0; i < src.whichVectors_.size(); i++)
            whichVectors(i)=src.whichVectors_[i];
          Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::view_type::t_dev,typename MVS::view_type::t_dev,DD,true,false>
                                                    (dst.getLocalView().template view<DD>(),
                                                     src.getLocalView().template view<DD>(),
                                                     whichVectors,whichVectors));
        } else {
          Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectors("MultiVector::createCopy::WhichVectors",src.whichVectors_.size());
          for(int i = 0; i < src.whichVectors_.size(); i++)
            whichVectors(i)=src.whichVectors_[i];
          Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::view_type::t_host,typename MVS::view_type::t_host,typename DD::host_mirror_device_type,true,false>
                                                    (dst.getLocalView().template view<typename DD::host_mirror_device_type>(),
                                                     src.getLocalView().template view<typename DD::host_mirror_device_type>(),
                                                     whichVectors,whichVectors));
        }
      } else {
        if(src.isConstantStride()) {
          if(src.getLocalView().modified_device>=src.getLocalView().modified_host) {
            Kokkos::View<int*,DD> whichVectors("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectors(i)=dst.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::view_type::t_dev,typename MVS::view_type::t_dev,DD,false,true>
                                                      (dst.getLocalView().template view<DD>(),
                                                       src.getLocalView().template view<DD>(),
                                                       whichVectors,whichVectors));
          } else {
            Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectors("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectors(i)=dst.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::view_type::t_host,typename MVS::view_type::t_host,typename DD::host_mirror_device_type,false,true>
                                                      (dst.getLocalView().template view<typename DD::host_mirror_device_type>(),
                                                       src.getLocalView().template view<typename DD::host_mirror_device_type>(),
                                                       whichVectors,whichVectors));
          }
        } else {
          if(src.getLocalView().modified_device>=src.getLocalView().modified_host) {
            Kokkos::View<int*,DD> whichVectorsDst("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsDst(i)=dst.whichVectors_[i];
            Kokkos::View<int*,DD> whichVectorsSrc("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsSrc(i)=src.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::view_type::t_dev,typename MVS::view_type::t_dev,DD,false,false>
                                                      (dst.getLocalView().template view<DD>(),
                                                       src.getLocalView().template view<DD>(),
                                                       whichVectorsDst,whichVectorsSrc));
          } else {
            Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectorsDst("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsDst(i)=dst.whichVectors_[i];
            Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectorsSrc("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsSrc(i)=src.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::view_type::t_host,typename MVS::view_type::t_host,typename DD::host_mirror_device_type,false,false>
                                                      (dst.getLocalView().template view<typename DD::host_mirror_device_type>(),
                                                       src.getLocalView().template view<typename DD::host_mirror_device_type>(),
                                                       whichVectorsDst,whichVectorsSrc));
          }
        }
      }
    }
  }
} // namespace Tpetra


#endif // TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP
