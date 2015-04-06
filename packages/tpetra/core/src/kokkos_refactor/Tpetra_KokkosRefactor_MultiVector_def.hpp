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

#include <Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp>

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_KokkosRefactor_MultiVector_decl.hpp"
#endif

#include <KokkosCompat_View.hpp>
#include <Kokkos_MV.hpp>
#include <Kokkos_MV_GEMM.hpp>
#include <Kokkos_Blas1_MV.hpp>
#include <Kokkos_Random.hpp>

namespace { // (anonymous)

  /// \brief Allocate and return a 2-D Kokkos::DualView for Tpetra::MultiVector.
  ///
  /// This function takes the same first three template parameters as
  /// Tpetra::MultiVector.  The fourth template parameter is the
  /// Kokkos "device type" (same as the "DeviceType" template
  /// parameter below).
  ///
  /// \param lclNumRows [in] Number of rows in the DualView.
  ///   "Local" means "local to the calling MPI process."
  /// \param numCols [in] Number of columns in the DualView.
  /// \param zeroOut [in] Whether to initialize all the entries of the
  ///   DualView to zero.  Kokkos does first-touch initialization.
  ///
  /// \return The allocated Kokkos::DualView.
  template<class S, class LO, class GO, class D>
  typename Tpetra::MultiVector<S, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<D>, false>::dual_view_type
  allocDualView (const size_t lclNumRows, const size_t numCols, const bool zeroOut = true)
  {
    typedef typename Tpetra::MultiVector<S, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<D>, false>::dual_view_type dual_view_type;
    const char* label = "MV::DualView";

    (void) zeroOut;
    return dual_view_type (label, lclNumRows, numCols);

    // if (zeroOut) {
    //   return dual_view_type (label, lclNumRows, numCols);
    // } else {
    //   // FIXME (mfh 18 Feb 2015) This is just a hack, until
    //   // Kokkos::DualView accepts an AllocationProperties initial
    //   // argument, just like Kokkos::View does.  However, the hack is
    //   // harmless, since it does what the (currently nonexistent)
    //   // equivalent DualView constructor would have done anyway.
    //   typename dual_view_type::t_dev d_view (Kokkos::ViewAllocateWithoutInitializing (label), lclNumRows, numCols);
    //   typename dual_view_type::t_host h_view = Kokkos::create_mirror_view (d_view);
    //   return dual_view_type (d_view, h_view);
    // }
  }

  // Convert 1-D Teuchos::ArrayView to an unmanaged 1-D host Kokkos::View.
  //
  // T: The type of the entries of the View.
  // ExecSpace: The execution space (corresponds to Tpetra's DeviceType).
  template<class T, class ExecSpace>
  struct MakeUnmanagedView {
    // The 'false' part of the branch carefully ensures that this
    // won't attempt to use a host execution space that hasn't been
    // initialized.  For example, if Kokkos::OpenMP is disabled and
    // Kokkos::Threads is enabled, the latter is always the default
    // execution space of Kokkos::HostSpace, even when ExecSpace is
    // Kokkos::Serial.  That's why we go through the trouble of asking
    // Kokkos::DualView what _its_ space is.  That seems to work
    // around this default execution space issue.
    typedef typename Kokkos::Impl::if_c<
      Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<ExecSpace, Kokkos::HostSpace>::value,
      typename ExecSpace::device_type,
      typename Kokkos::DualView<T*, ExecSpace>::host_mirror_space>::type host_exec_space;
    typedef Kokkos::LayoutLeft array_layout;
    typedef Kokkos::View<T*, array_layout, host_exec_space,
                         Kokkos::MemoryUnmanaged> view_type;

    static view_type getView (const Teuchos::ArrayView<T>& x_in)
    {
      const size_t numEnt = static_cast<size_t> (x_in.size ());
      if (numEnt == 0) {
        return view_type ();
      } else {
        return view_type (x_in.getRawPtr (), numEnt);
      }
    }
  };

} // namespace (anonymous)


namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  vectorIndexOutOfRange (const size_t VectorIndex) const {
    return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector () :
    base_type (Teuchos::rcp (new map_type ()))
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const size_t numVecs,
               const bool zeroOut) : /* default is true */
    base_type (map)
  {
   TEUCHOS_TEST_FOR_EXCEPTION(
     numVecs < 1, std::invalid_argument, "Tpetra::MultiVector::MultiVector"
     "(map,numVecs,zeroOut): numVecs = " << numVecs << " < 1.");
    const size_t lclNumRows = this->getLocalLength ();
    view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, DeviceType> (lclNumRows, numVecs, zeroOut);
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& source) :
    base_type (source),
    view_ (source.view_),
    origView_ (source.origView_),
    whichVectors_ (source.whichVectors_)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& source,
               const Teuchos::DataAccess copyOrView) :
    base_type (source),
    view_ (source.view_),
    origView_ (source.origView_),
    whichVectors_ (source.whichVectors_)
  {
    const char tfecfFuncName[] = "MultiVector(const MultiVector&, "
      "const Teuchos::DataAccess): ";

    if (copyOrView == Teuchos::Copy) {
      // Reuse the conveniently already existing function that creates
      // a deep copy.
      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> cpy =
        createCopy (source);
      this->view_ = cpy.view_;
      this->origView_ = cpy.origView_;
      this->whichVectors_ = cpy.whichVectors_;
    }
    else if (copyOrView == Teuchos::View) {
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        true, std::invalid_argument, "Second argument 'copyOrView' has an "
        "invalid value " << copyOrView << ".  Valid values include "
        "Teuchos::Copy = " << Teuchos::Copy << " and Teuchos::View = "
        << Teuchos::View << ".");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const dual_view_type& view) :
    base_type (map),
    view_ (view),
    origView_ (view)
  {
    const char tfecfFuncName[] = "MultiVector(map,view): ";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.dimension_1 () > 1) ? stride[1] :
      origView_.dimension_0 ();
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::DualView's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input view's first dimension (number "
      "of rows = " << view.dimension_0 () << ") does not not match the number "
      "of entries in the Map on the calling process.");
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const typename dual_view_type::t_dev& d_view) :
    base_type (map)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "MultiVector(map,d_view): ";

    // Get stride of view: if second dimension is 0, the stride might
    // be 0, so take view_dimension instead.
    size_t stride[8];
    d_view.stride (stride);
    const size_t LDA = (d_view.dimension_1 () > 1) ? stride[1] :
      d_view.dimension_0 ();
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::View's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input view's first dimension (number "
      "of rows = " << d_view.dimension_0 () << ") does not not match the "
      "number of entries in the Map on the calling process.");

    // The difference between create_mirror and create_mirror_view, is
    // that the latter copies to host.  We don't necessarily want to
    // do that; we just want to allocate the memory.
    view_ = dual_view_type (d_view, Kokkos::create_mirror (d_view));
    // The user gave us a device view.  We take it as canonical, which
    // means we mark it as "modified," so that the next sync will
    // synchronize it with the host view.
    view_.template modify<execution_space> ();
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const dual_view_type& view,
               const dual_view_type& origView) :
    base_type (map),
    view_ (view),
    origView_ (origView)
  {
    const char tfecfFuncName[] = "MultiVector(map,view,origView): ";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.dimension_1 () > 1) ? stride[1] :
      origView_.dimension_0 ();
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::DualView's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input origView's first dimension (number "
      "of rows = " << origView.dimension_0 () << ") does not not match the number "
      "of entries in the Map on the calling process.");
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const dual_view_type& view,
               const Teuchos::ArrayView<const size_t>& whichVectors) :
    base_type (map),
    view_ (view),
    origView_ (view),
    whichVectors_ (whichVectors.begin (), whichVectors.end ())
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    const char tfecfFuncName[] = "MultiVector(map,view,whichVectors): ";

    const size_t lclNumRows = map.is_null () ? size_t (0) :
      map->getNodeNumElements ();
    // Check dimensions of the input DualView.  We accept that Kokkos
    // might not allow construction of a 0 x m (Dual)View with m > 0,
    // so we only require the number of rows to match if the
    // (Dual)View has more than zero columns.  Likewise, we only
    // require the number of columns to match if the (Dual)View has
    // more than zero rows.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      view.dimension_1 () != 0 && static_cast<size_t> (view.dimension_0 ()) < lclNumRows,
      std::invalid_argument, "view.dimension_0() = " << view.dimension_0 ()
      << " < map->getNodeNumElements() = " << lclNumRows << ".");
    if (whichVectors.size () != 0) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        view.dimension_1 () != 0 && view.dimension_1 () == 0,
        std::invalid_argument, "view.dimension_1() = 0, but whichVectors.size()"
        " = " << whichVectors.size () << " > 0.");
      size_t maxColInd = 0;
      typedef Teuchos::ArrayView<const size_t>::size_type size_type;
      for (size_type k = 0; k < whichVectors.size (); ++k) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          whichVectors[k] == Teuchos::OrdinalTraits<size_t>::invalid (),
          std::invalid_argument, "whichVectors[" << k << "] = "
          "Teuchos::OrdinalTraits<size_t>::invalid().");
        maxColInd = std::max (maxColInd, whichVectors[k]);
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        view.dimension_1 () != 0 && static_cast<size_t> (view.dimension_1 ()) <= maxColInd,
        std::invalid_argument, "view.dimension_1() = " << view.dimension_1 ()
        << " <= max(whichVectors) = " << maxColInd << ".");
    }

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.dimension_1 () > 1) ? stride[1] :
      origView_.dimension_0 ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument,
      "LDA = " << LDA << " < this->getLocalLength() = " << lclNumRows << ".");

    if (whichVectors.size () == 1) {
      // If whichVectors has only one entry, we don't need to bother
      // with nonconstant stride.  Just shift the view over so it
      // points to the desired column.
      //
      // NOTE (mfh 10 May 2014) This is a special case where we set
      // origView_ just to view that one column, not all of the
      // original columns.  This ensures that the use of origView_ in
      // offsetView works correctly.
      const std::pair<size_t, size_t> colRng (whichVectors[0],
                                              whichVectors[0] + 1);
      view_ = subview (view_, ALL (), colRng);
      origView_ = subview (origView_, ALL (), colRng);
      // whichVectors_.size() == 0 means "constant stride."
      whichVectors_.clear ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,node_type> >& map,
               const dual_view_type& view,
               const dual_view_type& origView,
               const Teuchos::ArrayView<const size_t>& whichVectors) :
    base_type (map),
    view_ (view),
    origView_ (origView),
    whichVectors_ (whichVectors.begin (), whichVectors.end ())
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    const char tfecfFuncName[] = "MultiVector(map,view,origView,whichVectors): ";

    const size_t lclNumRows = this->getLocalLength ();
    // Check dimensions of the input DualView.  We accept that Kokkos
    // might not allow construction of a 0 x m (Dual)View with m > 0,
    // so we only require the number of rows to match if the
    // (Dual)View has more than zero columns.  Likewise, we only
    // require the number of columns to match if the (Dual)View has
    // more than zero rows.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      view.dimension_1 () != 0 && static_cast<size_t> (view.dimension_0 ()) < lclNumRows,
      std::invalid_argument, "view.dimension_0() = " << view.dimension_0 ()
      << " < map->getNodeNumElements() = " << lclNumRows << ".");
    if (whichVectors.size () != 0) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        view.dimension_1 () != 0 && view.dimension_1 () == 0,
        std::invalid_argument, "view.dimension_1() = 0, but whichVectors.size()"
        " = " << whichVectors.size () << " > 0.");
      size_t maxColInd = 0;
      typedef Teuchos::ArrayView<const size_t>::size_type size_type;
      for (size_type k = 0; k < whichVectors.size (); ++k) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          whichVectors[k] == Teuchos::OrdinalTraits<size_t>::invalid (),
          std::invalid_argument, "whichVectors[" << k << "] = "
          "Teuchos::OrdinalTraits<size_t>::invalid().");
        maxColInd = std::max (maxColInd, whichVectors[k]);
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        view.dimension_1 () != 0 && static_cast<size_t> (view.dimension_1 ()) <= maxColInd,
        std::invalid_argument, "view.dimension_1() = " << view.dimension_1 ()
        << " <= max(whichVectors) = " << maxColInd << ".");
    }
    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.dimension_1 () > 1) ? stride[1] :
      origView_.dimension_0 ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "Input DualView's column stride"
      " = " << LDA << " < this->getLocalLength() = " << lclNumRows << ".");

    if (whichVectors.size () == 1) {
      // If whichVectors has only one entry, we don't need to bother
      // with nonconstant stride.  Just shift the view over so it
      // points to the desired column.
      //
      // NOTE (mfh 10 May 2014) This is a special case where we set
      // origView_ just to view that one column, not all of the
      // original columns.  This ensures that the use of origView_ in
      // offsetView works correctly.
      const std::pair<size_t, size_t> colRng (whichVectors[0],
                                              whichVectors[0] + 1);
      view_ = subview (view_, ALL (), colRng);
      origView_ = subview (origView_, ALL (), colRng);
      // whichVectors_.size() == 0 means "constant stride."
      whichVectors_.clear ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,node_type> >& map,
               const Teuchos::ArrayView<const Scalar>& data,
               const size_t LDA,
               const size_t numVecs) :
    base_type (map)
  {
    using Kokkos::subview;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename dual_view_type::host_mirror_space HMS;
    typedef MakeUnmanagedView<const IST, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type in_view_type;
    typedef Kokkos::View<IST*, Kokkos::LayoutLeft, HMS> out_view_type;
    const char tfecfFuncName[] = "MultiVector(map,data,LDA,numVecs): ";

    // Deep copy constructor, constant stride (NO whichVectors_).
    // There is no need for a deep copy constructor with nonconstant stride.

    const size_t lclNumRows =
      map.is_null () ? size_t (0) : map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < lclNumRows, std::runtime_error,
      "LDA = " << LDA << " < numRows = " << lclNumRows << ".");

    view_ = allocDualView<Scalar, LO, GO, DeviceType> (lclNumRows, numVecs);
    view_.template modify<HMS> ();

    ArrayView<const IST> X_in_av = av_reinterpret_cast<const IST> (data);
    in_view_type X_in = view_getter_type::getView (X_in_av);
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    for (size_t j = 0; j < numVecs; ++j) {
      const std::pair<size_t, size_t> rng (j*LDA, j*LDA + lclNumRows);
      in_view_type X_j_in = subview (X_in, rng);
      out_view_type X_j_out = subview (view_.h_view, rowRng, j);
      Kokkos::deep_copy (X_j_out, X_j_in);
    }
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,node_type> >& map,
               const Teuchos::ArrayView<const ArrayView<const Scalar> >& ArrayOfPtrs,
               const size_t numVecs) :
    base_type (map)
  {
    using Kokkos::subview;
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename dual_view_type::host_mirror_space HMS;
    typedef MakeUnmanagedView<const IST, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type in_view_type;
    typedef Kokkos::View<IST*, Kokkos::LayoutLeft, HMS> out_view_type;
    const char tfecfFuncName[] = "MultiVector(map,ArrayOfPtrs,numVecs): ";

    const size_t lclNumRows =
      map.is_null () ? size_t (0) : map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs < 1 || numVecs != static_cast<size_t> (ArrayOfPtrs.size ()),
      std::runtime_error,
      "ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayView<const Scalar> X_j_av = ArrayOfPtrs[j];
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        static_cast<size_t> (X_j_av.size ()) < lclNumRows,
        std::invalid_argument, "ArrayOfPtrs[" << j << "].size() = "
        << X_j_av.size () << " < map->getNodeNumElements() = " << lclNumRows
        << ".");
    }

    view_ = allocDualView<Scalar, LO, GO, DeviceType> (lclNumRows, numVecs);
    view_.template modify<HMS> ();

    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayView<const IST> X_j_av =
        av_reinterpret_cast<const IST> (ArrayOfPtrs[j]);
      in_view_type X_j_in (X_j_av.getRawPtr (), lclNumRows);
      out_view_type X_j_out = subview (view_.h_view, rowRng, j);
      Kokkos::deep_copy (X_j_out, X_j_in);
    }
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  ~MultiVector () {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  isConstantStride () const {
    return whichVectors_.empty ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalLength () const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getNodeNumElements ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getGlobalLength () const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getGlobalNumElements ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getStride () const
  {
    if (isConstantStride ()) {
      // Get stride of view: if second dimension is 0, the
      // stride might be 0, so take view_dimension instead.
      size_t stride[8];
      origView_.stride (stride);
      const size_t LDA = (origView_.dimension_1 () > 1) ? stride[1] : origView_.dimension_0 ();
      return LDA;
    }
    else {
      return static_cast<size_t> (0);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  checkSizes (const SrcDistObject& sourceObj)
  {
    // Check whether the source object is a MultiVector.  If not, then
    // we can't even compare sizes, so it's definitely not OK to
    // Import or Export from it.
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> this_type;
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
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  constantNumberOfPackets () const {
    return this->getNumVectors ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  copyAndPermuteNew (
    const SrcDistObject& sourceObj,
    size_t numSameIDs,
    const Kokkos::View<const LocalOrdinal*, execution_space> &permuteToLIDs,
    const Kokkos::View<const LocalOrdinal*, execution_space> &permuteFromLIDs)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Kokkos::subview;
    typedef Kokkos::DualView<impl_scalar_type*,
      typename dual_view_type::array_layout,
      execution_space> col_dual_view_type;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;
    //typedef typename ArrayView<const LocalOrdinal>::size_type size_type; // unused
    const char tfecfFuncName[] = "copyAndPermute";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      ": permuteToLIDs and permuteFromLIDs must have the same size."
      << std::endl << "permuteToLIDs.size() = " << permuteToLIDs.size ()
      << " != permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);
    const size_t numCols = this->getNumVectors ();

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
    //
    // Note (ETP 2 Jul 2014)  We need to always copy one column at a
    // time, even when both multivectors are constant-stride, since
    // deep_copy between strided subviews with more than one column
    // doesn't currently work.
    if (numSameIDs > 0) {
      const std::pair<size_t, size_t> rows (0, numSameIDs);
      for (size_t j = 0; j < numCols; ++j) {
        const size_t dstCol = isConstantStride () ? j : whichVectors_[j];
        const size_t srcCol =
          sourceMV.isConstantStride () ? j : sourceMV.whichVectors_[j];
        col_dual_view_type dst_j =
          subview (view_, rows, dstCol);
        col_dual_view_type src_j =
          subview (sourceMV.view_, rows, srcCol);
        Kokkos::deep_copy (dst_j, src_j); // Copy src_j into dst_j
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
        getKokkosViewDeepCopy<execution_space> (whichVectors_ ()),
        getKokkosViewDeepCopy<execution_space> (sourceMV.whichVectors_ ()),
        numCols);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  packAndPrepareNew (const SrcDistObject& sourceObj,
                     const Kokkos::View<const local_ordinal_type*, execution_space> &exportLIDs,
                     Kokkos::View<impl_scalar_type*, execution_space> &exports,
                     const Kokkos::View<size_t*, execution_space> &numExportPacketsPerLID,
                     size_t& constantNumPackets,
                     Distributor & /* distor */ )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;
    //typedef Array<size_t>::size_type size_type; // unused

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
          getKokkosViewDeepCopy<execution_space> (sourceMV.whichVectors_ ()),
          numCols);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  unpackAndCombineNew (const Kokkos::View<const local_ordinal_type*, execution_space> &importLIDs,
                       const Kokkos::View<const impl_scalar_type*, execution_space> &imports,
                       const Kokkos::View<size_t*, execution_space> &numPacketsPerLID,
                       size_t constantNumPackets,
                       Distributor & /* distor */,
                       CombineMode CM)
  {
    using Teuchos::ArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    const char tfecfFuncName[] = "unpackAndCombine";

    // If we have no imports, there is nothing to do
    if (importLIDs.size () == 0) {
      return;
    }

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

      // NOTE (mfh 10 Mar 2012, 24 Mar 2014) If you want to implement
      // custom combine modes, start editing here.  Also, if you trust
      // inlining, it would be nice to condense this code by using a
      // binary function object f in the pack functors.
      if (CM == INSERT || CM == REPLACE) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::InsertOp(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<execution_space>(whichVectors_ ()),
            KokkosRefactor::Details::InsertOp(),
            numVecs);
        }
      }
      else if (CM == ADD) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::AddOp(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<execution_space>(whichVectors_ ()),
            KokkosRefactor::Details::AddOp(),
            numVecs);
        }
      }
      else if (CM == ABSMAX) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::AbsMaxOp(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<execution_space>(whichVectors_ ()),
            KokkosRefactor::Details::AbsMaxOp(),
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
  size_t
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getNumVectors () const
  {
    if (isConstantStride ()) {
      return static_cast<size_t> (view_.dimension_1 ());
    } else {
      return static_cast<size_t> (whichVectors_.size ());
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& A,
       const Teuchos::ArrayView<dot_type>& dots) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    // using Kokkos::LayoutLeft;
    // using Kokkos::LayoutRight;
    // using Kokkos::LayoutStride;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    // using Kokkos::Impl::if_c;
    // using Kokkos::Impl::is_same;
    // View of a MultiVector's local data (all columns).
    typedef typename dual_view_type::t_dev mv_view_type;
    //typedef typename mv_view_type::array_layout mv_array_layout;
    //
    // View of a single column of a MultiVector's local data.
    //
    // If mv_view_type has LayoutLeft, then vec_view_type's layout
    // should be LayoutLeft also; if mv_view_type has LayoutRight,
    // then vec_view_type's layout should be LayoutStride.
    //
    // FIXME (mfh 06 Jan 2015) What if the layout is neither?  The
    // code commented out below below sets the column layout to "void"
    // in that case, which is wrong because it would make the View's
    // layout the default.
    // typedef typename if_c<is_same<mv_array_layout, LayoutLeft>::value,
    //                       LayoutLeft,
    //                       typename if_c<is_same<mv_array_layout, LayoutRight>::value,
    //                                     LayoutStride,
    //                                     void>::type>::type col_array_layout;
    typedef Kokkos::LayoutLeft col_array_layout;
    typedef Kokkos::View<impl_scalar_type*, col_array_layout, execution_space> vec_view_type;
    typedef typename dual_view_type::host_mirror_space host_mirror_space;
    // View of all the dot product results.
    typedef MakeUnmanagedView<dot_type, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type host_dots_view_type;
    typedef Kokkos::View<dot_type*, Kokkos::LayoutLeft,
      host_mirror_space> host_dots_managed_view_type;
    const char tfecfFuncName[] = "Tpetra::MultiVector::dot: ";

#ifdef HAVE_TPETRA_DEBUG
    {
      const bool compat = this->getMap ()->isCompatible (* (A.getMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! compat, std::invalid_argument, "This MultiVector is not compatible "
        "with the input A.  We only test for this in a debug build.");
    }
#endif //  HAVE_TPETRA_DEBUG

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<dot_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<dot_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::runtime_error,
      ": MultiVectors do not have the same local length.  "
      "this->getLocalLength() = " << getLocalLength () << " != "
      "A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNumVectors () != A.getNumVectors (), std::runtime_error,
      ": MultiVectors must have the same number of columns (vectors).  "
      "this->getNumVectors() = " << getNumVectors () << " != "
      "A.getNumVectors() = " << A.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        static_cast<size_t> (dots.size () ) != getNumVectors (), std::runtime_error, ": The output "
      "array 'dots' must have the same number of entries as the number of "
      "columns (vectors) in *this and A.  dots.size() = " << dots.size ()
      << " != this->getNumVectors() = " << getNumVectors () << ".");

    // const size_t numVecs = getNumVectors (); // not used
    const size_t lclNumRows = getLocalLength ();
    const size_t numDots = static_cast<size_t> (dots.size ());

    // We're computing using the device's data, so we need to make
    // sure first that the device is in sync with the host.
    A.view_.template sync<DeviceType> ();
    view_.template sync<DeviceType> ();

    // In case the input dimensions don't match, make sure that we
    // don't overwrite memory that doesn't belong to us, by using
    // subset views with the minimum dimensions over all input.
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numDots);
    mv_view_type X = subview (view_.d_view, rowRng, colRng);
    mv_view_type Y = subview (A.view_.d_view, rowRng, colRng);

    if (numDots == 1) {
      // Special case 1: Both MultiVectors only have a single column.
      // The single-vector dot product kernel may be more efficient.
      const size_t ZERO = static_cast<size_t> (0);
      vec_view_type X_k = subview (X, ALL (), ZERO);
      vec_view_type Y_k = subview (Y, ALL (), ZERO);

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_k.dimension_0 () != lclNumRows, std::logic_error, "Tpetra::Multi"
        "Vector::dot: In the special single-vector case, X_k.dimension_0() = "
        << X_k.dimension_0 () << " != lclNumRows = " << lclNumRows
        << ".  Please report this bug to the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        Y_k.dimension_0 () != lclNumRows, std::logic_error, "Tpetra::Multi"
        "Vector::dot: In the special single-vector case, Y_k.dimension_0() = "
        << Y_k.dimension_0 () << " != lclNumRows = " << lclNumRows
        << ".  Please report this bug to the Tpetra developers.");

      dots[0] = Kokkos::V_Dot (X_k, Y_k, lclNumRows);
    }
    else if (isConstantStride () && A.isConstantStride ()) {
      // Special case 2: Both MultiVectors have constant stride.
      (void) Kokkos::MV_Dot (dots.getRawPtr (), X, Y, lclNumRows);
    }
    else {
      // FIXME (mfh 14 Jul 2014) This does a kernel launch for every
      // column.  It might be better to have a kernel that does the
      // work all at once.  On the other hand, we don't prioritize
      // performance of MultiVector views of noncontiguous columns.
      for (size_t k = 0; k < numDots; ++k) {
        const size_t X_col = isConstantStride () ? k : whichVectors_[k];
        const size_t Y_col = A.isConstantStride () ? k : A.whichVectors_[k];
        vec_view_type X_k = subview (X, ALL (), X_col);
        vec_view_type Y_k = subview (Y, ALL (), Y_col);
        dots[k] = Kokkos::V_Dot (X_k, Y_k, lclNumRows);
      }
    }

    // If the MultiVectors are distributed over multiple processes,
    // sum the results across processes.  We assume that the MPI
    // implementation can read from and write to device memory.
    //
    // replaceMap() may have removed some processes.  Those processes
    // have a null Map.  They must not participate in any collective
    // operations.  We ask first whether the Map is null, because
    // isDistributed() defers that question to the Map.  We still
    // compute and return local dot products for processes not
    // participating in collective operations; those probably don't
    // make any sense, but it doesn't hurt to do them, since it's
    // illegal to call dot() on those processes anyway.
    if (! this->getMap ().is_null () && this->isDistributed ()) {
      host_dots_view_type theDots (dots.getRawPtr (), numDots);
      // MPI doesn't allow aliasing of arguments, so we have to make a
      // copy of the local sum.
      host_dots_managed_view_type lclDots ("MV::dot lcl", numDots);

      TEUCHOS_TEST_FOR_EXCEPTION(
        lclDots.dimension_0 () != theDots.dimension_0 (), std::logic_error,
        "Tpetra::MultiVector::dot: lclDots and theDots have different sizes.  "
        "lclDots.dimension_0 () = " << lclDots.dimension_0 () << " != "
        "theDots.dimension_0 () = " << theDots.dimension_0 () << ".  "
        "Please report this bug to the Tpetra developers.");

      Kokkos::deep_copy (lclDots, theDots);
      const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
      const dot_type* const lclSum = lclDots.ptr_on_device ();
      dot_type* const gblSum = theDots.ptr_on_device ();
      reduceAll<int, dot_type> (comm, REDUCE_SUM, static_cast<int> (numDots),
                                lclSum, gblSum);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& A,
       const Kokkos::View<dot_type*, execution_space>& dots) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    // View of a MultiVector's local data (all columns).
    typedef typename dual_view_type::t_dev mv_view_type;
    // View of a single column of a MultiVector's local data.
    //
    // FIXME (mfh 14 Jul 2014) It would be better to get this typedef
    // from mv_view_type itself, in case the layout changes.
    typedef Kokkos::View<impl_scalar_type*, Kokkos::LayoutLeft, execution_space> vec_view_type;
    // View of all the dot product results.
    typedef Kokkos::View<dot_type*, execution_space> dots_view_type;
    // Scalar view; view of a single dot product result.
    typedef Kokkos::View<dot_type, execution_space> dot_view_type;
    const char tfecfFuncName[] = "Tpetra::MultiVector::dot";

#ifdef HAVE_TPETRA_DEBUG
    {
      const bool compat = this->getMap ()->isCompatible (* (A.getMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! compat, std::invalid_argument, "Tpetra::MultiVector::dot: *this is "
        "not compatible with the input MultiVector A.  We only test for this "
        "in a debug build.");
    }
#endif //  HAVE_TPETRA_DEBUG

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<dot_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<dot_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::runtime_error,
      ": MultiVectors do not have the same local length.  "
      "this->getLocalLength() = " << getLocalLength () << " != "
      "A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNumVectors () != A.getNumVectors (), std::runtime_error,
      ": MultiVectors must have the same number of columns (vectors).  "
      "this->getNumVectors() = " << getNumVectors () << " != "
      "A.getNumVectors() = " << A.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      dots.dimension_0 () != getNumVectors (), std::runtime_error,
      ": The output array 'dots' must have the same number of entries as the "
      "number of columns (vectors) in *this and A.  dots.dimension_0() = " <<
      dots.dimension_0 () << " != this->getNumVectors() = " << getNumVectors ()
      << ".");

    // We're computing using the device's data, so we need to make
    // sure first that the device is in sync with the host.
    A.view_.template sync<DeviceType> ();
    view_.template sync<DeviceType> ();

    // const size_t numVecs = getNumVectors (); // not used
    const size_t lclNumRows = getLocalLength ();
    const size_t numDots = dots.dimension_0 ();

    // In case the input dimensions don't match, make sure that we
    // don't overwrite memory that doesn't belong to us, by using
    // subset views with the minimum dimensions over all input.
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numDots);
    dots_view_type theDots = subview (dots, colRng);
    mv_view_type X = subview (view_.d_view, rowRng, colRng);
    mv_view_type Y = subview (A.view_.d_view, rowRng, colRng);

    // FIXME (mfh 14 Jul 2014) How come ALL() works as the first
    // argument, but not a row range?  The first line below doesn't
    // compile, but the second line does.  See
    // kokkos/core/unit_test/TestViewAPI.hpp, in particular
    // run_test_vector(), for an example of allowed subview arguments.
    //
    //vec_view_type X_0 = subview (X, rowRng, static_cast<size_t> (0));
    //vec_view_type X_0 = subview (X, ALL (), static_cast<size_t> (0));

    if (numDots == 1) {
      // Special case 1: Both MultiVectors only have a single column.
      // The single-vector dot product kernel may be more efficient.
      const size_t ZERO = static_cast<size_t> (0);
      vec_view_type X_k = subview (X, ALL (), ZERO);
      vec_view_type Y_k = subview (Y, ALL (), ZERO);
      dot_view_type dot_k = subview (theDots, ZERO);

      TEUCHOS_TEST_FOR_EXCEPTION(
        X_k.dimension_0 () != lclNumRows, std::logic_error, "Tpetra::Multi"
        "Vector::dot: In the special single-vector case, X_k.dimension_0() = "
        << X_k.dimension_0 () << " != lclNumRows = " << lclNumRows
        << ".  Please report this bug to the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        Y_k.dimension_0 () != lclNumRows, std::logic_error, "Tpetra::Multi"
        "Vector::dot: In the special single-vector case, Y_k.dimension_0() = "
        << Y_k.dimension_0 () << " != lclNumRows = " << lclNumRows
        << ".  Please report this bug to the Tpetra developers.");

      Kokkos::VecDotFunctor<vec_view_type> f (X_k, Y_k, dot_k);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else if (isConstantStride () && A.isConstantStride ()) {
      // Special case 2: Both MultiVectors have constant stride.
      Kokkos::MultiVecDotFunctor<mv_view_type> f (X, Y, theDots);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else {
      // FIXME (mfh 14 Jul 2014) This does a kernel launch for every
      // column.  It might be better to have a kernel that does the
      // work all at once.  On the other hand, we don't prioritize
      // performance of MultiVector views of noncontiguous columns.
      for (size_t k = 0; k < numDots; ++k) {
        const size_t X_col = isConstantStride () ? k : whichVectors_[k];
        const size_t Y_col = A.isConstantStride () ? k : A.whichVectors_[k];
        vec_view_type X_k = subview (X, ALL (), X_col);
        vec_view_type Y_k = subview (Y, ALL (), Y_col);
        dot_view_type dot_k = subview (theDots, k);
        Kokkos::VecDotFunctor<vec_view_type> f (X_k, Y_k, dot_k);
        Kokkos::parallel_reduce (lclNumRows, f);
      }
    }

    // If the MultiVectors are distributed over multiple processes,
    // sum the results across processes.  We assume that the MPI
    // implementation can read from and write to device memory.
    //
    // replaceMap() may have removed some processes.  Those processes
    // have a null Map.  They must not participate in any collective
    // operations.  We ask first whether the Map is null, because
    // isDistributed() defers that question to the Map.  We still
    // compute and return local dot products for processes not
    // participating in collective operations; those probably don't
    // make any sense, but it doesn't hurt to do them, since it's
    // illegal to call dot() on those processes anyway.
    if (! this->getMap ().is_null () && this->isDistributed ()) {
      // MPI doesn't allow aliasing of arguments, so we have to make a
      // copy of the local sum.
      dots_view_type lclDots ("MV::dot lcl", numDots);

      TEUCHOS_TEST_FOR_EXCEPTION(
        lclDots.dimension_0 () != theDots.dimension_0 (), std::logic_error,
        "Tpetra::MultiVector::dot: lclDots and theDots have different sizes.  "
        "lclDots.dimension_0 () = " << lclDots.dimension_0 () << " != "
        "theDots.dimension_0 () = " << theDots.dimension_0 () << ".  "
        "Please report this bug to the Tpetra developers.");

      Kokkos::deep_copy (lclDots, theDots);
      const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
      const dot_type* const lclSum = lclDots.ptr_on_device ();
      dot_type* const gblSum = theDots.ptr_on_device ();
      reduceAll<int, dot_type> (comm, REDUCE_SUM, static_cast<int> (numDots),
                                lclSum, gblSum);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  norm2 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef Kokkos::View<mag_type*, execution_space> dev_norms_view_type;
    typedef MakeUnmanagedView<mag_type, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type host_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    dev_norms_view_type normsDevView ("MV::norm2 tmp", numNorms);
    this->norm2 (normsDevView); // Do the computation on the device.
    Kokkos::deep_copy (normsHostView, normsDevView); // Bring back result to host
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  norm2 (const Kokkos::View<mag_type*, execution_space>& norms) const
  {
    this->normImpl (norms, NORM_TWO);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  normWeighted (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& weights,
                const Teuchos::ArrayView<mag_type> &norms) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;
    typedef Kokkos::Details::ArithTraits<mag_type> ATM;
    typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;
    const char tfecfFuncName[] = "normWeighted";

    const size_t lclNumRows = this->getLocalLength ();
    const size_t numVecs = this->getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (norms.size ()) != numVecs, std::runtime_error,
      ": norms.size() must be as large as the number of vectors in *this.");

    const bool OneW = (weights.getNumVectors () == 1);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! OneW && weights.getNumVectors () != numVecs, std::runtime_error,
      ": The input MultiVector of weights must contain either one column, "
      "or must have the same number of columns as *this.  "
      "weights.getNumVectors() = " << weights.getNumVectors ()
      << " and this->getNumVectors() = " << numVecs << ".");

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! this->getMap ()->isCompatible (*weights.getMap ()), std::runtime_error,
      ": MultiVectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << *this->getMap()
      << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != weights.getLocalLength (), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif

    // Unfortunately, we have to use a dot product function for
    // intermediate results.  This means that intermediate results
    // have type dot_type, not mag_type.  Thus, we need a temporary
    // array.  Once we finish the local part of the computation, we
    // can take magnitudes and do the MPI all-reduce over mag_type.
    typename Kokkos::View<dot_type*, DeviceType>::HostMirror lclDots ("lclDots", numVecs);

    view_.template sync<DeviceType> ();
    weights.view_.template sync<DeviceType> ();
    if (isConstantStride ()) {
      if (OneW) {
        view_type weights_0 =
          subview (weights.view_.d_view, ALL (), 0);
        Kokkos::MV_DotWeighted (lclDots.ptr_on_device (), weights_0,
                                view_.d_view, lclNumRows);
      } else
        Kokkos::MV_DotWeighted (lclDots.ptr_on_device (), weights.view_.d_view,
                                view_.d_view, lclNumRows);
    }
    else {
      // FIXME (mfh 11 Mar 2014) Once we have strided Views, we won't
      // have to write the explicit for loop over columns any more.
      if (OneW) {
        view_type weights_0 =
          subview (weights.view_.d_view, ALL (), 0);
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t curCol = whichVectors_[k];
          view_type vector_k = subview (view_.d_view, ALL (), curCol);
          lclDots(k) = Kokkos::V_DotWeighted (weights_0, vector_k, lclNumRows);
        }
      } else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t curCol = whichVectors_[k];
          view_type weights_k =
            subview (weights.view_.d_view, ALL (), curCol);
          view_type vector_k = subview (view_.d_view, ALL (), curCol);
          lclDots(k) = Kokkos::V_DotWeighted (weights_k, vector_k, lclNumRows);
        }
      }
    }

    const mag_type OneOverN =
      ATM::one () / static_cast<mag_type> (getGlobalLength ());

    if (this->isDistributed ()) {
      typename Kokkos::View<mag_type*, DeviceType>::HostMirror lclNorms ("lclNorms", numVecs);
      for (size_t k = 0; k < numVecs; ++k) {
        lclNorms(k) = ATS::abs (lclDots(k));
      }

      RCP<const Comm<int> > comm =
        this->getMap ().is_null () ? null : this->getMap ()->getComm ();
      if (comm.is_null ()) {
        for (size_t k = 0; k < numVecs; ++k) {
          norms[k] = ATM::zero ();
        }
      } else {
        reduceAll<int, mag_type> (*comm, REDUCE_SUM, static_cast<int> (numVecs),
                                  lclNorms.ptr_on_device (), norms.getRawPtr ());
        for (size_t k = 0; k < numVecs; ++k) {
          norms[k] = ATM::sqrt (norms[k] * OneOverN);
        }
      }
    }
    else {
      for (size_t k = 0; k < numVecs; ++k) {
        norms[k] = ATM::sqrt (ATS::magnitude (lclDots(k)) * OneOverN);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  norm1 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef Kokkos::View<mag_type*, execution_space> dev_norms_view_type;
    typedef MakeUnmanagedView<mag_type, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type host_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    dev_norms_view_type normsDevView ("MV::norm1 tmp", numNorms);
    this->norm1 (normsDevView); // Do the computation on the device.
    Kokkos::deep_copy (normsHostView, normsDevView); // Bring back result to host
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  norm1 (const Kokkos::View<mag_type*, execution_space>& norms) const
  {
    this->normImpl (norms, NORM_ONE);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  normInf (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef Kokkos::View<mag_type*, execution_space> dev_norms_view_type;
    typedef MakeUnmanagedView<mag_type, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type host_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    dev_norms_view_type normsDevView ("MV::normInf tmp", numNorms);
    this->normInf (normsDevView); // Do the computation on the device.
    Kokkos::deep_copy (normsHostView, normsDevView); // Bring back result to host
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  normInf (const Kokkos::View<mag_type*, execution_space>& norms) const
  {
    this->normImpl (norms, NORM_INF);
  }

  namespace { // (anonymous)

    //! Input argument for localNormImpl() (which see).
    enum EWhichNormImpl {
      IMPL_NORM_ONE, //<! Use the one-norm
      IMPL_NORM_TWO, //<! Use the two-norm
      IMPL_NORM_INF  //<! Use the infinity-norm
    };

    template<class RV, class XMV>
    void
    lclNormImpl (const RV& normsOut,
                 const XMV& X_lcl,
                 const size_t lclNumRows,
                 const size_t numVecs,
                 const Teuchos::ArrayView<const size_t>& whichVecs,
                 const bool constantStride,
                 const EWhichNormImpl whichNorm)
    {
      using Kokkos::ALL;
      using Kokkos::subview;
      typedef typename RV::non_const_value_type mag_type;

#ifdef KOKKOS_HAVE_CXX11
      static_assert (Kokkos::Impl::is_view<RV>::value,
                     "Tpetra::MultiVector::lclNormImpl: "
                     "The first argument RV is not a Kokkos::View.");
      static_assert (RV::rank == 1, "Tpetra::MultiVector::lclNormImpl: "
                     "The first argument normsOut must have rank 1.");
      static_assert (Kokkos::Impl::is_view<XMV>::value,
                     "Tpetra::MultiVector::lclNormImpl: "
                     "The second argument X_lcl is not a Kokkos::View.");
      static_assert (XMV::rank == 2, "Tpetra::MultiVector::lclNormImpl: "
                     "The second argument X_lcl must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

      // In case the input dimensions don't match, make sure that we
      // don't overwrite memory that doesn't belong to us, by using
      // subset views with the minimum dimensions over all input.
      const std::pair<size_t, size_t> rowRng (0, lclNumRows);
      const std::pair<size_t, size_t> colRng (0, numVecs);
      RV theNorms = subview (normsOut, colRng);
      XMV X = subview (X_lcl, rowRng, colRng);

      // mfh 10 Mar 2015: Kokkos::(Dual)View subviews don't quite
      // behave how you think when they have zero rows.  In that case,
      // it returns a 0 x 0 (Dual)View.
      TEUCHOS_TEST_FOR_EXCEPTION(
        lclNumRows != 0 && (X.dimension_0 () != lclNumRows || X.dimension_1 () != numVecs),
        std::logic_error, "X's dimensions are " << X.dimension_0 () << " x "
        << X.dimension_1 () << ", which differ from the local dimensions "
        << lclNumRows << " x " << numVecs << ".  Please report this bug to "
        "the Tpetra developers.");

      if (lclNumRows == 0) {
        const mag_type zeroMag = Kokkos::Details::ArithTraits<mag_type>::zero ();
        Kokkos::Impl::ViewFill<RV> (theNorms, zeroMag);
      }
      else { // lclNumRows != 0
        if (constantStride) {
          if (whichNorm == IMPL_NORM_INF) {
            KokkosBlas::nrmInf (theNorms, X);
          }
          else if (whichNorm == IMPL_NORM_ONE) {
            KokkosBlas::nrm1 (theNorms, X);
          }
          else if (whichNorm == IMPL_NORM_TWO) {
            KokkosBlas::nrm2_squared (theNorms, X);
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
          }
        }
        else { // not constant stride
          // NOTE (mfh 15 Jul 2014) This does a kernel launch for
          // every column.  It might be better to have a kernel that
          // does the work all at once.  On the other hand, we don't
          // prioritize performance of MultiVector views of
          // noncontiguous columns.
          for (size_t k = 0; k < numVecs; ++k) {
            const size_t X_col = constantStride ? k : whichVecs[k];
            if (whichNorm == IMPL_NORM_INF) {
              KokkosBlas::nrmInf (theNorms, k, X, X_col);
            }
            else if (whichNorm == IMPL_NORM_ONE) {
              KokkosBlas::nrm1 (theNorms, k, X, X_col);
            }
            else if (whichNorm == IMPL_NORM_TWO) {
              KokkosBlas::nrm2_squared (theNorms, k, X, X_col);
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
            }
          } // for each column
        } // constantStride
      } // lclNumRows != 0
    }

    template<class NormsViewType>
    void
    gblNormImpl (const NormsViewType& normsOut,
                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const bool distributed,
                 const EWhichNormImpl whichNorm)
    {
      using Teuchos::REDUCE_MAX;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      typedef typename NormsViewType::non_const_value_type mag_type;

      const size_t numVecs = normsOut.dimension_0 ();

      // If the MultiVector is distributed over multiple processes, do
      // the distributed (interprocess) part of the norm.  We assume
      // that the MPI implementation can read from and write to device
      // memory.
      //
      // replaceMap() may have removed some processes.  Those processes
      // have a null Map.  They must not participate in any collective
      // operations.  We ask first whether the Map is null, because
      // isDistributed() defers that question to the Map.  We still
      // compute and return local norms for processes not participating
      // in collective operations; those probably don't make any sense,
      // but it doesn't hurt to do them, since it's illegal to call
      // normInf() on those processes anyway.
      if (distributed && ! comm.is_null ()) {
        // The calling process only participates in the collective if
        // both the Map and its Comm on that process are nonnull.
        //
        // MPI doesn't allow aliasing of arguments, so we have to make
        // a copy of the local sum.
        NormsViewType lclNorms ("MV::normImpl lcl", numVecs);
        Kokkos::deep_copy (lclNorms, normsOut);
        const mag_type* const lclSum = lclNorms.ptr_on_device ();
        mag_type* const gblSum = normsOut.ptr_on_device ();
        const int nv = static_cast<int> (numVecs);
        if (whichNorm == IMPL_NORM_INF) {
          reduceAll<int, mag_type> (*comm, REDUCE_MAX, nv, lclSum, gblSum);
        } else {
          reduceAll<int, mag_type> (*comm, REDUCE_SUM, nv, lclSum, gblSum);
        }
      }

      if (whichNorm == IMPL_NORM_TWO) {
        // Replace the norm-squared results with their square roots in
        // place, to get the final output.  If the device memory and
        // the host memory are the same, it probably doesn't pay to
        // launch a parallel kernel for that, since there isn't enough
        // parallelism for the typical MultiVector case.
        const bool inHostMemory =
          Kokkos::Impl::is_same<typename NormsViewType::memory_space,
            typename NormsViewType::host_mirror_space::memory_space>::value;
        if (inHostMemory) {
          for (size_t j = 0; j < numVecs; ++j) {
            normsOut(j) = Kokkos::Details::ArithTraits<mag_type>::sqrt (normsOut(j));
          }
        }
        else {
          // There's not as much parallelism now, but that's OK.  The
          // point of doing parallel dispatch here is to keep the norm
          // results on the device, thus avoiding a copy to the host and
          // back again.
          Kokkos::SquareRootFunctor<NormsViewType> f (normsOut);
          Kokkos::parallel_for (numVecs, f);
        }
      }
    }

  } // namespace (anonymous)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  normImpl (const Kokkos::View<mag_type*, execution_space>& norms,
            const EWhichNorm whichNorm) const
  {
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    // View of all the norm results.
    typedef Kokkos::View<mag_type*, execution_space> RV;

    const size_t numVecs = this->getNumVectors ();
    if (numVecs == 0) {
      return; // nothing to do
    }
    const size_t lclNumRows = this->getLocalLength ();
    const size_t numNorms = static_cast<size_t> (norms.dimension_0 ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numNorms < numVecs, std::runtime_error, "Tpetra::MultiVector::normImpl: "
      "'norms' must have at least as many entries as the number of vectors in "
      "*this.  norms.dimension_0() = " << numNorms << " < this->getNumVectors()"
      " = " << numVecs << ".");

    const std::pair<size_t, size_t> colRng (0, numVecs);
    RV normsOut = subview (norms, colRng);

    EWhichNormImpl lclNormType;
    if (whichNorm == NORM_ONE) {
      lclNormType = IMPL_NORM_ONE;
    } else if (whichNorm == NORM_TWO) {
      lclNormType = IMPL_NORM_TWO;
    } else {
      lclNormType = IMPL_NORM_INF;
    }

    RCP<const Comm<int> > comm = this->getMap ().is_null () ? null :
      this->getMap ()->getComm ();

    // FIXME (mfh 05 Mar 2015) DualView flags are not indicative when
    // the two memory spaces are the same, so we check the latter.
    const bool oneMemorySpace =
      Kokkos::Impl::is_same<typename dual_view_type::t_dev::memory_space,
                            typename dual_view_type::t_host::memory_space>::value;
    if (! oneMemorySpace && view_.modified_host >= view_.modified_device) {
      // DualView was last modified on host, so run the local kernel there.
      // This means we need a host mirror of the array of norms too.
      typedef typename dual_view_type::t_host XMV;
      lclNormImpl<RV, XMV> (normsOut, view_.h_view, lclNumRows, numVecs,
                            this->whichVectors_, this->isConstantStride (),
                            lclNormType);
      typename RV::HostMirror normsOutHost =
        Kokkos::create_mirror_view (normsOut);
      gblNormImpl<typename RV::HostMirror> (normsOutHost, comm,
                                            this->isDistributed (),
                                            lclNormType);
      Kokkos::deep_copy (normsOut, normsOutHost);
    }
    else {
      // DualView was last modified on device, so run the local kernel there.
      typedef typename dual_view_type::t_dev XMV;
      lclNormImpl<RV, XMV> (normsOut, view_.d_view, lclNumRows, numVecs,
                            this->whichVectors_, this->isConstantStride (),
                            lclNormType);
      gblNormImpl<RV> (normsOut, comm, this->isDistributed (), lclNormType);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  meanValue (const Teuchos::ArrayView<impl_scalar_type>& means) const
  {
    // KR FIXME Overload this method to take a View.

    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Array;
    using Teuchos::arcp_const_cast;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;

    const size_t lclNumRows = getLocalLength ();
    const size_t numVecs = getNumVectors ();
    const size_t numMeans = static_cast<size_t> (means.size ());

    TEUCHOS_TEST_FOR_EXCEPTION(
      numMeans != numVecs, std::runtime_error,
      "Tpetra::MultiVector::meanValue: means.size() = " << numMeans
      << " != this->getNumVectors() = " << numVecs << ".");

    // compute local components of the means
    // sum these across all nodes
    view_.template sync<DeviceType> ();
    if (isConstantStride ()) {
      Kokkos::MV_Sum (means.getRawPtr (), view_.d_view, lclNumRows);
    }
    else {
      const std::pair<size_t, size_t> rowRng (0, lclNumRows);
      for (size_t j = 0; j < numVecs; ++j) {
        typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;
        const size_t col = whichVectors_[j];
        view_type X_col = subview (view_.d_view, rowRng, col);
        means[j] = Kokkos::V_Sum (X_col);
      }
    }
    if (this->isDistributed ()) {
      Teuchos::Array<impl_scalar_type> lmeans (means);
      // only combine if we are a distributed MV
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM,
                 static_cast<int> (numVecs), lmeans.getRawPtr (),
                 means.getRawPtr ());
    }
    // mfh 12 Apr 2012: Don't take out the cast from the ordinal type
    // to the magnitude type, since operator/ (std::complex<T>, int)
    // isn't necessarily defined.
    const impl_scalar_type OneOverN =
      ATS::one () / static_cast<mag_type> (getGlobalLength ());
    for (size_t k = 0; k < numMeans; ++k) {
      means[k] = means[k] * OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  randomize ()
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef impl_scalar_type IST;
    typedef Kokkos::Details::ArithTraits<IST> ATS;
    typedef Kokkos::Random_XorShift64_Pool<DeviceType> pool_type;
    typedef typename pool_type::generator_type generator_type;

    // Seed the pseudorandom number generator using the calling
    // process' rank.  This helps decorrelate different process'
    // pseudorandom streams.  It's not perfect but it's effective and
    // doesn't require MPI communication.  The seed also includes bits
    // from the standard library's rand().
    //
    // FIXME (mfh 07 Jan 2015) Should we save the seed for later use?
    // The code below just makes a new seed each time.

    const uint64_t myRank =
      static_cast<uint64_t> (this->getMap ()->getComm ()->getRank ());
    uint64_t seed64 = static_cast<uint64_t> (std::rand ()) + myRank + 17311uLL;
    unsigned int seed = static_cast<unsigned int> (seed64&0xffffffff);

    pool_type rand_pool (seed);
    const IST max = Kokkos::rand<generator_type, IST>::max ();
    const IST min = ATS::is_signed ? IST (-max) : ATS::zero ();

    if (isConstantStride ()) {
      Kokkos::fill_random (view_.d_view, rand_pool, min, max);
      view_.template modify<DeviceType> ();
    }
    else {
      const size_t numVecs = getNumVectors ();
      view_.template sync<DeviceType> ();
      typedef Kokkos::View<IST*, DeviceType> view_type;
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t col = whichVectors_[k];
        view_type X_k = subview (view_.d_view, ALL (), col);
        Kokkos::fill_random (X_k, rand_pool, min, max);
      }
      view_.template modify<DeviceType> ();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  putScalar (const Scalar& alpha)
  {
    using Kokkos::ALL;
    using Kokkos::Impl::ViewFill;
    using Kokkos::subview;
    typedef typename dual_view_type::t_dev::device_type DMS;
    typedef typename dual_view_type::t_host::device_type HMS;

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const size_t lclNumRows = getLocalLength ();
    const size_t numVecs = getNumVectors ();
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // Modify the most recently updated version of the data.  This
    // avoids sync'ing, which could violate users' expectations.
    if (view_.modified_device >= view_.modified_host) {
      //
      // Last modified in device memory, so modify data there.
      //
      // Type of the device memory View of the MultiVector's data.
      typedef typename dual_view_type::t_dev mv_view_type;
      // Type of a View of a single column of the MultiVector's data.
      typedef Kokkos::View<impl_scalar_type*,
        typename mv_view_type::array_layout, DMS> vec_view_type;

      this->template modify<DMS> (); // we are about to modify on the device
      mv_view_type X =
        subview (this->getDualView ().template view<DMS> (),
                               rowRng, colRng);
      if (numVecs == 1) {
        vec_view_type X_0 =
          subview (X, ALL (), static_cast<size_t> (0));
        // The constructor of ViewFill invokes the functor in
        // parallel, so we don't have to call parallel_for ourselves.
        ViewFill<vec_view_type> vf (X_0, theAlpha);
      }
      else if (isConstantStride ()) {
        ViewFill<mv_view_type> vf (X, theAlpha);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t col = whichVectors_[k];
          vec_view_type X_k = subview (X, ALL (), col);
          ViewFill<vec_view_type> vf (X_k, theAlpha);
        }
      }
    }
    else { // last modified in host memory, so modify data there.
      typedef typename dual_view_type::t_host mv_view_type;
      typedef Kokkos::View<impl_scalar_type*,
        typename mv_view_type::array_layout, HMS> vec_view_type;

      this->template modify<HMS> (); // we are about to modify on the host
      mv_view_type X =
        subview (this->getDualView ().template view<HMS> (),
                               rowRng, colRng);
      if (numVecs == 1) {
        vec_view_type X_0 =
          subview (X, ALL (), static_cast<size_t> (0));
        // The constructor of ViewFill invokes the functor in
        // parallel, so we don't have to call parallel_for ourselves.
        ViewFill<vec_view_type> vf (X_0, theAlpha);
      }
      else if (isConstantStride ()) {
        ViewFill<mv_view_type> vf (X, theAlpha);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t col = whichVectors_[k];
          vec_view_type X_k = subview (X, ALL (), col);
          ViewFill<vec_view_type> vf (X_k, theAlpha);
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceMap (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,node_type> >& newMap)
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
    // local DualView if it does not have the right number of rows.
    // If the new number of rows is nonzero, we'll fill the newly
    // allocated local data with zeros, as befits a projection
    // operation.
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
      // If this->getMap() is null, that means that this MultiVector
      // has already had replaceMap happen to it.  In that case, just
      // reallocate the DualView with the right size.

      TEUCHOS_TEST_FOR_EXCEPTION(
        newMap.is_null (), std::invalid_argument,
        "Tpetra::MultiVector::replaceMap: both current and new Maps are null.  "
        "This probably means that the input Map is incorrect.");

      // Case 3: current Map is null, new Map is nonnull.
      // Reallocate the DualView with the right dimensions.
      const size_t newNumRows = newMap->getNodeNumElements ();
      const size_t origNumRows = view_.dimension_0 ();
      const size_t numCols = this->getNumVectors ();

      if (origNumRows != newNumRows || view_.dimension_1 () != numCols) {
        view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, DeviceType> (newNumRows, numCols);
      }
    }
    else if (newMap.is_null ()) { // Case 2: current Map is nonnull, new Map is null
      // I am an excluded process.  Reinitialize my data so that I
      // have 0 rows.  Keep the number of columns as before.
      const size_t newNumRows = static_cast<size_t> (0);
      const size_t numCols = this->getNumVectors ();
      view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, DeviceType> (newNumRows, numCols);
    }

    this->map_ = newMap;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  scale (const Scalar& alpha)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;
    typedef typename dual_view_type::t_dev::device_type DMS;
    typedef typename dual_view_type::t_host::device_type HMS;

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const size_t lclNumRows = getLocalLength ();
    const size_t numVecs = getNumVectors ();
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // NOTE: can't substitute putScalar(0.0) for scale(0.0), because
    //       the former will overwrite NaNs present in the
    //       MultiVector, while the semantics of this call require
    //       multiplying them by 0, which IEEE requires to be NaN

    if (theAlpha == ATS::one ()) {
      return; // do nothing
    }

    // Modify the most recently updated version of the data.  This
    // avoids sync'ing, which could violate users' expectations.
    if (view_.modified_device >= view_.modified_host) {
      //
      // Last modified in device memory, so modify data there.
      //
      // Type of the device memory View of the MultiVector's data.
      typedef typename dual_view_type::t_dev mv_view_type;
      // Type of a View of a single column of the MultiVector's data.
      typedef Kokkos::View<impl_scalar_type*,
        typename mv_view_type::array_layout, DMS> vec_view_type;

      this->template modify<DMS> (); // we are about to modify on the device
      mv_view_type X =
        subview (this->getDualView ().template view<DMS> (),
                               rowRng, colRng);
      if (numVecs == 1) {
        vec_view_type X_0 =
          subview (X, ALL (), static_cast<size_t> (0));
        Kokkos::V_MulScalar (X_0, theAlpha, X_0);
      }
      else if (isConstantStride ()) {
        Kokkos::MV_MulScalar (X, theAlpha, X);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t col = whichVectors_[k];
          vec_view_type X_col = subview (X, ALL (), col);
          Kokkos::V_MulScalar (X_col, theAlpha, X_col);
        }
      }
    }
    else { // last modified in host memory, so modify data there.
      typedef typename dual_view_type::t_host mv_view_type;
      typedef Kokkos::View<impl_scalar_type*,
        typename mv_view_type::array_layout, HMS> vec_view_type;

      this->template modify<HMS> (); // we are about to modify on the host
      mv_view_type X =
        subview (this->getDualView ().template view<HMS> (),
                               rowRng, colRng);
      if (numVecs == 1) {
        vec_view_type X_0 =
          subview (X, ALL (), static_cast<size_t> (0));
        Kokkos::V_MulScalar (X_0, theAlpha, X_0);
      }
      else if (isConstantStride ()) {
        Kokkos::MV_MulScalar (X, theAlpha, X);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t col = whichVectors_[k];
          vec_view_type X_col = subview (X, ALL (), col);
          Kokkos::V_MulScalar (X_col, theAlpha, X_col);
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  scale (Teuchos::ArrayView<const Scalar> alphas)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const size_t numVecs = this->getNumVectors ();
    const size_t numAlphas = static_cast<size_t> (alphas.size ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numAlphas != numVecs, std::invalid_argument, "Tpetra::MultiVector::scale("
      "alphas): alphas.size() = " << numAlphas << " != this->getNumVectors() = "
      << numVecs << ".");

    const size_t myLen = view_.dimension_0 ();
    if (myLen == 0) {
      return;
    }

    if (isConstantStride ()) {
      // Use a DualView to copy the scaling constants onto the device.
      typedef Kokkos::DualView<impl_scalar_type*, execution_space> k_alphas_type ;
      k_alphas_type k_alphas ("alphas::tmp", numAlphas);
      k_alphas.template modify<typename k_alphas_type::host_mirror_space> ();
      for (size_t i=0; i < numAlphas; ++i) {
        k_alphas.h_view(i) = static_cast<impl_scalar_type> (alphas[i]);
      }
      k_alphas.template sync<execution_space> ();

      // Modify the MultiVector on the device.
      view_.template sync<DeviceType> ();
      view_.template modify<DeviceType> ();
      Kokkos::MV_MulScalar (view_.d_view, k_alphas.d_view, view_.d_view);
    }
    else {
      typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t col = isConstantStride () ? k : whichVectors_[k];
        const impl_scalar_type alpha_k = static_cast<impl_scalar_type> (alphas[k]);
        view_type vector_k = subview (view_.d_view, ALL (), col);
        Kokkos::V_MulScalar (vector_k, alpha_k, vector_k);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  scale (const Kokkos::View<const impl_scalar_type*, execution_space> alphas)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<impl_scalar_type*, execution_space> view_type;
    typedef typename view_type::HostMirror host_view_type;

    const size_t numVecs = this->getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (alphas.dimension_0 ()) != numVecs,
      std::invalid_argument, "Tpetra::MultiVector::scale(alphas): "
      "alphas.dimension_0() = " << alphas.dimension_0 ()
      << " != this->getNumVectors () = " << numVecs << ".");
    if (this->getLocalLength () == 0) {
      return;
    }

    if (this->isConstantStride ()) {
      view_.template sync<execution_space> ();
      view_.template modify<execution_space> ();
      Kokkos::MV_MulScalar (view_.d_view, alphas, view_.d_view);
    }
    else {
      host_view_type h_alphas = Kokkos::create_mirror_view (alphas);
      Kokkos::deep_copy (h_alphas, alphas);
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t curCol = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = subview (view_.d_view, ALL (), curCol);
        Kokkos::V_MulScalar (vector_k, alphas(k), vector_k);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  scale (const Scalar& alpha,
         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A)
  {
    const char tfecfFuncName[] = "scale(alpha,A): ";
    const size_t numVecs = getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
       getLocalLength () != A.getLocalLength (), std::runtime_error,
       "MultiVectors do not have the same local length.  "
       "this->getLocalLength() = " << getLocalLength ()
       << " != A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != numVecs, std::runtime_error,
      "MultiVectors do not have the same number of columns (vectors).  "
       "this->getNumVectors() = " << getNumVectors ()
       << " != A.getNumVectors() = " << A.getNumVectors () << ".");

    // FIXME (mfh 07 Jan 2015) We shouldn't call modify() or sync() on
    // the input MultiVector, because it's const.  We should instead
    // move _this_ MultiVector's data to the space where A's data have
    // been most recently updated.

    if (isConstantStride () && A.isConstantStride ()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_MulScalar (view_.d_view, alpha, A.view_.d_view);
    }
    else {
      using Kokkos::ALL;
      using Kokkos::subview;
      typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;

      view_.template sync<DeviceType> ();
      view_.template modify<DeviceType> ();
      A.view_.template sync<DeviceType> ();
      A.view_.template modify<DeviceType> ();
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = subview (view_.d_view, ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        view_type vector_Ak = subview (A.view_.d_view, ALL (), A_col);
        Kokkos::V_MulScalar (vector_k, alpha, vector_Ak);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  reciprocal (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A)
  {
    const char tfecfFuncName[] = "reciprocal";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
       getLocalLength () != A.getLocalLength (), std::runtime_error,
       ": MultiVectors do not have the same local length.  "
       "this->getLocalLength() = " << getLocalLength ()
       << " != A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != this->getNumVectors (), std::runtime_error,
      ": MultiVectors do not have the same number of columns (vectors).  "
       "this->getNumVectors() = " << getNumVectors ()
       << " != A.getNumVectors() = " << A.getNumVectors () << ".");

    // FIXME (mfh 07 Jan 2015) See note on two-argument scale() above.

    const size_t numVecs = getNumVectors ();
    try {
      if (isConstantStride () && A.isConstantStride ()) {
        view_.template sync<DeviceType> ();
        view_.template modify<DeviceType> ();
        Kokkos::MV_Reciprocal (view_.d_view, A.view_.d_view);
      }
      else {
        using Kokkos::ALL;
        using Kokkos::subview;
        typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;

        view_.template sync<DeviceType> ();
        view_.template modify<DeviceType> ();

        // FIXME (mfh 23 Jul 2014) I'm not sure if it should be our
        // responsibility to sync A.
        A.view_.template sync<DeviceType> ();
        A.view_.template modify<DeviceType> ();

        for (size_t k = 0; k < numVecs; ++k) {
          const size_t this_col = isConstantStride () ? k : whichVectors_[k];
          view_type vector_k = subview (view_.d_view, ALL (), this_col);
          const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
          view_type vector_Ak = subview (A.view_.d_view, ALL (), A_col);
          Kokkos::V_Reciprocal(vector_k, vector_Ak);
        }
      }
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error,
        ": Caught exception from Kokkos: " << e.what () << std::endl);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  abs (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A)
  {
    const char tfecfFuncName[] = "abs";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
       getLocalLength () != A.getLocalLength (), std::runtime_error,
       ": MultiVectors do not have the same local length.  "
       "this->getLocalLength() = " << getLocalLength ()
       << " != A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != this->getNumVectors (), std::runtime_error,
      ": MultiVectors do not have the same number of columns (vectors).  "
       "this->getNumVectors() = " << getNumVectors ()
       << " != A.getNumVectors() = " << A.getNumVectors () << ".");

    // FIXME (mfh 07 Jan 2015) See note on two-argument scale() above.

    const size_t numVecs = getNumVectors ();
    if (isConstantStride () && A.isConstantStride ()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_Abs(view_.d_view,A.view_.d_view);
    }
    else {
      using Kokkos::ALL;
      using Kokkos::subview;
      typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;

      view_.template sync<DeviceType> ();
      view_.template modify<DeviceType> ();
      A.view_.template sync<DeviceType> ();
      A.view_.template modify<DeviceType> ();

      for (size_t k=0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = subview (view_.d_view, ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        view_type vector_Ak = subview (A.view_.d_view, ALL (), A_col);
        Kokkos::V_Abs(vector_k, vector_Ak);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  update (const Scalar& alpha,
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
          const Scalar& beta)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<impl_scalar_type*, DeviceType> view_type;
    const char tfecfFuncName[] = "update: ";

    // FIXME (mfh 07 Jan 2015) See note on two-argument scale() above.

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::invalid_argument,
      "The input MultiVector A has " << A.getLocalLength () << " local "
      "row(s), but this MultiVector has " << getLocalLength () << " local "
      "row(s).");
    const size_t numVecs = getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != numVecs, std::invalid_argument,
      "The input MultiVector A has " << A.getNumVectors () << " column(s), "
      "but this MultiVector has " << numVecs << " column(s).");

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const impl_scalar_type theBeta = static_cast<impl_scalar_type> (beta);
    const size_t lclNumRows = getLocalLength ();
    if (isConstantStride () && A.isConstantStride ()) {
      Kokkos::MV_Add (view_.d_view, theAlpha, A.view_.d_view, theBeta,
                      view_.d_view, lclNumRows);
    }
    else {
      // Make sure that Kokkos only uses the local length for add.
      const std::pair<size_t, size_t> rowRng (0, lclNumRows);
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        const size_t A_col = A.isConstantStride () ? k : A.whichVectors_[k];
        view_type this_colView =
          subview (view_.d_view, rowRng, this_col);
        view_type A_colView =
          subview (A.view_.d_view, rowRng, A_col);
        Kokkos::V_Add (this_colView, theAlpha, A_colView,
                       theBeta, this_colView);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  update (const Scalar& alpha,
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
          const Scalar& beta,
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& B,
          const Scalar& gamma)
  {
    using Kokkos::ALL;
    using Kokkos::V_Add;
    using Kokkos::subview;
    const char tfecfFuncName[] = "update(alpha,A,beta,B,gamma): ";

    // FIXME (mfh 07 Jan 2015) See note on two-argument scale() above.

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::invalid_argument,
      "The input MultiVector A has " << A.getLocalLength () << " local "
      "row(s), but this MultiVector has " << getLocalLength () << " local "
      "row(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != B.getLocalLength (), std::invalid_argument,
      "The input MultiVector B has " << B.getLocalLength () << " local "
      "row(s), but this MultiVector has " << getLocalLength () << " local "
      "row(s).");
    const size_t numVecs = getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != numVecs, std::invalid_argument,
      "The input MultiVector A has " << A.getNumVectors () << " column(s), "
      "but this MultiVector has " << numVecs << " column(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      B.getNumVectors () != numVecs, std::invalid_argument,
      "The input MultiVector B has " << B.getNumVectors () << " column(s), "
      "but this MultiVector has " << numVecs << " column(s).");

    const impl_scalar_type zero = Kokkos::Details::ArithTraits<impl_scalar_type>::zero ();
    const impl_scalar_type one = Kokkos::Details::ArithTraits<impl_scalar_type>::one ();
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const impl_scalar_type theBeta = static_cast<impl_scalar_type> (beta);
    const impl_scalar_type theGamma = static_cast<impl_scalar_type> (gamma);

    if (isConstantStride() && A.isConstantStride () && B.isConstantStride ()) {
      if (theGamma == zero) {
        Kokkos::MV_Add (view_.d_view, theAlpha, A.view_.d_view, theBeta,
                        B.view_.d_view);
      } else {
        Kokkos::MV_Add (view_.d_view, theAlpha, A.view_.d_view, theGamma,
                        view_.d_view);
        Kokkos::MV_Add (view_.d_view, theBeta, B.view_.d_view, one,
                        view_.d_view);
      }
    } else {
      // Some input (or *this) is not constant stride,
      // so perform the update one column at a time.

      // Make sure that Kokkos only uses the local length for add.
      const size_t lclNumRows = this->getLocalLength ();
      const std::pair<size_t, size_t> rowRng (0, lclNumRows);

      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        const size_t A_col = A.isConstantStride () ? k : A.whichVectors_[k];
        const size_t B_col = B.isConstantStride () ? k : B.whichVectors_[k];
        if (theGamma == zero) {
          // TODO: make sure it only uses LocalLength for add.
          V_Add (subview (view_.d_view, rowRng, this_col),
                 theAlpha,
                 subview (A.view_.d_view, rowRng, A_col),
                 theBeta,
                 subview (B.view_.d_view, rowRng, B_col));
        } else {
          V_Add (subview (view_.d_view, rowRng, this_col),
                 theAlpha,
                 subview (A.view_.d_view, rowRng, A_col),
                 theGamma,
                 subview (view_.d_view, rowRng, this_col));
          V_Add (subview (view_.d_view, rowRng, this_col),
                 theBeta,
                 subview (B.view_.d_view, rowRng, B_col),
                 one,
                 subview (view_.d_view, rowRng, this_col));
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getData (size_t j) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef typename dual_view_type::host_mirror_space host_type;
    typedef typename dual_view_type::t_host host_view_type;

    // Any MultiVector method that called the (classic) Kokkos Node's
    // viewBuffer or viewBufferNonConst methods always implied a
    // device->host synchronization.  Thus, we synchronize here as
    // well.
    view_.template sync<host_type> ();

    // Get a host view of the entire MultiVector's data.
    host_view_type hostView = view_.template view<host_type> ();
    // Get a subview of column j.
    host_view_type hostView_j;
    if (isConstantStride ()) {
      hostView_j = subview (hostView, ALL (), Kokkos::pair<int,int>(j,j+1));
    } else {
      hostView_j = subview (hostView, ALL (), Kokkos::pair<int,int>(whichVectors_[j],whichVectors_[j]+1));
    }

    // Wrap up the subview of column j in an ArrayRCP<const impl_scalar_type>.
    Teuchos::ArrayRCP<const impl_scalar_type> dataAsArcp =
      Kokkos::Compat::persistingView (hostView_j, 0, getLocalLength ());

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t>(hostView_j.dimension_0 ()) < static_cast<size_t>(dataAsArcp.size ()), std::logic_error,
      "Tpetra::MultiVector::getData: hostView_j.dimension_0() = "
      << hostView_j.dimension_0 () << " < dataAsArcp.size() = "
      << dataAsArcp.size () << ".  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    return Teuchos::arcp_reinterpret_cast<const Scalar> (dataAsArcp);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,Kokkos
    ::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getDataNonConst (size_t j)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef typename dual_view_type::host_mirror_space host_type;
    typedef typename dual_view_type::t_host host_view_type;

    // Any MultiVector method that called the (classic) Kokkos Node's
    // viewBuffer or viewBufferNonConst methods always implied a
    // device->host synchronization.  Thus, we synchronize here as
    // well.
    view_.template sync<host_type> ();

    // Get a host view of the entire MultiVector's data.
    host_view_type hostView = view_.template view<host_type> ();
    // Get a subview of column j.
    host_view_type hostView_j;
    if (isConstantStride ()) {
      hostView_j = subview (hostView, ALL (), Kokkos::pair<int,int>(j,j+1));
    } else {
      hostView_j = subview (hostView, ALL (), Kokkos::pair<int,int>(whichVectors_[j],whichVectors_[j]+1));
    }

    // Calling getDataNonConst() implies that the user plans to modify
    // the values in the MultiVector, so we call modify on the view
    // here.
    view_.template modify<host_type> ();

    // Wrap up the subview of column j in an ArrayRCP<const impl_scalar_type>.
    Teuchos::ArrayRCP<impl_scalar_type> dataAsArcp =
      Kokkos::Compat::persistingView (hostView_j, 0, getLocalLength ());

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t>(hostView_j.dimension_0 ()) < static_cast<size_t>(dataAsArcp.size ()), std::logic_error,
      "Tpetra::MultiVector::getDataNonConst: hostView_j.dimension_0() = "
      << hostView_j.dimension_0 () << " < dataAsArcp.size() = "
      << dataAsArcp.size () << ".  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    return Teuchos::arcp_reinterpret_cast<Scalar> (dataAsArcp);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>&
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  operator= (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type, false>& source)
  {
    if (this != &source) {
      base_type::operator= (source);
      //
      // operator= implements view semantics (shallow copy).
      //

      // Kokkos::View operator= also implements view semantics.
      view_ = source.view_;
      origView_ = source.origView_;

      // NOTE (mfh 24 Mar 2014) Christian wrote here that assigning
      // whichVectors_ is "probably not ok" (probably constitutes deep
      // copy).  I would say that it's OK, because whichVectors_ is
      // immutable (from the user's perspective); it's analogous to
      // the dimensions or stride.  Once we make whichVectors_ a
      // Kokkos::View instead of a Teuchos::Array, all debate will go
      // away and we will unquestionably have view semantics.
      whichVectors_ = source.whichVectors_;
    }
    return *this;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  subCopy (const Teuchos::ArrayView<const size_t>& cols) const
  {
    using Teuchos::RCP;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;

    // Check whether the index set in cols is contiguous.  If it is,
    // use the more efficient Range1D version of subCopy.
    bool contiguous = true;
    const size_t numCopyVecs = static_cast<size_t> (cols.size ());
    for (size_t j = 1; j < numCopyVecs; ++j) {
      if (cols[j] != cols[j-1] + static_cast<size_t> (1)) {
        contiguous = false;
        break;
      }
    }
    if (contiguous && numCopyVecs > 0) {
      return this->subCopy (Teuchos::Range1D (cols[0], cols[numCopyVecs-1]));
    }
    else {
      RCP<const MV> X_sub = this->subView (cols);
      RCP<MV> Y (new MV (this->getMap (), numCopyVecs, false));
      Y->assign (*X_sub);
      return Y;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  subCopy (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::RCP;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;

    RCP<const MV> X_sub = this->subView (colRng);
    RCP<MV> Y (new MV (this->getMap (), static_cast<size_t> (colRng.size ()), false));
    Y->assign (*X_sub);
    return Y;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getOrigNumLocalRows () const {
    return origView_.dimension_0 ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getOrigNumLocalCols () const {
    return origView_.dimension_1 ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    const size_t newNumRows = subMap->getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > this->getOrigNumLocalRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > this->getLocalLength (), std::runtime_error,
        "Tpetra::MultiVector::offsetView(NonConst): Invalid input Map.  The "
        "input Map owns " << newNumRows << " entries on process " << myRank <<
        ".  offset = " << offset << ".  Yet, the MultiVector contains only "
        << this->getOrigNumLocalRows () << " rows on this process.");
    }

#ifdef HAVE_TPETRA_DEBUG
    const size_t strideBefore = this->isConstantStride () ?
      this->getStride () :
      static_cast<size_t> (0);
    const size_t lclNumRowsBefore = this->getLocalLength ();
    const size_t numColsBefore = this->getNumVectors ();
    const impl_scalar_type* hostPtrBefore =
      this->getDualView ().h_view.ptr_on_device ();
#endif // HAVE_TPETRA_DEBUG

    const std::pair<size_t, size_t> rowRng (offset, offset + newNumRows);
    // FIXME (mfh 10 May 2014) Use of origView_ instead of view_ for
    // the second argument may be wrong, if view_ resulted from a
    // previous call to offsetView with offset != 0.
    dual_view_type newView =
      subview (origView_, rowRng, ALL ());
    // NOTE (mfh 06 Jan 2015) Work-around to deal with Kokkos not
    // handling subviews of degenerate Views quite so well.  For some
    // reason, the ([0,0], [0,2]) subview of a 0 x 2 DualView is 0 x
    // 0.  We work around by creating a new empty DualView of the
    // desired (degenerate) dimensions.
    if (newView.dimension_0 () == 0 &&
        newView.dimension_1 () != view_.dimension_1 ()) {
      newView = allocDualView<Scalar,
                              LocalOrdinal,
                              GlobalOrdinal,
                              DeviceType> (size_t (0),
                                           this->getNumVectors ());
    }
    RCP<const MV> subViewMV;
    if (isConstantStride ()) {
      subViewMV = rcp (new MV (subMap, newView, origView_));
    } else {
      subViewMV = rcp (new MV (subMap, newView, origView_, whichVectors_ ()));
    }

#ifdef HAVE_TPETRA_DEBUG
    const size_t strideAfter = this->isConstantStride () ?
      this->getStride () :
      static_cast<size_t> (0);
    const size_t lclNumRowsAfter = this->getLocalLength ();
    const size_t numColsAfter = this->getNumVectors ();
    const impl_scalar_type* hostPtrAfter =
      this->getDualView ().h_view.ptr_on_device ();

    const size_t strideRet = subViewMV->isConstantStride () ?
      subViewMV->getStride () :
      static_cast<size_t> (0);
    const size_t lclNumRowsRet = subViewMV->getLocalLength ();
    const size_t numColsRet = subViewMV->getNumVectors ();

    const char prefix[] = "Tpetra::MultiVector::offsetView: ";
    const char suffix[] = ".  This should never happen.  Please report this "
      "bug to the Tpetra developers.";

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! subMap.is_null () && lclNumRowsRet != subMap->getNodeNumElements (),
      std::logic_error, prefix << "Returned MultiVector has a number of rows "
      "different than the number of local indices in the input Map.  "
      "lclNumRowsRet: " << lclNumRowsRet << ", subMap->getNodeNumElements(): "
      << subMap->getNodeNumElements () << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION(
      strideBefore != strideAfter || lclNumRowsBefore != lclNumRowsAfter ||
      numColsBefore != numColsAfter || hostPtrBefore != hostPtrAfter,
      std::logic_error, prefix << "Original MultiVector changed dimensions, "
      "stride, or host pointer after taking offset view.  strideBefore: " <<
      strideBefore << ", strideAfter: " << strideAfter << ", lclNumRowsBefore: "
      << lclNumRowsBefore << ", lclNumRowsAfter: " << lclNumRowsAfter <<
      ", numColsBefore: " << numColsBefore << ", numColsAfter: " <<
      numColsAfter << ", hostPtrBefore: " << hostPtrBefore << ", hostPtrAfter: "
      << hostPtrAfter << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION(
      strideBefore != strideRet, std::logic_error, prefix << "Returned "
      "MultiVector has different stride than original MultiVector.  "
      "strideBefore: " << strideBefore << ", strideRet: " << strideRet <<
      ", numColsBefore: " << numColsBefore << ", numColsRet: " << numColsRet
      << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION(
      numColsBefore != numColsRet, std::logic_error,
      prefix << "Returned MultiVector has a different number of columns than "
      "original MultiVector.  numColsBefore: " << numColsBefore << ", "
      "numColsRet: " << numColsRet << suffix);
#endif // HAVE_TPETRA_DEBUG

    return subViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    return Teuchos::rcp_const_cast<MV> (this->offsetView (subMap, offset));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  subView (const Teuchos::ArrayView<const size_t>& cols) const
  {
    using Teuchos::Array;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;

    const size_t numViewCols = static_cast<size_t> (cols.size ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numViewCols < 1, std::runtime_error, "Tpetra::MultiVector::subView"
      "(const Teuchos::ArrayView<const size_t>&): The input array cols must "
      "contain at least one entry, but cols.size() = " << cols.size ()
      << " == 0.");

    // Check whether the index set in cols is contiguous.  If it is,
    // use the more efficient Range1D version of subView.
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
        return this->subView (Teuchos::Range1D (cols[0], cols[numViewCols-1]));
      }
    }

    if (isConstantStride ()) {
      return rcp (new MV (this->getMap (), view_, origView_, cols));
    }
    else {
      Array<size_t> newcols (cols.size ());
      for (size_t j = 0; j < numViewCols; ++j) {
        newcols[j] = whichVectors_[cols[j]];
      }
      return rcp (new MV (this->getMap (), view_, origView_, newcols ()));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  subView (const Teuchos::Range1D& colRng) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Array;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;
    const char tfecfFuncName[] = "subView(Range1D): ";

    const size_t numVecs = this->getNumVectors ();
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   colRng.size() == 0, std::runtime_error, prefix << "Range must include "
    //   "at least one vector.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (colRng.size ()) > numVecs, std::runtime_error,
      "colRng.size() = " << colRng.size () << " > this->getNumVectors() = "
      << numVecs << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs != 0 && colRng.size () != 0 &&
      (colRng.lbound () < static_cast<Teuchos::Ordinal> (0) ||
       static_cast<size_t> (colRng.ubound ()) >= numVecs),
      std::invalid_argument, "Nonempty input range [" << colRng.lbound () <<
      "," << colRng.ubound () << "] exceeds the valid range of column indices "
      "[0, " << numVecs << "].");

    // FIXME (mfh 07 Mar 2015) The commented-out subview() invocations
    // using ALL() don't work when some processes have zero rows, but
    // the enabled code does work.
    const std::pair<size_t, size_t> rows (0, this->getLocalLength ());
    if (colRng.size () == 0) {
      const std::pair<size_t, size_t> cols (0, 0); // empty range
      //dual_view_type X_sub = subview<dual_view_type> (view_, ALL (), cols);
      dual_view_type X_sub = subview (view_, rows, cols);
      return rcp (new MV (this->getMap (), X_sub, origView_));
    }
    else {
      // resulting MultiVector is constant stride only if *this is
      if (isConstantStride ()) {
        const std::pair<size_t, size_t> cols (colRng.lbound (), colRng.ubound () + 1);
        //dual_view_type X_sub = subview<dual_view_type> (view_, ALL (), cols);
        dual_view_type X_sub = subview (view_, rows, cols);
        return rcp (new MV (this->getMap (), X_sub, origView_));
      }
      else {
        if (colRng.size () == 1) {
          // We're only asking for one column, so the result does have
          // constant stride, even though this MultiVector does not.
          const std::pair<size_t, size_t> col (whichVectors_[0] + colRng.lbound (),
                                               whichVectors_[0] + colRng.ubound () + 1);
          //dual_view_type X_sub = subview<dual_view_type> (view_, ALL (), col);
          dual_view_type X_sub = subview (view_, rows, col);
          return rcp (new MV (this->getMap (), X_sub, origView_));
        }
        else {
          Array<size_t> which (whichVectors_.begin () + colRng.lbound (),
                               whichVectors_.begin () + colRng.ubound () + 1);
          return rcp (new MV (this->getMap (), view_, origView_, which));
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  subViewNonConst (const ArrayView<const size_t> &cols)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    return Teuchos::rcp_const_cast<MV> (this->subView (cols));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  subViewNonConst (const Teuchos::Range1D &colRng)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    return Teuchos::rcp_const_cast<MV> (this->subView (colRng));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getVector (const size_t j) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type, false> V;

#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getVector(NonConst): ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      this->vectorIndexOutOfRange (j), std::runtime_error, "Input index j (== "
      << j << ") exceeds valid range [0, " << this->getNumVectors ()
      << " - 1].");
#endif // HAVE_TPETRA_DEBUG
    const size_t jj = this->isConstantStride () ?
      static_cast<size_t> (j) :
      static_cast<size_t> (this->whichVectors_[j]);
    const std::pair<size_t, size_t> rng (jj, jj+1);
    if (view_.dimension_0 () > 0) {
      return rcp (new V (this->getMap (),
                         subview (view_, ALL (), rng),
                         origView_));
    } else {
      // FIXME (mfh 04 Mar 2015) Doesn't this need to know about origView_?
      return rcp (new V (this->getMap ()));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getVectorNonConst (const size_t j)
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;
    return Teuchos::rcp_const_cast<V> (this->getVector (j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get1dCopy (const Teuchos::ArrayView<Scalar>& A, const size_t LDA) const
  {
    using Kokkos::subview;
    typedef impl_scalar_type IST;
    typedef MakeUnmanagedView<IST, execution_space> view_getter_type;
    typedef typename view_getter_type::view_type input_col_type;
    // Types of views of this MultiVector's data.
    typedef typename dual_view_type::t_host host_view_type;
    typedef typename dual_view_type::t_dev dev_view_type;
    typedef Kokkos::View<IST*,
      typename host_view_type::array_layout,
      typename host_view_type::memory_space> host_col_type;
    typedef Kokkos::View<IST*,
      typename dev_view_type::array_layout,
      typename dev_view_type::memory_space> dev_col_type;
    const char tfecfFuncName[] = "get1dCopy: ";

    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const std::pair<size_t, size_t> rowRange (0, numRows);
    const std::pair<size_t, size_t> colRange (0, numCols);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < numRows, std::runtime_error,
      "LDA = " << LDA << " < numRows = " << numRows << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numRows > static_cast<size_t> (0) &&
      numCols > static_cast<size_t> (0) &&
      static_cast<size_t> (A.size ()) < LDA * (numCols - 1) + numRows,
      std::runtime_error,
      "A.size() = " << A.size () << ", but its size must be at least "
      << (LDA * (numCols - 1) + numRows) << " to hold all the entries.");

    // FIXME (mfh 22 Jul 2014, 10 Dec 2014) Currently, it doesn't work
    // to do a 2-D copy, even if this MultiVector has constant stride.
    // This is because Kokkos can't currently tell the difference
    // between padding (which permits a single deep_copy for the whole
    // 2-D View) and stride > numRows (which does NOT permit a single
    // deep_copy for the whole 2-D View).  Carter is working on this,
    // but for now, the temporary fix is to copy one column at a time.

    for (size_t j = 0; j < numCols; ++j) {
      const size_t srcCol =
        this->isConstantStride () ? j : this->whichVectors_[j];
      const size_t dstCol = j;
      IST* const dstColRaw =
        reinterpret_cast<IST*> (A.getRawPtr () + LDA * dstCol);
      input_col_type dstColView (dstColRaw, numRows);
      // Use the most recently updated version of this MultiVector's
      // data.  This avoids sync'ing, which could violate users'
      // expectations.
      if (view_.modified_host >= view_.modified_device) {
        host_col_type srcColView =
          subview (view_.h_view, rowRange, srcCol);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          dstColView.dimension_0 () != srcColView.dimension_0 (),
          std::logic_error, ": srcColView and dstColView have different "
          "dimensions.  Please report this bug to the Tpetra developers.");
        Kokkos::deep_copy (dstColView, srcColView);
      }
      else {
        dev_col_type srcColView =
          subview (view_.d_view, rowRange, srcCol);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          dstColView.dimension_0 () != srcColView.dimension_0 (),
          std::logic_error, ": srcColView and dstColView have different "
          "dimensions.  Please report this bug to the Tpetra developers.");
        Kokkos::deep_copy (dstColView, srcColView);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get2dCopy (const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> >& ArrayOfPtrs) const
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> V;
    const char tfecfFuncName[] = "get2dCopy: ";
    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (ArrayOfPtrs.size ()) != numCols,
      std::runtime_error, "Input array of pointers must contain as many "
      "entries (arrays) as the MultiVector has columns.  ArrayOfPtrs.size() = "
      << ArrayOfPtrs.size () << " != getNumVectors() = " << numCols << ".");

    if (numRows != 0 && numCols != 0) {
      // No side effects until we've validated the input.
      for (size_t j = 0; j < numCols; ++j) {
        const size_t dstLen = static_cast<size_t> (ArrayOfPtrs[j].size ());
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          dstLen < numRows, std::invalid_argument, "Array j = " << j << " of "
          "the input array of arrays is not long enough to fit all entries in "
          "that column of the MultiVector.  ArrayOfPtrs[j].size() = " << dstLen
          << " < getLocalLength() = " << numRows << ".");
      }

      // We've validated the input, so it's safe to start copying.
      for (size_t j = 0; j < numCols; ++j) {
        RCP<const V> X_j = this->getVector (j);
        const size_t LDA = static_cast<size_t> (ArrayOfPtrs[j].size ());
        X_j->get1dCopy (ArrayOfPtrs[j], LDA);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get1dView () const
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
      // NOTE (mfh 09 2014) get1dView() and get1dViewNonConst() have
      // always been device->host synchronization points.  We might
      // want to change this in the future.
      typedef typename dual_view_type::host_mirror_space host_type;
      view_.template sync<host_type> ();
      // Both get1dView() and get1dViewNonConst() return a host view
      // of the data.
      Teuchos::ArrayRCP<const impl_scalar_type> dataAsArcp =
        Kokkos::Compat::persistingView (view_.template view<host_type> ());
      return Teuchos::arcp_reinterpret_cast<const Scalar> (dataAsArcp);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get1dViewNonConst ()
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
      // NOTE (mfh 09 May 2014) get1dView() and get1dViewNonConst()
      // have always been device->host synchronization points.  We
      // might want to change this in the future.
      typedef typename dual_view_type::host_mirror_space host_type;
      view_.template sync<host_type> ();
      // Both get1dView() and get1dViewNonConst() return a host view
      // of the data.
      Teuchos::ArrayRCP<impl_scalar_type> dataAsArcp =
        Kokkos::Compat::persistingView (view_.template view<host_type> ());
      return Teuchos::arcp_reinterpret_cast<Scalar> (dataAsArcp);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get2dViewNonConst ()
  {
    using Teuchos::ArrayRCP;
    typedef Kokkos::DualView<impl_scalar_type*,
      typename dual_view_type::array_layout, execution_space> col_dual_view_type;

    const size_t numCols = getNumVectors ();
    ArrayRCP<ArrayRCP<Scalar> > views (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      const size_t col = isConstantStride () ? j : whichVectors_[j];
      col_dual_view_type X_col =
        Kokkos::subview (view_, Kokkos::ALL (), col);
      ArrayRCP<impl_scalar_type> X_col_arcp =
        Kokkos::Compat::persistingView (X_col.d_view);
      views[j] = Teuchos::arcp_reinterpret_cast<Scalar> (X_col_arcp);
    }
    return views;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  get2dView () const
  {
    using Teuchos::ArrayRCP;
    typedef Kokkos::DualView<const impl_scalar_type*,
      typename dual_view_type::array_layout, execution_space> col_dual_view_type;

    const size_t numCols = getNumVectors ();
    ArrayRCP<ArrayRCP<const Scalar> > views (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      const size_t col = isConstantStride () ? j : whichVectors_[j];
      col_dual_view_type X_col =
        Kokkos::subview (view_, Kokkos::ALL (), col);
      ArrayRCP<const impl_scalar_type> X_col_arcp =
        Kokkos::Compat::persistingView (X_col.d_view);
      views[j] = Teuchos::arcp_reinterpret_cast<const Scalar> (X_col_arcp);
    }
    return views;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  multiply (Teuchos::ETransp transA,
            Teuchos::ETransp transB,
            const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& A,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& B,
            const Scalar& beta)
  {
    using Teuchos::CONJ_TRANS;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;
    const char errPrefix[] = "Tpetra::MultiVector::multiply: ";

    // This routine performs a variety of matrix-matrix multiply
    // operations, interpreting the MultiVector (this-aka C , A and B)
    // as 2D matrices.  Variations are due to the fact that A, B and C
    // can be local replicated or global distributed MultiVectors and
    // that we may or may not operate with the transpose of A and B.
    // Possible cases are:
    //
    //     Operations                          # Cases  Notes
    //  1) C(local) = A^X(local) * B^X(local)  4        X=Trans or Not, no comm needed
    //  2) C(local) = A^T(distr) * B  (distr)  1        2-D dot product, replicate C
    //  3) C(distr) = A  (distr) * B^X(local)  2        2-D vector update, no comm needed
    //
    // The following operations are not meaningful for 1-D
    // distributions:
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
    // Total number of cases: 32 (= 2^5).

    TEUCHOS_TEST_FOR_EXCEPTION(
      ATS::is_complex && (transA == TRANS || transB == TRANS),
      std::invalid_argument, errPrefix << "Transpose without conjugation "
      "(transA == TRANS || transB == TRANS) is not supported for complex Scalar "
      "types.");

    transA = (transA == NO_TRANS ? NO_TRANS : CONJ_TRANS);
    transB = (transB == NO_TRANS ? NO_TRANS : CONJ_TRANS);

    // Compute effective dimensions, w.r.t. transpose operations on
    size_t A_nrows = (transA==CONJ_TRANS) ? A.getNumVectors() : A.getLocalLength();
    size_t A_ncols = (transA==CONJ_TRANS) ? A.getLocalLength() : A.getNumVectors();
    size_t B_nrows = (transB==CONJ_TRANS) ? B.getNumVectors() : B.getLocalLength();
    size_t B_ncols = (transB==CONJ_TRANS) ? B.getLocalLength() : B.getNumVectors();

    impl_scalar_type beta_local = beta; // local copy of beta; might be reassigned below

    TEUCHOS_TEST_FOR_EXCEPTION(
      getLocalLength () != A_nrows || getNumVectors () != B_ncols ||
      A_ncols != B_nrows, std::runtime_error, errPrefix << "Dimensions of "
      "*this, op(A), and op(B) must be consistent.  Local part of *this is "
      << getLocalLength() << " x " << getNumVectors()
      << ", A is " << A_nrows << " x " << A_ncols
      << ", and B is " << B_nrows << " x " << B_ncols << ".");

    const bool A_is_local = ! A.isDistributed ();
    const bool B_is_local = ! B.isDistributed ();
    const bool C_is_local = ! this->isDistributed ();
    // Case 1: C(local) = A^X(local) * B^X(local)
    const bool Case1 = C_is_local && A_is_local && B_is_local;
    // Case 2: C(local) = A^T(distr) * B  (distr)
    const bool Case2 = C_is_local && ! A_is_local && ! B_is_local &&
      transA == CONJ_TRANS && transB == NO_TRANS;
    // Case 3: C(distr) = A  (distr) * B^X(local)
    const bool Case3 = ! C_is_local && ! A_is_local && B_is_local &&
      transA == NO_TRANS;

    // Test that we are considering a meaningful case
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! Case1 && ! Case2 && ! Case3, std::runtime_error, errPrefix
      << "Multiplication of op(A) and op(B) into *this is not a "
      "supported use case.");

    if (beta != STS::zero () && Case2) {
      // If Case2, then C is local and contributions must be summed
      // across all processes.  However, if beta != 0, then accumulate
      // beta*C into the sum.  When summing across all processes, we
      // only want to accumulate this once, so set beta == 0 on all
      // processes except Process 0.
      const int myRank = this->getMap ()->getComm ()->getRank ();
      if (myRank != 0) {
        beta_local = STS::zero ();
      }
    }

    // We only know how to do matrix-matrix multiplies if all the
    // MultiVectors have constant stride.  If not, we have to make
    // temporary copies of those MultiVectors (including possibly
    // *this) that don't have constant stride.
    RCP<MV> C_tmp;
    if (! isConstantStride ()) {
      C_tmp = rcp (new MV (*this, Teuchos::Copy)); // deep copy
    } else {
      C_tmp = rcp (this, false);
    }

    RCP<const MV> A_tmp;
    if (! A.isConstantStride ()) {
      A_tmp = rcp (new MV (A, Teuchos::Copy)); // deep copy
    } else {
      A_tmp = rcpFromRef (A);
    }

    RCP<const MV> B_tmp;
    if (! B.isConstantStride ()) {
      B_tmp = rcp (new MV (B, Teuchos::Copy)); // deep copy
    } else {
      B_tmp = rcpFromRef (B);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! C_tmp->isConstantStride () || ! B_tmp->isConstantStride () ||
      ! A_tmp->isConstantStride (), std::logic_error, errPrefix
      << "Failed to make temporary constant-stride copies of MultiVectors.");

    typedef Kokkos::DeviceGEMM<impl_scalar_type, execution_space> gemm_type;

    gemm_type::GEMM (transA, transB, alpha,
                     A_tmp->getDualView ().d_view, B_tmp->getDualView ().d_view,
                     beta_local, C_tmp->getDualView ().d_view);
    if (! isConstantStride ()) {
      deep_copy (*this, *C_tmp); // Copy the result back into *this.
    }

    // Dispose of (possibly) extra copies of A and B.
    A_tmp = Teuchos::null;
    B_tmp = Teuchos::null;

    // If Case 2 then sum up *this and distribute it to all processes.
    if (Case2) {
      this->reduce ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  elementWiseMultiply (Scalar scalarAB,
                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal, node_type, false>& A,
                       const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, node_type>& B,
                       Scalar scalarThis)
  {
    using Kokkos::ALL;
    using Kokkos::subview;

    typedef typename dual_view_type::t_dev view_2d_type;
    typedef Kokkos::View<impl_scalar_type*,
      typename view_2d_type::array_layout,
      typename view_2d_type::device_type,
      typename view_2d_type::memory_traits> view_1d_type;

    const char tfecfFuncName[] = "elementWiseMultiply: ";
    const size_t numVecs = this->getNumVectors ();

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength() != A.getLocalLength() ||
      getLocalLength() != B.getLocalLength(), std::runtime_error,
      "MultiVectors do not have the same local length.");
#endif // HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs != B.getNumVectors (), std::runtime_error, "this->getNumVectors"
      "() = " << numVecs << " != B.getNumVectors() = " << B.getNumVectors ()
      << ".");

    if (isConstantStride () && A.isConstantStride ()) {
      // FIXME (mfh 02 Oct 2014) Shouldn't it be asking if B has
      // constant stride?  A is just a Vector; it only has one column,
      // so it always has constant stride.
      //
      // If both *this and A have constant stride, we can do an
      // element-wise multiply on all columns at once.
      view_.template sync<DeviceType> ();
      view_.template modify<DeviceType> ();
      A.view_.template sync<DeviceType> ();
      A.view_.template modify<DeviceType> ();
      B.view_.template sync<DeviceType> ();
      B.view_.template modify<DeviceType> ();
      view_1d_type vector_A = subview (A.view_.d_view, ALL (), 0);
      Kokkos::MV_ElementWiseMultiply (scalarThis, view_.d_view,
                                      scalarAB, vector_A, B.view_.d_view);
    }
    else {
      view_.template sync<DeviceType> ();
      view_.template modify<DeviceType> ();
      A.view_.template sync<DeviceType> ();
      A.view_.template modify<DeviceType> ();
      B.view_.template sync<DeviceType> ();
      B.view_.template modify<DeviceType> ();
      view_1d_type vector_A = subview (A.view_.d_view, ALL (), 0);
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_1d_type vector_k =
          subview (view_.d_view, ALL (), this_col);
        const size_t B_col = isConstantStride () ? k : B.whichVectors_[k];
        view_1d_type vector_Bk =
          subview (B.view_.d_view, ALL (), B_col);
        Kokkos::V_ElementWiseMultiply (scalarThis, vector_k, scalarAB,
                                       vector_A, vector_Bk);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  reduce ()
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef typename dual_view_type::t_dev device_view_type;
    typedef typename dual_view_type::host_mirror_space host_mirror_space;

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->isDistributed (), std::runtime_error,
      "Tpetra::MultiVector::reduce() should only be called with locally "
      "replicated or otherwise not distributed MultiVector objects.");
    const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
    if (comm.getSize () == 1) {
      return;
    }

    const size_t numLclRows = getLocalLength ();
    const size_t numCols = getNumVectors ();

    // FIXME (mfh 16 June 2014) This exception will cause deadlock if
    // it triggers on only some processes.  We don't have a good way
    // to pack this result into the all-reduce below, but this would
    // be a good reason to set a "local error flag" and find other
    // opportunities to let it propagate.
    TEUCHOS_TEST_FOR_EXCEPTION(
      numLclRows > static_cast<size_t> (std::numeric_limits<int>::max ()),
      std::runtime_error, "Tpetra::MultiVector::reduce: On Process " <<
      comm.getRank () << ", the number of local rows " << numLclRows <<
      " does not fit in int.");

    //
    // Use MPI to sum the entries across all local blocks.
    //
    // If this MultiVector's local data are stored contiguously, we
    // can use the local View as the source buffer in the
    // MPI_Allreduce.  Otherwise, we have to allocate a temporary
    // source buffer and pack.
    const bool contig = isConstantStride () && getStride () == numLclRows;
    device_view_type srcBuf;
    if (contig) {
      srcBuf = view_.d_view;
    }
    else {
      srcBuf = device_view_type ("srcBuf", numLclRows, numCols);
      Kokkos::deep_copy (srcBuf, view_.d_view);
    }

    // MPI requires that the send and receive buffers don't alias one
    // another, so we have to copy temporary storage for the result.
    //
    // We expect that MPI implementations will know how to read device
    // pointers.
    device_view_type tgtBuf ("tgtBuf", numLclRows, numCols);

    const int reduceCount = static_cast<int> (numLclRows * numCols);
    reduceAll<int, impl_scalar_type> (comm, REDUCE_SUM, reduceCount,
                                 srcBuf.ptr_on_device (),
                                 tgtBuf.ptr_on_device ());

    // Tell the DualView that we plan to modify the device data.
    view_.template modify<execution_space> ();

    const std::pair<size_t, size_t> lclRowRange (0, numLclRows);
    device_view_type d_view =
      subview (view_.d_view, lclRowRange, ALL ());

    if (contig || isConstantStride ()) {
      Kokkos::deep_copy (d_view, tgtBuf);
    }
    else {
      for (size_t j = 0; j < numCols; ++j) {
        device_view_type d_view_j =
          subview (d_view, ALL (), std::pair<int,int>(j,j+1));
        device_view_type tgtBuf_j =
          subview (tgtBuf, ALL (), std::pair<int,int>(j,j+1));
        Kokkos::deep_copy (d_view_j, tgtBuf_j);
      }
    }

    // Synchronize the host with changes on the device.
    //
    // FIXME (mfh 16 June 2014) This raises the question of whether we
    // want to synchronize always.  Users will find it reassuring if
    // MultiVector methods always leave the MultiVector in a
    // synchronized state, but it seems silly to synchronize to host
    // if they hardly ever need host data.
    view_.template sync<host_mirror_space> ();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const impl_scalar_type& ScalarValue)
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
    view_.h_view (MyRow, colInd) = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sumIntoLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const impl_scalar_type& ScalarValue)
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
    view_.h_view (MyRow, colInd) += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  replaceGlobalValue (GlobalOrdinal GlobalRow,
                      size_t VectorIndex,
                      const impl_scalar_type& ScalarValue)
  {
    const LocalOrdinal MyRow = this->getMap ()->getLocalElement (GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid (),
      std::runtime_error,
      "Tpetra::MultiVector::replaceGlobalValue: Global row index " << GlobalRow
      << "is not present on this process "
      << this->getMap ()->getComm ()->getRank () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange (VectorIndex), std::runtime_error,
      "Tpetra::MultiVector::replaceGlobalValue: Vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif // HAVE_TPETRA_DEBUG
    replaceLocalValue (MyRow, VectorIndex, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  sumIntoGlobalValue (GlobalOrdinal GlobalRow,
                      size_t VectorIndex,
                      const impl_scalar_type& ScalarValue)
  {
    const LocalOrdinal MyRow = this->getMap ()->getLocalElement (GlobalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid (),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: Global row index " << GlobalRow
      << "is not present on this process "
      << this->getMap ()->getComm ()->getRank () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: Vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    sumIntoLocalValue (MyRow, VectorIndex, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class T>
  Teuchos::ArrayRCP<T>
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getSubArrayRCP (Teuchos::ArrayRCP<T> arr,
                  size_t j) const
  {
    typedef Kokkos::DualView<impl_scalar_type*,
      typename dual_view_type::array_layout,
      execution_space> col_dual_view_type;
    const size_t col = isConstantStride () ? j : whichVectors_[j];
    col_dual_view_type X_col =
      Kokkos::subview (view_, Kokkos::ALL (), col);
    return Kokkos::Compat::persistingView (X_col.d_view);
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  KokkosClassic::MultiVector<
    Scalar, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getLocalMV () const
  {
    using Teuchos::ArrayRCP;
    typedef KokkosClassic::MultiVector<Scalar, node_type> KMV;

    // This method creates the KokkosClassic object on the fly.  Thus,
    // it's OK to leave this in place for backwards compatibility.
    ArrayRCP<Scalar> data;
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      data = Teuchos::null;
    }
    else {
      ArrayRCP<impl_scalar_type> dataTmp = (getLocalLength() > 0) ?
        Kokkos::Compat::persistingView (view_.d_view) :
        Teuchos::null;
      data = Teuchos::arcp_reinterpret_cast<Scalar> (dataTmp);
    }
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA =
      origView_.dimension_1 () > 1 ? stride[1] : origView_.dimension_0 ();

    KMV kmv (this->getMap ()->getNode ());
    kmv.initializeValues (getLocalLength (), getNumVectors (),
                          data, LDA, getOrigNumLocalRows (),
                          getOrigNumLocalCols ());
    return kmv;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class DeviceType>
  typename MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::dual_view_type
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  getDualView () const {
    return view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  description () const
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
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
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
              << "Node: " << Teuchos::TypeNameTraits<node_type>::name () << endl;
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
            out << "localNumRows: " << getLocalLength() << endl
                << "isConstantStride: " << isConstantStride () << endl;
            if (vl != VERB_MEDIUM) {
              // >= VERB_HIGH: print isConstantStride() and getStride()
              if (isConstantStride()) {
                out << "columnStride: " << getStride() << endl;
              }
              if (vl == VERB_EXTREME) {
                // VERB_EXTREME: print all the values in the multivector.
                out << "values: " << endl;
                typename dual_view_type::t_host X = this->getDualView ().h_view;
                out << "[";
                for (size_t i = 0; i < getLocalLength (); ++i) {
                  for (size_t j = 0; j < getNumVectors (); ++j) {
                    const size_t col = isConstantStride () ? j : whichVectors_[j];
                    out << X(i,col);
                    if (j + 1 < getNumVectors ()) {
                      out << ", ";
                    }
                  } // for each column
                  if (i + 1 < getLocalLength ()) {
                    out << "; ";
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
        comm->barrier ();
      } // for each process in the communicator
    } // if vl != VERB_NONE
  }

#if TPETRA_USE_KOKKOS_DISTOBJECT
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  createViews () const
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  createViewsNonConst (KokkosClassic::ReadWriteOption rwo)
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  releaseViews () const
  {
    // Do nothing in Kokkos::View implementation
  }

#else // NOT TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  createViews () const
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  createViewsNonConst (KokkosClassic::ReadWriteOption /* rwo */ )
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  releaseViews () const
  {}

#endif // TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, node_type> >& newMap)
  {
    replaceMap (newMap);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>, false>::
  assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& src)
  {
    using Kokkos::parallel_for;
    typedef LocalOrdinal LO;
    typedef typename DeviceType::device_type DT;
    typedef typename dual_view_type::host_mirror_space::device_type HMDT;
    const bool debug = false;

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getGlobalLength () != src.getGlobalLength () ||
      this->getNumVectors () != src.getNumVectors (), std::invalid_argument,
      "Tpetra::deep_copy: Global dimensions of the two Tpetra::MultiVector "
      "objects do not match.  src has dimensions [" << src.getGlobalLength ()
      << "," << src.getNumVectors () << "], and *this has dimensions ["
      << this->getGlobalLength () << "," << this->getNumVectors () << "].");
    // FIXME (mfh 28 Jul 2014) Don't throw; just set a local error flag.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getLocalLength () != src.getLocalLength (), std::invalid_argument,
      "Tpetra::deep_copy: The local row counts of the two Tpetra::MultiVector "
      "objects do not match.  src has " << src.getLocalLength () << " row(s) "
      << " and *this has " << this->getLocalLength () << " row(s).");

    if (debug && this->getMap ()->getComm ()->getRank () == 0) {
      std::cout << "*** MultiVector::assign: ";
    }

    if (src.isConstantStride () && this->isConstantStride ()) {
      if (debug && this->getMap ()->getComm ()->getRank () == 0) {
        std::cout << "Both *this and src have constant stride" << std::endl;
      }

      if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
        // Device memory has the most recent version of src.
        this->template modify<DT> (); // We are about to modify dst on device.
        // Copy from src to dst on device.
        Details::localDeepCopyConstStride (this->getDualView ().template view<DT> (),
                                  src.getDualView ().template view<DT> ());
        this->template sync<HMDT> (); // Sync dst from device to host.
      }
      else { // Host memory has the most recent version of src.
        this->template modify<HMDT> (); // We are about to modify dst on host.
        // Copy from src to dst on host.
        Details::localDeepCopyConstStride (this->getDualView ().template view<HMDT> (),
                                  src.getDualView ().template view<HMDT> ());
        this->template sync<DT> (); // Sync dst from host to device.
      }
    }
    else {
      if (this->isConstantStride ()) {
        if (debug && this->getMap ()->getComm ()->getRank () == 0) {
          std::cout << "Only *this has constant stride";
        }

        const LO numWhichVecs = static_cast<LO> (src.whichVectors_.size ());
        const std::string whichVecsLabel ("MV::deep_copy::whichVecs");

        // We can't sync src, since it is only an input argument.
        // Thus, we have to use the most recently modified version of
        // src, device or host.
        if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
          if (debug && this->getMap ()->getComm ()->getRank () == 0) {
            std::cout << "; Copy from device version of src" << std::endl;
          }
          // Copy from the device version of src.
          //
          // whichVecs tells the kernel which vectors (columns) of src
          // to copy.  Fill whichVecs on the host, and sync to device.
          typedef Kokkos::DualView<LO*, DT> whichvecs_type;
          whichvecs_type srcWhichVecs (whichVecsLabel, numWhichVecs);
          srcWhichVecs.template modify<HMDT> ();
          for (LO i = 0; i < numWhichVecs; ++i) {
            srcWhichVecs.h_view(i) = static_cast<LO> (src.whichVectors_[i]);
          }
          // Sync the host version of srcWhichVecs to the device.
          srcWhichVecs.template sync<DT> ();

          // Mark the device version of dst's DualView as modified.
          this->template modify<DT> ();

          // Copy from the selected vectors of src to dst, on the
          // device.  The function ignores its dstWhichVecs argument
          // in this case.
          Details::localDeepCopy (this->getDualView ().template view<DT> (),
                         src.getDualView ().template view<DT> (),
                         true, false, srcWhichVecs.d_view, srcWhichVecs.d_view);
          // Sync *this' DualView to the host.  This is cheaper than
          // repeating the above copy from src to *this on the host.
          this->template sync<HMDT> ();
        }
        else { // host version of src was the most recently modified
          if (debug && this->getMap ()->getComm ()->getRank () == 0) {
            std::cout << "; Copy from host version of src" << std::endl;
          }
          // Copy from the host version of src.
          //
          // whichVecs tells the kernel which vectors (columns) of src
          // to copy.  Fill whichVecs on the host, and use it there.
          typedef Kokkos::View<LO*, HMDT> whichvecs_type;
          whichvecs_type srcWhichVecs (whichVecsLabel, numWhichVecs);
          for (LO i = 0; i < numWhichVecs; ++i) {
            srcWhichVecs(i) = static_cast<LO> (src.whichVectors_[i]);
          }
          // Copy from the selected vectors of src to dst, on the
          // host.  The function ignores its dstWhichVecs argument in
          // this case.
          Details::localDeepCopy (this->getDualView ().template view<HMDT> (),
                         src.getDualView ().template view<HMDT> (),
                         true, false, srcWhichVecs, srcWhichVecs);
          // Sync dst back to the device, since we only copied on the host.
          this->template sync<DT> ();
        }
      }
      else { // dst is NOT constant stride
        if (src.isConstantStride ()) {
          if (debug && this->getMap ()->getComm ()->getRank () == 0) {
            std::cout << "Only src has constant stride" << std::endl;
          }

          if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
            // Copy from the device version of src.
            //
            // whichVecs tells the kernel which vectors (columns) of dst
            // to copy.  Fill whichVecs on the host, and sync to device.
            typedef Kokkos::DualView<LO*, DT> whichvecs_type;
            const std::string whichVecsLabel ("MV::deep_copy::whichVecs");
            const LO numWhichVecs = static_cast<LO> (this->whichVectors_.size ());
            whichvecs_type whichVecs (whichVecsLabel, numWhichVecs);
            whichVecs.template modify<HMDT> ();
            for (LO i = 0; i < numWhichVecs; ++i) {
              whichVecs.h_view(i) = this->whichVectors_[i];
            }
            // Sync the host version of whichVecs to the device.
            whichVecs.template sync<DT> ();

            // Copy src to the selected vectors of dst, on the device.
            Details::localDeepCopy (this->getDualView ().template view<DT> (),
                           src.getDualView ().template view<DT> (),
                           this->isConstantStride (), src.isConstantStride (),
                           whichVecs.d_view, whichVecs.d_view);
            // We can't sync src and repeat the above copy on the
            // host, so sync dst back to the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            this->template sync<HMDT> ();
          }
          else { // host version of src was the most recently modified
            // Copy from the host version of src.
            //
            // whichVecs tells the kernel which vectors (columns) of src
            // to copy.  Fill whichVecs on the host, and use it there.
            typedef Kokkos::View<LO*, HMDT> whichvecs_type;
            const LO numWhichVecs = static_cast<LO> (this->whichVectors_.size ());
            whichvecs_type whichVecs ("MV::deep_copy::whichVecs", numWhichVecs);
            for (LO i = 0; i < numWhichVecs; ++i) {
              whichVecs(i) = static_cast<LO> (this->whichVectors_[i]);
            }
            // Copy from src to the selected vectors of dst, on the
            // host.  The functor ignores its 4th arg in this case.
            Details::localDeepCopy (this->getDualView ().template view<HMDT> (),
                           src.getDualView ().template view<HMDT> (),
                           this->isConstantStride (), src.isConstantStride (),
                           whichVecs, whichVecs);
            // Sync dst back to the device, since we only copied on the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            this->template sync<DT> ();
          }
        }
        else { // neither src nor dst have constant stride
          if (debug && this->getMap ()->getComm ()->getRank () == 0) {
            std::cout << "Neither *this nor src has constant stride" << std::endl;
          }

          if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
            // Copy from the device version of src.
            //
            // whichVectorsDst tells the kernel which vectors
            // (columns) of dst to copy.  Fill it on the host, and
            // sync to device.
            const LO dstNumWhichVecs = static_cast<LO> (this->whichVectors_.size ());
            Kokkos::DualView<LO*, DT> whichVecsDst ("MV::deep_copy::whichVecsDst",
                                                    dstNumWhichVecs);
            whichVecsDst.template modify<HMDT> ();
            for (LO i = 0; i < dstNumWhichVecs; ++i) {
              whichVecsDst.h_view(i) = static_cast<LO> (this->whichVectors_[i]);
            }
            // Sync the host version of whichVecsDst to the device.
            whichVecsDst.template sync<DT> ();

            // whichVectorsSrc tells the kernel which vectors
            // (columns) of src to copy.  Fill it on the host, and
            // sync to device.  Use the destination MultiVector's
            // LocalOrdinal type here.
            const LO srcNumWhichVecs = static_cast<LO> (src.whichVectors_.size ());
            Kokkos::DualView<LO*, DT> whichVecsSrc ("MV::deep_copy::whichVecsSrc",
                                                    srcNumWhichVecs);
            whichVecsSrc.template modify<HMDT> ();
            for (LO i = 0; i < srcNumWhichVecs; ++i) {
              whichVecsSrc.h_view(i) = static_cast<LO> (src.whichVectors_[i]);
            }
            // Sync the host version of whichVecsSrc to the device.
            whichVecsSrc.template sync<DT> ();

            // Copy from the selected vectors of src to the selected
            // vectors of dst, on the device.
            Details::localDeepCopy (this->getDualView ().template view<DT> (),
                           src.getDualView ().template view<DT> (),
                           this->isConstantStride (), src.isConstantStride (),
                           whichVecsDst.d_view, whichVecsSrc.d_view);
          }
          else {
            const LO dstNumWhichVecs = static_cast<LO> (this->whichVectors_.size ());
            Kokkos::View<LO*, HMDT> whichVectorsDst ("dstWhichVecs", dstNumWhichVecs);
            for (LO i = 0; i < dstNumWhichVecs; ++i) {
              whichVectorsDst(i) = this->whichVectors_[i];
            }

            // Use the destination MultiVector's LocalOrdinal type here.
            const LO srcNumWhichVecs = static_cast<LO> (src.whichVectors_.size ());
            Kokkos::View<LO*, HMDT> whichVectorsSrc ("srcWhichVecs", srcNumWhichVecs);
            for (LO i = 0; i < srcNumWhichVecs; ++i) {
              whichVectorsSrc(i) = src.whichVectors_[i];
            }

            // Copy from the selected vectors of src to the selected
            // vectors of dst, on the host.
            Details::localDeepCopy (this->getDualView ().template view<HMDT> (),
                           src.getDualView ().template view<HMDT> (),
                           this->isConstantStride (), src.isConstantStride (),
                           whichVectorsDst, whichVectorsSrc);

            // We can't sync src and repeat the above copy on the
            // host, so sync dst back to the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            this->template sync<HMDT> ();
          }
        }
      }
    }
  }


  template <class Scalar, class LO, class GO, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  createMultiVector (const Teuchos::RCP<const Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >& map,
                     size_t numVectors)
  {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;
    typedef MultiVector<Scalar, LO, GO, node_type> MV;
    return Teuchos::rcp (new MV (map, numVectors));
  }

  template <class ST, class LO, class GO, class DeviceType>
  MultiVector<ST, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  createCopy (const MultiVector<ST, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& src)
  {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;
    typedef MultiVector<ST, LO, GO, node_type> MV;

    MV cpy (src.getMap (), src.getNumVectors (), false);
    cpy.assign (src);
    return cpy;
  }

} // namespace Tpetra


#endif // TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP
