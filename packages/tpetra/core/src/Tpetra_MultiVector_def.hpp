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

#ifndef TPETRA_MULTIVECTOR_DEF_HPP
#define TPETRA_MULTIVECTOR_DEF_HPP

/// \file Tpetra_MultiVector_def.hpp
/// \brief Definition of the Tpetra::MultiVector class
///
/// If you want to use Tpetra::MultiVector, include
/// "Tpetra_MultiVector.hpp" (a file which CMake generates and
/// installs for you).  If you only want the declaration of
/// Tpetra::MultiVector, include "Tpetra_MultiVector_decl.hpp".

#include "Tpetra_Util.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_allReduceView.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_checkView.hpp"
#include "Tpetra_Details_fill.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_isInterComm.hpp"
#include "Tpetra_Details_lclDot.hpp"
#include "Tpetra_Details_normImpl.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"
#ifdef HAVE_TPETRACORE_TEUCHOSNUMERICS
#  include "Teuchos_SerialDenseMatrix.hpp"
#endif // HAVE_TPETRACORE_TEUCHOSNUMERICS
#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp"
#include "KokkosCompat_View.hpp"
#include "KokkosBlas.hpp"
#include "KokkosKernels_Utils.hpp"
#include "Kokkos_Random.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <memory>
#include <sstream>

#ifdef HAVE_TPETRA_INST_FLOAT128
namespace Kokkos {
  // FIXME (mfh 04 Sep 2015) Just a stub for now!
  template<class Generator>
  struct rand<Generator, __float128> {
    static KOKKOS_INLINE_FUNCTION __float128 max ()
    {
      return static_cast<__float128> (1.0);
    }
    static KOKKOS_INLINE_FUNCTION __float128
    draw (Generator& gen)
    {
      // Half the smallest normalized double, is the scaling factor of
      // the lower-order term in the double-double representation.
      const __float128 scalingFactor =
        static_cast<__float128> (std::numeric_limits<double>::min ()) /
        static_cast<__float128> (2.0);
      const __float128 higherOrderTerm = static_cast<__float128> (gen.drand ());
      const __float128 lowerOrderTerm =
        static_cast<__float128> (gen.drand ()) * scalingFactor;
      return higherOrderTerm + lowerOrderTerm;
    }
    static KOKKOS_INLINE_FUNCTION __float128
    draw (Generator& gen, const __float128& range)
    {
      // FIXME (mfh 05 Sep 2015) Not sure if this is right.
      const __float128 scalingFactor =
        static_cast<__float128> (std::numeric_limits<double>::min ()) /
        static_cast<__float128> (2.0);
      const __float128 higherOrderTerm =
        static_cast<__float128> (gen.drand (range));
      const __float128 lowerOrderTerm =
        static_cast<__float128> (gen.drand (range)) * scalingFactor;
      return higherOrderTerm + lowerOrderTerm;
    }
    static KOKKOS_INLINE_FUNCTION __float128
    draw (Generator& gen, const __float128& start, const __float128& end)
    {
      // FIXME (mfh 05 Sep 2015) Not sure if this is right.
      const __float128 scalingFactor =
        static_cast<__float128> (std::numeric_limits<double>::min ()) /
        static_cast<__float128> (2.0);
      const __float128 higherOrderTerm =
        static_cast<__float128> (gen.drand (start, end));
      const __float128 lowerOrderTerm =
        static_cast<__float128> (gen.drand (start, end)) * scalingFactor;
      return higherOrderTerm + lowerOrderTerm;
    }
  };
} // namespace Kokkos
#endif // HAVE_TPETRA_INST_FLOAT128

namespace { // (anonymous)

  /// \brief Allocate and return a 2-D Kokkos::DualView for Tpetra::MultiVector.
  ///
  /// This function takes the same first four template parameters as
  /// Tpetra::MultiVector.
  ///
  /// \param lclNumRows [in] Number of rows in the DualView.
  ///   "Local" means "local to the calling MPI process."
  /// \param numCols [in] Number of columns in the DualView.
  /// \param zeroOut [in] Whether to initialize all the entries of the
  ///   DualView to zero.  Kokkos does first-touch initialization.
  /// \param allowPadding [in] Whether to give Kokkos the option to
  ///   pad Views for alignment.
  ///
  /// \return The allocated Kokkos::DualView.
  template<class ST, class LO, class GO, class NT>
  typename Tpetra::MultiVector<ST, LO, GO, NT>::dual_view_type
  allocDualView (const size_t lclNumRows,
                 const size_t numCols,
                 const bool zeroOut = true,
                 const bool allowPadding = false)
  {
    using ::Tpetra::Details::Behavior;
    using Kokkos::AllowPadding;
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    typedef typename Tpetra::MultiVector<ST, LO, GO, NT>::dual_view_type dual_view_type;
    typedef typename dual_view_type::t_dev dev_view_type;
    // This needs to be a string and not a char*, if given as an
    // argument to Kokkos::view_alloc.  This is because view_alloc
    // also allows a raw pointer as its first argument.  See
    // https://github.com/kokkos/kokkos/issues/434.
    const std::string label ("MV::DualView");
    const bool debug = Behavior::debug ();

    // NOTE (mfh 18 Feb 2015, 12 Apr 2015, 22 Sep 2016) Our separate
    // creation of the DualView's Views works around
    // Kokkos::DualView's current inability to accept an
    // AllocationProperties initial argument (as Kokkos::View does).
    // However, the work-around is harmless, since it does what the
    // (currently nonexistent) equivalent DualView constructor would
    // have done anyway.

    dev_view_type d_view;
    if (zeroOut) {
      if (allowPadding) {
        d_view = dev_view_type (view_alloc (label, AllowPadding),
                                lclNumRows, numCols);
      }
      else {
        d_view = dev_view_type (view_alloc (label),
                                lclNumRows, numCols);
      }
    }
    else {
      if (allowPadding) {
        d_view = dev_view_type (view_alloc (label,
                                            WithoutInitializing,
                                            AllowPadding),
                                lclNumRows, numCols);
      }
      else {
        d_view = dev_view_type (view_alloc (label, WithoutInitializing),
                                lclNumRows, numCols);
      }
      if (debug) {
        // Filling with NaN is a cheap and effective way to tell if
        // downstream code is trying to use a MultiVector's data
        // without them having been initialized.  ArithTraits lets us
        // call nan() even if the scalar type doesn't define it; it
        // just returns some undefined value in the latter case.  This
        // won't hurt anything because by setting zeroOut=false, users
        // already agreed that they don't care about the contents of
        // the MultiVector.
        const ST nan = Kokkos::Details::ArithTraits<ST>::nan ();
        KokkosBlas::fill (d_view, nan);
      }
    }
    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (d_view.extent (0)) != lclNumRows ||
         static_cast<size_t> (d_view.extent (1)) != numCols, std::logic_error,
         "allocDualView: d_view's dimensions actual dimensions do not match "
         "requested dimensions.  d_view is " << d_view.extent (0) << " x " <<
         d_view.extent (1) << "; requested " << lclNumRows << " x " << numCols
         << ".  Please report this bug to the Tpetra developers.");
    }

    dual_view_type dv (d_view, Kokkos::create_mirror_view (d_view));
    // Whether or not the user cares about the initial contents of the
    // MultiVector, the device and host views are out of sync.  We
    // prefer to work in device memory.  The way to ensure this
    // happens is to mark the device view as modified.
    dv.modify_device ();

    return dv;
  }

  // Convert 1-D Teuchos::ArrayView to an unmanaged 1-D host Kokkos::View.
  //
  // T: The type of the entries of the View.
  // ExecSpace: The Kokkos execution space.
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
    //
    // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
    // a memory space, rather than an execution space, as the first
    // argument of VerifyExecutionCanAccessMemorySpace.
    typedef typename Kokkos::Impl::if_c<
      Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
        typename ExecSpace::memory_space,
        Kokkos::HostSpace>::value,
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

  // mfh 14 Apr 2015: Work-around for bug in Kokkos::subview, where
  // taking a subview of a 0 x N DualView incorrectly always results
  // in a 0 x 0 DualView.
  template<class DualViewType>
  DualViewType
  takeSubview (const DualViewType& X,
               const Kokkos::Impl::ALL_t&,
               const std::pair<size_t, size_t>& colRng)
  {
    if (X.extent (0) == 0 && X.extent (1) != 0) {
      return DualViewType ("MV::DualView", 0, colRng.second - colRng.first);
    }
    else {
      return subview (X, Kokkos::ALL (), colRng);
    }
  }

  // mfh 14 Apr 2015: Work-around for bug in Kokkos::subview, where
  // taking a subview of a 0 x N DualView incorrectly always results
  // in a 0 x 0 DualView.
  template<class DualViewType>
  DualViewType
  takeSubview (const DualViewType& X,
               const std::pair<size_t, size_t>& rowRng,
               const std::pair<size_t, size_t>& colRng)
  {
    if (X.extent (0) == 0 && X.extent (1) != 0) {
      return DualViewType ("MV::DualView", 0, colRng.second - colRng.first);
    }
    else {
      return subview (X, rowRng, colRng);
    }
  }

  template<class DualViewType>
  size_t
  getDualViewStride (const DualViewType& dv)
  {
    // FIXME (mfh 15 Mar 2019) DualView doesn't have a stride
    // method yet, but its Views do.
    const size_t LDA = dv.d_view.stride (1);
    const size_t numRows = dv.extent (0);

    if (LDA == 0) {
      return (numRows == 0) ? size_t (1) : numRows;
    }
    else {
      return LDA;
    }
  }

  template<class ViewType>
  size_t
  getViewStride (const ViewType& view)
  {
    const size_t LDA = view.stride (1);
    const size_t numRows = view.extent (0);

    if (LDA == 0) {
      return (numRows == 0) ? size_t (1) : numRows;
    }
    else {
      return LDA;
    }
  }

  template <class impl_scalar_type, class buffer_device_type>
  bool
  runKernelOnHost ( Kokkos::DualView<impl_scalar_type*, buffer_device_type> imports )
  {
    if (! imports.need_sync_device ()) {
      return false; // most up-to-date on device
    }
    else { // most up-to-date on host
      size_t localLengthThreshold = Tpetra::Details::Behavior::multivectorKernelLocationThreshold();
      return imports.extent(0) <= localLengthThreshold;
    }
  }


  template <class SC, class LO, class GO, class NT>
  bool
  runKernelOnHost (const ::Tpetra::MultiVector<SC, LO, GO, NT>& X)
  {
    if (! X.need_sync_device ()) {
      return false; // most up-to-date on device
    }
    else { // most up-to-date on host
      constexpr size_t localLengthThreshold = 10000;
      return X.getLocalLength () <= localLengthThreshold;
    }
  }

  template <class SC, class LO, class GO, class NT>
  void
  multiVectorNormImpl (typename ::Tpetra::MultiVector<SC, LO, GO, NT>::mag_type norms[],
                       ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
                       const ::Tpetra::Details::EWhichNorm whichNorm)
  {
    using ::Tpetra::Details::normImpl;
    using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
    using val_type = typename MV::impl_scalar_type;
    using mag_type = typename MV::mag_type;
    using dual_view_type = typename MV::dual_view_type;

    auto map = X.getMap ();
    auto comm = map.is_null () ? nullptr : map->getComm ().getRawPtr ();
    auto whichVecs = getMultiVectorWhichVectors (X);
    const bool isConstantStride = X.isConstantStride ();
    const bool isDistributed = X.isDistributed ();

    const bool runOnHost = runKernelOnHost (X);
    if (runOnHost) {
      using view_type = typename dual_view_type::t_host;
      using array_layout = typename view_type::array_layout;
      using device_type = typename view_type::device_type;

      if (X.need_sync_host ()) {
        X.sync_host ();
      }
      view_type X_lcl = X.getLocalViewHost ();
      normImpl<val_type, array_layout, device_type,
        mag_type> (norms, X_lcl, whichNorm, whichVecs,
                   isConstantStride, isDistributed, comm);
    }
    else {
      using view_type = typename dual_view_type::t_dev;
      using array_layout = typename view_type::array_layout;
      using device_type = typename view_type::device_type;

      if (X.need_sync_device ()) {
        X.sync_device ();
      }
      view_type X_lcl = X.getLocalViewDevice ();
      normImpl<val_type, array_layout, device_type,
        mag_type> (norms, X_lcl, whichNorm, whichVecs,
                   isConstantStride, isDistributed, comm);
    }
  }
} // namespace (anonymous)


namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  vectorIndexOutOfRange (const size_t VectorIndex) const {
    return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector () :
    base_type (Teuchos::rcp (new map_type ()))
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const size_t numVecs,
               const bool zeroOut) : /* default is true */
    base_type (map)
  {
    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV ctor (map,numVecs,zeroOut)");

    const size_t lclNumRows = this->getLocalLength ();
    view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, Node> (lclNumRows, numVecs, zeroOut);
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
               const Teuchos::DataAccess copyOrView) :
    base_type (source),
    view_ (source.view_),
    origView_ (source.origView_),
    whichVectors_ (source.whichVectors_)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "MultiVector(const MultiVector&, "
      "const Teuchos::DataAccess): ";

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV 2-arg \"copy\" ctor");

    if (copyOrView == Teuchos::Copy) {
      // Reuse the conveniently already existing function that creates
      // a deep copy.
      MV cpy = createCopy (source);
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const dual_view_type& view) :
    base_type (map),
    view_ (view),
    origView_ (view)
  {
    const char tfecfFuncName[] = "MultiVector(Map,DualView): ";
    const size_t lclNumRows_map = map.is_null () ? size_t (0) :
      map->getNodeNumElements ();
    const size_t lclNumRows_view = view.extent (0);
    const size_t LDA = getDualViewStride (view);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (LDA < lclNumRows_map || lclNumRows_map != lclNumRows_view,
       std::invalid_argument, "Kokkos::DualView does not match Map. "
       "map->getNodeNumElements() = " << lclNumRows_map
       << ", view.extent(0) = " << lclNumRows_view
       << ", and getStride() = " << LDA << ".");

    using ::Tpetra::Details::Behavior;
    const bool debug = Behavior::debug ();
    if (debug) {
      using ::Tpetra::Details::checkGlobalDualViewValidity;
      std::ostringstream gblErrStrm;
      const bool verbose = Behavior::verbose ();
      const auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const bool gblValid =
        checkGlobalDualViewValidity (&gblErrStrm, view, verbose,
                                     comm.getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! gblValid, std::runtime_error, gblErrStrm.str ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const typename dual_view_type::t_dev& d_view) :
    base_type (map)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "MultiVector(map,d_view): ";

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV ctor (map,d_view)");

    const size_t LDA = getViewStride (d_view);
    const size_t lclNumRows = map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (LDA < lclNumRows, std::invalid_argument, "Map does not match "
       "Kokkos::View.  map->getNodeNumElements() = " << lclNumRows
       << ", View's column stride = " << LDA
       << ", and View's extent(0) = " << d_view.extent (0) << ".");

    auto h_view = Kokkos::create_mirror_view (d_view);
    view_ = dual_view_type (d_view, h_view);

    using ::Tpetra::Details::Behavior;
    const bool debug = Behavior::debug ();
    if (debug) {
      using ::Tpetra::Details::checkGlobalDualViewValidity;
      std::ostringstream gblErrStrm;
      const bool verbose = Behavior::verbose ();
      const auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const bool gblValid =
        checkGlobalDualViewValidity (&gblErrStrm, view_, verbose,
                                     comm.getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! gblValid, std::runtime_error, gblErrStrm.str ());
    }
    // The user gave us a device view.  In order to respect its
    // initial contents, we mark the DualView as "modified on device."
    // That way, the next sync will synchronize it with the host view.
    this->modify_device ();
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const dual_view_type& view,
               const dual_view_type& origView) :
    base_type (map),
    view_ (view),
    origView_ (origView)
  {
    const char tfecfFuncName[] = "MultiVector(map,view,origView): ";

    const size_t LDA = getDualViewStride (origView);
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::DualView's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input origView's first dimension (number "
      "of rows = " << origView.extent (0) << ") does not not match the number "
      "of entries in the Map on the calling process.");

    using ::Tpetra::Details::Behavior;
    const bool debug = Behavior::debug ();
    if (debug) {
      using ::Tpetra::Details::checkGlobalDualViewValidity;
      std::ostringstream gblErrStrm;
      const bool verbose = Behavior::verbose ();
      const auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const bool gblValid_0 =
        checkGlobalDualViewValidity (&gblErrStrm, view, verbose,
                                     comm.getRawPtr ());
      const bool gblValid_1 =
        checkGlobalDualViewValidity (&gblErrStrm, origView, verbose,
                                     comm.getRawPtr ());
      const bool gblValid = gblValid_0 && gblValid_1;
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! gblValid, std::runtime_error, gblErrStrm.str ());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
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

    using ::Tpetra::Details::Behavior;
    const bool debug = Behavior::debug ();
    if (debug) {
      using ::Tpetra::Details::checkGlobalDualViewValidity;
      std::ostringstream gblErrStrm;
      const bool verbose = Behavior::verbose ();
      const auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const bool gblValid =
        checkGlobalDualViewValidity (&gblErrStrm, view, verbose,
                                     comm.getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! gblValid, std::runtime_error, gblErrStrm.str ());
    }

    const size_t lclNumRows = map.is_null () ? size_t (0) :
      map->getNodeNumElements ();
    // Check dimensions of the input DualView.  We accept that Kokkos
    // might not allow construction of a 0 x m (Dual)View with m > 0,
    // so we only require the number of rows to match if the
    // (Dual)View has more than zero columns.  Likewise, we only
    // require the number of columns to match if the (Dual)View has
    // more than zero rows.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      view.extent (1) != 0 && static_cast<size_t> (view.extent (0)) < lclNumRows,
      std::invalid_argument, "view.extent(0) = " << view.extent (0)
      << " < map->getNodeNumElements() = " << lclNumRows << ".");
    if (whichVectors.size () != 0) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        view.extent (1) != 0 && view.extent (1) == 0,
        std::invalid_argument, "view.extent(1) = 0, but whichVectors.size()"
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
        view.extent (1) != 0 && static_cast<size_t> (view.extent (1)) <= maxColInd,
        std::invalid_argument, "view.extent(1) = " << view.extent (1)
        << " <= max(whichVectors) = " << maxColInd << ".");
    }

    // If extent(1) is 0, the stride might be 0.  BLAS doesn't like
    // zero strides, so modify in that case.
    const size_t LDA = getDualViewStride (view);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (LDA < lclNumRows, std::invalid_argument,
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
      view_ = takeSubview (view_, ALL (), colRng);
      origView_ = takeSubview (origView_, ALL (), colRng);
      // whichVectors_.size() == 0 means "constant stride."
      whichVectors_.clear ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
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

    using ::Tpetra::Details::Behavior;
    const bool debug = Behavior::debug ();
    if (debug) {
      using ::Tpetra::Details::checkGlobalDualViewValidity;
      std::ostringstream gblErrStrm;
      const bool verbose = Behavior::verbose ();
      const auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const bool gblValid_0 =
        checkGlobalDualViewValidity (&gblErrStrm, view, verbose,
                                     comm.getRawPtr ());
      const bool gblValid_1 =
        checkGlobalDualViewValidity (&gblErrStrm, origView, verbose,
                                     comm.getRawPtr ());
      const bool gblValid = gblValid_0 && gblValid_1;
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! gblValid, std::runtime_error, gblErrStrm.str ());
    }

    const size_t lclNumRows = this->getLocalLength ();
    // Check dimensions of the input DualView.  We accept that Kokkos
    // might not allow construction of a 0 x m (Dual)View with m > 0,
    // so we only require the number of rows to match if the
    // (Dual)View has more than zero columns.  Likewise, we only
    // require the number of columns to match if the (Dual)View has
    // more than zero rows.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      view.extent (1) != 0 && static_cast<size_t> (view.extent (0)) < lclNumRows,
      std::invalid_argument, "view.extent(0) = " << view.extent (0)
      << " < map->getNodeNumElements() = " << lclNumRows << ".");
    if (whichVectors.size () != 0) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        view.extent (1) != 0 && view.extent (1) == 0,
        std::invalid_argument, "view.extent(1) = 0, but whichVectors.size()"
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
        view.extent (1) != 0 && static_cast<size_t> (view.extent (1)) <= maxColInd,
        std::invalid_argument, "view.extent(1) = " << view.extent (1)
        << " <= max(whichVectors) = " << maxColInd << ".");
    }

    // If extent(1) is 0, the stride might be 0.  BLAS doesn't like
    // zero strides, so modify in that case.
    const size_t LDA = getDualViewStride (origView);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (LDA < lclNumRows, std::invalid_argument, "Map and DualView origView "
       "do not match.  LDA = " << LDA << " < this->getLocalLength() = " <<
       lclNumRows << ".  origView.extent(0) = " << origView.extent(0)
       << ", origView.stride(1) = " << origView.d_view.stride(1) << ".");

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
      view_ = takeSubview (view_, ALL (), colRng);
      origView_ = takeSubview (origView_, ALL (), colRng);
      // whichVectors_.size() == 0 means "constant stride."
      whichVectors_.clear ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const Teuchos::ArrayView<const Scalar>& data,
               const size_t LDA,
               const size_t numVecs) :
    base_type (map)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "MultiVector(map,data,LDA,numVecs): ";
    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV ctor (map,Teuchos::ArrayView,LDA,numVecs)");

    // Deep copy constructor, constant stride (NO whichVectors_).
    // There is no need for a deep copy constructor with nonconstant stride.

    const size_t lclNumRows =
      map.is_null () ? size_t (0) : map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (LDA < lclNumRows, std::invalid_argument, "LDA = " << LDA << " < "
       "map->getNodeNumElements() = " << lclNumRows << ".");
    if (numVecs != 0) {
      const size_t minNumEntries = LDA * (numVecs - 1) + lclNumRows;
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (data.size ()) < minNumEntries,
         std::invalid_argument, "Input Teuchos::ArrayView does not have enough "
         "entries, given the input Map and number of vectors in the MultiVector."
         "  data.size() = " << data.size () << " < (LDA*(numVecs-1)) + "
         "map->getNodeNumElements () = " << minNumEntries << ".");
    }

    this->view_ = allocDualView<Scalar, LO, GO, Node> (lclNumRows, numVecs);
    this->modify_device ();
    auto X_out = this->getLocalViewDevice ();
    origView_ = view_;

    // Make an unmanaged host Kokkos::View of the input data.  First
    // create a View (X_in_orig) with the original stride.  Then,
    // create a subview (X_in) with the right number of columns.
    const impl_scalar_type* const X_in_raw =
      reinterpret_cast<const impl_scalar_type*> (data.getRawPtr ());
    Kokkos::View<const impl_scalar_type**,
      Kokkos::LayoutLeft,
      Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged> X_in_orig (X_in_raw, LDA, numVecs);
    const Kokkos::pair<size_t, size_t> rowRng (0, lclNumRows);
    auto X_in = Kokkos::subview (X_in_orig, rowRng, Kokkos::ALL ());

    // If LDA != X_out's column stride, then we need to copy one
    // column at a time; Kokkos::deep_copy refuses to work in that
    // case.
    const size_t outStride =
      X_out.extent (1) == 0 ? size_t (1) : X_out.stride (1);
    if (LDA == outStride) { // strides are the same; deep_copy once
      // This only works because MultiVector uses LayoutLeft.
      // We would need a custom copy functor otherwise.
      Kokkos::deep_copy (X_out, X_in);
    }
    else { // strides differ; copy one column at a time
      typedef decltype (Kokkos::subview (X_out, Kokkos::ALL (), 0))
        out_col_view_type;
      typedef decltype (Kokkos::subview (X_in, Kokkos::ALL (), 0))
        in_col_view_type;
      for (size_t j = 0; j < numVecs; ++j) {
        out_col_view_type X_out_j = Kokkos::subview (X_out, Kokkos::ALL (), j);
        in_col_view_type X_in_j = Kokkos::subview (X_in, Kokkos::ALL (), j);
        Kokkos::deep_copy (X_out_j, X_in_j);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const Teuchos::RCP<const map_type>& map,
               const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >& ArrayOfPtrs,
               const size_t numVecs) :
    base_type (map)
  {
    typedef impl_scalar_type IST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "MultiVector(map,ArrayOfPtrs,numVecs): ";
    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV ctor (map,Teuchos::ArrayView of ArrayView,numVecs)");

    const size_t lclNumRows =
      map.is_null () ? size_t (0) : map->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numVecs < 1 || numVecs != static_cast<size_t> (ArrayOfPtrs.size ()),
       std::runtime_error, "Either numVecs (= " << numVecs << ") < 1, or "
       "ArrayOfPtrs.size() (= " << ArrayOfPtrs.size () << ") != numVecs.");
    for (size_t j = 0; j < numVecs; ++j) {
      Teuchos::ArrayView<const Scalar> X_j_av = ArrayOfPtrs[j];
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        static_cast<size_t> (X_j_av.size ()) < lclNumRows,
        std::invalid_argument, "ArrayOfPtrs[" << j << "].size() = "
        << X_j_av.size () << " < map->getNodeNumElements() = " << lclNumRows
        << ".");
    }

    view_ = allocDualView<Scalar, LO, GO, Node> (lclNumRows, numVecs);
    this->modify_device ();
    auto X_out = this->getLocalViewDevice ();

    // Make sure that the type of a single input column has the same
    // array layout as each output column does, so we can deep_copy.
    using array_layout = typename decltype (X_out)::array_layout;
    using input_col_view_type = typename Kokkos::View<const IST*,
      array_layout,
      Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged>;

    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    for (size_t j = 0; j < numVecs; ++j) {
      Teuchos::ArrayView<const IST> X_j_av =
        Teuchos::av_reinterpret_cast<const IST> (ArrayOfPtrs[j]);
      input_col_view_type X_j_in (X_j_av.getRawPtr (), lclNumRows);
      auto X_j_out = Kokkos::subview (X_out, rowRng, j);
      Kokkos::deep_copy (X_j_out, X_j_in);
    }
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isConstantStride () const {
    return whichVectors_.empty ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalLength () const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getNodeNumElements ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalLength () const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getGlobalNumElements ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getStride () const
  {
    return isConstantStride () ? getDualViewStride (origView_) : size_t (0);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  checkSizes (const SrcDistObject& sourceObj)
  {
    // Check whether the source object is a MultiVector.  If not, then
    // we can't even compare sizes, so it's definitely not OK to
    // Import or Export from it.
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
    const this_type* src = dynamic_cast<const this_type*> (&sourceObj);
    if (src == nullptr) {
      return false;
    }
    else {
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
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  constantNumberOfPackets () const {
    return this->getNumVectors ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermute
  (const SrcDistObject& sourceObj,
   const size_t numSameIDs,
   const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteToLIDs,
   const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteFromLIDs)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::getDualViewCopyFromArrayView;
    using ::Tpetra::Details::ProfilingRegion;
    using std::endl;
    using KokkosRefactor::Details::permute_array_multi_column;
    using KokkosRefactor::Details::permute_array_multi_column_variable_stride;
    using Kokkos::Compat::create_const_view;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char tfecfFuncName[] = "copyAndPermute: ";
    ProfilingRegion regionCAP ("Tpetra::MultiVector::copyAndPermute");

    const bool verbose = Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      auto map = this->getMap ();
      auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": MV::copyAndPermute: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "Start" << endl;
      std::cerr << os.str ();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (permuteToLIDs.extent (0) != permuteFromLIDs.extent (0),
       std::logic_error, "permuteToLIDs.extent(0) = "
       << permuteToLIDs.extent (0) << " != permuteFromLIDs.extent(0) = "
       << permuteFromLIDs.extent (0) << ".");

    // We've already called checkSizes(), so this cast must succeed.
    MV& sourceMV = const_cast<MV &>(dynamic_cast<const MV&> (sourceObj));
    const size_t numCols = this->getNumVectors ();

    // sourceMV doesn't belong to us, so we can't sync it.  Do the
    // copying where it's currently sync'd.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (sourceMV.need_sync_device () && sourceMV.need_sync_host (),
       std::logic_error, "Input MultiVector needs sync to both host "
       "and device.");
    const bool copyOnHost = runKernelOnHost(sourceMV);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "copyOnHost=" << (copyOnHost ? "true" : "false") << endl;
      std::cerr << os.str ();
    }

    if (copyOnHost) {
      sourceMV.sync_host();
      this->sync_host ();
      this->modify_host ();
    }
    else {
      sourceMV.sync_device();
      if (this->need_sync_device ()) {
        this->sync_device ();
      }
      this->modify_device ();
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Copy" << endl;
      std::cerr << os.str ();
    }

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
    // Note (ETP 2 Jul 2014)  We need to always copy one column at a
    // time, even when both multivectors are constant-stride, since
    // deep_copy between strided subviews with more than one column
    // doesn't currently work.

    // FIXME (mfh 04 Feb 2019) Need to optimize for the case where
    // both source and target are constant stride and have multiple
    // columns.
    if (numSameIDs > 0) {
      const std::pair<size_t, size_t> rows (0, numSameIDs);
      if (copyOnHost) {
        auto tgt_h = this->getLocalViewHost ();
        auto src_h = create_const_view (sourceMV.getLocalViewHost ());

        for (size_t j = 0; j < numCols; ++j) {
          const size_t tgtCol = isConstantStride () ? j : whichVectors_[j];
          const size_t srcCol =
            sourceMV.isConstantStride () ? j : sourceMV.whichVectors_[j];

          auto tgt_j = Kokkos::subview (tgt_h, rows, tgtCol);
          auto src_j = Kokkos::subview (src_h, rows, srcCol);
          Kokkos::deep_copy (tgt_j, src_j); // Copy src_j into tgt_j
        }
      }
      else { // copy on device
        auto tgt_d = this->getLocalViewDevice ();
        auto src_d = create_const_view (sourceMV.getLocalViewDevice ());

        for (size_t j = 0; j < numCols; ++j) {
          const size_t tgtCol = isConstantStride () ? j : whichVectors_[j];
          const size_t srcCol =
            sourceMV.isConstantStride () ? j : sourceMV.whichVectors_[j];

          auto tgt_j = Kokkos::subview (tgt_d, rows, tgtCol);
          auto src_j = Kokkos::subview (src_d, rows, srcCol);
          Kokkos::deep_copy (tgt_j, src_j); // Copy src_j into tgt_j
        }
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
    if (permuteFromLIDs.extent (0) == 0 ||
        permuteToLIDs.extent (0) == 0) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "No permutations. Done!" << endl;
        std::cerr << os.str ();
      }
      return;
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Permute" << endl;
      std::cerr << os.str ();
    }

    // We could in theory optimize for the case where exactly one of
    // them is constant stride, but we don't currently do that.
    const bool nonConstStride =
      ! this->isConstantStride () || ! sourceMV.isConstantStride ();

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "nonConstStride="
         << (nonConstStride ? "true" : "false") << endl;
      std::cerr << os.str ();
    }

    // We only need the "which vectors" arrays if either the source or
    // target MV is not constant stride.  Since we only have one
    // kernel that must do double-duty, we have to create a "which
    // vectors" array for the MV that _is_ constant stride.
    Kokkos::DualView<const size_t*, device_type> srcWhichVecs;
    Kokkos::DualView<const size_t*, device_type> tgtWhichVecs;
    if (nonConstStride) {
      if (this->whichVectors_.size () == 0) {
        Kokkos::DualView<size_t*, device_type> tmpTgt ("tgtWhichVecs", numCols);
        tmpTgt.modify_host ();
        for (size_t j = 0; j < numCols; ++j) {
          tmpTgt.h_view(j) = j;
        }
        if (! copyOnHost) {
          tmpTgt.sync_device ();
        }
        tgtWhichVecs = tmpTgt;
      }
      else {
        Teuchos::ArrayView<const size_t> tgtWhichVecsT = this->whichVectors_ ();
        tgtWhichVecs =
          getDualViewCopyFromArrayView<size_t, device_type> (tgtWhichVecsT,
                                                             "tgtWhichVecs",
                                                             copyOnHost);
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (tgtWhichVecs.extent (0)) !=
         this->getNumVectors (),
         std::logic_error, "tgtWhichVecs.extent(0) = " <<
         tgtWhichVecs.extent (0) << " != this->getNumVectors() = " <<
         this->getNumVectors () << ".");

      if (sourceMV.whichVectors_.size () == 0) {
        Kokkos::DualView<size_t*, device_type> tmpSrc ("srcWhichVecs", numCols);
        tmpSrc.modify_host ();
        for (size_t j = 0; j < numCols; ++j) {
          tmpSrc.h_view(j) = j;
        }
        if (! copyOnHost) {
          tmpSrc.sync_device ();
        }
        srcWhichVecs = tmpSrc;
      }
      else {
        Teuchos::ArrayView<const size_t> srcWhichVecsT =
          sourceMV.whichVectors_ ();
        srcWhichVecs =
          getDualViewCopyFromArrayView<size_t, device_type> (srcWhichVecsT,
                                                             "srcWhichVecs",
                                                             copyOnHost);
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (srcWhichVecs.extent (0)) !=
         sourceMV.getNumVectors (), std::logic_error,
         "srcWhichVecs.extent(0) = " << srcWhichVecs.extent (0)
         << " != sourceMV.getNumVectors() = " << sourceMV.getNumVectors ()
         << ".");
    }

    if (copyOnHost) { // permute on host too
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Get permute LIDs on host" << std::endl;
        std::cerr << os.str ();
      }
      auto tgt_h = this->getLocalViewHost ();
      auto src_h = create_const_view (sourceMV.getLocalViewHost ());

      TEUCHOS_ASSERT( ! permuteToLIDs.need_sync_host () );
      auto permuteToLIDs_h = create_const_view (permuteToLIDs.view_host ());
      TEUCHOS_ASSERT( ! permuteFromLIDs.need_sync_host () );
      auto permuteFromLIDs_h =
        create_const_view (permuteFromLIDs.view_host ());

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Permute on host" << endl;
        std::cerr << os.str ();
      }
      if (nonConstStride) {
        // No need to sync first, because copyOnHost argument to
        // getDualViewCopyFromArrayView puts them in the right place.
        auto tgtWhichVecs_h =
          create_const_view (tgtWhichVecs.view_host ());
        auto srcWhichVecs_h =
          create_const_view (srcWhichVecs.view_host ());
        permute_array_multi_column_variable_stride (tgt_h, src_h,
                                                    permuteToLIDs_h,
                                                    permuteFromLIDs_h,
                                                    tgtWhichVecs_h,
                                                    srcWhichVecs_h, numCols);
      }
      else {
        permute_array_multi_column (tgt_h, src_h, permuteToLIDs_h,
                                    permuteFromLIDs_h, numCols);
      }
    }
    else { // permute on device
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Get permute LIDs on device" << endl;
        std::cerr << os.str ();
      }
      auto tgt_d = this->getLocalViewDevice ();
      auto src_d = create_const_view (sourceMV.getLocalViewDevice ());

      TEUCHOS_ASSERT( ! permuteToLIDs.need_sync_device () );
      auto permuteToLIDs_d = create_const_view (permuteToLIDs.view_device ());
      TEUCHOS_ASSERT( ! permuteFromLIDs.need_sync_device () );
      auto permuteFromLIDs_d =
        create_const_view (permuteFromLIDs.view_device ());

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Permute on device" << endl;
        std::cerr << os.str ();
      }
      if (nonConstStride) {
        // No need to sync first, because copyOnHost argument to
        // getDualViewCopyFromArrayView puts them in the right place.
        auto tgtWhichVecs_d = create_const_view (tgtWhichVecs.view_device ());
        auto srcWhichVecs_d = create_const_view (srcWhichVecs.view_device ());
        permute_array_multi_column_variable_stride (tgt_d, src_d,
                                                    permuteToLIDs_d,
                                                    permuteFromLIDs_d,
                                                    tgtWhichVecs_d,
                                                    srcWhichVecs_d, numCols);
      }
      else {
        permute_array_multi_column (tgt_d, src_d, permuteToLIDs_d,
                                    permuteFromLIDs_d, numCols);
      }
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done!" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packAndPrepare
  (const SrcDistObject& sourceObj,
   const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
   Kokkos::DualView<impl_scalar_type*, buffer_device_type>& exports,
   Kokkos::DualView<size_t*, buffer_device_type> /* numExportPacketsPerLID */,
   size_t& constantNumPackets,
   Distributor & /* distor */ )
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::reallocDualViewIfNeeded;
    using Kokkos::Compat::create_const_view;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using std::endl;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char tfecfFuncName[] = "packAndPrepare: ";
    ProfilingRegion regionPAP ("Tpetra::MultiVector::packAndPrepare");

    // mfh 09 Sep 2016, 26 Sep 2017: The pack and unpack functions now
    // have the option to check indices.  We do so when Tpetra is in
    // debug mode.  It is in debug mode by default in a debug build,
    // but you may control this at run time, before launching the
    // executable, by setting the TPETRA_DEBUG environment variable to
    // "1" (or "TRUE").
    const bool debugCheckIndices = Behavior::debug ();
    // mfh 03 Aug 2017, 27 Sep 2017: Set the TPETRA_VERBOSE
    // environment variable to "1" (or "TRUE") for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool printDebugOutput = Behavior::verbose ();

    std::unique_ptr<std::string> prefix;
    if (printDebugOutput) {
      auto map = this->getMap ();
      auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": MV::packAndPrepare: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "Start" << endl;
      std::cerr << os.str ();
    }

    // We've already called checkSizes(), so this cast must succeed.
    MV& sourceMV = const_cast<MV&>(dynamic_cast<const MV&> (sourceObj));

    const size_t numCols = sourceMV.getNumVectors ();

    // This spares us from needing to fill numExportPacketsPerLID.
    // Setting constantNumPackets to a nonzero value signals that
    // all packets have the same number of entries.
    constantNumPackets = numCols;

    // If we have no exports, there is nothing to do.  Make sure this
    // goes _after_ setting constantNumPackets correctly.
    if (exportLIDs.extent (0) == 0) {
      if (printDebugOutput) {
        std::ostringstream os;
        os << *prefix << "No exports on this proc, DONE" << std::endl;
        std::cerr << os.str ();
      }
      return;
    }

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

    const size_t numExportLIDs = exportLIDs.extent (0);
    const size_t newExportsSize = numCols * numExportLIDs;
    if (printDebugOutput) {
      std::ostringstream os;
      os << *prefix << "realloc: "
         << "numExportLIDs: " << numExportLIDs
         << ", exports.extent(0): " << exports.extent (0)
         << ", newExportsSize: " << newExportsSize << std::endl;
      std::cerr << os.str ();
    }
    reallocDualViewIfNeeded (exports, newExportsSize, "exports");

    // mfh 04 Feb 2019: sourceMV doesn't belong to us, so we can't
    // sync it.  Pack it where it's currently sync'd.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (sourceMV.need_sync_device () && sourceMV.need_sync_host (),
       std::logic_error, "Input MultiVector needs sync to both host "
       "and device.");
    const bool packOnHost = runKernelOnHost(sourceMV);
    auto src_dev = sourceMV.getLocalViewHost ();
    auto src_host = sourceMV.getLocalViewDevice ();
    if (printDebugOutput) {
      std::ostringstream os;
      os << *prefix << "packOnHost=" << (packOnHost ? "true" : "false") << endl;
      std::cerr << os.str ();
    }

    // Mark 'exports' here, since we might have resized it above.
    // Resizing currently requires calling the constructor, which
    // clears out the 'modified' flags.
    if (packOnHost) {
      // nde 06 Feb 2020: If 'exports' does not require resize
      // when reallocDualViewIfNeeded is called, the modified flags 
      // are not cleared out. This can result in host and device views
      // being out-of-sync, resuling in an error in exports.modify_* calls.
      // Clearing the sync flags prevents this possible case.
      exports.clear_sync_state ();
      exports.modify_host ();
      sourceMV.sync_host();
    }
    else {
      // nde 06 Feb 2020: If 'exports' does not require resize
      // when reallocDualViewIfNeeded is called, the modified flags 
      // are not cleared out. This can result in host and device views
      // being out-of-sync, resuling in an error in exports.modify_* calls.
      // Clearing the sync flags prevents this possible case.
      exports.clear_sync_state ();
      exports.modify_device ();
      sourceMV.sync_device();
    }

    if (numCols == 1) { // special case for one column only
      // MultiVector always represents a single column with constant
      // stride, but it doesn't hurt to implement both cases anyway.
      //
      // ETP:  I'm not sure I agree with the above statement.  Can't a single-
      // column multivector be a subview of another multi-vector, in which case
      // sourceMV.whichVectors_[0] != 0 ?  I think we have to handle that case
      // separately.
      //
      // mfh 18 Jan 2016: In answer to ETP's comment above:
      // MultiVector treats single-column MultiVectors created using a
      // "nonconstant stride constructor" as a special case, and makes
      // them constant stride (by making whichVectors_ have length 0).
      if (sourceMV.isConstantStride ()) {
        using KokkosRefactor::Details::pack_array_single_column;
        if (printDebugOutput) {
          std::ostringstream os;
          os << *prefix << "Pack numCols=1 const stride" << endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_single_column (exports.view_host (),
                                    create_const_view (src_host),
                                    exportLIDs.view_host (),
                                    0,
                                    debugCheckIndices);
        }
        else { // pack on device
          pack_array_single_column (exports.view_device (),
                                    create_const_view (src_dev),
                                    exportLIDs.view_device (),
                                    0,
                                    debugCheckIndices);
        }
      }
      else {
        using KokkosRefactor::Details::pack_array_single_column;
        if (printDebugOutput) {
          std::ostringstream os;
          os << *prefix << "Pack numCols=1 nonconst stride" << endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_single_column (exports.view_host (),
                                    create_const_view (src_host),
                                    exportLIDs.view_host (),
                                    sourceMV.whichVectors_[0],
                                    debugCheckIndices);
        }
        else { // pack on device
          pack_array_single_column (exports.view_device (),
                                    create_const_view (src_dev),
                                    exportLIDs.view_device (),
                                    sourceMV.whichVectors_[0],
                                    debugCheckIndices);
        }
      }
    }
    else { // the source MultiVector has multiple columns
      if (sourceMV.isConstantStride ()) {
        using KokkosRefactor::Details::pack_array_multi_column;
        if (printDebugOutput) {
          std::ostringstream os;
          os << *prefix << "Pack numCols=" << numCols << " const stride" << endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_multi_column (exports.view_host (),
                                   create_const_view (src_host),
                                   exportLIDs.view_host (),
                                   numCols,
                                   debugCheckIndices);
        }
        else { // pack on device
          pack_array_multi_column (exports.view_device (),
                                   create_const_view (src_dev),
                                   exportLIDs.view_device (),
                                   numCols,
                                   debugCheckIndices);
        }
      }
      else {
        using KokkosRefactor::Details::pack_array_multi_column_variable_stride;
        if (printDebugOutput) {
          std::ostringstream os;
          os << *prefix << "Pack numCols=" << numCols << " nonconst stride"
             << endl;
          std::cerr << os.str ();
        }
        // FIXME (mfh 04 Feb 2019) Creating a Kokkos::View for
        // whichVectors_ can be expensive, but pack and unpack for
        // nonconstant-stride MultiVectors is slower anyway.
        using IST = impl_scalar_type;
        using DV = Kokkos::DualView<IST*, device_type>;
        using HES = typename DV::t_host::execution_space;
        using DES = typename DV::t_dev::execution_space;
        Teuchos::ArrayView<const size_t> whichVecs = sourceMV.whichVectors_ ();
        if (packOnHost) {
          pack_array_multi_column_variable_stride
            (exports.view_host (),
             create_const_view (src_host),
             exportLIDs.view_host (),
             getKokkosViewDeepCopy<HES> (whichVecs),
             numCols,
             debugCheckIndices);
        }
        else { // pack on device
          pack_array_multi_column_variable_stride
            (exports.view_device (),
             create_const_view (src_dev),
             exportLIDs.view_device (),
             getKokkosViewDeepCopy<DES> (whichVecs),
             numCols,
             debugCheckIndices);
        }
      }
    }

    if (printDebugOutput) {
      std::ostringstream os;
      os << *prefix << "Done!" << endl;
      std::cerr << os.str ();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombine
  (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& importLIDs,
   Kokkos::DualView<impl_scalar_type*, buffer_device_type> imports,
   Kokkos::DualView<size_t*, buffer_device_type> /* numPacketsPerLID */,
   const size_t constantNumPackets,
   Distributor& /* distor */,
   const CombineMode CM)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::ProfilingRegion;
    using KokkosRefactor::Details::unpack_array_multi_column;
    using KokkosRefactor::Details::unpack_array_multi_column_variable_stride;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using std::endl;
    using IST = impl_scalar_type;
    const char longFuncName[] = "Tpetra::MultiVector::unpackAndCombine";
    const char tfecfFuncName[] = "unpackAndCombine: ";
    ProfilingRegion regionUAC (longFuncName);

    // mfh 09 Sep 2016, 26 Sep 2017: The pack and unpack functions now
    // have the option to check indices.  We do so when Tpetra is in
    // debug mode.  It is in debug mode by default in a debug build,
    // but you may control this at run time, before launching the
    // executable, by setting the TPETRA_DEBUG environment variable to
    // "1" (or "TRUE").
    const bool debugCheckIndices = Behavior::debug ();

    const bool printDebugOutput = Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (printDebugOutput) {
      auto map = this->getMap ();
      auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const int myRank = comm.is_null () ? -1 : comm->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": " << longFuncName << ": ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << "Start" << endl;
      std::cerr << os.str ();
    }

    // If we have no imports, there is nothing to do
    if (importLIDs.extent (0) == 0) {
      if (printDebugOutput) {
        std::ostringstream os;
        os << *prefix << "No imports. Done!" << endl;
      }
      return;
    }

    const size_t numVecs = getNumVectors ();
    if (debugCheckIndices) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (imports.extent (0)) !=
         numVecs * importLIDs.extent (0),
         std::runtime_error,
         "imports.extent(0) = " << imports.extent (0)
         << " != getNumVectors() * importLIDs.extent(0) = " << numVecs
         << " * " << importLIDs.extent (0) << " = "
         << numVecs * importLIDs.extent (0) << ".");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (constantNumPackets == static_cast<size_t> (0), std::runtime_error,
         "constantNumPackets input argument must be nonzero.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (numVecs) !=
         static_cast<size_t> (constantNumPackets),
         std::runtime_error, "constantNumPackets must equal numVecs.");
    }

    // mfh 12 Apr 2016, 04 Feb 2019: Decide where to unpack based on
    // the memory space in which the imports buffer was last modified.
    // DistObject::doTransferNew gets to decide this.
    const bool unpackOnHost = runKernelOnHost(imports);

    if (printDebugOutput) {
      std::ostringstream os;
      os << *prefix << "unpackOnHost=" << (unpackOnHost ? "true" : "false")
         << endl;
      std::cerr << os.str ();
    }

    // We have to sync before modifying, because this method may read
    // as well as write (depending on the CombineMode).
    if (unpackOnHost) {
      imports.sync_host();
      this->sync_host ();
      this->modify_host ();
    }
    else {
      imports.sync_device();
      this->sync_device ();
      this->modify_device ();
    }
    auto X_d = this->getLocalViewDevice ();
    auto X_h = this->getLocalViewHost ();
    auto imports_d = imports.view_device ();
    auto imports_h = imports.view_host ();
    auto importLIDs_d = importLIDs.view_device ();
    auto importLIDs_h = importLIDs.view_host ();

    Kokkos::DualView<size_t*, device_type> whichVecs;
    if (! isConstantStride ()) {
      Kokkos::View<const size_t*, Kokkos::HostSpace,
        Kokkos::MemoryUnmanaged> whichVecsIn (whichVectors_.getRawPtr (),
                                              numVecs);
      whichVecs = Kokkos::DualView<size_t*, device_type> ("whichVecs", numVecs);
      if (unpackOnHost) {
        whichVecs.modify_host ();
        Kokkos::deep_copy (whichVecs.view_host (), whichVecsIn);
      }
      else {
        whichVecs.modify_device ();
        Kokkos::deep_copy (whichVecs.view_device (), whichVecsIn);
      }
    }
    auto whichVecs_d = whichVecs.view_device ();
    auto whichVecs_h = whichVecs.view_host ();

    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    if (numVecs > 0 && importLIDs.extent (0) > 0) {
      using dev_exec_space = typename dual_view_type::t_dev::execution_space;
      using host_exec_space = typename dual_view_type::t_host::execution_space;

      // This fixes GitHub Issue #4418.
      const bool use_atomic_updates = unpackOnHost ?
        host_exec_space::concurrency () != 1 :
        dev_exec_space::concurrency () != 1;

      if (printDebugOutput) {
        std::ostringstream os;
        os << *prefix << "Unpack: " << combineModeToString (CM) << endl;
        std::cerr << os.str ();
      }

      // NOTE (mfh 10 Mar 2012, 24 Mar 2014) If you want to implement
      // custom combine modes, start editing here.

      if (CM == INSERT || CM == REPLACE) {
        using op_type = KokkosRefactor::Details::InsertOp<IST>;
        if (isConstantStride ()) {
          if (unpackOnHost) {
            unpack_array_multi_column (host_exec_space (),
                                       X_h, imports_h, importLIDs_h,
                                       op_type (), numVecs,
                                       use_atomic_updates,
                                       debugCheckIndices);

          }
          else { // unpack on device
            unpack_array_multi_column (dev_exec_space (),
                                       X_d, imports_d, importLIDs_d,
                                       op_type (), numVecs,
                                       use_atomic_updates,
                                       debugCheckIndices);
          }
        }
        else { // not constant stride
          if (unpackOnHost) {
            unpack_array_multi_column_variable_stride (host_exec_space (),
                                                       X_h, imports_h,
                                                       importLIDs_h,
                                                       whichVecs_h,
                                                       op_type (),
                                                       numVecs,
                                                       use_atomic_updates,
                                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column_variable_stride (dev_exec_space (),
                                                       X_d, imports_d,
                                                       importLIDs_d,
                                                       whichVecs_d,
                                                       op_type (),
                                                       numVecs,
                                                       use_atomic_updates,
                                                       debugCheckIndices);
          }
        }
      }
      else if (CM == ADD) {
        using op_type = KokkosRefactor::Details::AddOp<IST>;
        if (isConstantStride ()) {
          if (unpackOnHost) {
            unpack_array_multi_column (host_exec_space (),
                                       X_h, imports_h, importLIDs_h,
                                       op_type (), numVecs,
                                       use_atomic_updates,
                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column (dev_exec_space (),
                                       X_d, imports_d, importLIDs_d,
                                       op_type (), numVecs,
                                       use_atomic_updates,
                                       debugCheckIndices);
          }
        }
        else { // not constant stride
          if (unpackOnHost) {
            unpack_array_multi_column_variable_stride (host_exec_space (),
                                                       X_h, imports_h,
                                                       importLIDs_h,
                                                       whichVecs_h,
                                                       op_type (),
                                                       numVecs,
                                                       use_atomic_updates,
                                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column_variable_stride (dev_exec_space (),
                                                       X_d, imports_d,
                                                       importLIDs_d,
                                                       whichVecs_d,
                                                       op_type (),
                                                       numVecs,
                                                       use_atomic_updates,
                                                       debugCheckIndices);
          }
        }
      }
      else if (CM == ABSMAX) {
        using op_type = KokkosRefactor::Details::AbsMaxOp<IST>;
        if (isConstantStride ()) {
          if (unpackOnHost) {
            unpack_array_multi_column (host_exec_space (),
                                       X_h, imports_h, importLIDs_h,
                                       op_type (), numVecs,
                                       use_atomic_updates,
                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column (dev_exec_space (),
                                       X_d, imports_d, importLIDs_d,
                                       op_type (), numVecs,
                                       use_atomic_updates,
                                       debugCheckIndices);
          }
        }
        else {
          if (unpackOnHost) {
            unpack_array_multi_column_variable_stride (host_exec_space (),
                                                       X_h, imports_h,
                                                       importLIDs_h,
                                                       whichVecs_h,
                                                       op_type (),
                                                       numVecs,
                                                       use_atomic_updates,
                                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column_variable_stride (dev_exec_space (),
                                                       X_d, imports_d,
                                                       importLIDs_d,
                                                       whichVecs_d,
                                                       op_type (),
                                                       numVecs,
                                                       use_atomic_updates,
                                                       debugCheckIndices);
          }
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::logic_error, "Invalid CombineMode");
      }
    }
    else {
      if (printDebugOutput) {
        std::ostringstream os;
        os << *prefix << "Nothing to unpack" << endl;
        std::cerr << os.str ();
      }
    }

    if (printDebugOutput) {
      std::ostringstream os;
      os << *prefix << "Done!" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getNumVectors () const
  {
    if (isConstantStride ()) {
      return static_cast<size_t> (view_.extent (1));
    } else {
      return static_cast<size_t> (whichVectors_.size ());
    }
  }

  namespace { // (anonymous)

    template<class RV>
    void
    gblDotImpl (const RV& dotsOut,
                const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                const bool distributed)
    {
      using Teuchos::REDUCE_MAX;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      typedef typename RV::non_const_value_type dot_type;

      const size_t numVecs = dotsOut.extent (0);

      // If the MultiVector is distributed over multiple processes, do
      // the distributed (interprocess) part of the dot product.  We
      // assume that the MPI implementation can read from and write to
      // device memory.
      //
      // replaceMap() may have removed some processes.  Those
      // processes have a null Map.  They must not participate in any
      // collective operations.  We ask first whether the Map is null,
      // because isDistributed() defers that question to the Map.  We
      // still compute and return local dots for processes not
      // participating in collective operations; those probably don't
      // make any sense, but it doesn't hurt to do them, since it's
      // illegal to call dot() on those processes anyway.
      if (distributed && ! comm.is_null ()) {
        // The calling process only participates in the collective if
        // both the Map and its Comm on that process are nonnull.
        const int nv = static_cast<int> (numVecs);
        const bool commIsInterComm = ::Tpetra::Details::isInterComm (*comm);

        if (commIsInterComm) {
          // If comm is an intercomm, then we may not alias input and
          // output buffers, so we have to make a copy of the local
          // sum.
          typename RV::non_const_type lclDots (Kokkos::ViewAllocateWithoutInitializing ("tmp"), numVecs);
          Kokkos::deep_copy (lclDots, dotsOut);
          const dot_type* const lclSum = lclDots.data ();
          dot_type* const gblSum = dotsOut.data ();
          reduceAll<int, dot_type> (*comm, REDUCE_SUM, nv, lclSum, gblSum);
        }
        else {
          dot_type* const inout = dotsOut.data ();
          reduceAll<int, dot_type> (*comm, REDUCE_SUM, nv, inout, inout);
        }
      }
    }
  } // namespace (anonymous)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
       const Kokkos::View<dot_type*, Kokkos::HostSpace>& dots) const
  {
    using ::Tpetra::Details::Behavior;
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    // View of all the dot product results.
    typedef Kokkos::View<dot_type*, Kokkos::HostSpace> RV;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    typedef typename dual_view_type::t_dev XMV;
    const char tfecfFuncName[] = "Tpetra::MultiVector::dot: ";

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::dot (Kokkos::View)");

    const size_t numVecs = this->getNumVectors ();
    if (numVecs == 0) {
      return; // nothing to do
    }
    const size_t lclNumRows = this->getLocalLength ();
    const size_t numDots = static_cast<size_t> (dots.extent (0));
    const bool debug = Behavior::debug ();

    if (debug) {
      const bool compat = this->getMap ()->isCompatible (* (A.getMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! compat, std::invalid_argument, "'*this' MultiVector is not "
        "compatible with the input MultiVector A.  We only test for this "
        "in debug mode.");
    }

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<dot_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<dot_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != A.getLocalLength (), std::runtime_error,
      "MultiVectors do not have the same local length.  "
      "this->getLocalLength() = " << lclNumRows << " != "
      "A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs != A.getNumVectors (), std::runtime_error,
      "MultiVectors must have the same number of columns (vectors).  "
      "this->getNumVectors() = " << numVecs << " != "
      "A.getNumVectors() = " << A.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numDots != numVecs, std::runtime_error,
      "The output array 'dots' must have the same number of entries as the "
      "number of columns (vectors) in *this and A.  dots.extent(0) = " <<
      numDots << " != this->getNumVectors() = " << numVecs << ".");

    const std::pair<size_t, size_t> colRng (0, numVecs);
    RV dotsOut = subview (dots, colRng);
    RCP<const Comm<int> > comm = this->getMap ().is_null () ? null :
      this->getMap ()->getComm ();

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) {
      const_cast<MV&>(*this).sync_device ();
    }
    if (A.need_sync_device ()) {
      const_cast<MV&>(A).sync_device ();
    }

    auto thisView = this->getLocalViewDevice ();
    auto A_view = A.getLocalViewDevice ();

    ::Tpetra::Details::lclDot<RV, XMV> (dotsOut, thisView, A_view, lclNumRows, numVecs,
                     this->whichVectors_.getRawPtr (),
                     A.whichVectors_.getRawPtr (),
                     this->isConstantStride (), A.isConstantStride ());
    gblDotImpl (dotsOut, comm, this->isDistributed ());
  }

  namespace { // (anonymous)
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type
    multiVectorSingleColumnDot (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x,
                                const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y)
    {
      using ::Tpetra::Details::ProfilingRegion;
      using MV = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using dot_type = typename MV::dot_type;
      ProfilingRegion region ("Tpetra::multiVectorSingleColumnDot");

      auto map = x.getMap ();
      Teuchos::RCP<const Teuchos::Comm<int> > comm =
        map.is_null () ? Teuchos::null : map->getComm ();
      if (comm.is_null ()) {
        return Kokkos::ArithTraits<dot_type>::zero ();
      }
      else {
        using LO = LocalOrdinal;
        // The min just ensures that we don't overwrite memory that
        // doesn't belong to us, in the erroneous input case where x
        // and y have different numbers of rows.
        const LO lclNumRows = static_cast<LO> (std::min (x.getLocalLength (),
                                                         y.getLocalLength ()));
        const Kokkos::pair<LO, LO> rowRng (0, lclNumRows);
        dot_type lclDot = Kokkos::ArithTraits<dot_type>::zero ();
        dot_type gblDot = Kokkos::ArithTraits<dot_type>::zero ();

        // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
        if (x.need_sync_device ()) {
          x.sync_device ();
        }
        if (y.need_sync_device ()) {
          const_cast<MV&>(y).sync_device ();
        }

        x.modify_device ();
        auto x_2d = x.getLocalViewDevice ();
        auto x_1d = Kokkos::subview (x_2d, rowRng, 0);
        auto y_2d = y.getLocalViewDevice ();
        auto y_1d = Kokkos::subview (y_2d, rowRng, 0);
        lclDot = KokkosBlas::dot (x_1d, y_1d);

        if (x.isDistributed ()) {
          using Teuchos::outArg;
          using Teuchos::REDUCE_SUM;
          using Teuchos::reduceAll;
          reduceAll<int, dot_type> (*comm, REDUCE_SUM, lclDot, outArg (gblDot));
        }
        else {
          gblDot = lclDot;
        }
        return gblDot;
      }
    }
  } // namespace (anonymous)



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
       const Teuchos::ArrayView<dot_type>& dots) const
  {
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "dot: ";
    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::dot (Teuchos::ArrayView)");

    const size_t numVecs = this->getNumVectors ();
    const size_t lclNumRows = this->getLocalLength ();
    const size_t numDots = static_cast<size_t> (dots.size ());

    // FIXME (mfh 11 Jul 2014, 31 May 2017) These exception tests may
    // not necessarily be thrown on all processes consistently.  We
    // keep them for now, because MultiVector's unit tests insist on
    // them.  In the future, we should instead pass along error state
    // with the inner product.  We could do this by setting an extra
    // slot to Kokkos::Details::ArithTraits<dot_type>::one() on error.
    // The final sum should be
    // Kokkos::Details::ArithTraits<dot_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclNumRows != A.getLocalLength (), std::runtime_error,
       "MultiVectors do not have the same local length.  "
       "this->getLocalLength() = " << lclNumRows << " != "
       "A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numVecs != A.getNumVectors (), std::runtime_error,
       "MultiVectors must have the same number of columns (vectors).  "
       "this->getNumVectors() = " << numVecs << " != "
       "A.getNumVectors() = " << A.getNumVectors () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numDots != numVecs, std::runtime_error,
       "The output array 'dots' must have the same number of entries as the "
       "number of columns (vectors) in *this and A.  dots.extent(0) = " <<
       numDots << " != this->getNumVectors() = " << numVecs << ".");

    if (numVecs == 1 && this->isConstantStride () && A.isConstantStride ()) {
      const dot_type gblDot = multiVectorSingleColumnDot (const_cast<MV&> (*this), A);
      dots[0] = gblDot;
    }
    else {
      this->dot (A, Kokkos::View<dot_type*, Kokkos::HostSpace>(dots.getRawPtr (), numDots));
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm2 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using ::Tpetra::Details::NORM_TWO;
    using ::Tpetra::Details::ProfilingRegion;
    ProfilingRegion region ("Tpetra::MV::norm2 (host output)");

    // The function needs to be able to sync X.
    MV& X = const_cast<MV&> (*this);
    multiVectorNormImpl (norms.getRawPtr (), X, NORM_TWO);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm2 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const
  {
    Teuchos::ArrayView<mag_type> norms_av (norms.data (), norms.extent (0));
    this->norm2 (norms_av);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm1 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using ::Tpetra::Details::NORM_ONE;
    using ::Tpetra::Details::ProfilingRegion;
    ProfilingRegion region ("Tpetra::MV::norm1 (host output)");

    // The function needs to be able to sync X.
    MV& X = const_cast<MV&> (*this);
    multiVectorNormImpl (norms.data (), X, NORM_ONE);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm1 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const
  {
    Teuchos::ArrayView<mag_type> norms_av (norms.data (), norms.extent (0));
    this->norm1 (norms_av);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normInf (const Teuchos::ArrayView<mag_type>& norms) const
  {
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using ::Tpetra::Details::NORM_INF;
    using ::Tpetra::Details::ProfilingRegion;
    ProfilingRegion region ("Tpetra::MV::normInf (host output)");

    // The function needs to be able to sync X.
    MV& X = const_cast<MV&> (*this);
    multiVectorNormImpl (norms.getRawPtr (), X, NORM_INF);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normInf (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const
  {
    Teuchos::ArrayView<mag_type> norms_av (norms.data (), norms.extent (0));
    this->normInf (norms_av);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  meanValue (const Teuchos::ArrayView<impl_scalar_type>& means) const
  {
    // KR FIXME Overload this method to take a View.
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;

    const size_t lclNumRows = this->getLocalLength ();
    const size_t numVecs = this->getNumVectors ();
    const size_t numMeans = static_cast<size_t> (means.size ());

    TEUCHOS_TEST_FOR_EXCEPTION(
      numMeans != numVecs, std::runtime_error,
      "Tpetra::MultiVector::meanValue: means.size() = " << numMeans
      << " != this->getNumVectors() = " << numVecs << ".");

    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // Make sure that the final output view has the same layout as the
    // temporary view's HostMirror.  Left or Right doesn't matter for
    // a 1-D array anyway; this is just to placate the compiler.
    typedef Kokkos::View<impl_scalar_type*, device_type> local_view_type;
    typedef Kokkos::View<impl_scalar_type*,
      typename local_view_type::HostMirror::array_layout,
      Kokkos::HostSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > host_local_view_type;
    host_local_view_type meansOut (means.getRawPtr (), numMeans);

    RCP<const Comm<int> > comm = this->getMap ().is_null () ? Teuchos::null :
      this->getMap ()->getComm ();

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->need_sync_device ();
    if (useHostVersion) {
      // DualView was last modified on host, so run the local kernel there.
      auto X_lcl = subview (getLocalViewHost (),
                            rowRng, Kokkos::ALL ());
      // Compute the local sum of each column.
      Kokkos::View<impl_scalar_type*, Kokkos::HostSpace> lclSums ("MV::meanValue tmp", numVecs);
      if (isConstantStride ()) {
        KokkosBlas::sum (lclSums, X_lcl);
      }
      else {
        for (size_t j = 0; j < numVecs; ++j) {
          const size_t col = whichVectors_[j];
          KokkosBlas::sum (subview (lclSums, j), subview (X_lcl, ALL (), col));
        }
      }

      // If there are multiple MPI processes, the all-reduce reads
      // from lclSums, and writes to meansOut.  Otherwise, we just
      // copy lclSums into meansOut.
      if (! comm.is_null () && this->isDistributed ()) {
        reduceAll (*comm, REDUCE_SUM, static_cast<int> (numVecs),
                   lclSums.data (), meansOut.data ());
      }
      else {
        Kokkos::deep_copy (meansOut, lclSums);
      }
    }
    else {
      // DualView was last modified on device, so run the local kernel there.
      auto X_lcl = subview (this->getLocalViewDevice (),
                            rowRng, Kokkos::ALL ());

      // Compute the local sum of each column.
      Kokkos::View<impl_scalar_type*, Kokkos::HostSpace> lclSums ("MV::meanValue tmp", numVecs);
      if (isConstantStride ()) {
        KokkosBlas::sum (lclSums, X_lcl);
      }
      else {
        for (size_t j = 0; j < numVecs; ++j) {
          const size_t col = whichVectors_[j];
          KokkosBlas::sum (subview (lclSums, j), subview (X_lcl, ALL (), col));
        }
      }

      // If there are multiple MPI processes, the all-reduce reads
      // from lclSums, and writes to meansOut.  (We assume that MPI
      // can read device memory.)  Otherwise, we just copy lclSums
      // into meansOut.
      if (! comm.is_null () && this->isDistributed ()) {
        reduceAll (*comm, REDUCE_SUM, static_cast<int> (numVecs),
                   lclSums.data (), meansOut.data ());
      }
      else {
        Kokkos::deep_copy (meansOut, lclSums);
      }
    }

    // mfh 12 Apr 2012: Don't take out the cast from the ordinal type
    // to the magnitude type, since operator/ (std::complex<T>, int)
    // isn't necessarily defined.
    const impl_scalar_type OneOverN =
      ATS::one () / static_cast<mag_type> (this->getGlobalLength ());
    for (size_t k = 0; k < numMeans; ++k) {
      meansOut(k) = meansOut(k) * OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  randomize ()
  {
    typedef impl_scalar_type IST;
    typedef Kokkos::Details::ArithTraits<IST> ATS;
    typedef Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> pool_type;
    typedef typename pool_type::generator_type generator_type;

    const IST max = Kokkos::rand<generator_type, IST>::max ();
    const IST min = ATS::is_signed ? IST (-max) : ATS::zero ();

    this->randomize (static_cast<Scalar> (min), static_cast<Scalar> (max));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  randomize (const Scalar& minVal, const Scalar& maxVal)
  {
    typedef impl_scalar_type IST;
    typedef Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> pool_type;

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
    const IST max = static_cast<IST> (maxVal);
    const IST min = static_cast<IST> (minVal);

    // See #1510.  In case diag has already been marked modified on
    // host or device, we need to clear those flags, since the code
    // below works on device.
    this->view_.clear_sync_state();

    this->modify_device ();
    auto thisView = this->getLocalViewDevice ();

    if (isConstantStride ()) {
      Kokkos::fill_random (thisView, rand_pool, min, max);
    }
    else {
      const size_t numVecs = getNumVectors ();
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t col = whichVectors_[k];
        auto X_k = Kokkos::subview (thisView, Kokkos::ALL (), col);
        Kokkos::fill_random (X_k, rand_pool, min, max);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  putScalar (const Scalar& alpha)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::Blas::fill;
    using DES = typename dual_view_type::t_dev::execution_space;
    using HES = typename dual_view_type::t_host::execution_space;
    using LO = LocalOrdinal;
    ProfilingRegion region ("Tpetra::MultiVector::putScalar");

    // We need this cast for cases like Scalar = std::complex<T> but
    // impl_scalar_type = Kokkos::complex<T>.
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const LO lclNumRows = static_cast<LO> (this->getLocalLength ());
    const LO numVecs = static_cast<LO> (this->getNumVectors ());

    // Modify the most recently updated version of the data.  This
    // avoids sync'ing, which could violate users' expectations.
    //
    // If we need sync to device, then host has the most recent version.
    const bool runOnHost = runKernelOnHost(*this);

    this->clear_sync_state();
    if (! runOnHost) {
      this->modify_device ();
      auto X = this->getLocalViewDevice ();
      if (this->isConstantStride ()) {
        fill (DES (), X, theAlpha, lclNumRows, numVecs);
      }
      else {
        fill (DES (), X, theAlpha, lclNumRows, numVecs,
              this->whichVectors_.getRawPtr ());
      }
    }
    else { // last modified in host memory, so modify data there.
      this->modify_host ();
      auto X = this->getLocalViewHost ();
      if (this->isConstantStride ()) {
        fill (HES (), X, theAlpha, lclNumRows, numVecs);
      }
      else {
        fill (HES (), X, theAlpha, lclNumRows, numVecs,
              this->whichVectors_.getRawPtr ());
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceMap (const Teuchos::RCP<const map_type>& newMap)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using ST = Scalar;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;

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
      const size_t origNumRows = view_.extent (0);
      const size_t numCols = this->getNumVectors ();

      if (origNumRows != newNumRows || view_.extent (1) != numCols) {
        view_ = allocDualView<ST, LO, GO, Node> (newNumRows, numCols);
      }
    }
    else if (newMap.is_null ()) { // Case 2: current Map is nonnull, new Map is null
      // I am an excluded process.  Reinitialize my data so that I
      // have 0 rows.  Keep the number of columns as before.
      const size_t newNumRows = static_cast<size_t> (0);
      const size_t numCols = this->getNumVectors ();
      view_ = allocDualView<ST, LO, GO, Node> (newNumRows, numCols);
    }

    this->map_ = newMap;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Scalar& alpha)
  {
    using Kokkos::ALL;
    using IST = impl_scalar_type;

    const IST theAlpha = static_cast<IST> (alpha);
    if (theAlpha == Kokkos::ArithTraits<IST>::one ()) {
      return; // do nothing
    }
    const size_t lclNumRows = getLocalLength ();
    const size_t numVecs = getNumVectors ();
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // We can't substitute putScalar(0.0) for scale(0.0), because the
    // former will overwrite NaNs present in the MultiVector.  The
    // semantics of this call require multiplying them by 0, which
    // IEEE 754 requires to be NaN.

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = need_sync_device ();
    if (useHostVersion) {
      auto Y_lcl = Kokkos::subview (getLocalViewHost (), rowRng, ALL ());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, theAlpha, Y_lcl);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = isConstantStride () ? k : whichVectors_[k];
          auto Y_k = Kokkos::subview (Y_lcl, ALL (), Y_col);
          KokkosBlas::scal (Y_k, theAlpha, Y_k);
        }
      }
    }
    else { // work on device
      auto Y_lcl = Kokkos::subview (getLocalViewDevice (), rowRng, ALL ());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, theAlpha, Y_lcl);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = isConstantStride () ? k : whichVectors_[k];
          auto Y_k = Kokkos::subview (Y_lcl, ALL (), Y_col);
          KokkosBlas::scal (Y_k, theAlpha, Y_k);
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Teuchos::ArrayView<const Scalar>& alphas)
  {
    const size_t numVecs = this->getNumVectors ();
    const size_t numAlphas = static_cast<size_t> (alphas.size ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      numAlphas != numVecs, std::invalid_argument, "Tpetra::MultiVector::"
      "scale: alphas.size() = " << numAlphas << " != this->getNumVectors() = "
      << numVecs << ".");

    // Use a DualView to copy the scaling constants onto the device.
    using k_alphas_type = Kokkos::DualView<impl_scalar_type*, device_type>;
    k_alphas_type k_alphas ("alphas::tmp", numAlphas);
    k_alphas.modify_host ();
    for (size_t i = 0; i < numAlphas; ++i) {
      k_alphas.h_view(i) = static_cast<impl_scalar_type> (alphas[i]);
    }
    k_alphas.sync_device ();
    // Invoke the scale() overload that takes a device View of coefficients.
    this->scale (k_alphas.view_device ());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Kokkos::View<const impl_scalar_type*, device_type>& alphas)
  {
    using Kokkos::ALL;
    using Kokkos::subview;

    const size_t lclNumRows = this->getLocalLength ();
    const size_t numVecs = this->getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (alphas.extent (0)) != numVecs,
      std::invalid_argument, "Tpetra::MultiVector::scale(alphas): "
      "alphas.extent(0) = " << alphas.extent (0)
      << " != this->getNumVectors () = " << numVecs << ".");
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // NOTE (mfh 08 Apr 2015) We prefer to let the compiler deduce the
    // type of the return value of subview.  This is because if we
    // switch the array layout from LayoutLeft to LayoutRight
    // (preferred for performance of block operations), the types
    // below won't be valid.  (A view of a column of a LayoutRight
    // multivector has LayoutStride, not LayoutLeft.)

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->need_sync_device ();
    if (useHostVersion) {
      // Work in host memory.  This means we need to create a host
      // mirror of the input View of coefficients.
      auto alphas_h = Kokkos::create_mirror_view (alphas);
      Kokkos::deep_copy (alphas_h, alphas);

      auto Y_lcl = subview (this->getLocalViewHost (), rowRng, ALL ());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, alphas_h, Y_lcl);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = this->isConstantStride () ? k :
            this->whichVectors_[k];
          auto Y_k = subview (Y_lcl, ALL (), Y_col);
          // We don't have to use the entire 1-D View here; we can use
          // the version that takes a scalar coefficient.
          KokkosBlas::scal (Y_k, alphas_h(k), Y_k);
        }
      }
    }
    else { // Work in device memory, using the input View 'alphas' directly.
      auto Y_lcl = subview (this->getLocalViewDevice (), rowRng, ALL ());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, alphas, Y_lcl);
      }
      else {
        // FIXME (mfh 15 Mar 2019) We need one coefficient at a time,
        // as values on host, so copy them to host.  Another approach
        // would be to fix scal() so that it takes a 0-D View as the
        // second argument.
        auto alphas_h = Kokkos::create_mirror_view (alphas);
        Kokkos::deep_copy (alphas_h, alphas);

        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = this->isConstantStride () ? k :
            this->whichVectors_[k];
          auto Y_k = subview (Y_lcl, ALL (), Y_col);
          KokkosBlas::scal (Y_k, alphas_h(k), Y_k);
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Scalar& alpha,
         const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "scale: ";

    const size_t lclNumRows = getLocalLength ();
    const size_t numVecs = getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != A.getLocalLength (), std::invalid_argument,
      "this->getLocalLength() = " << lclNumRows << " != A.getLocalLength() = "
      << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs != A.getNumVectors (), std::invalid_argument,
      "this->getNumVectors() = " << numVecs << " != A.getNumVectors() = "
      << A.getNumVectors () << ".");

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) {
      this->sync_device ();
    }
    if (A.need_sync_device ()) {
      const_cast<MV&>(A).sync_device ();
    }

    this->modify_device ();
    auto Y_lcl_orig = this->getLocalViewDevice ();
    auto X_lcl_orig = A.getLocalViewDevice ();
    auto Y_lcl = subview (Y_lcl_orig, rowRng, ALL ());
    auto X_lcl = subview (X_lcl_orig, rowRng, ALL ());

    if (isConstantStride () && A.isConstantStride ()) {
      KokkosBlas::scal (Y_lcl, theAlpha, X_lcl);
    }
    else {
      // Make sure that Kokkos only uses the local length for add.
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t Y_col = this->isConstantStride () ? k : this->whichVectors_[k];
        const size_t X_col = A.isConstantStride () ? k : A.whichVectors_[k];
        auto Y_k = subview (Y_lcl, ALL (), Y_col);
          auto X_k = subview (X_lcl, ALL (), X_col);

          KokkosBlas::scal (Y_k, theAlpha, X_k);
      }
    }
  }



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  reciprocal (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A)
  {
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char tfecfFuncName[] = "reciprocal: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
       getLocalLength () != A.getLocalLength (), std::runtime_error,
       "MultiVectors do not have the same local length.  "
       "this->getLocalLength() = " << getLocalLength ()
       << " != A.getLocalLength() = " << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != this->getNumVectors (), std::runtime_error,
      ": MultiVectors do not have the same number of columns (vectors).  "
       "this->getNumVectors() = " << getNumVectors ()
       << " != A.getNumVectors() = " << A.getNumVectors () << ".");

    const size_t numVecs = getNumVectors ();

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) {
      this->sync_device ();
    }
    if (A.need_sync_device ()) {
      const_cast<MV&>(A).sync_device ();
    }
    this->modify_device ();

    auto this_view_dev = this->getLocalViewDevice ();
    auto A_view_dev = A.getLocalViewDevice ();

    if (isConstantStride () && A.isConstantStride ()) {
      KokkosBlas::reciprocal (this_view_dev, A_view_dev);
    }
    else {
      using Kokkos::ALL;
      using Kokkos::subview;
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        auto vector_k = subview (this_view_dev, ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        auto vector_Ak = subview (A_view_dev, ALL (), A_col);
        KokkosBlas::reciprocal (vector_k, vector_Ak);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  abs (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A)
  {
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
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
    const size_t numVecs = getNumVectors ();

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) {
      this->sync_device ();
    }
    if (A.need_sync_device ()) {
      const_cast<MV&>(A).sync_device ();
    }
    this->modify_device ();

    auto this_view_dev = this->getLocalViewDevice ();
    auto A_view_dev = A.getLocalViewDevice ();

    if (isConstantStride () && A.isConstantStride ()) {
      KokkosBlas::abs (this_view_dev, A_view_dev);
    }
    else {
      using Kokkos::ALL;
      using Kokkos::subview;

      for (size_t k=0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        auto vector_k = subview (this_view_dev, ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        auto vector_Ak = subview (A_view_dev, ALL (), A_col);
        KokkosBlas::abs (vector_k, vector_Ak);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  update (const Scalar& alpha,
          const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
          const Scalar& beta)
  {
    const char tfecfFuncName[] = "update: ";
    using Kokkos::subview;
    using Kokkos::ALL;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::update(alpha,A,beta)");

    const size_t lclNumRows = getLocalLength ();
    const size_t numVecs = getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != A.getLocalLength (), std::invalid_argument,
      "this->getLocalLength() = " << lclNumRows << " != A.getLocalLength() = "
      << A.getLocalLength () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs != A.getNumVectors (), std::invalid_argument,
      "this->getNumVectors() = " << numVecs << " != A.getNumVectors() = "
      << A.getNumVectors () << ".");

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) {
      this->sync_device ();
    }
    if (A.need_sync_device ()) {
      const_cast<MV&>(A).sync_device ();
    }

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const impl_scalar_type theBeta = static_cast<impl_scalar_type> (beta);
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    auto Y_lcl_orig = this->getLocalViewDevice ();
    auto Y_lcl = subview (Y_lcl_orig, rowRng, Kokkos::ALL ());
    auto X_lcl_orig = A.getLocalViewDevice ();
    auto X_lcl = subview (X_lcl_orig, rowRng, Kokkos::ALL ());

    // The device memory of *this is about to be modified
    this->modify_device ();
    if (isConstantStride () && A.isConstantStride ()) {
      KokkosBlas::axpby (theAlpha, X_lcl, theBeta, Y_lcl);
    }
    else {
      // Make sure that Kokkos only uses the local length for add.
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t Y_col = this->isConstantStride () ? k : this->whichVectors_[k];
        const size_t X_col = A.isConstantStride () ? k : A.whichVectors_[k];
        auto Y_k = subview (Y_lcl, ALL (), Y_col);
        auto X_k = subview (X_lcl, ALL (), X_col);

        KokkosBlas::axpby (theAlpha, X_k, theBeta, Y_k);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  update (const Scalar& alpha,
          const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
          const Scalar& beta,
          const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
          const Scalar& gamma)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    const char tfecfFuncName[] = "update(alpha,A,beta,B,gamma): ";

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::update(alpha,A,beta,B,gamma)");

    const size_t lclNumRows = this->getLocalLength ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != A.getLocalLength (), std::invalid_argument,
      "The input MultiVector A has " << A.getLocalLength () << " local "
      "row(s), but this MultiVector has " << lclNumRows << " local "
      "row(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      lclNumRows != B.getLocalLength (), std::invalid_argument,
      "The input MultiVector B has " << B.getLocalLength () << " local "
      "row(s), but this MultiVector has " << lclNumRows << " local "
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

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const impl_scalar_type theBeta = static_cast<impl_scalar_type> (beta);
    const impl_scalar_type theGamma = static_cast<impl_scalar_type> (gamma);

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) this->sync_device ();
    if (A.need_sync_device ())     const_cast<MV&>(A).sync_device ();
    if (B.need_sync_device ())     const_cast<MV&>(B).sync_device ();

    // This method modifies *this.
    this->modify_device ();

    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // Prefer 'auto' over specifying the type explicitly.  This avoids
    // issues with a subview possibly having a different type than the
    // original view.
    auto C_lcl = subview (this->getLocalViewDevice (), rowRng, ALL ());
    auto A_lcl = subview (A.getLocalViewDevice (), rowRng, ALL ());
    auto B_lcl = subview (B.getLocalViewDevice (), rowRng, ALL ());

    if (isConstantStride () && A.isConstantStride () && B.isConstantStride ()) {
      KokkosBlas::update (theAlpha, A_lcl, theBeta, B_lcl, theGamma, C_lcl);
    }
    else {
      // Some input (or *this) is not constant stride,
      // so perform the update one column at a time.
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        const size_t A_col = A.isConstantStride () ? k : A.whichVectors_[k];
        const size_t B_col = B.isConstantStride () ? k : B.whichVectors_[k];
        KokkosBlas::update (theAlpha, subview (A_lcl, rowRng, A_col),
                            theBeta, subview (B_lcl, rowRng, B_col),
                            theGamma, subview (C_lcl, rowRng, this_col));
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getData (size_t j) const
  {
    using Kokkos::ALL;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using IST = impl_scalar_type;
    const char tfecfFuncName[] = "getData: ";

    // Any MultiVector method that called the (classic) Kokkos Node's
    // viewBuffer or viewBufferNonConst methods always implied a
    // device->host synchronization.  Thus, we synchronize here as
    // well.
    const_cast<MV&> (*this).sync_host ();

    auto hostView = getLocalViewHost ();
    const size_t col = isConstantStride () ? j : whichVectors_[j];
    auto hostView_j = Kokkos::subview (hostView, ALL (), col);
    Teuchos::ArrayRCP<const IST> dataAsArcp =
      Kokkos::Compat::persistingView (hostView_j, 0, getLocalLength ());

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (hostView_j.extent (0)) <
       static_cast<size_t> (dataAsArcp.size ()), std::logic_error,
       "hostView_j.extent(0) = " << hostView_j.extent (0)
       << " < dataAsArcp.size() = " << dataAsArcp.size () << ".  "
       "Please report this bug to the Tpetra developers.");

    return Teuchos::arcp_reinterpret_cast<const Scalar> (dataAsArcp);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getDataNonConst (size_t j)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using IST = impl_scalar_type;
    const char tfecfFuncName[] = "getDataNonConst: ";

    // Any MultiVector method that called the (classic) Kokkos Node's
    // viewBuffer or viewBufferNonConst methods always implied a
    // device->host synchronization.  Thus, we synchronize here as
    // well.
    const_cast<MV*> (this)->sync_host ();
    // Calling getDataNonConst() implies that the user plans to modify
    // the values in the MultiVector, so we mark the host data as modified.
    modify_host ();

    auto hostView = getLocalViewHost ();
    const size_t col = isConstantStride () ? j : whichVectors_[j];
    auto hostView_j = subview (hostView, ALL (), col);
    Teuchos::ArrayRCP<IST> dataAsArcp =
      Kokkos::Compat::persistingView (hostView_j, 0, getLocalLength ());

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (hostView_j.extent (0)) <
       static_cast<size_t> (dataAsArcp.size ()), std::logic_error,
       "hostView_j.extent(0) = " << hostView_j.extent (0)
       << " < dataAsArcp.size() = " << dataAsArcp.size () << ".  "
       "Please report this bug to the Tpetra developers.");

    return Teuchos::arcp_reinterpret_cast<Scalar> (dataAsArcp);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  subCopy (const Teuchos::ArrayView<const size_t>& cols) const
  {
    using Teuchos::RCP;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;

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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  subCopy (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::RCP;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;

    RCP<const MV> X_sub = this->subView (colRng);
    RCP<MV> Y (new MV (this->getMap (), static_cast<size_t> (colRng.size ()), false));
    Y->assign (*X_sub);
    return Y;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getOrigNumLocalRows () const {
    return origView_.extent (0);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getOrigNumLocalCols () const {
    return origView_.extent (1);
  }

  template <class Scalar, class LO, class GO, class Node>
  MultiVector<Scalar, LO, GO, Node>::
  MultiVector (const MultiVector<Scalar, LO, GO, Node>& X,
               const Teuchos::RCP<const map_type>& subMap,
               const local_ordinal_type rowOffset) :
    base_type (subMap)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_MIN;
    using std::endl;
    using MV = MultiVector<Scalar, LO, GO, Node>;
    const char prefix[] = "Tpetra::MultiVector constructor (offsetView): ";
    const char suffix[] = "Please report this bug to the Tpetra developers.";
    int lclGood = 1;
    int gblGood = 1;
    std::unique_ptr<std::ostringstream> errStrm;
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();

    // Be careful to use the input Map's communicator, not X's.  The
    // idea is that, on return, *this is a subview of X, using the
    // input Map.
    const auto comm = subMap->getComm ();
    TEUCHOS_ASSERT( ! comm.is_null () );
    const int myRank = comm->getRank ();

    const LO lclNumRowsBefore = static_cast<LO> (X.getLocalLength ());
    const LO numCols = static_cast<LO> (X.getNumVectors ());
    const LO newNumRows = static_cast<LO> (subMap->getNodeNumElements ());
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": " << prefix
         << "X: {lclNumRows: " << lclNumRowsBefore
         << ", origLclNumRows: " << X.getOrigNumLocalRows ()
         << ", numCols: " << numCols << "}, "
         << "subMap: {lclNumRows: " << newNumRows << "}" << endl;
      std::cerr << os.str ();
    }
    // We ask for the _original_ number of rows in X, because X could
    // be a shorter (fewer rows) view of a longer MV.  (For example, X
    // could be a domain Map view of a column Map MV.)
    const bool tooManyElts =
      newNumRows + rowOffset > static_cast<LO> (X.getOrigNumLocalRows ());
    if (tooManyElts) {
      errStrm = std::unique_ptr<std::ostringstream> (new std::ostringstream);
      *errStrm << "  Proc " << myRank << ": subMap->getNodeNumElements() (="
               << newNumRows << ") + rowOffset (=" << rowOffset
               << ") > X.getOrigNumLocalRows() (=" << X.getOrigNumLocalRows ()
               << ")." << endl;
      lclGood = 0;
      TEUCHOS_TEST_FOR_EXCEPTION
        (! debug && tooManyElts, std::invalid_argument,
         prefix << errStrm->str () << suffix);
    }

    if (debug) {
      reduceAll<int, int> (*comm, REDUCE_MIN, lclGood, outArg (gblGood));
      if (gblGood != 1) {
        std::ostringstream gblErrStrm;
        const std::string myErrStr =
          errStrm.get () != nullptr ? errStrm->str () : std::string ("");
        ::Tpetra::Details::gathervPrint (gblErrStrm, myErrStr, *comm);
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, gblErrStrm.str ());
      }
    }

    using range_type = std::pair<LO, LO>;
    const range_type origRowRng
      (rowOffset, static_cast<LO> (X.origView_.extent (0)));
    const range_type rowRng
      (rowOffset, rowOffset + newNumRows);

    dual_view_type newOrigView = subview (X.origView_, origRowRng, ALL ());
    // FIXME (mfh 29 Sep 2016) If we just use X.view_ here, it breaks
    // CrsMatrix's Gauss-Seidel implementation (which assumes the
    // ability to create domain Map views of column Map MultiVectors,
    // and then get the original column Map MultiVector out again).
    // If we just use X.origView_ here, it breaks the fix for #46.
    // The test for rowOffset == 0 is a hack that makes both tests
    // pass, but doesn't actually fix the more general issue.  In
    // particular, the right way to fix Gauss-Seidel would be to fix
    // #385; that would make "getting the original column Map
    // MultiVector out again" unnecessary.
    dual_view_type newView =
      subview (rowOffset == 0 ? X.origView_ : X.view_, rowRng, ALL ());

    // NOTE (mfh 06 Jan 2015) Work-around to deal with Kokkos not
    // handling subviews of degenerate Views quite so well.  For some
    // reason, the ([0,0], [0,2]) subview of a 0 x 2 DualView is 0 x
    // 0.  We work around by creating a new empty DualView of the
    // desired (degenerate) dimensions.
    if (newOrigView.extent (0) == 0 &&
        newOrigView.extent (1) != X.origView_.extent (1)) {
      newOrigView =
        allocDualView<Scalar, LO, GO, Node> (0, X.getNumVectors ());
    }
    if (newView.extent (0) == 0 &&
        newView.extent (1) != X.view_.extent (1)) {
      newView =
        allocDualView<Scalar, LO, GO, Node> (0, X.getNumVectors ());
    }

    MV subViewMV = X.isConstantStride () ?
      MV (subMap, newView, newOrigView) :
      MV (subMap, newView, newOrigView, X.whichVectors_ ());

    if (debug) {
      const LO lclNumRowsRet = static_cast<LO> (subViewMV.getLocalLength ());
      const LO numColsRet = static_cast<LO> (subViewMV.getNumVectors ());
      if (newNumRows != lclNumRowsRet || numCols != numColsRet) {
        lclGood = 0;
        if (errStrm.get () == nullptr) {
          errStrm = std::unique_ptr<std::ostringstream> (new std::ostringstream);
        }
        *errStrm << "  Proc " << myRank <<
          ": subMap.getNodeNumElements(): " << newNumRows <<
          ", subViewMV.getLocalLength(): " << lclNumRowsRet <<
          ", X.getNumVectors(): " << numCols <<
          ", subViewMV.getNumVectors(): " << numColsRet << endl;
      }
      reduceAll<int, int> (*comm, REDUCE_MIN, lclGood, outArg (gblGood));
      if (gblGood != 1) {
        std::ostringstream gblErrStrm;
        if (myRank == 0) {
          gblErrStrm << prefix << "Returned MultiVector has the wrong local "
            "dimensions on one or more processes:" << endl;
        }
        const std::string myErrStr =
          errStrm.get () != nullptr ? errStrm->str () : std::string ("");
        ::Tpetra::Details::gathervPrint (gblErrStrm, myErrStr, *comm);
        gblErrStrm << suffix << endl;
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, gblErrStrm.str ());
      }
    }

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": " << prefix << "Call op=" << endl;
      std::cerr << os.str ();
    }

    *this = subViewMV; // shallow copy

    if (verbose) {
      std::ostringstream os;
      os << "Proc " << myRank << ": " << prefix << "Done!" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Scalar, class LO, class GO, class Node>
  MultiVector<Scalar, LO, GO, Node>::
  MultiVector (const MultiVector<Scalar, LO, GO, Node>& X,
               const map_type& subMap,
               const size_t rowOffset) :
    MultiVector (X, Teuchos::RCP<const map_type> (new map_type (subMap)),
                 static_cast<local_ordinal_type> (rowOffset))
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  offsetView (const Teuchos::RCP<const map_type>& subMap,
              const size_t offset) const
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    return Teuchos::rcp (new MV (*this, *subMap, offset));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                      const size_t offset)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    return Teuchos::rcp (new MV (*this, *subMap, offset));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  subView (const Teuchos::ArrayView<const size_t>& cols) const
  {
    using Teuchos::Array;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;

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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  subView (const Teuchos::Range1D& colRng) const
  {
    using ::Tpetra::Details::Behavior;
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "subView(Range1D): ";

    const size_t lclNumRows = this->getLocalLength ();
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

    RCP<const MV> X_ret; // the MultiVector subview to return

    // FIXME (mfh 14 Apr 2015) Apparently subview on DualView is still
    // broken for the case of views with zero rows.  I will brutally
    // enforce that the subview has the correct dimensions.  In
    // particular, in the case of zero rows, I will, if necessary,
    // create a new dual_view_type with zero rows and the correct
    // number of columns.  In a debug build, I will use an all-reduce
    // to ensure that it has the correct dimensions on all processes.

    const std::pair<size_t, size_t> rows (0, lclNumRows);
    if (colRng.size () == 0) {
      const std::pair<size_t, size_t> cols (0, 0); // empty range
      dual_view_type X_sub = takeSubview (this->view_, ALL (), cols);
      X_ret = rcp (new MV (this->getMap (), X_sub, origView_));
    }
    else {
      // Returned MultiVector is constant stride only if *this is.
      if (isConstantStride ()) {
        const std::pair<size_t, size_t> cols (colRng.lbound (),
                                              colRng.ubound () + 1);
        dual_view_type X_sub = takeSubview (this->view_, ALL (), cols);
        X_ret = rcp (new MV (this->getMap (), X_sub, origView_));
      }
      else {
        if (static_cast<size_t> (colRng.size ()) == static_cast<size_t> (1)) {
          // We're only asking for one column, so the result does have
          // constant stride, even though this MultiVector does not.
          const std::pair<size_t, size_t> col (whichVectors_[0] + colRng.lbound (),
                                               whichVectors_[0] + colRng.ubound () + 1);
          dual_view_type X_sub = takeSubview (view_, ALL (), col);
          X_ret = rcp (new MV (this->getMap (), X_sub, origView_));
        }
        else {
          Array<size_t> which (whichVectors_.begin () + colRng.lbound (),
                               whichVectors_.begin () + colRng.ubound () + 1);
          X_ret = rcp (new MV (this->getMap (), view_, origView_, which));
        }
      }
    }

    const bool debug = Behavior::debug ();
    if (debug) {
      using Teuchos::Comm;
      using Teuchos::outArg;
      using Teuchos::REDUCE_MIN;
      using Teuchos::reduceAll;

      RCP<const Comm<int> > comm = this->getMap ().is_null () ?
        Teuchos::null : this->getMap ()->getComm ();
      if (! comm.is_null ()) {
        int lclSuccess = 1;
        int gblSuccess = 1;

        if (X_ret.is_null ()) {
          lclSuccess = 0;
        }
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (lclSuccess != 1, std::logic_error, "X_ret (the subview of this "
           "MultiVector; the return value of this method) is null on some MPI "
           "process in this MultiVector's communicator.  This should never "
           "happen.  Please report this bug to the Tpetra developers.");
        if (! X_ret.is_null () &&
            X_ret->getNumVectors () != static_cast<size_t> (colRng.size ())) {
          lclSuccess = 0;
        }
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                             outArg (gblSuccess));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (lclSuccess != 1, std::logic_error, "X_ret->getNumVectors() != "
           "colRng.size(), on at least one MPI process in this MultiVector's "
           "communicator.  This should never happen.  "
           "Please report this bug to the Tpetra developers.");
      }
    }
    return X_ret;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  subViewNonConst (const Teuchos::ArrayView<const size_t> &cols)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    return Teuchos::rcp_const_cast<MV> (this->subView (cols));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  subViewNonConst (const Teuchos::Range1D &colRng)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    return Teuchos::rcp_const_cast<MV> (this->subView (colRng));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
               const size_t j)
    : base_type (X.getMap ())
  {
    using Kokkos::subview;
    typedef std::pair<size_t, size_t> range_type;
    const char tfecfFuncName[] = "MultiVector(const MultiVector&, const size_t): ";

    const size_t numCols = X.getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (j >= numCols, std::invalid_argument, "Input index j (== " << j
       << ") exceeds valid column index range [0, " << numCols << " - 1].");
    const size_t jj = X.isConstantStride () ?
      static_cast<size_t> (j) :
      static_cast<size_t> (X.whichVectors_[j]);
    this->view_ = takeSubview (X.view_, Kokkos::ALL (), range_type (jj, jj+1));
    this->origView_ = X.origView_;

    // mfh 31 Jul 2017: It would be unwise to execute concurrent
    // Export or Import operations with different subviews of a
    // MultiVector.  Thus, it is safe to reuse communication buffers.
    // See #1560 discussion.
    //
    // We only need one column's worth of buffer for imports_ and
    // exports_.  Taking subviews now ensures that their lengths will
    // be exactly what we need, so we won't have to resize them later.
    {
      const size_t newSize = X.imports_.extent (0) / numCols;
      const size_t offset = jj*newSize;
      auto newImports = X.imports_;
      newImports.d_view = subview (X.imports_.d_view,
                                   range_type (offset, offset+newSize));
      newImports.h_view = subview (X.imports_.h_view,
                                   range_type (offset, offset+newSize));
      this->imports_ = newImports;
    }
    {
      const size_t newSize = X.exports_.extent (0) / numCols;
      const size_t offset = jj*newSize;
      auto newExports = X.exports_;
      newExports.d_view = subview (X.exports_.d_view,
                                   range_type (offset, offset+newSize));
      newExports.h_view = subview (X.exports_.h_view,
                                   range_type (offset, offset+newSize));
      this->exports_ = newExports;
    }
    // These two DualViews already either have the right number of
    // entries, or zero entries.  This means that we don't need to
    // resize them.
    this->numImportPacketsPerLID_ = X.numImportPacketsPerLID_;
    this->numExportPacketsPerLID_ = X.numExportPacketsPerLID_;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getVector (const size_t j) const
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;
    return Teuchos::rcp (new V (*this, j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getVectorNonConst (const size_t j)
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;
    return Teuchos::rcp (new V (*this, j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get1dCopy (const Teuchos::ArrayView<Scalar>& A, const size_t LDA) const
  {
    using dev_view_type = typename dual_view_type::t_dev;
    using host_view_type = typename dual_view_type::t_host;
    using IST = impl_scalar_type;
    using input_view_type = Kokkos::View<IST**, Kokkos::LayoutLeft,
                                         Kokkos::HostSpace,
                                         Kokkos::MemoryUnmanaged>;
    const char tfecfFuncName[] = "get1dCopy: ";

    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (LDA < numRows, std::runtime_error,
       "LDA = " << LDA << " < numRows = " << numRows << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numRows > size_t (0) && numCols > size_t (0) &&
       size_t (A.size ()) < LDA * (numCols - 1) + numRows,
       std::runtime_error,
       "A.size() = " << A.size () << ", but its size must be at least "
       << (LDA * (numCols - 1) + numRows) << " to hold all the entries.");

    const std::pair<size_t, size_t> rowRange (0, numRows);
    const std::pair<size_t, size_t> colRange (0, numCols);

    input_view_type A_view_orig (reinterpret_cast<IST*> (A.getRawPtr ()),
                                 LDA, numCols);
    auto A_view = Kokkos::subview (A_view_orig, rowRange, colRange);

    // Use the most recently updated version of this MultiVector's
    // data.  This avoids sync'ing, which could violate users'
    // expectations.
    //
    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->need_sync_device ();

    dev_view_type srcView_dev;
    host_view_type srcView_host;
    if (useHostVersion) {
      srcView_host = this->getLocalViewHost ();
    }
    else {
      srcView_dev = this->getLocalViewDevice ();
    }

    if (this->isConstantStride ()) {
      if (useHostVersion) {
        Kokkos::deep_copy (A_view, srcView_host);
      }
      else {
        Kokkos::deep_copy (A_view, srcView_dev);
      }
    }
    else {
      for (size_t j = 0; j < numCols; ++j) {
        const size_t srcCol = this->whichVectors_[j];
        auto dstColView = Kokkos::subview (A_view, rowRange, j);

        if (useHostVersion) {
          auto srcColView_host = Kokkos::subview (srcView_host, rowRange, srcCol);
          Kokkos::deep_copy (dstColView, srcColView_host);
        }
        else {
          auto srcColView_dev = Kokkos::subview (srcView_dev, rowRange, srcCol);
          Kokkos::deep_copy (dstColView, srcColView_dev);
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get2dCopy (const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> >& ArrayOfPtrs) const
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;
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
        Teuchos::RCP<const V> X_j = this->getVector (j);
        const size_t LDA = static_cast<size_t> (ArrayOfPtrs[j].size ());
        X_j->get1dCopy (ArrayOfPtrs[j], LDA);
      }
    }
  }

  namespace { // (anonymous)
    template <class SC, class LO, class GO, class NT>
    typename MultiVector<SC, LO, GO, NT>::dual_view_type::t_host
    syncMVToHostIfNeededAndGetHostView (MultiVector<SC, LO, GO, NT>& X,
                                        const bool markModified)
    {
      // NOTE (mfh 16 May 2016) get?dView() and get?dViewNonConst()
      // (replace ? with 1 or 2) have always been device->host
      // synchronization points, since <= 2012.  We retain this
      // behavior for backwards compatibility.
      if (X.need_sync_host ()) {
        X.sync_host ();
      }
      if (markModified) {
        X.modify_host ();
      }
      return X.getLocalViewHost ();
    }
  } // namespace (anonymous)


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
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
      // Since get1dView() is and was always marked const, I have to
      // cast away const here in order not to break backwards
      // compatibility.
      constexpr bool markModified = false;
      using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      auto X_lcl = syncMVToHostIfNeededAndGetHostView (const_cast<MV&> (*this),
                                                       markModified);
      Teuchos::ArrayRCP<const impl_scalar_type> dataAsArcp =
        Kokkos::Compat::persistingView (X_lcl);
      return Teuchos::arcp_reinterpret_cast<const Scalar> (dataAsArcp);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get1dViewNonConst ()
  {
    if (this->getLocalLength () == 0 || this->getNumVectors () == 0) {
      return Teuchos::null;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
         "get1dViewNonConst: This MultiVector does not have constant stride, "
         "so it is not possible to view its data as a single array.  You may "
         "check whether a MultiVector has constant stride by calling "
         "isConstantStride().");
      constexpr bool markModified = true;
      auto X_lcl = syncMVToHostIfNeededAndGetHostView (*this, markModified);
      Teuchos::ArrayRCP<impl_scalar_type> dataAsArcp =
        Kokkos::Compat::persistingView (X_lcl);
      return Teuchos::arcp_reinterpret_cast<Scalar> (dataAsArcp);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get2dViewNonConst ()
  {
    constexpr bool markModified = true;
    auto X_lcl = syncMVToHostIfNeededAndGetHostView (*this, markModified);

    // Don't use the row range here on the outside, in order to avoid
    // a strided return type (in case Kokkos::subview is conservative
    // about that).  Instead, use the row range for the column views
    // in the loop.
    const size_t myNumRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const Kokkos::pair<size_t, size_t> rowRange (0, myNumRows);

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > views (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      const size_t col = this->isConstantStride () ? j : this->whichVectors_[j];
      auto X_lcl_j = Kokkos::subview (X_lcl, rowRange, col);
      Teuchos::ArrayRCP<impl_scalar_type> X_lcl_j_arcp =
        Kokkos::Compat::persistingView (X_lcl_j);
      views[j] = Teuchos::arcp_reinterpret_cast<Scalar> (X_lcl_j_arcp);
    }
    return views;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get2dView () const
  {
    // Since get2dView() is and was always marked const, I have to
    // cast away const here in order not to break backwards
    // compatibility.
    constexpr bool markModified = false;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    auto X_lcl = syncMVToHostIfNeededAndGetHostView (const_cast<MV&> (*this),
                                                     markModified);
    // Don't use the row range here on the outside, in order to avoid
    // a strided return type (in case Kokkos::subview is conservative
    // about that).  Instead, use the row range for the column views
    // in the loop.
    const size_t myNumRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const Kokkos::pair<size_t, size_t> rowRange (0, myNumRows);

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > views (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      const size_t col = this->isConstantStride () ? j : this->whichVectors_[j];
      auto X_lcl_j = Kokkos::subview (X_lcl, rowRange, col);
      Teuchos::ArrayRCP<const impl_scalar_type> X_lcl_j_arcp =
        Kokkos::Compat::persistingView (X_lcl_j);
      views[j] = Teuchos::arcp_reinterpret_cast<const Scalar> (X_lcl_j_arcp);
    }
    return views;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  multiply (Teuchos::ETransp transA,
            Teuchos::ETransp transB,
            const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
            const Scalar& beta)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::CONJ_TRANS;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using std::endl;
    using ATS = Kokkos::ArithTraits<impl_scalar_type>;
    using LO = local_ordinal_type;
    using STS = Teuchos::ScalarTraits<Scalar>;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char tfecfFuncName[] = "multiply: ";
    ProfilingRegion region ("Tpetra::MV::multiply");

    // This routine performs a variety of matrix-matrix multiply
    // operations, interpreting the MultiVector (this-aka C , A and B)
    // as 2D matrices.  Variations are due to the fact that A, B and C
    // can be locally replicated or globally distributed MultiVectors
    // and that we may or may not operate with the transpose of A and
    // B.  Possible cases are:
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

    impl_scalar_type beta_local = beta; // local copy of beta; might be reassigned below

    const bool A_is_local = ! A.isDistributed ();
    const bool B_is_local = ! B.isDistributed ();
    const bool C_is_local = ! this->isDistributed ();

    // In debug mode, check compatibility of local dimensions.  We
    // only do this in debug mode, since it requires an all-reduce
    // to ensure correctness on all processses.  It's entirely
    // possible that only some processes may have incompatible local
    // dimensions.  Throwing an exception only on those processes
    // could cause this method to hang.
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if (debug) {
      auto myMap = this->getMap ();
      if (! myMap.is_null () && ! myMap->getComm ().is_null ()) {
        using Teuchos::REDUCE_MIN;
        using Teuchos::reduceAll;
        using Teuchos::outArg;

        auto comm = myMap->getComm ();
        const size_t A_nrows =
          (transA != NO_TRANS) ? A.getNumVectors () : A.getLocalLength ();
        const size_t A_ncols =
          (transA != NO_TRANS) ? A.getLocalLength () : A.getNumVectors ();
        const size_t B_nrows =
          (transB != NO_TRANS) ? B.getNumVectors () : B.getLocalLength ();
        const size_t B_ncols =
          (transB != NO_TRANS) ? B.getLocalLength () : B.getNumVectors ();

        const bool lclBad = this->getLocalLength () != A_nrows ||
          this->getNumVectors () != B_ncols || A_ncols != B_nrows;

        const int myRank = comm->getRank ();
        std::ostringstream errStrm;
        if (this->getLocalLength () != A_nrows) {
          errStrm << "Proc " << myRank << ": this->getLocalLength()="
            << this->getLocalLength () << " != A_nrows=" << A_nrows
            << "." << std::endl;
        }
        if (this->getNumVectors () != B_ncols) {
          errStrm << "Proc " << myRank << ": this->getNumVectors()="
            << this->getNumVectors () << " != B_ncols=" << B_ncols
            << "." << std::endl;
        }
        if (A_ncols != B_nrows) {
          errStrm << "Proc " << myRank << ": A_ncols="
            << A_ncols << " != B_nrows=" << B_nrows
            << "." << std::endl;
        }

        const int lclGood = lclBad ? 0 : 1;
        int gblGood = 0;
        reduceAll<int, int> (*comm, REDUCE_MIN, lclGood, outArg (gblGood));
        if (gblGood != 1) {
          std::ostringstream os;
          ::Tpetra::Details::gathervPrint (os, errStrm.str (), *comm);

          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (true, std::runtime_error, "Inconsistent local dimensions on at "
             "least one process in this object's communicator." << std::endl
             << "Operation: "
             << "C(" << (C_is_local ? "local" : "distr") << ") = "
             << alpha << "*A"
             << (transA == Teuchos::TRANS ? "^T" :
                 (transA == Teuchos::CONJ_TRANS ? "^H" : ""))
             << "(" << (A_is_local ? "local" : "distr") << ") * "
             << beta << "*B"
             << (transA == Teuchos::TRANS ? "^T" :
                 (transA == Teuchos::CONJ_TRANS ? "^H" : ""))
             << "(" << (B_is_local ? "local" : "distr") << ")." << std::endl
             << "Global dimensions: C(" << this->getGlobalLength () << ", "
             << this->getNumVectors () << "), A(" << A.getGlobalLength ()
             << ", " << A.getNumVectors () << "), B(" << B.getGlobalLength ()
             << ", " << B.getNumVectors () << ")." << std::endl
             << os.str ());
        }
      }
    }

    // Case 1: C(local) = A^X(local) * B^X(local)
    const bool Case1 = C_is_local && A_is_local && B_is_local;
    // Case 2: C(local) = A^T(distr) * B  (distr)
    const bool Case2 = C_is_local && ! A_is_local && ! B_is_local &&
      transA != NO_TRANS &&
      transB == NO_TRANS;
    // Case 3: C(distr) = A  (distr) * B^X(local)
    const bool Case3 = ! C_is_local && ! A_is_local && B_is_local &&
      transA == NO_TRANS;

    // Test that we are considering a meaningful case
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! Case1 && ! Case2 && ! Case3, std::runtime_error,
       "Multiplication of op(A) and op(B) into *this is not a "
       "supported use case.");

    if (beta != STS::zero () && Case2) {
      // If Case2, then C is local and contributions must be summed
      // across all processes.  However, if beta != 0, then accumulate
      // beta*C into the sum.  When summing across all processes, we
      // only want to accumulate this once, so set beta == 0 on all
      // processes except Process 0.
      const int myRank = this->getMap ()->getComm ()->getRank ();
      if (myRank != 0) {
        beta_local = ATS::zero ();
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

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! C_tmp->isConstantStride () || ! B_tmp->isConstantStride () ||
       ! A_tmp->isConstantStride (), std::logic_error,
       "Failed to make temporary constant-stride copies of MultiVectors.");

    {
      if (A_tmp->need_sync_device ()) {
        const_cast<MV&> (*A_tmp).sync_device (); // const is a lie
      }
      const LO A_lclNumRows = A_tmp->getLocalLength ();
      const LO A_numVecs = A_tmp->getNumVectors ();
      auto A_lcl = A_tmp->getLocalViewDevice ();
      auto A_sub = Kokkos::subview (A_lcl,
                                    std::make_pair (LO (0), A_lclNumRows),
                                    std::make_pair (LO (0), A_numVecs));

      if (B_tmp->need_sync_device ()) {
        const_cast<MV&> (*B_tmp).sync_device (); // const is a lie
      }
      const LO B_lclNumRows = B_tmp->getLocalLength ();
      const LO B_numVecs = B_tmp->getNumVectors ();
      auto B_lcl = B_tmp->getLocalViewDevice ();
      auto B_sub = Kokkos::subview (B_lcl,
                                    std::make_pair (LO (0), B_lclNumRows),
                                    std::make_pair (LO (0), B_numVecs));

      if (C_tmp->need_sync_device ()) {
        const_cast<MV&> (*C_tmp).sync_device (); // const is a lie
      }
      const LO C_lclNumRows = C_tmp->getLocalLength ();
      const LO C_numVecs = C_tmp->getNumVectors ();
      auto C_lcl = C_tmp->getLocalViewDevice ();
      auto C_sub = Kokkos::subview (C_lcl,
                                    std::make_pair (LO (0), C_lclNumRows),
                                    std::make_pair (LO (0), C_numVecs));
      const char ctransA = (transA == Teuchos::NO_TRANS ? 'N' :
                            (transA == Teuchos::TRANS ? 'T' : 'C'));
      const char ctransB = (transB == Teuchos::NO_TRANS ? 'N' :
                            (transB == Teuchos::TRANS ? 'T' : 'C'));
      const impl_scalar_type alpha_IST (alpha);

      ProfilingRegion regionGemm ("Tpetra::MV::multiply-call-gemm");

      this->modify_device ();

      KokkosBlas::gemm (&ctransA, &ctransB, alpha_IST, A_sub, B_sub,
                        beta_local, C_sub);
    }

    if (! isConstantStride ()) {
      ::Tpetra::deep_copy (*this, *C_tmp); // Copy the result back into *this.
    }

    // Dispose of (possibly) extra copies of A and B.
    A_tmp = Teuchos::null;
    B_tmp = Teuchos::null;

    // If Case 2 then sum up *this and distribute it to all processes.
    if (Case2) {
      this->reduce ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  elementWiseMultiply (Scalar scalarAB,
                       const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                       const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                       Scalar scalarThis)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using V = Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const char tfecfFuncName[] = "elementWiseMultiply: ";

    const size_t lclNumRows = this->getLocalLength ();
    const size_t numVecs = this->getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclNumRows != A.getLocalLength () || lclNumRows != B.getLocalLength (),
       std::runtime_error, "MultiVectors do not have the same local length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numVecs != B.getNumVectors (), std::runtime_error, "this->getNumVectors"
      "() = " << numVecs << " != B.getNumVectors() = " << B.getNumVectors ()
      << ".");

    // All non-unary kernels are executed on the device as per Tpetra policy.  Sync to device if needed.
    if (this->need_sync_device ()) {
      this->sync_device ();
    }
    if (A.need_sync_device ()) {
      const_cast<V&>(A).sync_device ();
    }
    if (B.need_sync_device ()) {
      const_cast<MV&>(B).sync_device ();
    }
    this->modify_device ();

    auto this_view = this->getLocalViewDevice ();
    auto A_view = A.getLocalViewDevice ();
    auto B_view = B.getLocalViewDevice ();

    if (isConstantStride () && B.isConstantStride ()) {
      // A is just a Vector; it only has one column, so it always has
      // constant stride.
      //
      // If both *this and B have constant stride, we can do an
      // element-wise multiply on all columns at once.
      KokkosBlas::mult (scalarThis,
                        this_view,
                        scalarAB,
                        subview (A_view, ALL (), 0),
                        B_view);
    }
    else {
      for (size_t j = 0; j < numVecs; ++j) {
        const size_t C_col = isConstantStride () ? j : whichVectors_[j];
        const size_t B_col = B.isConstantStride () ? j : B.whichVectors_[j];
        KokkosBlas::mult (scalarThis,
                          subview (this_view, ALL (), C_col),
                          scalarAB,
                          subview (A_view, ALL (), 0),
                          subview (B_view, ALL (), B_col));
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  reduce ()
  {
    using ::Tpetra::Details::allReduceView;
    using ::Tpetra::Details::ProfilingRegion;
    ProfilingRegion region ("Tpetra::MV::reduce");

    const auto map = this->getMap ();
    if (map.get () == nullptr) {
      return;
    }
    const auto comm = map->getComm ();
    if (comm.get () == nullptr) {
      return;
    }

    // Avoid giving device buffers to MPI.  Even if MPI handles them
    // correctly, doing so may not perform well.
    const bool changed_on_device = this->need_sync_host ();
    if (changed_on_device) {
      // NOTE (mfh 17 Mar 2019) If we ever get rid of UVM, then device
      // and host will be separate allocations.  In that case, it may
      // pay to do the all-reduce from device to host.
      Kokkos::fence(); // for UVM getLocalViewDevice is UVM which can be read as host by allReduceView, so we must not read until device is fenced
      this->modify_device ();
      auto X_lcl = this->getLocalViewDevice ();
      allReduceView (X_lcl, X_lcl, *comm);
    }
    else {
      this->modify_host ();
      auto X_lcl = this->getLocalViewHost ();
      allReduceView (X_lcl, X_lcl, *comm);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceLocalValue (const LocalOrdinal lclRow,
                     const size_t col,
                     const impl_scalar_type& ScalarValue) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const LocalOrdinal minLocalIndex = this->getMap()->getMinLocalIndex();
    const LocalOrdinal maxLocalIndex = this->getMap()->getMaxLocalIndex();
    TEUCHOS_TEST_FOR_EXCEPTION(
      lclRow < minLocalIndex || lclRow > maxLocalIndex,
      std::runtime_error,
      "Tpetra::MultiVector::replaceLocalValue: row index " << lclRow
      << " is invalid.  The range of valid row indices on this process "
      << this->getMap()->getComm()->getRank() << " is [" << minLocalIndex
      << ", " << maxLocalIndex << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(col),
      std::runtime_error,
      "Tpetra::MultiVector::replaceLocalValue: vector index " << col
      << " of the multivector is invalid.");
#endif
    const size_t colInd = isConstantStride () ? col : whichVectors_[col];
    view_.h_view (lclRow, colInd) = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoLocalValue (const LocalOrdinal lclRow,
                     const size_t col,
                     const impl_scalar_type& value,
                     const bool atomic) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const LocalOrdinal minLocalIndex = this->getMap()->getMinLocalIndex();
    const LocalOrdinal maxLocalIndex = this->getMap()->getMaxLocalIndex();
    TEUCHOS_TEST_FOR_EXCEPTION(
      lclRow < minLocalIndex || lclRow > maxLocalIndex,
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoLocalValue: row index " << lclRow
      << " is invalid.  The range of valid row indices on this process "
      << this->getMap()->getComm()->getRank() << " is [" << minLocalIndex
      << ", " << maxLocalIndex << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(col),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoLocalValue: vector index " << col
      << " of the multivector is invalid.");
#endif
    const size_t colInd = isConstantStride () ? col : whichVectors_[col];
    if (atomic) {
      Kokkos::atomic_add (& (view_.h_view(lclRow, colInd)), value);
    }
    else {
      view_.h_view (lclRow, colInd) += value;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  replaceGlobalValue (const GlobalOrdinal gblRow,
                      const size_t col,
                      const impl_scalar_type& ScalarValue) const
  {
    // mfh 23 Nov 2015: Use map_ and not getMap(), because the latter
    // touches the RCP's reference count, which isn't thread safe.
    const LocalOrdinal lclRow = this->map_->getLocalElement (gblRow);
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "replaceGlobalValue: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (lclRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid (),
       std::runtime_error,
       "Global row index " << gblRow << " is not present on this process "
       << this->getMap ()->getComm ()->getRank () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->vectorIndexOutOfRange (col), std::runtime_error,
       "Vector index " << col << " of the MultiVector is invalid.");
#endif // HAVE_TPETRA_DEBUG
    this->replaceLocalValue (lclRow, col, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sumIntoGlobalValue (const GlobalOrdinal globalRow,
                      const size_t col,
                      const impl_scalar_type& value,
                      const bool atomic) const
  {
    // mfh 23 Nov 2015: Use map_ and not getMap(), because the latter
    // touches the RCP's reference count, which isn't thread safe.
    const LocalOrdinal lclRow = this->map_->getLocalElement (globalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      lclRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid (),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: Global row index " << globalRow
      << " is not present on this process "
      << this->getMap ()->getComm ()->getRank () << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(col),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: Vector index " << col
      << " of the multivector is invalid.");
#endif
    this->sumIntoLocalValue (lclRow, col, value, atomic);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class T>
  Teuchos::ArrayRCP<T>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
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


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  clear_sync_state () {
    view_.clear_sync_state ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sync_host () {
    view_.sync_host ();

    // This fence was motivated by the following specific situation:
    // For transform Y to X:
    //  Y.putScalar()    // acts on device
    //  Y.sync_host()    // now need_sync_host() and need_sync_device() are false
    //  transform (on device)
    //  Y.sync_host()    // no modifications so no fence - this usually will be a fence
    //  read Y           // crashes
    // The expectation is that Tpetra developers would not normally be using sync_host
    // so this fence should not be an issue for internal performance.
    execution_space ().fence ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  sync_device () {
    view_.sync_device ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  need_sync_host () const {
    return view_.need_sync_host ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  need_sync_device () const {
    return view_.need_sync_device ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  modify_device () {
    view_.modify_device ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  modify_host () {
    view_.modify_host ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalViewDevice () const {
    return view_.view_device ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_host
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getLocalViewHost () const {
    return view_.view_host ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  descriptionImpl (const std::string& className) const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream out;
    out << "\"" << className << "\": {";
    out << "Template parameters: {Scalar: " << TypeNameTraits<Scalar>::name ()
        << ", LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name ()
        << ", GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name ()
        << ", Node" << Node::name ()
        << "}, ";
    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\", ";
    }
    out << ", numRows: " << this->getGlobalLength ();
    if (className != "Tpetra::Vector") {
      out << ", numCols: " << this->getNumVectors ()
          << ", isConstantStride: " << this->isConstantStride ();
    }
    if (this->isConstantStride ()) {
      out << ", columnStride: " << this->getStride ();
    }
    out << "}";

    return out.str ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  description () const
  {
    return this->descriptionImpl ("Tpetra::MultiVector");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  localDescribeToString (const Teuchos::EVerbosityLevel vl) const
  {
    typedef LocalOrdinal LO;
    using std::endl;

    if (vl <= Teuchos::VERB_LOW) {
      return std::string ();
    }
    auto map = this->getMap ();
    if (map.is_null ()) {
      return std::string ();
    }
    auto outStringP = Teuchos::rcp (new std::ostringstream ());
    auto outp = Teuchos::getFancyOStream (outStringP);
    Teuchos::FancyOStream& out = *outp;
    auto comm = map->getComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    out << "Process " << myRank << " of " << numProcs << ":" << endl;
    Teuchos::OSTab tab1 (out);

    // VERB_MEDIUM and higher prints getLocalLength()
    const LO lclNumRows = static_cast<LO> (this->getLocalLength ());
    out << "Local number of rows: " << lclNumRows << endl;

    if (vl > Teuchos::VERB_MEDIUM) {
      // VERB_HIGH and higher prints isConstantStride() and getStride().
      // The first is only relevant if the Vector has multiple columns.
      if (this->getNumVectors () != static_cast<size_t> (1)) {
        out << "isConstantStride: " << this->isConstantStride () << endl;
      }
      if (this->isConstantStride ()) {
        out << "Column stride: " << this->getStride () << endl;
      }

      if (vl > Teuchos::VERB_HIGH && lclNumRows > 0) {
        // VERB_EXTREME prints values.  Get a host View of the
        // Vector's local data, so we can print it.  (Don't assume
        // that we can access device data directly in host code.)
        typename dual_view_type::t_host X_host;
        if (need_sync_host ()) {
          // Device memory has the latest version.  Don't actually
          // sync to host; that changes the Vector's state, and may
          // change future computations (that use the data's current
          // place to decide where to compute).  Instead, create a
          // temporary host copy and print that.
          auto X_dev = getLocalViewDevice ();
          auto X_host_copy = Kokkos::create_mirror_view (X_dev);
          Kokkos::deep_copy (X_host_copy, X_dev);
          X_host = X_host_copy;
        }
        else {
          // Either host and device are in sync, or host has the
          // latest version of the Vector's data.  Thus, we can use
          // the host version directly.
          X_host = getLocalViewHost ();
        }
        // The square braces [] and their contents are in Matlab
        // format, so users may copy and paste directly into Matlab.
        out << "Values: " << endl
            << "[";
        const LO numCols = static_cast<LO> (this->getNumVectors ());
        if (numCols == 1) {
          for (LO i = 0; i < lclNumRows; ++i) {
            out << X_host(i,0);
            if (i + 1 < lclNumRows) {
              out << "; ";
            }
          }
        }
        else {
          for (LO i = 0; i < lclNumRows; ++i) {
            for (LO j = 0; j < numCols; ++j) {
              out << X_host(i,j);
              if (j + 1 < numCols) {
                out << ", ";
              }
            }
            if (i + 1 < lclNumRows) {
              out << "; ";
            }
          }
        }
        out << "]" << endl;
      }
    }

    out.flush (); // make sure the ostringstream got everything
    return outStringP->str ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  describeImpl (Teuchos::FancyOStream& out,
                const std::string& className,
                const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using std::endl;
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // don't print anything
    }
    // If this Vector's Comm is null, then the Vector does not
    // participate in collective operations with the other processes.
    // In that case, it is not even legal to call this method.  The
    // reasonable thing to do in that case is nothing.
    auto map = this->getMap ();
    if (map.is_null ()) {
      return;
    }
    auto comm = map->getComm ();
    if (comm.is_null ()) {
      return;
    }

    const int myRank = comm->getRank ();

    // Only Process 0 should touch the output stream, but this method
    // in general may need to do communication.  Thus, we may need to
    // preserve the current tab level across multiple "if (myRank ==
    // 0) { ... }" inner scopes.  This is why we sometimes create
    // OSTab instances by pointer, instead of by value.  We only need
    // to create them by pointer if the tab level must persist through
    // multiple inner scopes.
    Teuchos::RCP<Teuchos::OSTab> tab0, tab1;

    // VERB_LOW and higher prints the equivalent of description().
    if (myRank == 0) {
      tab0 = Teuchos::rcp (new Teuchos::OSTab (out));
      out << "\"" << className << "\":" << endl;
      tab1 = Teuchos::rcp (new Teuchos::OSTab (out));
      {
        out << "Template parameters:" << endl;
        Teuchos::OSTab tab2 (out);
        out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
            << "LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name () << endl
            << "GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name () << endl
            << "Node: " << Node::name () << endl;
      }
      if (this->getObjectLabel () != "") {
        out << "Label: \"" << this->getObjectLabel () << "\", ";
      }
      out << "Global number of rows: " << this->getGlobalLength () << endl;
      if (className != "Tpetra::Vector") {
        out << "Number of columns: " << this->getNumVectors () << endl;
      }
      // getStride() may differ on different processes, so it (and
      // isConstantStride()) properly belong to per-process data.
    }

    // This is collective over the Map's communicator.
    if (vl > VERB_LOW) { // VERB_MEDIUM, VERB_HIGH, or VERB_EXTREME
      const std::string lclStr = this->localDescribeToString (vl);
      ::Tpetra::Details::gathervPrint (out, lclStr, *comm);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    this->describeImpl (out, "Tpetra::MultiVector", verbLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    replaceMap (newMap);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src)
  {
    using ::Tpetra::Details::localDeepCopy;
    const char prefix[] = "Tpetra::MultiVector::assign: ";

    TEUCHOS_TEST_FOR_EXCEPTION
      (this->getGlobalLength () != src.getGlobalLength () ||
       this->getNumVectors () != src.getNumVectors (), std::invalid_argument,
       prefix << "Global dimensions of the two Tpetra::MultiVector "
       "objects do not match.  src has dimensions [" << src.getGlobalLength ()
       << "," << src.getNumVectors () << "], and *this has dimensions ["
       << this->getGlobalLength () << "," << this->getNumVectors () << "].");
    // FIXME (mfh 28 Jul 2014) Don't throw; just set a local error flag.
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->getLocalLength () != src.getLocalLength (), std::invalid_argument,
       prefix << "The local row counts of the two Tpetra::MultiVector "
       "objects do not match.  src has " << src.getLocalLength () << " row(s) "
       << " and *this has " << this->getLocalLength () << " row(s).");

    // See #1510.  We're writing to *this, so we don't care about its
    // contents in either memory space, and we don't want
    // DualView::modify to complain about "concurrent modification" of
    // host and device Views.
    this->clear_sync_state();
    this->modify_device ();

    // If need sync to device, then host has most recent version.
    const bool src_last_updated_on_host = src.need_sync_device ();

    if (src_last_updated_on_host) {
      localDeepCopy (this->getLocalViewDevice (),
                     src.getLocalViewHost (),
                     this->isConstantStride (),
                     src.isConstantStride (),
                     this->whichVectors_ (),
                     src.whichVectors_ ());
    }
    else {
      localDeepCopy (this->getLocalViewDevice (),
                     src.getLocalViewDevice (),
                     this->isConstantStride (),
                     src.isConstantStride (),
                     this->whichVectors_ (),
                     src.whichVectors_ ());
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isSameSize (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& vec) const
  {
    using ::Tpetra::Details::PackTraits;
    using ST = impl_scalar_type;

    const size_t l1 = this->getLocalLength();
    const size_t l2 = vec.getLocalLength();
    if ((l1!=l2) || (this->getNumVectors() != vec.getNumVectors())) {
      return false;
    }
    if (l1==0) {
      return true;
    }

    auto v1 = this->getLocalViewHost ();
    auto v2 = vec.getLocalViewHost ();
    if (PackTraits<ST>::packValueCount (v1(0,0)) !=
        PackTraits<ST>::packValueCount (v2(0,0))) {
      return false;
    }

    return true;
  }

  template <class ST, class LO, class GO, class NT>
  void MultiVector<ST, LO, GO, NT>::
  swap(MultiVector<ST, LO, GO, NT> & mv) {
    std::swap(mv.map_, this->map_);
    std::swap(mv.view_, this->view_);
    std::swap(mv.origView_, this->origView_);
    std::swap(mv.whichVectors_, this->whichVectors_);
  }

#ifdef HAVE_TPETRACORE_TEUCHOSNUMERICS
  template <class ST, class LO, class GO, class NT>
  void
  deep_copy (MultiVector<ST, LO, GO, NT>& dst,
             const Teuchos::SerialDenseMatrix<int, ST>& src)
  {
    using ::Tpetra::Details::localDeepCopy;
    using MV = MultiVector<ST, LO, GO, NT>;
    using IST = typename MV::impl_scalar_type;
    using input_view_type =
      Kokkos::View<const IST**, Kokkos::LayoutLeft,
        Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
    using pair_type = std::pair<LO, LO>;

    const LO numRows = static_cast<LO> (src.numRows ());
    const LO numCols = static_cast<LO> (src.numCols ());

    TEUCHOS_TEST_FOR_EXCEPTION
      (numRows != static_cast<LO> (dst.getLocalLength ()) ||
       numCols != static_cast<LO> (dst.getNumVectors ()),
       std::invalid_argument, "Tpetra::deep_copy: On Process "
       << dst.getMap ()->getComm ()->getRank () << ", dst is "
       << dst.getLocalLength () << " x " << dst.getNumVectors ()
       << ", but src is " << numRows << " x " << numCols << ".");

    const IST* src_raw = reinterpret_cast<const IST*> (src.values ());
    input_view_type src_orig (src_raw, src.stride (), numCols);
    input_view_type src_view =
      Kokkos::subview (src_orig, pair_type (0, numRows), Kokkos::ALL ());

    dst.clear_sync_state ();
    dst.modify_device ();
    constexpr bool src_isConstantStride = true;
    Teuchos::ArrayView<const size_t> srcWhichVectors (nullptr, 0);
    localDeepCopy (dst.getLocalViewDevice (),
                   src_view,
                   dst.isConstantStride (),
                   src_isConstantStride,
                   getMultiVectorWhichVectors (dst),
                   srcWhichVectors);
  }

  template <class ST, class LO, class GO, class NT>
  void
  deep_copy (Teuchos::SerialDenseMatrix<int, ST>& dst,
             const MultiVector<ST, LO, GO, NT>& src)
  {
    using ::Tpetra::Details::localDeepCopy;
    using MV = MultiVector<ST, LO, GO, NT>;
    using IST = typename MV::impl_scalar_type;
    using output_view_type =
      Kokkos::View<IST**, Kokkos::LayoutLeft,
        Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
    using pair_type = std::pair<LO, LO>;

    const LO numRows = static_cast<LO> (dst.numRows ());
    const LO numCols = static_cast<LO> (dst.numCols ());

    TEUCHOS_TEST_FOR_EXCEPTION
      (numRows != static_cast<LO> (src.getLocalLength ()) ||
       numCols != static_cast<LO> (src.getNumVectors ()),
       std::invalid_argument, "Tpetra::deep_copy: On Process "
       << src.getMap ()->getComm ()->getRank () << ", src is "
       << src.getLocalLength () << " x " << src.getNumVectors ()
       << ", but dst is " << numRows << " x " << numCols << ".");

    IST* dst_raw = reinterpret_cast<IST*> (dst.values ());
    output_view_type dst_orig (dst_raw, dst.stride (), numCols);
    auto dst_view =
      Kokkos::subview (dst_orig, pair_type (0, numRows), Kokkos::ALL ());

    constexpr bool dst_isConstantStride = true;
    Teuchos::ArrayView<const size_t> dstWhichVectors (nullptr, 0);

    // Prefer the host version of src's data.
    if (src.need_sync_host ()) { // last modified on device
      localDeepCopy (dst_view,
                     src.getLocalViewDevice (),
                     dst_isConstantStride,
                     src.isConstantStride (),
                     dstWhichVectors,
                     getMultiVectorWhichVectors (src));
    }
    else {
      localDeepCopy (dst_view,
                     src.getLocalViewHost (),
                     dst_isConstantStride,
                     src.isConstantStride (),
                     dstWhichVectors,
                     getMultiVectorWhichVectors (src));
    }
  }
#endif // HAVE_TPETRACORE_TEUCHOSNUMERICS

  template <class Scalar, class LO, class GO, class NT>
  Teuchos::RCP<MultiVector<Scalar, LO, GO, NT> >
  createMultiVector (const Teuchos::RCP<const Map<LO, GO, NT> >& map,
                     size_t numVectors)
  {
    typedef MultiVector<Scalar, LO, GO, NT> MV;
    return Teuchos::rcp (new MV (map, numVectors));
  }

  template <class ST, class LO, class GO, class NT>
  MultiVector<ST, LO, GO, NT>
  createCopy (const MultiVector<ST, LO, GO, NT>& src)
  {
    typedef MultiVector<ST, LO, GO, NT> MV;
    MV cpy (src.getMap (), src.getNumVectors (), false);
    cpy.assign (src);
    return cpy;
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#ifdef HAVE_TPETRACORE_TEUCHOSNUMERICS
#  define TPETRA_MULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class MultiVector< SCALAR , LO , GO , NODE >; \
  template MultiVector< SCALAR , LO , GO , NODE > createCopy( const MultiVector< SCALAR , LO , GO , NODE >& src); \
  template Teuchos::RCP<MultiVector< SCALAR , LO , GO , NODE > > createMultiVector (const Teuchos::RCP<const Map<LO, GO, NODE> >& map, size_t numVectors); \
  template void deep_copy (MultiVector<SCALAR, LO, GO, NODE>& dst, const Teuchos::SerialDenseMatrix<int, SCALAR>& src); \
  template void deep_copy (Teuchos::SerialDenseMatrix<int, SCALAR>& dst, const MultiVector<SCALAR, LO, GO, NODE>& src);

#else
#  define TPETRA_MULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class MultiVector< SCALAR , LO , GO , NODE >; \
  template MultiVector< SCALAR , LO , GO , NODE > createCopy( const MultiVector< SCALAR , LO , GO , NODE >& src);

#endif // HAVE_TPETRACORE_TEUCHOSNUMERICS

#endif // TPETRA_MULTIVECTOR_DEF_HPP
