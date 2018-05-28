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

/// \file Tpetra_MultiVector_def.hpp
/// \brief Definition of the Tpetra::MultiVector class
///
/// If you want to use Tpetra::MultiVector, include
/// "Tpetra_MultiVector.hpp" (a file which CMake generates and
/// installs for you).  If you only want the declaration of
/// Tpetra::MultiVector, include "Tpetra_MultiVector_decl.hpp".

#include "Tpetra_Util.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_castAwayConstDualView.hpp"
#include "Tpetra_Details_fill.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_gemm.hpp"
#include "Tpetra_Details_isInterComm.hpp"
#include "Tpetra_Details_lclDot.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp"



#include "KokkosCompat_View.hpp"
#include "KokkosBlas.hpp"
#include "KokkosKernels_Utils.hpp"
#include "Kokkos_Random.hpp"

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

    typename dual_view_type::t_host h_view = Kokkos::create_mirror_view (d_view);

    dual_view_type dv (d_view, h_view);
    // Whether or not the user cares about the initial contents of the
    // MultiVector, the device and host views are out of sync.  We
    // prefer to work in device memory.  The way to ensure this
    // happens is to mark the device view as modified.
    dv.template modify<typename dev_view_type::memory_space> ();

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
  MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source) :
    base_type (source),
    view_ (source.view_),
    origView_ (source.origView_),
    whichVectors_ (source.whichVectors_)
  {}

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
    const char tfecfFuncName[] = "MultiVector(map,view): ";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.extent (1) > 1) ? stride[1] :
      origView_.extent (0);
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::DualView's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input view's first dimension (number "
      "of rows = " << view.extent (0) << ") does not not match the number "
      "of entries in the Map on the calling process.");
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

    // Get stride of view: if second dimension is 0, the stride might
    // be 0, so take view_dimension instead.
    size_t stride[8];
    d_view.stride (stride);
    const size_t LDA = (d_view.extent (1) > 1) ? stride[1] :
      d_view.extent (0);
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::View's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input view's first dimension (number "
      "of rows = " << d_view.extent (0) << ") does not not match the "
      "number of entries in the Map on the calling process.");

    // The difference between create_mirror and create_mirror_view, is
    // that the latter copies to host.  We don't necessarily want to
    // do that; we just want to allocate the memory.
    view_ = dual_view_type (d_view, Kokkos::create_mirror (d_view));
    // The user gave us a device view.  We take it as canonical, which
    // means we mark it as "modified," so that the next sync will
    // synchronize it with the host view.
    this->template modify<device_type> ();
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

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.extent (1) > 1) ? stride[1] :
      origView_.extent (0);
    const size_t lclNumRows = this->getLocalLength (); // comes from the Map
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < lclNumRows, std::invalid_argument, "The input Kokkos::DualView's "
      "column stride LDA = " << LDA << " < getLocalLength() = " << lclNumRows
      << ".  This may also mean that the input origView's first dimension (number "
      "of rows = " << origView.extent (0) << ") does not not match the number "
      "of entries in the Map on the calling process.");
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

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.extent (1) > 1) ? stride[1] :
      origView_.extent (0);
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
    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    origView_.stride (stride);
    const size_t LDA = (origView_.extent (1) > 1) ? stride[1] :
      origView_.extent (0);
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
    this->template modify<device_type> ();
    auto X_out = this->template getLocalView<device_type> ();
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
    size_t outStrides[8];
    X_out.stride (outStrides);
    const size_t outStride = (X_out.extent (1) > 1) ? outStrides[1] :
      X_out.extent (0);
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
    this->template modify<device_type> ();
    auto X_out = this->getLocalView<device_type> ();

    // Make sure that the type of a single input column has the same
    // array layout as each output column does, so we can deep_copy.
    typedef typename decltype (X_out)::array_layout array_layout;
    typedef typename Kokkos::View<const IST*,
      array_layout,
      Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged> input_col_view_type;

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
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ~MultiVector () {}

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
    if (isConstantStride ()) {
      // Get stride of view: if second dimension is 0, the
      // stride might be 0, so take view_dimension instead.
      size_t stride[8];
      origView_.stride (stride);
      const size_t LDA = (origView_.extent (1) > 1) ? stride[1] :
        origView_.extent (0);
      return LDA;
    }
    else {
      return static_cast<size_t> (0);
    }
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
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  constantNumberOfPackets () const {
    return this->getNumVectors ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermuteNew (const SrcDistObject& sourceObj,
                     const size_t numSameIDs,
                     const Kokkos::DualView<const LocalOrdinal*, device_type>& permuteToLIDs,
                     const Kokkos::DualView<const LocalOrdinal*, device_type>& permuteFromLIDs)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::castAwayConstDualView;
    using ::Tpetra::Details::getDualViewCopyFromArrayView;
    using ::Tpetra::Details::ProfilingRegion;
    using KokkosRefactor::Details::permute_array_multi_column;
    using KokkosRefactor::Details::permute_array_multi_column_variable_stride;
    using Kokkos::Compat::create_const_view;
    typedef typename device_type::memory_space DMS;
    typedef Kokkos::HostSpace HMS;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "copyAndPermuteNew: ";
    ProfilingRegion regionCAP ("Tpetra::MultiVector::copyAndPermute");

    // mfh 03 Aug 2017, 27 Sep 2017: Set the TPETRA_VERBOSE
    // environment variable to "1" (or "TRUE") for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool debug = Behavior::verbose ();
    int myRank = 0;
    if (debug && ! this->getMap ().is_null () &&
        ! this->getMap ()->getComm ().is_null ()) {
      myRank = this->getMap ()->getComm ()->getRank ();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (permuteToLIDs.extent (0) != permuteFromLIDs.extent (0),
       std::runtime_error, "permuteToLIDs.extent(0) = "
       << permuteToLIDs.extent (0) << " != permuteFromLIDs.extent(0) = "
       << permuteFromLIDs.extent (0) << ".");

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);
    const size_t numCols = this->getNumVectors ();

    // The input sourceObj controls whether we copy on host or device.
    // This is because this method is not allowed to modify sourceObj,
    // so it may not sync sourceObj; it must use the most recently
    // updated version (host or device) of sourceObj's data.
    //
    // If we need sync to device, then host has the most recent version.
    const bool copyOnHost = sourceMV.template need_sync<DMS> ();
    if (debug) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
        "copyOnHost=" << (copyOnHost ? "true" : "false") << std::endl;
      std::cerr << os.str ();
    }

    if (copyOnHost) {
      this->template sync<HMS> ();
      this->template modify<HMS> ();
    }
    else {
      this->template sync<DMS> ();
      this->template modify<DMS> ();
    }

    if (debug) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
        "Copy source to target" << std::endl;
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

    if (numSameIDs > 0) {
      const std::pair<size_t, size_t> rows (0, numSameIDs);
      if (copyOnHost) {
        auto tgt_h = this->template getLocalView<HMS> ();
        auto src_h = create_const_view (sourceMV.template getLocalView<HMS> ());

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
        auto tgt_d = this->template getLocalView<DMS> ();
        auto src_d = create_const_view (sourceMV.template getLocalView<DMS> ());

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
      if (debug) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
          "No permutations; done." << std::endl;
        std::cerr << os.str ();
      }
      return;
    }

    if (debug) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
        "Permute source to target" << std::endl;
      std::cerr << os.str ();
    }

    // mfh 29 Aug 2017: Kokkos::DualView annoyingly forbids (at run
    // time, no less!) sync'ing a DualView<const T, D>.  We work
    // around this by managing the copying and flags ourselves --
    // essentially reimplementing sync and modify.
    Kokkos::DualView<LocalOrdinal*, device_type> permuteToLIDs_nc =
      castAwayConstDualView (permuteToLIDs);
    Kokkos::DualView<LocalOrdinal*, device_type> permuteFromLIDs_nc =
      castAwayConstDualView (permuteFromLIDs);

    // We could in theory optimize for the case where exactly one of
    // them is constant stride, but we don't currently do that.
    const bool nonConstStride =
      ! this->isConstantStride () || ! sourceMV.isConstantStride ();

    if (debug) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
        "nonConstStride=" << (nonConstStride ? "true" : "false") << std::endl;
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
        tmpTgt.template modify<HMS> ();
        for (size_t j = 0; j < numCols; ++j) {
          tmpTgt.h_view(j) = j;
        }
        if (! copyOnHost) {
          tmpTgt.template sync<DMS> ();
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
        tmpSrc.template modify<HMS> ();
        for (size_t j = 0; j < numCols; ++j) {
          tmpSrc.h_view(j) = j;
        }
        if (! copyOnHost) {
          tmpSrc.template sync<DMS> ();
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
      if (debug) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
          "Set up permutation arrays on host" << std::endl;
        std::cerr << os.str ();
      }
      auto tgt_h = this->template getLocalView<HMS> ();
      auto src_h = create_const_view (sourceMV.template getLocalView<HMS> ());
      permuteToLIDs_nc.template sync<HMS> ();
      auto permuteToLIDs_h =
        create_const_view (permuteToLIDs_nc.template view<HMS> ());
      permuteFromLIDs_nc.template sync<HMS> ();
      auto permuteFromLIDs_h =
        create_const_view (permuteFromLIDs_nc.template view<HMS> ());

      if (debug) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
          "Call permute kernel on host" << std::endl;
        std::cerr << os.str ();
      }
      if (nonConstStride) {
        // No need to sync first, because copyOnHost argument to
        // getDualViewCopyFromArrayView puts them in the right place.
        auto tgtWhichVecs_h =
          create_const_view (tgtWhichVecs.template view<HMS> ());
        auto srcWhichVecs_h =
          create_const_view (srcWhichVecs.template view<HMS> ());
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
      if (debug) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
          "Set up permutation arrays on device" << std::endl;
        std::cerr << os.str ();
      }
      auto tgt_d = this->template getLocalView<DMS> ();
      auto src_d = create_const_view (sourceMV.template getLocalView<DMS> ());
      permuteToLIDs_nc.template sync<DMS> ();
      auto permuteToLIDs_d =
        create_const_view (permuteToLIDs_nc.template view<DMS> ());
      permuteFromLIDs_nc.template sync<DMS> ();
      auto permuteFromLIDs_d =
        create_const_view (permuteFromLIDs_nc.template view<DMS> ());

      if (debug) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::copyAndPermuteNew: "
          "Call permute kernel on device" << std::endl;
        std::cerr << os.str ();
      }
      if (nonConstStride) {
        // No need to sync first, because copyOnHost argument to
        // getDualViewCopyFromArrayView puts them in the right place.
        auto tgtWhichVecs_d =
          create_const_view (tgtWhichVecs.template view<DMS> ());
        auto srcWhichVecs_d =
          create_const_view (srcWhichVecs.template view<DMS> ());
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

    if (debug) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::copyAndPermuteNew: Done" << std::endl;
      std::cerr << os.str ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packAndPrepareNew (const SrcDistObject& sourceObj,
                     const Kokkos::DualView<const local_ordinal_type*, device_type>& exportLIDs,
                     Kokkos::DualView<impl_scalar_type*, buffer_device_type>& exports,
                     const Kokkos::DualView<size_t*, buffer_device_type>& /* numExportPacketsPerLID */,
                     size_t& constantNumPackets,
                     Distributor & /* distor */ )
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::ProfilingRegion;
    using Kokkos::Compat::create_const_view;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    typedef impl_scalar_type IST;
    typedef Kokkos::HostSpace host_memory_space;
    typedef typename Kokkos::DualView<IST*, device_type>::t_dev::memory_space
      dev_memory_space;
    typedef typename Kokkos::DualView<IST*, device_type>::t_host::execution_space
      host_execution_space;
    typedef typename Kokkos::DualView<IST*, device_type>::t_dev::execution_space
      dev_execution_space;
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

    int myRank = 0;
    if (printDebugOutput && ! this->getMap ().is_null () &&
        ! this->getMap ()->getComm ().is_null ()) {
      myRank = this->getMap ()->getComm ()->getRank ();
    }

    if (printDebugOutput) {
      std::ostringstream os;
      os << "Proc " << myRank << ": MV::packAndPrepareNew" << std::endl;
      std::cerr << os.str ();
    }
    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    // mfh 12 Apr 2016: packAndPrepareNew decides where to pack based
    // on the memory space in which exportLIDs was last modified.
    // (DistObject::doTransferNew decides this currently.)
    //
    // We unfortunately can't change the source object sourceMV.
    // Thus, we can't sync it to the memory space where we want to
    // pack communication buffers.  As a result, for example, if
    // exportLIDs wants us to pack on host, but sourceMV was last
    // modified on device, we have to allocate a temporary host
    // version and copy back to host before we can pack.  Similarly,
    // if exportLIDs wants us to pack on device, but sourceMV was last
    // modified on host, we have to allocate a temporary device
    // version and copy back to device before we can pack.
    const bool packOnHost =
      exportLIDs.modified_host () > exportLIDs.modified_device ();
    auto src_dev = sourceMV.template getLocalView<dev_memory_space> ();
    auto src_host = sourceMV.template getLocalView<host_memory_space> ();
    if (packOnHost) {
      if (sourceMV.template need_sync<Kokkos::HostSpace> ()) {
        if (printDebugOutput) {
          std::ostringstream os;
          os << "Proc " << myRank << ": MV::packAndPrepareNew: "
            "Pack on host, need sync to host" << std::endl;
          std::cerr << os.str ();
        }
        // sourceMV was most recently updated on device; copy to host.
        // Allocate a new host mirror.  We'll use it for packing below.
        src_host = decltype (src_host) ("MV::DualView::h_view",
                                        src_dev.extent (0),
                                        src_dev.extent (1));
        Kokkos::deep_copy (src_host, src_dev);
      }
    }
    else { // pack on device
      if (sourceMV.template need_sync<device_type> ()) {
        if (printDebugOutput) {
          std::ostringstream os;
          os << "Proc " << myRank << ": MV::packAndPrepareNew: "
            "Pack on device, need sync to device" << std::endl;
          std::cerr << os.str ();
        }
        // sourceMV was most recently updated on host; copy to device.
        // Allocate a new "device mirror."  We'll use it for packing below.
        src_dev = decltype (src_dev) ("MV::DualView::d_view",
                                      src_host.extent (0),
                                      src_host.extent (1));
        Kokkos::deep_copy (src_dev, src_host);
      }
    }

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
        os << "Proc " << myRank << ": MV::packAndPrepareNew: "
          "No exports on this proc, DONE" << std::endl;
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
      os << "Proc " << myRank << ": MV::packAndPrepareNew: realloc: "
         << "numExportLIDs = " << numExportLIDs
         << ", exports.extent(0) = " << exports.extent (0)
         << ", newExportsSize = " << newExportsSize << std::endl;
      std::cerr << os.str ();
    }
    Details::reallocDualViewIfNeeded (exports, newExportsSize, "exports");

    // 'exports' may have different memory spaces than device_type
    // would indicate.  See GitHub issue #1088.  Abbreviations:
    // "exports host memory space" and "exports device memory spaces."
    typedef typename std::decay<decltype (exports) >::type::t_host::memory_space EHMS;
    typedef typename std::decay<decltype (exports) >::type::t_dev::memory_space EDMS;

    // Mark 'exports' here, since we might have resized it above.
    // Resizing currently requires calling the constructor, which
    // clears out the 'modified' flags.
    if (packOnHost) {
      exports.template modify<EHMS> ();
    }
    else {
      exports.template modify<EDMS> ();
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
          os << "Proc " << myRank << ": MV::packAndPrepareNew: "
            "pack numCols=1 const stride" << std::endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_single_column (exports.template view<EHMS> (),
                                    create_const_view (src_host),
                                    exportLIDs.template view<host_memory_space> (),
                                    0,
                                    debugCheckIndices);
        }
        else { // pack on device
          pack_array_single_column (exports.template view<EDMS> (),
                                    create_const_view (src_dev),
                                    exportLIDs.template view<dev_memory_space> (),
                                    0,
                                    debugCheckIndices);
        }
      }
      else {
        using KokkosRefactor::Details::pack_array_single_column;
        if (printDebugOutput) {
          std::ostringstream os;
          os << "Proc " << myRank << ": MV::packAndPrepareNew: "
            "pack numCols=1 nonconst stride" << std::endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_single_column (exports.template view<EHMS> (),
                                    create_const_view (src_host),
                                    exportLIDs.template view<host_memory_space> (),
                                    sourceMV.whichVectors_[0],
                                    debugCheckIndices);
        }
        else { // pack on device
          pack_array_single_column (exports.template view<EDMS> (),
                                    create_const_view (src_dev),
                                    exportLIDs.template view<dev_memory_space> (),
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
          os << "Proc " << myRank << ": MV::packAndPrepareNew: "
            "pack numCols>1 const stride" << std::endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_multi_column (exports.template view<EHMS> (),
                                   create_const_view (src_host),
                                   exportLIDs.template view<host_memory_space> (),
                                   numCols,
                                   debugCheckIndices);
        }
        else { // pack on device
          pack_array_multi_column (exports.template view<EDMS> (),
                                   create_const_view (src_dev),
                                   exportLIDs.template view<dev_memory_space> (),
                                   numCols,
                                   debugCheckIndices);
        }
      }
      else {
        using KokkosRefactor::Details::pack_array_multi_column_variable_stride;
        if (printDebugOutput) {
          std::ostringstream os;
          os << "Proc " << myRank << ": MV::packAndPrepareNew: "
            "pack numCols>1 nonconst stride" << std::endl;
          std::cerr << os.str ();
        }
        if (packOnHost) {
          pack_array_multi_column_variable_stride (exports.template view<EHMS> (),
                                                   create_const_view (src_host),
                                                   exportLIDs.template view<host_memory_space> (),
                                                   getKokkosViewDeepCopy<host_execution_space> (sourceMV.whichVectors_ ()),
                                                   numCols,
                                                   debugCheckIndices);
        }
        else { // pack on device
          pack_array_multi_column_variable_stride (exports.template view<EDMS> (),
                                                   create_const_view (src_dev),
                                                   exportLIDs.template view<dev_memory_space> (),
                                                   getKokkosViewDeepCopy<dev_execution_space> (sourceMV.whichVectors_ ()),
                                                   numCols,
                                                   debugCheckIndices);
        }
      }
    }

    if (printDebugOutput) {
      std::ostringstream os;
      os << "Proc " << myRank << ": MV::packAndPrepareNew: DONE" << std::endl;
      std::cerr << os.str ();
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombineNew (const Kokkos::DualView<const local_ordinal_type*, device_type>& importLIDs,
                       const Kokkos::DualView<const impl_scalar_type*, buffer_device_type>& imports,
                       const Kokkos::DualView<const size_t*, buffer_device_type>& /* numPacketsPerLID */,
                       const size_t constantNumPackets,
                       Distributor& /* distor */,
                       const CombineMode CM)
  {
    using ::Tpetra::Details::Behavior;
    using ::Tpetra::Details::ProfilingRegion;
    using ::Tpetra::Details::castAwayConstDualView;
    using KokkosRefactor::Details::unpack_array_multi_column;
    using KokkosRefactor::Details::unpack_array_multi_column_variable_stride;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef impl_scalar_type IST;
    typedef typename Kokkos::DualView<IST*,
      device_type>::t_dev::memory_space DMS;
    // For correct UVM use, make the "host memory space" (template
    // parameter of sync and modify) different than the "device memory
    // space."  Otherwise, sync() won't fence (indeed, it won't do
    // anything).
    typedef Kokkos::HostSpace HMS;
    const char tfecfFuncName[] = "unpackAndCombineNew: ";
    ProfilingRegion regionUAC ("Tpetra::MultiVector::unpackAndCombine");

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
    int myRank = 0;
    if (printDebugOutput && ! this->getMap ().is_null () &&
        ! this->getMap ()->getComm ().is_null ()) {
      myRank = this->getMap ()->getComm ()->getRank ();
    }

    // If we have no imports, there is nothing to do
    if (importLIDs.extent (0) == 0) {
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

    // mfh 12 Apr 2016: Decide where to unpack based on the memory
    // space in which the imports buffer was last modified.
    // DistObject::doTransferNew gets to decide this.  We currently
    // require importLIDs to match (its most recent version must be in
    // the same memory space as imports' most recent version).
    const bool unpackOnHost =
      imports.modified_host () != 0 &&
      imports.modified_host () > imports.modified_device ();

    if (printDebugOutput) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::unpackAndCombine: unpackOnHost: "
         << (unpackOnHost ? "true" : "false") << std::endl;
      std::cerr << os.str ();
    }

    // Kokkos::DualView<const T*, ...> forbids sync, so "cast away
    // const."  It still makes sense for the importLIDs to be const,
    // because this method does not _modify_ them.  We just need them
    // to be in the right place.
    auto theImportLIDs = castAwayConstDualView (importLIDs);
    if (unpackOnHost && importLIDs.modified_host () < importLIDs.modified_device ()) {
      if (printDebugOutput) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::unpackAndCombine: sync importLIDs "
          "to host" << std::endl;
        std::cerr << os.str ();
      }
      // 'imports' was last modified on host, but importLIDs was last
      // modified on device.  Sync importLIDs to host and do the work
      // there, since imports usually at least as much data (and often
      // has more) than importLIDs.
      theImportLIDs.template sync<HMS> ();
    }
    else if (! unpackOnHost && importLIDs.modified_host () > importLIDs.modified_device ()) {
      if (printDebugOutput) {
        std::ostringstream os;
        os << "(Proc " << myRank << ") MV::unpackAndCombine: sync importLIDs "
          "to device" << std::endl;
        std::cerr << os.str ();
      }
      // 'imports' was last modified on device, but importLIDs was
      // last modified on host.  Sync importLIDs to device and do the
      // work there, since imports usually at least as much data (and
      // often has more) than importLIDs.
      theImportLIDs.template sync<DMS> ();
    }
    else if (printDebugOutput) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::unpackAndCombine: no need to sync "
        "importLIDs" << std::endl;
      std::cerr << os.str ();
    }

    // We have to sync before modifying, because this method may read
    // as well as write (depending on the CombineMode).  This matters
    // because copyAndPermute may have modified *this in the other
    // memory space.
    if (unpackOnHost) {
      this->template sync<HMS> ();
      this->template modify<HMS> ();
    }
    else { // unpack on device
      this->template sync<DMS> ();
      this->template modify<DMS> ();
    }
    auto X_d = this->template getLocalView<DMS> ();
    auto X_h = this->template getLocalView<HMS> ();
    // 'imports' may have a different device memory space (see #1088).
    auto imports_d =
      imports.template view<typename buffer_device_type::memory_space> ();
    auto imports_h = imports.template view<Kokkos::HostSpace> ();
    auto importLIDs_d = importLIDs.template view<DMS> ();
    auto importLIDs_h = importLIDs.template view<HMS> ();

    Kokkos::DualView<size_t*, device_type> whichVecs;
    if (! isConstantStride ()) {
      Kokkos::View<const size_t*, HMS,
        Kokkos::MemoryUnmanaged> whichVecsIn (whichVectors_.getRawPtr (),
                                              numVecs);
      whichVecs = Kokkos::DualView<size_t*, device_type> ("whichVecs", numVecs);
      if (unpackOnHost) {
        whichVecs.template modify<HMS> ();
        Kokkos::deep_copy (whichVecs.template view<HMS> (), whichVecsIn);
      }
      else {
        whichVecs.template modify<DMS> ();
        Kokkos::deep_copy (whichVecs.template view<DMS> (), whichVecsIn);
      }
    }
    auto whichVecs_d = whichVecs.template view<DMS> ();
    auto whichVecs_h = whichVecs.template view<HMS> ();

    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    if (printDebugOutput) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::unpackAndCombine: unpacking"
         << std::endl;
      std::cerr << os.str ();
    }

    if (numVecs > 0 && importLIDs.extent (0) > 0) {
      typedef typename Kokkos::DualView<IST*,
        device_type>::t_dev::execution_space dev_exec_space;
      typedef typename Kokkos::DualView<IST*,
        device_type>::t_host::execution_space host_exec_space;

      // NOTE (mfh 10 Mar 2012, 24 Mar 2014) If you want to implement
      // custom combine modes, start editing here.  Also, if you trust
      // inlining, it would be nice to condense this code by using a
      // binary function object f in the pack functors.
      if (CM == INSERT || CM == REPLACE) {
        KokkosRefactor::Details::InsertOp<execution_space> op;

        if (isConstantStride ()) {
          if (unpackOnHost) {
            // FIXME (mfh 29 Aug 2017) Problem: X_h's execution space
            // controls where unpacking takes place, but the
            // HostMirror of a CudaUVMSpace View has execution_space
            // Cuda!  This is bad if imports_h is actually a host
            // View.  This means that we need to control the execution
            // space in which unpack_array_multi_column works.
            unpack_array_multi_column (host_exec_space (),
                                       X_h, imports_h, importLIDs_h, op,
                                       numVecs, debugCheckIndices);

          }
          else { // unpack on device
            unpack_array_multi_column (dev_exec_space (),
                                       X_d, imports_d, importLIDs_d, op,
                                       numVecs, debugCheckIndices);
          }
        }
        else { // not constant stride
          if (unpackOnHost) {
            // FIXME (mfh 29 Aug 2017) Problem: X_h's execution space
            // controls where unpacking takes place, but the
            // HostMirror of a CudaUVMSpace View has execution_space
            // Cuda!  This is bad if imports_h is actually a host
            // View.  This means that we need to control the execution
            // space in which
            // unpack_array_multi_column_variable_stride works.
            unpack_array_multi_column_variable_stride (host_exec_space (),
                                                       X_h, imports_h,
                                                       importLIDs_h,
                                                       whichVecs_h, op,
                                                       numVecs,
                                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column_variable_stride (dev_exec_space (),
                                                       X_d, imports_d,
                                                       importLIDs_d,
                                                       whichVecs_d, op,
                                                       numVecs,
                                                       debugCheckIndices);
          }
        }
      }
      else if (CM == ADD) {
        KokkosRefactor::Details::AddOp<execution_space> op;

        if (isConstantStride ()) {
          if (unpackOnHost) {
            unpack_array_multi_column (host_exec_space (),
                                       X_h, imports_h, importLIDs_h, op,
                                       numVecs, debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column (dev_exec_space (),
                                       X_d, imports_d, importLIDs_d, op,
                                       numVecs, debugCheckIndices);
          }
        }
        else { // not constant stride
          if (unpackOnHost) {
            unpack_array_multi_column_variable_stride (host_exec_space (),
                                                       X_h, imports_h,
                                                       importLIDs_h,
                                                       whichVecs_h, op,
                                                       numVecs,
                                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column_variable_stride (dev_exec_space (),
                                                       X_d, imports_d,
                                                       importLIDs_d,
                                                       whichVecs_d, op,
                                                       numVecs,
                                                       debugCheckIndices);
          }
        }
      }
      else if (CM == ABSMAX) {
        KokkosRefactor::Details::AbsMaxOp<execution_space> op;

        if (isConstantStride ()) {
          if (unpackOnHost) {
            unpack_array_multi_column (host_exec_space (),
                                       X_h, imports_h, importLIDs_h, op,
                                       numVecs, debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column (dev_exec_space (),
                                       X_d, imports_d, importLIDs_d, op,
                                       numVecs, debugCheckIndices);
          }
        }
        else {
          if (unpackOnHost) {
            unpack_array_multi_column_variable_stride (host_exec_space (),
                                                       X_h, imports_h,
                                                       importLIDs_h,
                                                       whichVecs_h, op,
                                                       numVecs,
                                                       debugCheckIndices);
          }
          else { // unpack on device
            unpack_array_multi_column_variable_stride (dev_exec_space (),
                                                       X_d, imports_d,
                                                       importLIDs_d,
                                                       whichVecs_d, op,
                                                       numVecs,
                                                       debugCheckIndices);
          }
        }
      }
    }

    if (printDebugOutput) {
      std::ostringstream os;
      os << "(Proc " << myRank << ") MV::unpackAndCombine: Done!"
         << std::endl;
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
        const bool commIsInterComm = Details::isInterComm (*comm);

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
    using Kokkos::create_mirror_view;
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    // View of all the dot product results.
    typedef Kokkos::View<dot_type*, Kokkos::HostSpace> RV;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
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

    // If we need sync to device, then host has the most recent
    // version.  A is a guest of this method, so we should sync it.
    // Thus, let A control where execution happens.
    const bool useHostVersion = A.template need_sync<device_type> ();
    if (useHostVersion) {
      // A was last modified on host, so run the local kernel there.
      // This means we need a host mirror of the array of norms too.
      typedef typename dual_view_type::t_host XMV;
      typedef typename XMV::memory_space cur_memory_space;

      // I consider it more polite to sync *this, then to sync A.
      // A is a "guest" of this method, and is passed in const.
      const_cast<MV*> (this)->template sync<cur_memory_space> ();
      auto thisView = this->template getLocalView<cur_memory_space> ();
      auto A_view = A.template getLocalView<cur_memory_space> ();

      using Tpetra::Details::lclDot;
      lclDot<RV, XMV> (dotsOut, thisView, A_view, lclNumRows, numVecs,
                       this->whichVectors_.getRawPtr (),
                       A.whichVectors_.getRawPtr (),
                       this->isConstantStride (), A.isConstantStride ());
      gblDotImpl (dotsOut, comm, this->isDistributed ());
    }
    else {
      // A was last modified on device, so run the local kernel there.
      typedef typename dual_view_type::t_dev XMV;
      typedef typename XMV::memory_space cur_memory_space;

      // I consider it more polite to sync *this, then to sync A.
      // A is a "guest" of this method, and is passed in const.
      //
      // Yes, "const" is a lie.
      const_cast<MV*> (this)->template sync<cur_memory_space> ();
      auto thisView = this->template getLocalView<cur_memory_space> ();
      auto A_view = A.template getLocalView<cur_memory_space> ();

      using Tpetra::Details::lclDot;
      lclDot<RV, XMV> (dotsOut, thisView, A_view, lclNumRows, numVecs,
                       this->whichVectors_.getRawPtr (),
                       A.whichVectors_.getRawPtr (),
                       this->isConstantStride (), A.isConstantStride ());
      gblDotImpl (dotsOut, comm, this->isDistributed ());
    }
  }

  namespace { // (anonymous)
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type
    multiVectorSingleColumnDot (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x,
                                const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y)
    {
      using Teuchos::Comm;
      using Teuchos::RCP;
      typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
      typedef typename MV::dot_type dot_type;
      typedef typename MV::dual_view_type::t_dev::memory_space dev_memory_space;
      typedef typename MV::dual_view_type::t_host::memory_space host_memory_space;
      ::Tpetra::Details::ProfilingRegion region ("Tpetra::multiVectorSingleColumnDot");

      RCP<const typename MV::map_type> map = x.getMap ();
      RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
      if (comm.is_null ()) {
        return Kokkos::ArithTraits<dot_type>::zero ();
      }
      else {
        typedef LocalOrdinal LO;
        // The min just ensures that we don't overwrite memory that
        // doesn't belong to us, in the erroneous input case where x
        // and y have different numbers of rows.
        const LO lclNumRows = static_cast<LO> (std::min (x.getLocalLength (),
                                                         y.getLocalLength ()));
        const Kokkos::pair<LO, LO> rowRng (0, lclNumRows);
        dot_type lclDot = Kokkos::ArithTraits<dot_type>::zero ();
        dot_type gblDot = Kokkos::ArithTraits<dot_type>::zero ();

        // This function is meant to be called by x.  Therefore, we
        // can sync x, but we can't sync y.  y is a "guest" of this
        // method.  Thus, we let y control where execution happens.
        // If we need sync to device, then data were last modified on
        // host, so we should run on host.

        //const bool runOnHost = false;
        const bool runOnHost = y.template need_sync<dev_memory_space> ();
        // const_cast<MV&> (y).template sync<dev_memory_space> ();
        // const_cast<MV&> (x).template sync<dev_memory_space> ();

        if (runOnHost) {
          typedef host_memory_space cur_memory_space;
          // x is nonconst, so we may sync x where we need to sync it.
          x.template sync<cur_memory_space> ();
          auto x_2d = x.template getLocalView<cur_memory_space> ();
          auto x_1d = Kokkos::subview (x_2d, rowRng, 0);
          auto y_2d = y.template getLocalView<cur_memory_space> ();
          auto y_1d = Kokkos::subview (y_2d, rowRng, 0);
          lclDot = KokkosBlas::dot (x_1d, y_1d);
        }
        else { // run on device
          typedef dev_memory_space cur_memory_space;
          // x is nonconst, so we may sync x where we need to sync it.
          x.template sync<cur_memory_space> ();
          auto x_2d = x.template getLocalView<cur_memory_space> ();
          auto x_1d = Kokkos::subview (x_2d, rowRng, 0);
          auto y_2d = y.template getLocalView<cur_memory_space> ();
          auto y_1d = Kokkos::subview (y_2d, rowRng, 0);
          lclDot = KokkosBlas::dot (x_1d, y_1d);
        }

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
    typedef Kokkos::View<mag_type*, Kokkos::HostSpace> host_norms_view_type;

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::norm2 (Teuchos::ArrayView)");

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    this->norm2 (normsHostView); // Do the computation on the device.
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm2 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const
  {
    this->normImpl (norms, NORM_TWO);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void TPETRA_DEPRECATED
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normWeighted (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& weights,
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
    typedef Kokkos::View<mag_type*, Kokkos::HostSpace> norms_view_type;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "normWeighted: ";

    const size_t numVecs = this->getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (norms.size ()) != numVecs, std::runtime_error,
      "norms.size() = " << norms.size () << " != this->getNumVectors() = "
      << numVecs << ".");

    const bool OneW = (weights.getNumVectors () == 1);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! OneW && weights.getNumVectors () != numVecs, std::runtime_error,
      "The input MultiVector of weights must contain either one column, "
      "or must have the same number of columns as *this.  "
      "weights.getNumVectors() = " << weights.getNumVectors ()
      << " and this->getNumVectors() = " << numVecs << ".");

    const bool debug = ::Tpetra::Details::Behavior::debug ();

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! this->getMap ()->isCompatible (*weights.getMap ()), std::runtime_error,
         "MultiVectors do not have compatible Maps:" << std::endl
         << "this->getMap(): " << std::endl << *this->getMap()
         << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
    }
    else {
      const size_t lclNumRows = this->getLocalLength ();
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (lclNumRows != weights.getLocalLength (), std::runtime_error,
         "MultiVectors do not have the same local length.");
    }

    norms_view_type lclNrms ("Tpetra::MV::lclNrms", numVecs);

    // FIXME (mfh 18 May 2016) Yes, I know "const" is a lie.
    const_cast<MV*> (this)->template sync<device_type> ();
    const_cast<MV&> (weights).template sync<device_type> ();

    auto X_lcl = this->template getLocalView<device_type> ();
    auto W_lcl = weights.template getLocalView<device_type> ();

    if (isConstantStride () && ! OneW) {
      KokkosBlas::nrm2w_squared (lclNrms, X_lcl, W_lcl);
    }
    else {
      for (size_t j = 0; j < numVecs; ++j) {
        const size_t X_col = this->isConstantStride () ? j :
          this->whichVectors_[j];
        const size_t W_col = OneW ? static_cast<size_t> (0) :
          (weights.isConstantStride () ? j : weights.whichVectors_[j]);
        KokkosBlas::nrm2w_squared (subview (lclNrms, j),
                                   subview (X_lcl, ALL (), X_col),
                                   subview (W_lcl, ALL (), W_col));
      }
    }

    const mag_type OneOverN =
      ATM::one () / static_cast<mag_type> (this->getGlobalLength ());
    RCP<const Comm<int> > comm = this->getMap ().is_null () ?
      Teuchos::null : this->getMap ()->getComm ();

    if (! comm.is_null () && this->isDistributed ()) {
      // Assume that MPI can access device memory.
      reduceAll<int, mag_type> (*comm, REDUCE_SUM, static_cast<int> (numVecs),
                                lclNrms.data (), norms.getRawPtr ());
      for (size_t k = 0; k < numVecs; ++k) {
        norms[k] = ATM::sqrt (norms[k] * OneOverN);
      }
    }
    else {
      for (size_t k = 0; k < numVecs; ++k) {
        norms[k] = ATM::sqrt (ATS::magnitude (lclNrms(k)) * OneOverN);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm1 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef typename Kokkos::View<mag_type*, Kokkos::HostSpace> host_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    this->norm1 (normsHostView); // Do the computation on the device.
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  norm1 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const
  {
    this->normImpl (norms, NORM_ONE);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normInf (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef Kokkos::View<mag_type*, Kokkos::HostSpace> host_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    this->normInf (normsHostView); // Do the computation on the device.
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normInf (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const
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

      // In case the input dimensions don't match, make sure that we
      // don't overwrite memory that doesn't belong to us, by using
      // subset views with the minimum dimensions over all input.
      const std::pair<size_t, size_t> rowRng (0, lclNumRows);
      const std::pair<size_t, size_t> colRng (0, numVecs);
      RV theNorms = subview (normsOut, colRng);
      XMV X = subview (X_lcl, rowRng, Kokkos::ALL());

      // mfh 10 Mar 2015: Kokkos::(Dual)View subviews don't quite
      // behave how you think when they have zero rows.  In that case,
      // it returns a 0 x 0 (Dual)View.
      TEUCHOS_TEST_FOR_EXCEPTION(
        lclNumRows != 0 && constantStride && ( \
          ( X.extent (0) != lclNumRows ) ||
          ( X.extent (1) != numVecs    ) ),
        std::logic_error, "Constant Stride X's dimensions are " << X.extent (0) << " x "
        << X.extent (1) << ", which differ from the local dimensions "
        << lclNumRows << " x " << numVecs << ".  Please report this bug to "
        "the Tpetra developers.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        lclNumRows != 0 && !constantStride && ( \
          ( X.extent (0) != lclNumRows ) ||
          ( X.extent (1) < numVecs    ) ),
        std::logic_error, "Strided X's dimensions are " << X.extent (0) << " x "
        << X.extent (1) << ", which are incompatible with the local dimensions "
        << lclNumRows << " x " << numVecs << ".  Please report this bug to "
        "the Tpetra developers.");

      if (lclNumRows == 0) {
        const mag_type zeroMag = Kokkos::Details::ArithTraits<mag_type>::zero ();
        Kokkos::deep_copy(theNorms, zeroMag);
      }
      else { // lclNumRows != 0
        if (constantStride) {
          if (whichNorm == IMPL_NORM_INF) {
            KokkosBlas::nrminf (theNorms, X);
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
              KokkosBlas::nrminf (subview (theNorms, k),
                                  subview (X, ALL (), X_col));
            }
            else if (whichNorm == IMPL_NORM_ONE) {
              KokkosBlas::nrm1 (subview (theNorms, k),
                                subview (X, ALL (), X_col));
            }
            else if (whichNorm == IMPL_NORM_TWO) {
              KokkosBlas::nrm2_squared (subview (theNorms, k),
                                        subview (X, ALL (), X_col));
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
            }
          } // for each column
        } // constantStride
      } // lclNumRows != 0
    }

    // Kokkos::parallel_for functor for applying square root to each
    // entry of a 1-D Kokkos::View.
    template<class ViewType>
    class SquareRootFunctor {
    public:
      typedef typename ViewType::execution_space execution_space;
      typedef typename ViewType::size_type size_type;

      SquareRootFunctor (const ViewType& theView) :
        theView_ (theView)
      {}

      KOKKOS_INLINE_FUNCTION void
      operator() (const size_type& i) const
      {
        typedef typename ViewType::non_const_value_type value_type;
        typedef Kokkos::Details::ArithTraits<value_type> KAT;
        theView_(i) = KAT::sqrt (theView_(i));
      }

    private:
      ViewType theView_;
    };

    template<class RV>
    void
    gblNormImpl (const RV& normsOut,
                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const bool distributed,
                 const EWhichNormImpl whichNorm)
    {
      using Teuchos::REDUCE_MAX;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      typedef typename RV::non_const_value_type mag_type;

      const size_t numVecs = normsOut.extent (0);

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
      // norm*() on those processes anyway.
      if (distributed && ! comm.is_null ()) {
        // The calling process only participates in the collective if
        // both the Map and its Comm on that process are nonnull.
        //
        // MPI doesn't allow aliasing of arguments, so we have to make
        // a copy of the local sum.
        RV lclNorms ("MV::normImpl lcl", numVecs);
        Kokkos::deep_copy (lclNorms, normsOut);
        const mag_type* const lclSum = lclNorms.data ();
        mag_type* const gblSum = normsOut.data ();
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
          Kokkos::Impl::is_same<typename RV::memory_space,
            typename RV::host_mirror_space::memory_space>::value;
        if (inHostMemory) {
          for (size_t j = 0; j < numVecs; ++j) {
            normsOut(j) = Kokkos::Details::ArithTraits<mag_type>::sqrt (normsOut(j));
          }
        }
        else {
          // There's not as much parallelism now, but that's OK.  The
          // point of doing parallel dispatch here is to keep the norm
          // results on the device, thus avoiding a copy to the host
          // and back again.
          SquareRootFunctor<RV> f (normsOut);
          typedef typename RV::execution_space execution_space;
          typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
          Kokkos::parallel_for (range_type (0, numVecs), f);
        }
      }
    }

  } // namespace (anonymous)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  normImpl (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms,
            const EWhichNorm whichNorm) const
  {
    using Kokkos::create_mirror_view;
    using Kokkos::subview;
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::RCP;
    // View of all the norm results.
    typedef Kokkos::View<mag_type*, Kokkos::HostSpace> RV;

    const size_t numVecs = this->getNumVectors ();
    if (numVecs == 0) {
      return; // nothing to do
    }
    const size_t lclNumRows = this->getLocalLength ();
    const size_t numNorms = static_cast<size_t> (norms.extent (0));
    TEUCHOS_TEST_FOR_EXCEPTION(
      numNorms < numVecs, std::runtime_error, "Tpetra::MultiVector::normImpl: "
      "'norms' must have at least as many entries as the number of vectors in "
      "*this.  norms.extent(0) = " << numNorms << " < this->getNumVectors()"
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

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->template need_sync<device_type> ();
    if (useHostVersion) {
      // DualView was last modified on host, so run the local kernel there.
      // This means we need a host mirror of the array of norms too.
      typedef typename dual_view_type::t_host XMV;
      typedef typename XMV::memory_space cur_memory_space;

      auto thisView = this->template getLocalView<cur_memory_space> ();
      lclNormImpl<RV, XMV> (normsOut, thisView, lclNumRows, numVecs,
                            this->whichVectors_, this->isConstantStride (),
                            lclNormType);
      gblNormImpl (normsOut, comm, this->isDistributed (), lclNormType);
    }
    else {
      // DualView was last modified on device, so run the local kernel there.
      typedef typename dual_view_type::t_dev XMV;
      typedef typename XMV::memory_space cur_memory_space;

      auto thisView = this->template getLocalView<cur_memory_space> ();
      lclNormImpl<RV, XMV> (normsOut, thisView, lclNumRows, numVecs,
                            this->whichVectors_, this->isConstantStride (),
                            lclNormType);
      gblNormImpl (normsOut, comm, this->isDistributed (), lclNormType);
    }
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
    const bool useHostVersion = this->template need_sync<device_type> ();
    if (useHostVersion) {
      // DualView was last modified on host, so run the local kernel there.
      auto X_lcl = subview (this->template getLocalView<Kokkos::HostSpace> (),
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
      auto X_lcl = subview (this->template getLocalView<device_type> (),
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
    this->view_.modified_device () = 0;
    this->view_.modified_host () = 0;

    this->template modify<device_type> ();
    auto thisView = this->getLocalView<device_type> ();

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
    typedef typename dual_view_type::t_dev::memory_space DMS;
    //typedef typename dual_view_type::t_host::memory_space HMS;
    typedef Kokkos::HostSpace HMS; // avoid CudaUVMSpace issues
    typedef typename dual_view_type::t_dev::execution_space DES;
    typedef typename dual_view_type::t_host::execution_space HES;
    typedef LocalOrdinal LO;
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
    const bool useHostVersion = this->template need_sync<DMS> ();

    if (! useHostVersion) { // last modified in device memory
      this->template modify<DMS> (); // we are about to modify on the device
      auto X = this->template getLocalView<DMS> ();
      if (this->isConstantStride ()) {
        fill (DES (), X, theAlpha, lclNumRows, numVecs);
      }
      else {
        fill (DES (), X, theAlpha, lclNumRows, numVecs,
              this->whichVectors_.getRawPtr ());
      }
    }
    else { // last modified in host memory, so modify data there.
      this->template modify<HMS> (); // we are about to modify on the host
      auto X = this->template getLocalView<HMS> ();
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
        view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, Node> (newNumRows, numCols);
      }
    }
    else if (newMap.is_null ()) { // Case 2: current Map is nonnull, new Map is null
      // I am an excluded process.  Reinitialize my data so that I
      // have 0 rows.  Keep the number of columns as before.
      const size_t newNumRows = static_cast<size_t> (0);
      const size_t numCols = this->getNumVectors ();
      view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, Node> (newNumRows, numCols);
    }

    this->map_ = newMap;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  scale (const Scalar& alpha)
  {
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    if (theAlpha == Kokkos::Details::ArithTraits<impl_scalar_type>::one ()) {
      return; // do nothing
    }
    const size_t lclNumRows = this->getLocalLength ();
    const size_t numVecs = this->getNumVectors ();
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // We can't substitute putScalar(0.0) for scale(0.0), because the
    // former will overwrite NaNs present in the MultiVector.  The
    // semantics of this call require multiplying them by 0, which
    // IEEE 754 requires to be NaN.

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->template need_sync<device_type> ();
    if (useHostVersion) {
      auto Y_lcl =
        Kokkos::subview (this->template getLocalView<Kokkos::HostSpace> (),
                         rowRng, Kokkos::ALL ());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, theAlpha, Y_lcl);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = this->isConstantStride () ? k :
            this->whichVectors_[k];
          auto Y_k = Kokkos::subview (Y_lcl, Kokkos::ALL (), Y_col);
          KokkosBlas::scal (Y_k, theAlpha, Y_k);
        }
      }
    }
    else { // work on device
      auto Y_lcl =
        Kokkos::subview (this->template getLocalView<device_type> (),
                         rowRng, Kokkos::ALL ());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, theAlpha, Y_lcl);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = this->isConstantStride () ? k :
            this->whichVectors_[k];
          auto Y_k = Kokkos::subview (Y_lcl, Kokkos::ALL (), Y_col);
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
    typedef Kokkos::DualView<impl_scalar_type*, device_type> k_alphas_type ;
    k_alphas_type k_alphas ("alphas::tmp", numAlphas);
    k_alphas.template modify<typename k_alphas_type::host_mirror_space> ();
    for (size_t i = 0; i < numAlphas; ++i) {
      k_alphas.h_view(i) = static_cast<impl_scalar_type> (alphas[i]);
    }
    k_alphas.template sync<typename k_alphas_type::memory_space> ();
    // Invoke the scale() overload that takes a device View of coefficients.
    this->scale (k_alphas.d_view);
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
    const bool useHostVersion = this->template need_sync<device_type> ();
    if (useHostVersion) {
      // Work in host memory.  This means we need to create a host
      // mirror of the input View of coefficients.
      auto alphas_h = Kokkos::create_mirror_view (alphas);
      typedef typename decltype (alphas_h)::memory_space cur_memory_space;
      Kokkos::deep_copy (alphas_h, alphas);

      auto Y_lcl =
        Kokkos::subview (this->template getLocalView<cur_memory_space> (),
                         rowRng, Kokkos::ALL ());
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
      auto Y_lcl =
        Kokkos::subview (this->template getLocalView<device_type> (),
                         rowRng, Kokkos::ALL());
      if (isConstantStride ()) {
        KokkosBlas::scal (Y_lcl, alphas, Y_lcl);
      }
      else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t Y_col = this->isConstantStride () ? k :
            this->whichVectors_[k];
          auto Y_k = subview (Y_lcl, ALL (), Y_col);
          //
          // FIXME (mfh 08 Apr 2015) This assumes UVM.  It would be
          // better to fix scal() so that it takes a 0-D View as the
          // second argument.
          //
          KokkosBlas::scal (Y_k, alphas(k), Y_k);
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

    typedef typename dual_view_type::t_dev dev_view_type;
    typedef typename dual_view_type::t_host host_view_type;

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->template need_sync<device_type> ();
    if (useHostVersion) {
      // Work on host, where A's data were most recently modified.  A
      // is a "guest" of this method, so it's more polite to sync
      // *this, than to sync A.
      typedef typename host_view_type::memory_space cur_memory_space;

      this->template sync<cur_memory_space> ();
      this->template modify<cur_memory_space> ();
      auto Y_lcl_orig = this->template getLocalView<cur_memory_space> ();
      auto X_lcl_orig = A.template getLocalView<cur_memory_space> ();
      auto Y_lcl = subview (Y_lcl_orig, rowRng, Kokkos::ALL ());
      auto X_lcl = subview (X_lcl_orig, rowRng, Kokkos::ALL ());

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
    else { // work on device
      typedef typename dev_view_type::memory_space cur_memory_space;

      this->template sync<cur_memory_space> ();
      this->template modify<cur_memory_space> ();
      auto Y_lcl_orig = this->template getLocalView<cur_memory_space> ();
      auto X_lcl_orig = A.template getLocalView<cur_memory_space> ();
      auto Y_lcl = subview (Y_lcl_orig, rowRng, Kokkos::ALL ());
      auto X_lcl = subview (X_lcl_orig, rowRng, Kokkos::ALL ());

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
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  reciprocal (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
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

    // FIXME (mfh 07 Jan 2015) See note on two-argument scale() above.

    const size_t numVecs = getNumVectors ();
    try {
      this->template sync<device_type> ();
      this->template modify<device_type> ();
      // FIXME (mfh 23 Jul 2014) I'm not sure if it should be our
      // responsibility to sync A.
      //
      // FIXME (mfh 18 May 2016) It's rude to sync the input
      // argument, since it is marked const.
      const_cast<MV&> (A).template sync<device_type> ();

      auto this_view_dev = this->template getLocalView<device_type> ();
      auto A_view_dev = A.template getLocalView<device_type> ();

      if (isConstantStride () && A.isConstantStride ()) {
        KokkosBlas::reciprocal (this_view_dev, A_view_dev);
      }
      else {
        using Kokkos::ALL;
        using Kokkos::subview;
        typedef Kokkos::View<impl_scalar_type*, device_type> view_type;

        for (size_t k = 0; k < numVecs; ++k) {
          const size_t this_col = isConstantStride () ? k : whichVectors_[k];
          view_type vector_k = subview (this_view_dev, ALL (), this_col);
          const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
          view_type vector_Ak = subview (A_view_dev, ALL (), A_col);
          KokkosBlas::reciprocal (vector_k, vector_Ak);
        }
      }
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error,
        ": Caught exception from Kokkos: " << e.what () << std::endl);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  abs (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
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

    this->template sync<device_type> ();
    this->template modify<device_type> ();
    // FIXME (mfh 23 Jul 2014) I'm not sure if it should be our
    // responsibility to sync A.
    //
    // FIXME (mfh 18 May 2016) It's rude to sync the input
    // argument, since it is marked const.
    const_cast<MV&> (A).template sync<device_type> ();

    auto this_view_dev = this->template getLocalView<device_type> ();
    auto A_view_dev = A.template getLocalView<device_type> ();

    if (isConstantStride () && A.isConstantStride ()) {
      KokkosBlas::abs (this_view_dev, A_view_dev);
    }
    else {
      using Kokkos::ALL;
      using Kokkos::subview;
      typedef Kokkos::View<impl_scalar_type*, device_type> view_type;

      for (size_t k=0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = subview (this_view_dev, ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        view_type vector_Ak = subview (A_view_dev, ALL (), A_col);
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
    using Kokkos::ALL;
    using Kokkos::subview;
    const char tfecfFuncName[] = "update: ";

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

    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    const impl_scalar_type theBeta = static_cast<impl_scalar_type> (beta);
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    typedef typename dual_view_type::t_dev dev_view_type;
    typedef typename dual_view_type::t_host host_view_type;

    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->template need_sync<device_type> ();
    if (useHostVersion) {
      // Work on host, where A's data were most recently modified.  A
      // is a "guest" of this method, so it's more polite to sync
      // *this, than to sync A.
      typedef typename host_view_type::memory_space cur_memory_space;

      this->template sync<cur_memory_space> ();
      this->template modify<cur_memory_space> ();
      auto Y_lcl_orig = this->template getLocalView<cur_memory_space> ();
      auto Y_lcl = subview (Y_lcl_orig, rowRng, Kokkos::ALL ());
      auto X_lcl_orig = A.template getLocalView<cur_memory_space> ();
      auto X_lcl = subview (X_lcl_orig, rowRng, Kokkos::ALL ());

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
    else { // work on device
      // A is a "guest" of this method, so it's more polite to sync
      // *this, than to sync A.
      typedef typename dev_view_type::memory_space cur_memory_space;
      this->template sync<cur_memory_space> ();
      this->template modify<cur_memory_space> ();
      auto Y_lcl_orig = this->template getLocalView<cur_memory_space> ();
      auto Y_lcl = subview (Y_lcl_orig, rowRng, Kokkos::ALL ());
      auto X_lcl_orig = A.template getLocalView<cur_memory_space> ();
      auto X_lcl = subview (X_lcl_orig, rowRng, Kokkos::ALL ());

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
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
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

    // We're lucky if *this, A, and B are all sync'd to the same
    // memory space.  If not, we have to sync _something_.  Unlike
    // three-argument update() or (say) dot(), we may have to sync one
    // of the inputs.  For now, we just sync _everything_ to device.
    this->template sync<device_type> ();
    const_cast<MV&> (A).template sync<device_type> ();
    const_cast<MV&> (B).template sync<device_type> ();

    // This method modifies *this.
    this->template modify<device_type> ();

    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numVecs);

    // Prefer 'auto' over specifying the type explicitly.  This avoids
    // issues with a subview possibly having a different type than the
    // original view.
    auto C_lcl = subview (this->template getLocalView<device_type> (), rowRng, ALL ());
    auto A_lcl = subview (A.template getLocalView<device_type> (), rowRng, ALL ());
    auto B_lcl = subview (B.template getLocalView<device_type> (), rowRng, ALL ());

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
    using Kokkos::subview;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "getData: ";

    // Any MultiVector method that called the (classic) Kokkos Node's
    // viewBuffer or viewBufferNonConst methods always implied a
    // device->host synchronization.  Thus, we synchronize here as
    // well.
    const_cast<MV*> (this)->template sync<Kokkos::HostSpace> ();

    // Get a host view of the entire MultiVector's data.
    auto hostView = this->template getLocalView<Kokkos::HostSpace> ();

    // Get a subview of column j.
    const size_t col = isConstantStride () ? j : whichVectors_[j];
    auto hostView_j = subview (hostView, ALL (), col);

    // Wrap up the subview of column j in an ArrayRCP<const impl_scalar_type>.
    Teuchos::ArrayRCP<const impl_scalar_type> dataAsArcp =
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
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "getDataNonConst: ";

    // Any MultiVector method that called the (classic) Kokkos Node's
    // viewBuffer or viewBufferNonConst methods always implied a
    // device->host synchronization.  Thus, we synchronize here as
    // well.
    const_cast<MV*> (this)->template sync<Kokkos::HostSpace> ();

    // Calling getDataNonConst() implies that the user plans to modify
    // the values in the MultiVector, so we mark the host data as modified.
    this->template modify<Kokkos::HostSpace> ();

    // Get a host view of the entire MultiVector's data.
    auto hostView = this->template getLocalView<Kokkos::HostSpace> ();

    // Get a subview of column j.
    const size_t col = isConstantStride () ? j : whichVectors_[j];
    auto hostView_j = subview (hostView, ALL (), col);

    // Wrap up the subview of column j in an ArrayRCP<const impl_scalar_type>.
    Teuchos::ArrayRCP<impl_scalar_type> dataAsArcp =
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
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  operator= (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source)
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
               const map_type& subMap,
               const size_t offset) :
    base_type (Teuchos::null) // to be replaced below
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LO, GO, Node> MV;
    const char prefix[] = "Tpetra::MultiVector constructor (offsetView): ";

    const size_t newNumRows = subMap.getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > X.getOrigNumLocalRows ();
    if (tooManyElts) {
      const int myRank = X.getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > X.getLocalLength (), std::runtime_error,
        prefix << "Invalid input Map.  The input Map owns " << newNumRows <<
        " entries on process " << myRank << ".  offset = " << offset << ".  "
        "Yet, the MultiVector contains only " << X.getOrigNumLocalRows () <<
        " rows on this process.");
    }

    // mfh 27 Sep 2017: This debugging code depends on being able to
    // declare variables that we can use below, so it's too hard for
    // now to use run-time control to enable this code.
#ifdef HAVE_TPETRA_DEBUG
    const size_t strideBefore =
      X.isConstantStride () ? X.getStride () : static_cast<size_t> (0);
    const size_t lclNumRowsBefore = X.getLocalLength ();
    const size_t numColsBefore = X.getNumVectors ();
    const impl_scalar_type* hostPtrBefore =
      X.template getLocalView<Kokkos::HostSpace> ().data ();
#endif // HAVE_TPETRA_DEBUG

    const std::pair<size_t, size_t> origRowRng (offset, X.origView_.extent (0));
    const std::pair<size_t, size_t> rowRng (offset, offset + newNumRows);

    dual_view_type newOrigView = subview (X.origView_, origRowRng, ALL ());
    // FIXME (mfh 29 Sep 2016) If we just use X.view_ here, it breaks
    // CrsMatrix's Gauss-Seidel implementation (which assumes the
    // ability to create domain Map views of column Map MultiVectors,
    // and then get the original column Map MultiVector out again).
    // If we just use X.origView_ here, it breaks the fix for #46.
    // The test for offset == 0 is a hack that makes both tests pass,
    // but doesn't actually fix the more general issue.  In
    // particular, the right way to fix Gauss-Seidel would be to fix
    // #385; that would make "getting the original column Map
    // MultiVector out again" unnecessary.
    dual_view_type newView = subview (offset == 0 ? X.origView_ : X.view_, rowRng, ALL ());

    // NOTE (mfh 06 Jan 2015) Work-around to deal with Kokkos not
    // handling subviews of degenerate Views quite so well.  For some
    // reason, the ([0,0], [0,2]) subview of a 0 x 2 DualView is 0 x
    // 0.  We work around by creating a new empty DualView of the
    // desired (degenerate) dimensions.
    if (newOrigView.extent (0) == 0 &&
        newOrigView.extent (1) != X.origView_.extent (1)) {
      newOrigView = allocDualView<Scalar, LO, GO, Node> (size_t (0),
                                                         X.getNumVectors ());
    }
    if (newView.extent (0) == 0 &&
        newView.extent (1) != X.view_.extent (1)) {
      newView = allocDualView<Scalar, LO, GO, Node> (size_t (0),
                                                     X.getNumVectors ());
    }

    MV subViewMV = X.isConstantStride () ?
      MV (Teuchos::rcp (new map_type (subMap)), newView, newOrigView) :
      MV (Teuchos::rcp (new map_type (subMap)), newView, newOrigView, X.whichVectors_ ());

#ifdef HAVE_TPETRA_DEBUG
    const size_t strideAfter = X.isConstantStride () ?
      X.getStride () :
      static_cast<size_t> (0);
    const size_t lclNumRowsAfter = X.getLocalLength ();
    const size_t numColsAfter = X.getNumVectors ();
    const impl_scalar_type* hostPtrAfter =
      X.template getLocalView<Kokkos::HostSpace> ().data ();

    const size_t strideRet = subViewMV.isConstantStride () ?
      subViewMV.getStride () :
      static_cast<size_t> (0);
    const size_t lclNumRowsRet = subViewMV.getLocalLength ();
    const size_t numColsRet = subViewMV.getNumVectors ();

    const char suffix[] = ".  This should never happen.  Please report this "
      "bug to the Tpetra developers.";

    TEUCHOS_TEST_FOR_EXCEPTION(
      lclNumRowsRet != subMap.getNodeNumElements (),
      std::logic_error, prefix << "Returned MultiVector has a number of rows "
      "different than the number of local indices in the input Map.  "
      "lclNumRowsRet: " << lclNumRowsRet << ", subMap.getNodeNumElements(): "
      << subMap.getNodeNumElements () << suffix);
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

    *this = subViewMV; // shallow copy
  }

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
      auto newImports = X.imports_;
      newImports.d_view = subview (X.imports_.d_view, range_type (0, newSize));
      newImports.h_view = subview (X.imports_.h_view, range_type (0, newSize));
    }
    {
      const size_t newSize = X.exports_.extent (0) / numCols;
      auto newExports = X.exports_;
      newExports.d_view = subview (X.exports_.d_view, range_type (0, newSize));
      newExports.h_view = subview (X.exports_.h_view, range_type (0, newSize));
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
    typedef decltype (this->template getLocalView<device_type> ())
      dev_view_type;
    typedef decltype (this->template getLocalView<Kokkos::HostSpace> ())
      host_view_type;
    typedef impl_scalar_type IST;
    typedef Kokkos::View<IST*, typename host_view_type::array_layout,
      Kokkos::HostSpace, Kokkos::MemoryUnmanaged> input_col_type;
    const char tfecfFuncName[] = "get1dCopy: ";

    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const std::pair<size_t, size_t> rowRange (0, numRows);

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

    // Use the most recently updated version of this MultiVector's
    // data.  This avoids sync'ing, which could violate users'
    // expectations.
    //
    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->template need_sync<device_type> ();

    dev_view_type srcView_dev;
    host_view_type srcView_host;
    if (useHostVersion) {
      srcView_host = this->template getLocalView<Kokkos::HostSpace> ();
    }
    else {
      srcView_dev = this->template getLocalView<device_type> ();
    }

    for (size_t j = 0; j < numCols; ++j) {
      const size_t srcCol =
        this->isConstantStride () ? j : this->whichVectors_[j];
      const size_t dstCol = j;
      IST* const dstColRaw =
        reinterpret_cast<IST*> (A.getRawPtr () + LDA * dstCol);
      input_col_type dstColView (dstColRaw, numRows);

      if (useHostVersion) {
        auto srcColView_host = Kokkos::subview (srcView_host, rowRange, srcCol);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          dstColView.extent (0) != srcColView_host.extent (0),
          std::logic_error, ": srcColView and dstColView_host have different "
          "dimensions.  Please report this bug to the Tpetra developers.");
        Kokkos::deep_copy (dstColView, srcColView_host);
      }
      else {
        auto srcColView_dev = Kokkos::subview (srcView_dev, rowRange, srcCol);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          dstColView.extent (0) != srcColView_dev.extent (0),
          std::logic_error, ": srcColView and dstColView_dev have different "
          "dimensions.  Please report this bug to the Tpetra developers.");
        Kokkos::deep_copy (dstColView, srcColView_dev);
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  get1dView () const
  {
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      return Teuchos::null;
    } else {
      typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
        "get1dView: This MultiVector does not have constant stride, so it is "
        "not possible to view its data as a single array.  You may check "
        "whether a MultiVector has constant stride by calling "
        "isConstantStride().");
      // mfh 18 May 2016: get?dView() and get?dViewNonConst() (replace
      // ? with 1 or 2) have always been device->host synchronization
      // points, since <= 2012.  We retain this behavior for backwards
      // compatibility.
      //
      // Yes, "const" is a lie.
      const_cast<MV*> (this)->template sync<Kokkos::HostSpace> ();
      // Both get1dView() and get1dViewNonConst() return a host view
      // of the data.
      auto X_lcl = this->template getLocalView<Kokkos::HostSpace> ();
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
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      return Teuchos::null;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
         "get1dViewNonConst: This MultiVector does not have constant stride, "
         "so it is not possible to view its data as a single array.  You may "
         "check whether a MultiVector has constant stride by calling "
         "isConstantStride().");
      // get1dView() and get1dViewNonConst() have always been
      // device->host synchronization points, since <= 2012.  We
      // retain this behavior for backwards compatibility.
      this->template sync<Kokkos::HostSpace> ();
      // Both get1dView() and get1dViewNonConst() return a host view
      // of the data.
      auto X_lcl = this->template getLocalView<Kokkos::HostSpace> ();
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
    // NOTE (mfh 16 May 2016) get?dView() and get?dViewNonConst()
    // (replace ? with 1 or 2) have always been device->host
    // synchronization points, since <= 2012.  We retain this behavior
    // for backwards compatibility.
    this->sync<Kokkos::HostSpace> ();
    // When users call the NonConst variants, it implies that they
    // want to change the data.  Thus, it is appropriate to mark this
    // MultiVector as modified.
    this->modify<Kokkos::HostSpace> ();

    const size_t myNumRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const Kokkos::pair<size_t, size_t> rowRange (0, myNumRows);
    // Don't use the row range here on the outside, in order to avoid
    // a strided return type (in case Kokkos::subview is conservative
    // about that).  Instead, use the row range for the column views
    // in the loop.
    auto X_lcl = this->getLocalView<Kokkos::HostSpace> ();

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
    // NOTE (mfh 16 May 2016) get?dView() and get?dViewNonConst()
    // (replace ? with 1 or 2) have always been device->host
    // synchronization points, since <= 2012.  We retain this behavior
    // for backwards compatibility.
    //
    // Since get2dView() is and was (unfortunately) always marked
    // const, I have to cast away const here in order not to break
    // backwards compatibility.
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
    const_cast<this_type*> (this)->sync<Kokkos::HostSpace> ();

    const size_t myNumRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const Kokkos::pair<size_t, size_t> rowRange (0, myNumRows);
    // Don't use the row range here on the outside, in order to avoid
    // a strided return type (in case Kokkos::subview is conservative
    // about that).  Instead, use the row range for the column views
    // in the loop.
    auto X_lcl = this->getLocalView<Kokkos::HostSpace> ();

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
    using Teuchos::CONJ_TRANS;
    using Teuchos::NO_TRANS;
    using Teuchos::TRANS;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    typedef Kokkos::Details::ArithTraits<impl_scalar_type> ATS;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    const char tfecfFuncName[] = "multiply: ";
    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::multiply");

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

    impl_scalar_type beta_local = beta; // local copy of beta; might be reassigned below

    // In debug mode, check compatibility of local dimensions.  We
    // only do this in debug mode, since it requires an all-reduce
    // to ensure correctness on all processses.  It's entirely
    // possible that only some processes may have incompatible local
    // dimensions.  Throwing an exception only on those processes
    // could cause this method to hang.
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if (debug) {
      if (! this->getMap ().is_null () && ! this->getMap ()->getComm ().is_null ()) {
        using Teuchos::Comm;
        using Teuchos::RCP;
        using Teuchos::REDUCE_MIN;
        using Teuchos::reduceAll;
        using Teuchos::outArg;

        RCP<const Comm<int> > comm = this->getMap ()->getComm ();
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
        const int lclGood = lclBad ? 0 : 1;
        int gblGood = 0;
        reduceAll<int, int> (*comm, REDUCE_MIN, lclGood, outArg (gblGood));

        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (gblGood != 1, std::runtime_error, "Local dimensions of '*this', "
           "op(A), and op(B) are not consistent on at least one process "
           "in this object's communicator.");
      }
    }

    const bool A_is_local = ! A.isDistributed ();
    const bool B_is_local = ! B.isDistributed ();
    const bool C_is_local = ! this->isDistributed ();
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
      typedef LocalOrdinal LO;

      const LO A_lclNumRows = A_tmp->getLocalLength ();
      const LO A_numVecs = A_tmp->getNumVectors ();
      auto A_lcl = A_tmp->template getLocalView<device_type> ();
      auto A_sub = Kokkos::subview (A_lcl,
                                    std::make_pair (LO (0), A_lclNumRows),
                                    std::make_pair (LO (0), A_numVecs));
      const LO B_lclNumRows = B_tmp->getLocalLength ();
      const LO B_numVecs = B_tmp->getNumVectors ();
      auto B_lcl = B_tmp->template getLocalView<device_type> ();
      auto B_sub = Kokkos::subview (B_lcl,
                                    std::make_pair (LO (0), B_lclNumRows),
                                    std::make_pair (LO (0), B_numVecs));
      const LO C_lclNumRows = C_tmp->getLocalLength ();
      const LO C_numVecs = C_tmp->getNumVectors ();
      auto C_lcl = C_tmp->template getLocalView<device_type> ();
      auto C_sub = Kokkos::subview (C_lcl,
                                    std::make_pair (LO (0), C_lclNumRows),
                                    std::make_pair (LO (0), C_numVecs));
      const char ctransA = (transA == Teuchos::NO_TRANS ? 'N' :
                            (transA == Teuchos::TRANS ? 'T' : 'C'));
      const char ctransB = (transB == Teuchos::NO_TRANS ? 'N' :
                            (transB == Teuchos::TRANS ? 'T' : 'C'));
      const impl_scalar_type alpha_IST (alpha);

      ::Tpetra::Details::ProfilingRegion regionGemm ("Tpetra::MV::multiply-call-gemm");
      ::Tpetra::Details::Blas::gemm (ctransA, ctransB, alpha_IST, A_sub,
                                     B_sub, beta_local, C_sub);
    }

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
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> V;
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

    // FIXME (mfh 18 May 2016) It would be rude to sync A and B here,
    // because they are guests of this method.  Instead, the polite
    // thing to do would be to copy them (if necessary) so we get
    // their most recently updated version.  *this should always
    // control where execution happens.

    this->template sync<device_type> ();
    this->template modify<device_type> ();
    const_cast<V&> (A).template sync<device_type> ();
    const_cast<MV&> (B).template sync<device_type> ();
    auto this_view = this->template getLocalView<device_type> ();
    auto A_view = A.template getLocalView<device_type> ();
    auto B_view = B.template getLocalView<device_type> ();

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
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef decltype (this->template getLocalView<device_type> ())
      dev_view_type;
    typedef decltype (this->template getLocalView<Kokkos::HostSpace> ())
      host_view_type;

    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV::reduce");

    TEUCHOS_TEST_FOR_EXCEPTION
      (this->isDistributed (), std::runtime_error,
      "Tpetra::MultiVector::reduce should only be called with locally "
      "replicated or otherwise not distributed MultiVector objects.");
    const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
    if (comm.getSize () == 1) {
      return; // MultiVector is already "reduced" to one process
    }

    const size_t lclNumRows = getLocalLength ();
    const size_t numCols = getNumVectors ();
    const size_t totalAllocSize = lclNumRows * numCols;

    // FIXME (mfh 16 June 2014) This exception will cause deadlock if
    // it triggers on only some processes.  We don't have a good way
    // to pack this result into the all-reduce below, but this would
    // be a good reason to set a "local error flag" and find other
    // opportunities to let it propagate.
    TEUCHOS_TEST_FOR_EXCEPTION(
      lclNumRows > static_cast<size_t> (std::numeric_limits<int>::max ()),
      std::runtime_error, "Tpetra::MultiVector::reduce: On Process " <<
      comm.getRank () << ", the number of local rows " << lclNumRows <<
      " does not fit in int.");

    //
    // Use MPI to sum the entries across all local blocks.
    //

    // Get the most recently updated version of the data.
    //
    // If we need sync to device, then host has the most recent version.
    const bool useHostVersion = this->template need_sync<device_type> ();
    // Only one of these (the most recently updated one) is active.
    dev_view_type srcView_dev;
    host_view_type srcView_host;
    if (useHostVersion) {
      srcView_host = this->template getLocalView<Kokkos::HostSpace> ();
      if (lclNumRows != static_cast<size_t> (srcView_host.extent (0))) {
        // Make sure the number of rows is correct.  If not, take a subview.
        const Kokkos::pair<size_t, size_t> rowRng (0, lclNumRows);
        srcView_host = Kokkos::subview (srcView_host, rowRng, Kokkos::ALL ());
      }
    }
    else {
      srcView_dev = this->template getLocalView<device_type> ();
      if (lclNumRows != static_cast<size_t> (srcView_dev.extent (0))) {
        // Make sure the number of rows is correct.  If not, take a subview.
        const Kokkos::pair<size_t, size_t> rowRng (0, lclNumRows);
        srcView_dev = Kokkos::subview (srcView_dev, rowRng, Kokkos::ALL ());
      }
    }

    // If this MultiVector's local data are stored contiguously, we
    // can use the local View as the source buffer in the
    // MPI_Allreduce.  Otherwise, we have to allocate a temporary
    // source buffer and pack.
    const bool contig = isConstantStride () && getStride () == lclNumRows;
    dev_view_type srcBuf_dev;
    host_view_type srcBuf_host;
    if (useHostVersion) {
      if (contig) {
        srcBuf_host = srcView_host; // use the View as-is
      }
      else { // need to repack into a contiguous buffer
        srcBuf_host = decltype (srcBuf_host) ("srcBuf", lclNumRows, numCols);
        Kokkos::deep_copy (srcBuf_host, srcView_host);
      }
    }
    else { // use device version
      if (contig) {
        srcBuf_dev = srcView_dev; // use the View as-is
      }
      else { // need to repack into a contiguous buffer
        srcBuf_dev = decltype (srcBuf_dev) ("srcBuf", lclNumRows, numCols);
        Kokkos::deep_copy (srcBuf_dev, srcView_dev);
      }
    }

    // Check expected invariant of the above block of code.  At this
    // point, either the srcBuf of choice points to the srcView of
    // choice, or it has the right allocation size.
    {
      // Use >=, not ==, because if srcBuf just points to srcView,
      // then srcView may actually be bigger than what we need.
      const bool correct =
        (useHostVersion && static_cast<size_t> (srcBuf_host.size ()) >= totalAllocSize) ||
        (! useHostVersion && static_cast<size_t> (srcBuf_dev.size ()) >= totalAllocSize);
      TEUCHOS_TEST_FOR_EXCEPTION
        (! correct, std::logic_error, "Tpetra::MultiVector::reduce: Violated "
         "invariant of temporary source buffer construction.  Please report "
         "this bug to the Tpetra developers.");
    }

    // MPI requires that the send and receive buffers don't alias one
    // another, so we have to copy temporary storage for the result.
    //
    // We expect that MPI implementations will know how to read device
    // pointers.
    dev_view_type tgtBuf_dev;
    host_view_type tgtBuf_host;
    if (useHostVersion) {
      tgtBuf_host = host_view_type ("tgtBuf", lclNumRows, numCols);
    }
    else {
      tgtBuf_dev = dev_view_type ("tgtBuf", lclNumRows, numCols);
    }

    const int reduceCount = static_cast<int> (totalAllocSize);
    if (useHostVersion) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (tgtBuf_host.size ()) < totalAllocSize,
         std::logic_error, "Tpetra::MultiVector::reduce: tgtBuf_host.size() = "
         << tgtBuf_host.size () << " < lclNumRows*numCols = " << totalAllocSize
         << ".  Please report this bug to the Tpetra developers.");
      reduceAll<int, impl_scalar_type> (comm, REDUCE_SUM, reduceCount,
                                        srcBuf_host.data (),
                                        tgtBuf_host.data ());
    }
    else { // use device version
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (tgtBuf_dev.size ()) < totalAllocSize,
         std::logic_error, "Tpetra::MultiVector::reduce: tgtBuf_dev.size() = "
         << tgtBuf_dev.size () << " < lclNumRows*numCols = " << totalAllocSize
         << ".  Please report this bug to the Tpetra developers.");
      reduceAll<int, impl_scalar_type> (comm, REDUCE_SUM, reduceCount,
                                        srcBuf_dev.data (),
                                        tgtBuf_dev.data ());
    }

    // Write back the results to *this.
    if (useHostVersion) {
      this->template modify<Kokkos::HostSpace> ();
      if (contig || isConstantStride ()) {
        Kokkos::deep_copy (srcView_host, tgtBuf_host);
      }
      else { // need to copy one column at a time
        for (size_t j = 0; j < numCols; ++j) {
          auto X_j_out = Kokkos::subview (srcView_host, Kokkos::ALL (), j);
          auto X_j_in = Kokkos::subview (tgtBuf_host, Kokkos::ALL (), j);
          Kokkos::deep_copy (X_j_out, X_j_in);
        }
      }
    }
    else { // use device version
      this->template modify<device_type> ();
      if (contig || isConstantStride ()) {
        Kokkos::deep_copy (srcView_dev, tgtBuf_dev);
      }
      else { // need to copy one column at a time
        for (size_t j = 0; j < numCols; ++j) {
          auto X_j_out = Kokkos::subview (srcView_dev, Kokkos::ALL (), j);
          auto X_j_in = Kokkos::subview (tgtBuf_dev, Kokkos::ALL (), j);
          Kokkos::deep_copy (X_j_out, X_j_in);
        }
      }
    }

    // We leave *this unsynchronized.
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
       "Global row index " << gblRow << "is not present on this process "
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
      << "is not present on this process "
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
  typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  getDualView () const {
    return view_;
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
        auto X_dual = this->getDualView ();
        typename dual_view_type::t_host X_host;
        if (X_dual.template need_sync<Kokkos::HostSpace> ()) {
          // Device memory has the latest version.  Don't actually
          // sync to host; that changes the Vector's state, and may
          // change future computations (that use the data's current
          // place to decide where to compute).  Instead, create a
          // temporary host copy and print that.
          auto X_dev = this->template getLocalView<device_type> ();
          auto X_host_copy = Kokkos::create_mirror_view (X_dev);
          Kokkos::deep_copy (X_host_copy, X_dev);
          X_host = X_host_copy;
        }
        else {
          // Either host and device are in sync, or host has the
          // latest version of the Vector's data.  Thus, we can use
          // the host version directly.
          X_host = this->template getLocalView<Kokkos::HostSpace> ();
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
      Tpetra::Details::gathervPrint (out, lclStr, *comm);
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
    typedef LocalOrdinal LO;
    typedef device_type DT;
    typedef typename dual_view_type::host_mirror_space::device_type HMDT;
    const char prefix[] = "Tpetra::deep_copy (MultiVector): ";
    const bool debug = false;

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
    this->view_.modified_device () = 0;
    this->view_.modified_host () = 0;

    if (debug && this->getMap ()->getComm ()->getRank () == 0) {
      std::cout << "*** MultiVector::assign: ";
    }

    if (src.isConstantStride () && this->isConstantStride ()) {
      if (debug && this->getMap ()->getComm ()->getRank () == 0) {
        std::cout << "Both *this and src have constant stride" << std::endl;
      }

      // If we need sync to device, then host has the most recent version.
      const bool useHostVersion = src.template need_sync<device_type> ();

      if (useHostVersion) { // Host memory has the most recent version of src.
        this->template modify<HMDT> (); // We are about to modify dst on host.
        // Copy from src to dst on host.
        Details::localDeepCopyConstStride (this->template getLocalView<HMDT> (),
                                           src.template getLocalView<HMDT> ());
        this->template sync<DT> (); // Sync dst from host to device.
      }
      else { // Device memory has the most recent version of src.
        this->template modify<DT> (); // We are about to modify dst on device.
        // Copy from src to dst on device.
        Details::localDeepCopyConstStride (this->template getLocalView<DT> (),
                                           src.template getLocalView<DT> ());
        this->template sync<HMDT> (); // Sync dst from device to host.
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
        //
        // If we need sync to device, then host has the most recent version.
        const bool useHostVersion = src.template need_sync<device_type> ();
        if (useHostVersion) { // host version of src most recently modified
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
          Details::localDeepCopy (this->template getLocalView<HMDT> (),
                                  src.template getLocalView<HMDT> (),
                                  true, false, srcWhichVecs, srcWhichVecs);
          // Sync dst back to the device, since we only copied on the host.
          this->template sync<DT> ();
        }
        else { // copy from the device version of src
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
          Details::localDeepCopy (this->template getLocalView<DT> (),
                                  src.template getLocalView<DT> (),
                                  true, false, srcWhichVecs.d_view,
                                  srcWhichVecs.d_view);
          // Sync *this' DualView to the host.  This is cheaper than
          // repeating the above copy from src to *this on the host.
          this->template sync<HMDT> ();
        }

      }
      else { // dst is NOT constant stride
        if (src.isConstantStride ()) {
          if (debug && this->getMap ()->getComm ()->getRank () == 0) {
            std::cout << "Only src has constant stride" << std::endl;
          }

          // If we need sync to device, then host has the most recent version.
          const bool useHostVersion = src.template need_sync<device_type> ();
          if (useHostVersion) { // src most recently modified on host
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
            Details::localDeepCopy (this->template getLocalView<HMDT> (),
                                    src.template getLocalView<HMDT> (),
                                    this->isConstantStride (),
                                    src.isConstantStride (),
                                    whichVecs, whichVecs);
            // Sync dst back to the device, since we only copied on the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            this->template sync<DT> ();
          }
          else { // Copy from the device version of src.
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
            Details::localDeepCopy (this->template getLocalView<DT> (),
                                    src.template getLocalView<DT> (),
                                    this->isConstantStride (),
                                    src.isConstantStride (),
                                    whichVecs.d_view, whichVecs.d_view);
            // We can't sync src and repeat the above copy on the
            // host, so sync dst back to the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            this->template sync<HMDT> ();
          }
        }
        else { // neither src nor dst have constant stride
          if (debug && this->getMap ()->getComm ()->getRank () == 0) {
            std::cout << "Neither *this nor src has constant stride" << std::endl;
          }

          // If we need sync to device, then host has the most recent version.
          const bool useHostVersion = src.template need_sync<device_type> ();
          if (useHostVersion) { // copy from the host version of src
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
            Details::localDeepCopy (this->template getLocalView<HMDT> (),
                                    src.template getLocalView<HMDT> (),
                                    this->isConstantStride (),
                                    src.isConstantStride (),
                                    whichVectorsDst, whichVectorsSrc);

            // We can't sync src and repeat the above copy on the
            // host, so sync dst back to the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            this->template sync<HMDT> ();
          }
          else { // copy from the device version of src
            // whichVectorsDst tells the kernel which columns of dst
            // to copy.  Fill it on the host, and sync to device.
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
            Details::localDeepCopy (this->template getLocalView<DT> (),
                                    src.template getLocalView<DT> (),
                                    this->isConstantStride (),
                                    src.isConstantStride (),
                                    whichVecsDst.d_view, whichVecsSrc.d_view);
          }
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  isSameSize (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & vec) const {
    typedef impl_scalar_type ST;
    typedef typename Kokkos::View<int*, device_type>::HostMirror::execution_space HES;
    size_t l1 = this->getLocalLength();
    size_t l2 = vec.getLocalLength();
    if ((l1!=l2) || (this->getNumVectors() != vec.getNumVectors()))
      return false;
    if(l1==0)  return true;

    auto v1 = this->template getLocalView<HES>();
    auto v2 = vec.template getLocalView<HES>();
    if(Details::PackTraits<ST,HES>::packValueCount(v1(0,0)) != Details::PackTraits<ST,HES>::packValueCount(v2(0,0)))
      return false;

    return true;
  }


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

  template <class ST, class LO, class GO, class NT>
  void MultiVector<ST, LO, GO, NT>::
  swap(MultiVector<ST, LO, GO, NT> & mv) {
    // Cache maps & views
    Teuchos::RCP<const map_type> map = mv.map_;
    dual_view_type  view, origView;
    Teuchos::Array<size_t> whichVectors; // FIXME: This is a deep copy 
    view         = mv.view_;
    origView     = mv.origView_;
    whichVectors = mv.whichVectors_;

    // Swap this-> mv
    mv.map_          = this->map_;
    mv.view_         = this->view_;
    mv.origView_     = this->origView_;
    mv.whichVectors_ = this->whichVectors_;
    
    // Swap mv -> this
    this->map_          = map;
    this->view_         = view;
    this->origView_     = origView;
    this->whichVectors_ = whichVectors;  
  }


} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class MultiVector< SCALAR , LO , GO , NODE >; \
  template MultiVector< SCALAR , LO , GO , NODE > createCopy( const MultiVector< SCALAR , LO , GO , NODE >& src); \
  template Teuchos::RCP<MultiVector< SCALAR , LO , GO , NODE > > createMultiVector (const Teuchos::RCP<const Map<LO, GO, NODE> >& map, size_t numVectors);

#endif // TPETRA_MULTIVECTOR_DEF_HPP
