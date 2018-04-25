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

#ifndef TPETRA_FEMULTIVECTOR_DEF_HPP
#define TPETRA_FEMULTIVECTOR_DEF_HPP

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
  // TODO: @CMS:  This is copied verbatim from MultiVector to get a compile--
  //              should this function get moved somewhere so we aren't copying
  //              code like this?
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
        (static_cast<size_t> (d_view.dimension_0 ()) != lclNumRows ||
         static_cast<size_t> (d_view.dimension_1 ()) != numCols, std::logic_error,
         "allocDualView: d_view's dimensions actual dimensions do not match "
         "requested dimensions.  d_view is " << d_view.dimension_0 () << " x " <<
         d_view.dimension_1 () << "; requested " << lclNumRows << " x " << numCols
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

#if 0
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
    if (X.dimension_0 () == 0 && X.dimension_1 () != 0) {
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
    if (X.dimension_0 () == 0 && X.dimension_1 () != 0) {
      return DualViewType ("MV::DualView", 0, colRng.second - colRng.first);
    }
    else {
      return subview (X, rowRng, colRng);
    }
  }
  #endif

} // namespace (anonymous)



namespace Tpetra {



template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FEMultiVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > & map,
              const Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >& importer,
              const size_t numVecs,
              const bool zeroOut)
              {
                // TODO: @CMS this should probably do something..
              }



#if 0  // WCMCLEN - use something like this for the c'tor ^^^
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  FEMultiVector (const Teuchos::RCP<const map_type>& map,
                 const size_t numVecs,
                 const bool zeroOut) : /* default is true */
    base_type (map)
  {
    ::Tpetra::Details::ProfilingRegion region ("Tpetra::MV ctor (map,numVecs,zeroOut)");

    const size_t lclNumRows = this->getLocalLength ();
    view_ = allocDualView<Scalar, LocalOrdinal, GlobalOrdinal, Node> (lclNumRows, numVecs, zeroOut);
    origView_ = view_;
  }
#endif



template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceMap (const Teuchos::RCP<const map_type>& newMap) {
  throw std::runtime_error("Tpetra::FEMultiVecotr::replaceMap() is not implemented");
}  // replaceMap ()


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doTargetToSource(const CombineMode CM) {
  throw std::runtime_error("stub");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::switchActiveMultiVector() {
  throw std::runtime_error("stub");
}



} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_FEMULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class FEMultiVector< SCALAR , LO , GO , NODE >;

#endif // TPETRA_FEMULTIVECTOR_DEF_HPP


