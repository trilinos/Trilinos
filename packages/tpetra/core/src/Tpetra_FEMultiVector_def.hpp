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

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Details_Behavior.hpp"


namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FEMultiVector (const Teuchos::RCP<const map_type>& map,
               const Teuchos::RCP<const Import<local_ordinal_type, global_ordinal_type, node_type>>& importer,
               const size_t numVecs,
               const bool zeroOut) :
  base_type (importer.is_null () ? map : importer->getTargetMap (),
             numVecs, zeroOut),
  activeMultiVector_ (Teuchos::rcp (new FEWhichActive (FE_ACTIVE_OWNED_PLUS_SHARED))),
  importer_ (importer)
{
  const char tfecfFuncName[] = "FEMultiVector constructor: ";

  if (! importer_.is_null ()) {
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if (debug) {
      // Checking Map sameness may require an all-reduce, so we should
      // reserve it for debug mode.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! importer_->getSourceMap ()->isSameAs (*map),
         std::runtime_error,
         "If you provide a nonnull Import, then the input Map "
         "must be the same as the input Import's source Map.");

      // Checking whether one Map is locally fitted to another could be
      // expensive.
      const bool locallyFitted =
        importer->getTargetMap ()->isLocallyFitted (* (importer->getSourceMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! locallyFitted, std::runtime_error,
         "If you provide a nonnull Import, then its target Map must be "
         "locally fitted (see Map::isLocallyFitted documentation) to its "
         "source Map.");
    }

    using range_type = Kokkos::pair<size_t, size_t>;
    auto dv = Kokkos::subview (this->view_,
                               range_type (0, map->getNodeNumElements ()),
                               Kokkos::ALL ());
    // Memory aliasing is required for FEMultiVector
    inactiveMultiVector_ =
      Teuchos::rcp (new base_type (importer_->getSourceMap (), dv));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
beginFill ()
{
  // The FEMultiVector is in owned+shared mode on construction, so we
  // do not throw in that case.
  if (*activeMultiVector_ == FE_ACTIVE_OWNED) {
    switchActiveMultiVector ();
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
endFill ()
{
  const char tfecfFuncName[] = "endFill: ";

  if (*activeMultiVector_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    doOwnedPlusSharedToOwned (Tpetra::ADD);
    switchActiveMultiVector ();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (true, std::runtime_error, "Owned+Shared MultiVector already active; "
       "cannot call endFill.");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
globalAssemble ()
{
  endFill ();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
replaceMap (const Teuchos::RCP<const map_type>& /* newMap */)
{
  const char tfecfFuncName[] = "replaceMap: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (true, std::runtime_error, "This method is not implemented.");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
doOwnedPlusSharedToOwned (const CombineMode CM)
{
  if (! importer_.is_null () &&
      *activeMultiVector_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    inactiveMultiVector_->doExport (*this, *importer_, CM);
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
doOwnedToOwnedPlusShared (const CombineMode CM)
{
  if (! importer_.is_null () &&
      *activeMultiVector_ == FE_ACTIVE_OWNED) {
    inactiveMultiVector_->doImport (*this, *importer_, CM);
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
switchActiveMultiVector ()
{
  if (*activeMultiVector_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    *activeMultiVector_ = FE_ACTIVE_OWNED;
  }
  else {
    *activeMultiVector_ = FE_ACTIVE_OWNED_PLUS_SHARED;
  }

  if (importer_.is_null ()) {
    return;
  }

  // Use MultiVector's swap routine here
  this->swap (*inactiveMultiVector_);
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
