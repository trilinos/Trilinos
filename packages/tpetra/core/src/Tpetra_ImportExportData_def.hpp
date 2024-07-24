// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_IMPORTEXPORTDATA_DEF_HPP
#define TPETRA_IMPORTEXPORTDATA_DEF_HPP

#include "Tpetra_Map.hpp"
#include "Tpetra_Details_makeValidVerboseStream.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal, GlobalOrdinal, Node>::
  ImportExportData (const Teuchos::RCP<const map_type>& source,
                    const Teuchos::RCP<const map_type>& target,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    source_ (source), // NOT allowed to be null
    target_ (target), // allowed to be null
    out_ (::Tpetra::Details::makeValidVerboseStream (out)),
    numSameIDs_ (0), // Import/Export constructor may change this
    distributor_ (source->getComm (), out_, plist), // Im/Ex ctor will init
    isLocallyComplete_ (true) // Im/Ex ctor may change this
  {
    TEUCHOS_ASSERT( ! out_.is_null () );
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal, GlobalOrdinal, Node>::
  ImportExportData (const Teuchos::RCP<const map_type>& source,
                    const Teuchos::RCP<const map_type>& target) :
    ImportExportData (source, target, Teuchos::null, Teuchos::null)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal, GlobalOrdinal, Node>::
  ImportExportData (const Teuchos::RCP<const map_type>& source,
                    const Teuchos::RCP<const map_type>& target,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out) :
    ImportExportData (source, target, out, Teuchos::null)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal, GlobalOrdinal, Node>::
  ImportExportData (const Teuchos::RCP<const map_type>& source,
                    const Teuchos::RCP<const map_type>& target,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    ImportExportData (source, target, Teuchos::null, plist)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<ImportExportData<LocalOrdinal, GlobalOrdinal, Node> >
  ImportExportData<LocalOrdinal, GlobalOrdinal, Node>::
  reverseClone ()
  {
    using Teuchos::ArrayView;
    using data_type = ImportExportData<LocalOrdinal, GlobalOrdinal, Node>;

    auto tData = Teuchos::rcp (new data_type (target_, source_, out_));

    // Things that stay the same
    tData->numSameIDs_       = numSameIDs_;

    // Things that reverse
    tData->distributor_      = *distributor_.getReverse();
    tData->permuteToLIDs_    = permuteFromLIDs_;
    tData->permuteFromLIDs_  = permuteToLIDs_;

    // Remotes / exports (easy part)
    tData->exportLIDs_       = remoteLIDs_;
    tData->remoteLIDs_       = exportLIDs_;
    tData->exportPIDs_.resize (tData->exportLIDs_.extent (0));

    // Remotes / exports (hard part) - extract the exportPIDs from the remotes of my distributor
    const size_t NumReceives            = distributor_.getNumReceives();
    ArrayView<const int> ProcsFrom      = distributor_.getProcsFrom();
    ArrayView<const size_t> LengthsFrom = distributor_.getLengthsFrom();

    // isLocallyComplete is a local predicate.
    // It could be true in one direction but false in another.

    bool isLocallyComplete = true; // by default
    for (size_t i = 0, j = 0; i < NumReceives; ++i) {
      const int pid = ProcsFrom[i];
      if (pid == -1) {
        isLocallyComplete = false;
      }
      for (size_t k = 0; k < LengthsFrom[i]; ++k) {
        tData->exportPIDs_[j] = pid;
        ++j;
      }
    }
    tData->isLocallyComplete_ = isLocallyComplete;

    return tData;
  }

} // namespace Tpetra

// Explicit instantiation macro.
// Only invoke this when in the Tpetra namespace.
// Most users do not need to use this.
//
// LO: The local ordinal type.
// GO: The global ordinal type.
// NODE: The Kokkos Node type.
#define TPETRA_IMPORTEXPORTDATA_INSTANT(LO, GO, NODE) \
  template class ImportExportData< LO , GO , NODE >;

#endif // TPETRA_IMPORTEXPORTDATA_DEF_HPP
