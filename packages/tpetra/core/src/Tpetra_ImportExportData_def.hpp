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

#ifndef TPETRA_IMPORTEXPORTDATA_DEF_HPP
#define TPETRA_IMPORTEXPORTDATA_DEF_HPP

#include "Tpetra_Map.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"

namespace { // (anonymous)

  Teuchos::RCP<Teuchos::FancyOStream>
  makeValidVerboseStream (const Teuchos::RCP<Teuchos::FancyOStream>& out =
			  Teuchos::null)
  {
    if (out.is_null ()) {
      try {
	return Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
      }
      catch (std::exception& e) {
	TEUCHOS_TEST_FOR_EXCEPTION
	  (true, std::logic_error, "getFancyOStream threw: " << e.what ());
      }
      catch (...) {
	TEUCHOS_TEST_FOR_EXCEPTION
	  (true, std::logic_error, "getFancyOStream threw an exception "
	   "not a subclass of std::exception.");
      }
    }
    else {
      try {
	return out;
      }
      catch (std::exception& e) {
	TEUCHOS_TEST_FOR_EXCEPTION
	  (true, std::logic_error, "Returning 'out' threw: " << e.what ());
      }
      catch (...) {
	TEUCHOS_TEST_FOR_EXCEPTION
	  (true, std::logic_error, "Returning 'out' threw an exception "
	   "not a subclass of std::exception.");
      }
    }
  }

} // namespace (anonymous)

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal, GlobalOrdinal, Node>::
  ImportExportData (const Teuchos::RCP<const map_type>& source,
                    const Teuchos::RCP<const map_type>& target,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    source_ (source), // NOT allowed to be null
    target_ (target), // allowed to be null
    out_ (makeValidVerboseStream (out)),
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
