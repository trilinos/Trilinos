// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_IMPORT_UTILS_HPP_
#define PACKAGES_MUELU_IMPORT_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Map.hpp"  // definition of UnderlyingLib
#include "Xpetra_Import.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

#include <utility>

#include "Xpetra_TpetraImport.hpp"
#include "Tpetra_Import_Util.hpp"

namespace MueLu {

/*!
  @class ImportUtils
  @brief MueLu utility class for Import-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.

*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ImportUtils {
#undef MUELU_IMPORTUTILS_SHORT

 public:
  /// \brief For each GID in the TargetMap, find who owns the GID in the SourceMap.
  ///
  /// This only uses the Distributor and does not communicate.  It
  /// returns (as an output argument) an array of (PID,GID) pairs.
  /// If use_minus_one_for_local is true, any GIDs owned by this
  /// processor get -1 instead of their PID.
  void
  getPidGidPairs(const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>& Importer,
                 Teuchos::Array<std::pair<int, GlobalOrdinal> >& gpids,
                 bool use_minus_one_for_local) {
    Xpetra::UnderlyingLib lib = Importer.getSourceMap()->lib();
    if (lib == Xpetra::UseTpetra) {
      Tpetra::Import_Util::getPidGidPairs(Xpetra::toTpetra(Importer), gpids, use_minus_one_for_local);
    }
  }

  //! Like getPidGidPairs, but just gets the PIDs, ordered by the column Map.
  void
  getPids(const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>& Importer,
          Teuchos::Array<int>& pids,
          bool use_minus_one_for_local) {
    Xpetra::UnderlyingLib lib = Importer.getSourceMap()->lib();
    if (lib == Xpetra::UseTpetra) {
      Tpetra::Import_Util::getPids(Xpetra::toTpetra(Importer), pids, use_minus_one_for_local);
    }
  }

  //! Like getPidGidPairs, but just gets the PIDs, ordered by the column Map.
  // Like the above, but without the resize
  void
  getPids(const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>& Importer,
          Teuchos::ArrayView<int>& pids,
          bool use_minus_one_for_local) {
    Xpetra::UnderlyingLib lib = Importer.getSourceMap()->lib();
    if (lib == Xpetra::UseTpetra) {
      Tpetra::Import_Util::getPids(Xpetra::toTpetra(Importer), pids, use_minus_one_for_local);
    }
  }

  /// \brief Get a list of remote PIDs from an importer in the order
  ///   corresponding to the remote LIDs.
  void
  getRemotePIDs(const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>& Importer,
                Teuchos::Array<int>& RemotePIDs) {
    Xpetra::UnderlyingLib lib = Importer.getSourceMap()->lib();
    if (lib == Xpetra::UseTpetra) {
      Tpetra::Import_Util::getRemotePIDs(Xpetra::toTpetra(Importer), RemotePIDs);
    }
  }

};  // end class ImportUtils

}  // namespace MueLu

#define MUELU_IMPORTUTILS_SHORT

#endif  // PACKAGES_MUELU_IMPORT_UTILS_HPP_
