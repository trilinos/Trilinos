// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MAPEXTRACTORFACTORY_HPP_
#define XPETRA_MAPEXTRACTORFACTORY_HPP_

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Xpetra_MapExtractor.hpp>

namespace Xpetra {

// factory class
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class MapExtractorFactory {
#undef XPETRA_MAPEXTRACTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private construtor, since this is a static class
  MapExtractorFactory() = delete;

 public:
  /*!
    @brief Build MapExtrasctor from a given set of partial maps.

    The Maps indirectly specify the linear algebra library to use (Tpetra or Epetra).
  */
  static Teuchos::RCP<Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Map>& fullmap,
        const std::vector<Teuchos::RCP<const Map>>& maps,
        bool bThyraMode = false) {
    return rcp(new Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullmap, maps, bThyraMode));
  }

  /*!
    @brief Build MapExtrasctor from a given BlockedMap.

    The Maps indirectly specify the linear algebra library to use (Tpetra or Epetra).
  */
  static Teuchos::RCP<Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>& blockedMap) {
    return rcp(new MapExtractor(blockedMap));
  }
};

}  // namespace Xpetra

#define XPETRA_MAPEXTRACTORFACTORY_SHORT
#endif /* XPETRA_MAPEXTRACTORFACTORY_HPP_ */
