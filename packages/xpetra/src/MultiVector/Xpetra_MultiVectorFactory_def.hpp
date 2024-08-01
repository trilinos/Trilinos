// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MULTIVECTORFACTORY_DEF_HPP
#define XPETRA_MULTIVECTORFACTORY_DEF_HPP

#include "Xpetra_MultiVectorFactory_decl.hpp"

#include "Xpetra_BlockedMultiVector.hpp"

#include "Xpetra_BlockedMap.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
          size_t NumVectors,
          bool zeroOut) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

  RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

  if (!bmap.is_null()) {
    return rcp(new Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, NumVectors, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, NumVectors, zeroOut));
  }
#endif

  XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
  XPETRA_FACTORY_END;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>>& ArrayOfPtrs,
          size_t NumVectors) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, ArrayOfPtrs, NumVectors));
  }
#endif

  XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
  XPETRA_FACTORY_END;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& source,
          Teuchos::DataAccess copyOrView) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (source->getMap()->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*source, copyOrView));
  }
#endif

  XPETRA_FACTORY_ERROR_IF_EPETRA(source->getMap()->lib());
  XPETRA_FACTORY_END;
}

}  // namespace Xpetra

#endif
