// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_MultiVectorFactory_decl.hpp"
#include "Xpetra_BlockedMultiVector.hpp"

namespace Xpetra {

// we need the Epetra specialization only if Epetra is enabled
#if defined(HAVE_XPETRA_EPETRA)

#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

MultiVectorFactory<double, int, int, EpetraNode>::
    MultiVectorFactory() {
}

RCP<MultiVector<double, int, int, EpetraNode>>
MultiVectorFactory<double, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>& map, size_t NumVectors, bool zeroOut) {
  using BlockedMultiVector = Xpetra::BlockedMultiVector<double, int, int, EpetraNode>

      XPETRA_MONITOR("MultiVectorFactory::Build");

  RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, int, EpetraNode>>(map);

  if (!bmap.is_null()) {
    return rcp(new BlockedMultiVector(bmap, NumVectors, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<double, int, int, EpetraNode>(map, NumVectors, zeroOut));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraMultiVectorT<int, EpetraNode>(map, NumVectors, zeroOut));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<double, int, int, EpetraNode>>
MultiVectorFactory<double, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const double>>& ArrayOfPtrs,
          size_t NumVectors) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<double, int, int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraMultiVectorT<int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<double, int, int, EpetraNode>>
MultiVectorFactory<double, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const MultiVector<double, int, int, EpetraNode>>& source,
          Teuchos::DataAccess copyOrView) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (source->getMap()->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<double, int, int, EpetraNode>(*source, copyOrView));
  }
#endif

  if (source->getMap()->lib() == UseEpetra) {
    return rcp(new EpetraMultiVectorT<int, EpetraNode>(*source, copyOrView));
  }

  XPETRA_FACTORY_END;
}

MultiVectorFactory<int, int, int, EpetraNode>::
    MultiVectorFactory() {
}

RCP<MultiVector<int, int, int, EpetraNode>>
MultiVectorFactory<int, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>& map, size_t NumVectors, bool zeroOut) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

  RCP<const BlockedMap<int, int, EpetraNode>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, int, EpetraNode>>(map);

  if (!bmap.is_null()) {
    return rcp(new BlockedMultiVector<int, int, int, EpetraNode>(bmap, NumVectors, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<int, int, int, EpetraNode>(map, NumVectors, zeroOut));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraIntMultiVectorT<int, EpetraNode>(map, NumVectors, zeroOut));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<int, int, int, EpetraNode>>
MultiVectorFactory<int, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const int>>& ArrayOfPtrs,
          size_t NumVectors) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<int, int, int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraIntMultiVectorT<int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<int, int, int, EpetraNode>>
MultiVectorFactory<int, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const MultiVector<int, int, int, EpetraNode>>& source,
          Teuchos::DataAccess copyOrView) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (source->getMap()->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<int, int, int, EpetraNode>(*source, copyOrView));
  }
#endif

  if (source->getMap()->lib() == UseEpetra) {
    return rcp(new EpetraIntMultiVectorT<int, EpetraNode>(*source, copyOrView));
  }

  XPETRA_FACTORY_END;
}

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

MultiVectorFactory<double, int, long long, EpetraNode>::MultiVectorFactory() {}

RCP<MultiVector<double, int, long long, EpetraNode>>
MultiVectorFactory<double, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>>& map,
          size_t NumVectors,
          bool zeroOut) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

  RCP<const BlockedMap<int, long long, EpetraNode>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, long long, EpetraNode>>(map);

  if (!bmap.is_null()) {
    return rcp(new BlockedMultiVector<double, int, long long, EpetraNode>(bmap, NumVectors, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<double, int, long long, EpetraNode>(map, NumVectors, zeroOut));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraMultiVectorT<long long, EpetraNode>(map, NumVectors, zeroOut));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<double, int, long long, EpetraNode>>
MultiVectorFactory<double, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>>& ArrayOfPtrs,
          size_t NumVectors) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<double, int, long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraMultiVectorT<long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<double, int, long long, EpetraNode>>
MultiVectorFactory<double, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const MultiVector<double, int, long long, EpetraNode>>& source,
          Teuchos::DataAccess copyOrView) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (source->getMap()->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<double, int, long long, EpetraNode>(*source, copyOrView));
  }
#endif

  if (source->getMap()->lib() == UseEpetra) {
    return rcp(new EpetraMultiVectorT<long long, EpetraNode>(*source, copyOrView));
  }

  XPETRA_FACTORY_END;
}

MultiVectorFactory<int, int, long long, EpetraNode>::
    MultiVectorFactory() {
}

RCP<MultiVector<int, int, long long, EpetraNode>>
MultiVectorFactory<int, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>>& map,
          size_t NumVectors,
          bool zeroOut) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

  RCP<const BlockedMap<int, long long, EpetraNode>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, long long, EpetraNode>>(map);

  if (!bmap.is_null()) {
    return rcp(new BlockedMultiVector<int, int, long long, EpetraNode>(bmap, NumVectors, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<int, int, long long, EpetraNode>(map, NumVectors, zeroOut));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraIntMultiVectorT<long long, EpetraNode>(map, NumVectors, zeroOut));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<int, int, long long, EpetraNode>>
MultiVectorFactory<int, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, long long, Node>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const int>>& ArrayOfPtrs,
          size_t NumVectors) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<int, int, long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraIntMultiVectorT<long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<MultiVector<int, int, long long, EpetraNode>>
MultiVectorFactory<int, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const MultiVector<int, int, long long, EpetraNode>>& source,
          Teuchos::DataAccess copyOrView) {
  XPETRA_MONITOR("MultiVectorFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (source->getMap()->lib() == UseTpetra) {
    return rcp(new TpetraMultiVector<int, int, long long, EpetraNode>(*source, copyOrView));
  }
#endif

  if (source->getMap()->lib() == UseEpetra) {
    return rcp(new EpetraIntMultiVectorT<long long, EpetraNode>(*source, copyOrView));
  }

  XPETRA_FACTORY_END;
}

#endif  // END !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

#endif  // END !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

#endif  // END HAVE_XPETRA_EPETRA

}  // namespace Xpetra
