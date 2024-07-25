// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_BlockedVector.hpp"

namespace Xpetra {

#if defined(HAVE_XPETRA_EPETRA)

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<double, int, int, EpetraNode>>
VectorFactory<double, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Xpetra::Map<int, int, EpetraNode>>& map, bool zeroOut) {
  XPETRA_MONITOR("VectorFactory::Build");

  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
      bmap = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

  if (!bmap.is_null()) {
    return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraVector(map, zeroOut));
  }
#endif  // HAVE_XPETRA_TPETRA

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraVectorT<GlobalOrdinal, EpetraNode>(map, zeroOut));
  }

  XPETRA_FACTORY_END;
}

#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<double, int, long long, EpetraNode>>
VectorFactory<double, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Xpetra::Map<int, long long, EpetraNode>>& map, bool zeroOut) {
  XPETRA_MONITOR("VectorFactory::Build");

  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
  if (!bmap.is_null()) {
    return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraVector(map, zeroOut));
  }
#endif

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraVectorT<GlobalOrdinal, Node>(map, zeroOut));
  }

  XPETRA_FACTORY_END;
}

#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<int, int, int, EpetraNode>>
VectorFactory<int, int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Xpetra::Map<int, int, EpetraNode>>& map, bool zeroOut) {
  XPETRA_MONITOR("VectorFactory::Build");

  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
  if (!bmap.is_null()) {
    return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraVector(map, zeroOut));
  }
#endif  // HAVE_XPETRA_TPETRA

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraIntVectorT<GlobalOrdinal, Node>(map, zeroOut));
  }

  XPETRA_FACTORY_END;
}

#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<int, int, long long, EpetraNode>>
VectorFactory<int, int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Xpetra::Map<int, long long, EpetraNode>>& map, bool zeroOut) {
  XPETRA_MONITOR("VectorFactory::Build");

  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

  if (!bmap.is_null()) {
    return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraVector(map, zeroOut));
  }
#endif  // HAVE_XPETRA_TPETRA

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraIntVectorT<GlobalOrdinal, Node>(map, zeroOut));
  }

  XPETRA_FACTORY_END;
}

#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

#endif  // #if defined(HAVE_XPETRA_EPETRA)

}  // namespace Xpetra
