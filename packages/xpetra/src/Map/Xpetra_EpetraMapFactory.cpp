// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_MapFactory.hpp"

#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_EpetraMap.hpp"
#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#endif

namespace Xpetra {

#if defined(HAVE_XPETRA_EPETRA)

#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

MapFactory<int, int, EpetraNode>::
    MapFactory() {
}

RCP<Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
          LocalGlobal lg) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
    return rcp(new Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, indexBase, comm, lg));
#endif

  if (lib == UseEpetra)
    return rcp(new EpetraMapT<int, Node>(numGlobalElements, indexBase, comm, lg));

  XPETRA_FACTORY_END;
}

RCP<Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          size_t numLocalElements,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, numLocalElements, indexBase, comm));
#endif

  if (lib == UseEpetra)
    return rcp(new EpetraMapT<int, Node>(numGlobalElements, numLocalElements, indexBase, comm));

  XPETRA_FACTORY_END;
}

RCP<Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const int> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");
#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, elementList, indexBase, comm));
#endif  // HAVE_XPETRA_TPETRA

  if (lib == UseEpetra)
    return rcp(new EpetraMapT<int, Node>(numGlobalElements, elementList, indexBase, comm));

  XPETRA_FACTORY_END;
}

//! Map constructor transforming degrees of freedom
//! for numDofPerNode this acts like a deep copy
Teuchos::RCP<Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, int, EpetraNode>> &map,
          const int numDofPerNode, const int gidOffset) {
  XPETRA_MONITOR("MapFactory::Build");

  RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
  if (!bmap.is_null()) {
    TEUCHOS_TEST_FOR_EXCEPTION(numDofPerNode != 1, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::Build: When provided a BlockedMap numDofPerNode must set to be one. It is set to " << numDofPerNode << ".");
    return rcp(new BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(*bmap));
  }

  LocalOrdinal N                                      = Teuchos::as<LocalOrdinal>(map->getLocalNumElements());
  Teuchos::ArrayView<const GlobalOrdinal> oldElements = map->getLocalElementList();
  Teuchos::Array<GlobalOrdinal> newElements(map->getLocalNumElements() * numDofPerNode);
  for (LocalOrdinal i = 0; i < N; i++) {
    for (LocalOrdinal j = 0; j < numDofPerNode; j++) {
      newElements[i * numDofPerNode + j] = oldElements[i] * numDofPerNode + j + gidOffset;
    }
  }

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra) {
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(map->getGlobalNumElements() * numDofPerNode, newElements, map->getIndexBase(), map->getComm()));
  }
#endif  // HAVE_XPETRA_TPETRA

  if (map->lib() == UseEpetra) {
    return rcp(new EpetraMapT<int, Node>(map->getGlobalNumElements() * numDofPerNode, newElements, map->getIndexBase(), map->getComm()));
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    createLocalMap(UnderlyingLib lib,
                   size_t numElements,
                   const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createLocalMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<int, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<int, Node>((Xpetra::global_size_t)numElements,  // num elements, global and local
                                                 0,                                   // index base is zero
                                                 comm, LocallyReplicated));
    return map.getConst();
  }

  XPETRA_FACTORY_END;
}

// TODO remove this

#ifdef HAVE_XPETRA_TPETRA
Teuchos::RCP<Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Kokkos::View<const int *, typename Node::device_type> &indexList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, indexList, indexBase, comm));
  if (lib == UseEpetra) {
    Teuchos::ArrayView<const int> v(indexList.data(), indexList.size());
    return rcp(new EpetraMapT<int, Node>(numGlobalElements, v, indexBase, comm));
  }
  XPETRA_FACTORY_END;
}
#endif  // HAVE_XPETRA_TPETRA

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createLocalMapWithNode<int, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createLocalMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<int, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<int, Node>((Xpetra::global_size_t)numElements,  // num elements, global and local
                                                 0,                                   // index base is zero
                                                 comm, LocallyReplicated));
    return map.getConst();
  }

  XPETRA_FACTORY_END;
}

// TODO remove this

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    createUniformContigMapWithNode(UnderlyingLib lib,
                                   global_size_t numElements,
                                   const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createUniformContigMapWithNode<int, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createUniformContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<int, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<int, Node>(numElements,  // num elements, global and local
                                                 0,            // index base is zero
                                                 comm, GloballyDistributed));
    return map.getConst();
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    createUniformContigMap(UnderlyingLib lib,
                           global_size_t numElements,
                           const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createUniformContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createUniformContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<int, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<int, Node>(numElements,  // num elements, global and local
                                                 0,            // index base is zero
                                                 comm, GloballyDistributed));
    return map.getConst();
  }
  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    createContigMap(UnderlyingLib lib,
                    global_size_t numElements,
                    size_t localNumElements,
                    const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createContigMapWithNode<int, GlobalOrdinal, Node>(numElements, localNumElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createContigMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif

  if (lib == UseEpetra) {
    return MapFactory<int, GlobalOrdinal, Node>::createContigMapWithNode(lib, numElements, localNumElements, comm);
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createContigMapWithNode<int, GlobalOrdinal, Node>(numElements, localNumElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<int, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<int, Node>(numElements, localNumElements,
                                                 0,  // index base is zero
                                                 comm));
    return map.getConst();
  }
  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, int, EpetraNode>>
MapFactory<int, int, EpetraNode>::copyMapWithNewComm(const Teuchos::RCP<const Map<int, int, EpetraNode>> &oldmap,
                                                     const Teuchos::RCP<const Teuchos::Comm<int>> &newComm) {
  XPETRA_MONITOR("MapFactory::Build");
  using XMF             = Xpetra::MapFactory<int, int, EpetraNode>;
  global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();

  size_t Nlocal         = oldmap->getLocalNumElements();
  global_size_t Nglobal = oldmap->getGlobalNumElements();

  // Sanity check -- if there's no comm, we can't keep elements on the map  (vice versa is OK)
  TEUCHOS_TEST_FOR_EXCEPTION(Nlocal && newComm.is_null(),
                             std::logic_error, "MapFactory::copyMapWithNewComm needs the comm to match the map.");

  // We'll return null if we don't have a Comm on this rank
  RCP<const Map<int, int, Node>> newMap;
  if (!newComm.is_null()) {
    if (oldmap->isContiguous()) {
      newMap = XMF::Build(oldmap->lib(), INVALID, Nlocal, oldmap->getIndexBase(), newComm);
    } else {
      newMap = XMF::Build(oldmap->lib(), Nglobal, oldmap->getLocalElementList(), oldmap->getIndexBase(), newComm);
    }
  }

  return newMap;
  XPETRA_FACTORY_END;
}

#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

//! Private constructor. This is a static class.
MapFactory<int, long long, EpetraNode>::
    MapFactory() {
}

RCP<Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
          LocalGlobal lg) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, indexBase, comm, lg));
#endif

  if (lib == UseEpetra)
    return rcp(new EpetraMapT<long long, Node>(numGlobalElements, indexBase, comm, lg));

  XPETRA_FACTORY_END;
}

RCP<Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          size_t numLocalElements,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, numLocalElements, indexBase, comm));
#endif

  if (lib == UseEpetra)
    return rcp(new EpetraMapT<long long, Node>(numGlobalElements, numLocalElements, indexBase, comm));

  XPETRA_FACTORY_END;
}

RCP<Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const long long> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, elementList, indexBase, comm));
#endif

  if (lib == UseEpetra)
    return rcp(new EpetraMapT<long long, Node>(numGlobalElements, elementList, indexBase, comm));

  XPETRA_FACTORY_END;
}

//! Map constructor transforming degrees of freedom
//! for numDofPerNode this acts like a deep copy
Teuchos::RCP<Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>> &map,
          int numDofPerNode) {
  XPETRA_MONITOR("MapFactory::Build");

  RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
  if (!bmap.is_null()) {
    TEUCHOS_TEST_FOR_EXCEPTION(numDofPerNode != 1, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::Build: When provided a BlockedMap numDofPerNode must set to be one. It is set to " << numDofPerNode << ".");
    return rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(*bmap));
  }

  LocalOrdinal N                                      = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> oldElements = map->getLocalElementList();
  Teuchos::Array<GlobalOrdinal> newElements(map->getLocalNumElements() * numDofPerNode);
  for (LocalOrdinal i = 0; i < N; i++)
    for (LocalOrdinal j = 0; j < numDofPerNode; j++)
      newElements[i * numDofPerNode + j] = oldElements[i] * numDofPerNode + j;

#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(map->getGlobalNumElements() * numDofPerNode, newElements, map->getIndexBase(), map->getComm()));
#endif

  if (map->lib() == UseEpetra)
    return rcp(new EpetraMapT<long long, Node>(map->getGlobalNumElements() * numDofPerNode, newElements, map->getIndexBase(), map->getComm()));

  XPETRA_FACTORY_END;
}

#ifdef HAVE_XPETRA_TPETRA
Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
MapFactory<int, long long, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Kokkos::View<const long long *, typename Node::device_type> &indexList,
          long long indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");
  if (lib == UseTpetra)
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements, indexList, indexBase, comm));
  if (lib == UseEpetra) {
    Teuchos::ArrayView<const long long> v(indexList.data(), indexList.size());
    return rcp(new EpetraMapT<long long, Node>(numGlobalElements, v, indexBase, comm));
  }
  XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
  XPETRA_FACTORY_END;
}
#endif  // HAVE_XPETRA_TPETRA

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    createLocalMap(UnderlyingLib lib,
                   size_t numElements,
                   const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createLocalMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

  if (lib == UseEpetra)
    return MapFactory<int, GlobalOrdinal, Node>::createLocalMapWithNode(lib, numElements, comm);

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createLocalMapWithNode<int, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createLocalMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<long long, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<long long, Node>((Xpetra::global_size_t)numElements,  // num elements, global and local
                                                       0,                                   // index base is zero
                                                       comm, LocallyReplicated));
    return map.getConst();
  }
  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    createUniformContigMapWithNode(UnderlyingLib lib,
                                   global_size_t numElements,
                                   const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createUniformContigMapWithNode<int, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createUniformContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<long long, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<long long, Node>(numElements,  // num elements, global and local
                                                       0,            // index base is zero
                                                       comm, GloballyDistributed));
    return map.getConst();
  }
  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    createUniformContigMap(UnderlyingLib lib,
                           global_size_t numElements,
                           const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createUniformContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(numElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createUniformContigMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

  if (lib == UseEpetra)
    return MapFactory<int, GlobalOrdinal, Node>::createUniformContigMapWithNode(lib, numElements, comm);

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::createContigMap(UnderlyingLib lib,
                                                        global_size_t numElements,
                                                        size_t localNumElements,
                                                        const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(
        Tpetra::createContigMapWithNode<int, GlobalOrdinal, Node>(numElements, localNumElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true,
                               Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createContigMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on "
                               "EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

  if (lib == UseEpetra)
    return MapFactory<int, GlobalOrdinal, Node>::createContigMapWithNode(lib, numElements, localNumElements, comm);

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP<const Teuchos::Comm<int>> &comm) {
  XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(
        Tpetra::createContigMapWithNode<int, GlobalOrdinal, Node>(numElements, localNumElements, comm)));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true,
                               Xpetra::Exceptions::RuntimeError,
                               "Xpetra::MapFactory::createContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not "
                               "instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif  // HAVE_XPETRA_TPETRA

  if (lib == UseEpetra) {
    Teuchos::RCP<EpetraMapT<long long, Node>> map;
    map = Teuchos::rcp(new EpetraMapT<long long, Node>(numElements,
                                                       localNumElements,
                                                       0,  // index base is zero
                                                       comm));
    return map.getConst();
  }

  XPETRA_FACTORY_END;
}

Teuchos::RCP<const Map<int, long long, EpetraNode>>
MapFactory<int, long long, EpetraNode>::copyMapWithNewComm(const Teuchos::RCP<const Map<int, long long, EpetraNode>> &oldmap,
                                                           const Teuchos::RCP<const Teuchos::Comm<int>> &newComm) {
  XPETRA_MONITOR("MapFactory::Build");
  using XMF             = Xpetra::MapFactory<int, ilong long, EpetraNode>;
  global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();

  size_t Nlocal         = oldmap->getLocalNumElements();
  global_size_t Nglobal = oldmap->getGlobalNumElements();

  // Sanity check -- if there's no comm, we can't keep elements on the map  (vice versa is OK)
  TEUCHOS_TEST_FOR_EXCEPTION(Nlocal && newComm.is_null(),
                             std::logic_error, "MapFactory::copyMapWithNewComm needs the comm to match the map.");

  // We'll return null if we don't have a Comm on this rank
  RCP<const Map<int, long long, Node>> newMap;
  if (!newComm.is_null()) {
    if (oldmap->isContiguous()) {
      newMap = XMF::Build(oldmap->lib(), INVALID, Nlocal, oldmap->getIndexBase(), newComm);
    } else {
      newMap = XMF::Build(oldmap->lib(), Nglobal, oldmap->getLocalElementList(), oldmap->getIndexBase(), newComm);
    }
  }

  return newMap;
  XPETRA_FACTORY_END;
}

#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

#endif  // #if defined(HAVE_XPETRA_EPETRA)

}  // namespace Xpetra

// EOF
