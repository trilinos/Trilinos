// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDMAPFACTORY_DEF_HPP
#define XPETRA_STRIDEDMAPFACTORY_DEF_HPP

#include "Xpetra_StridedMapFactory_decl.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          GlobalOrdinal indexBase,
          std::vector<size_t>& stridingInfo,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          LocalOrdinal stridedBlockId,
          GlobalOrdinal offset,
          LocalGlobal lg) {
  return rcp(new Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>(lib, numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, lg));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          size_t numLocalElements,
          GlobalOrdinal indexBase,
          std::vector<size_t>& stridingInfo,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          LocalOrdinal stridedBlockId,
          GlobalOrdinal offset) {
  return rcp(new StridedMap(lib, numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(const RCP<const Map>& map,
          std::vector<size_t>& stridingInfo,
          LocalOrdinal stridedBlockId,
          GlobalOrdinal offset) {
  return rcp(new StridedMap(map, stridingInfo, map->getIndexBase(), stridedBlockId, offset));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(const RCP<const StridedMap>& map, LocalOrdinal stridedBlockId) {
  TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId < 0, Exceptions::RuntimeError,
                             "Xpetra::StridedMapFactory::Build: constructor expects stridedBlockId > -1.");

  TEUCHOS_TEST_FOR_EXCEPTION(map->getStridedBlockId() != -1, Exceptions::RuntimeError,
                             "Xpetra::StridedMapFactory::Build: constructor expects a full map (stridedBlockId == -1).");

  std::vector<size_t> stridingInfo = map->getStridingData();

  Teuchos::ArrayView<const GlobalOrdinal> dofGids = map->getLocalElementList();

  // determine nStridedOffset
  size_t nStridedOffset = 0;
  for (int j = 0; j < map->getStridedBlockId(); j++) {
    nStridedOffset += stridingInfo[j];
  }

  const size_t numMyBlockDofs = (stridingInfo[stridedBlockId] * map->getLocalNumElements()) / map->getFixedBlockSize();

  std::vector<GlobalOrdinal> subBlockDofGids(numMyBlockDofs);

  // TODO fill vector with dofs
  LocalOrdinal ind = 0;
  for (typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it = dofGids.begin(); it != dofGids.end(); ++it) {
    if (map->GID2StridingBlockId(*it) == Teuchos::as<size_t>(stridedBlockId)) {
      subBlockDofGids[ind++] = *it;
    }
  }

  const Teuchos::ArrayView<const GlobalOrdinal> subBlockDofGids_view(&subBlockDofGids[0], subBlockDofGids.size());

  return rcp(new StridedMap(map->lib(),
                            Teuchos::OrdinalTraits<global_size_t>::invalid(),
                            subBlockDofGids_view,
                            map->getIndexBase(),
                            stridingInfo,
                            map->getComm(),
                            stridedBlockId));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(const StridedMap& map) {
  XPETRA_MONITOR("MapFactory::Build");

  LocalOrdinal N                                      = map.getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> oldElements = map.getLocalElementList();
  Teuchos::Array<GlobalOrdinal> newElements(map.getLocalNumElements());

  for (LocalOrdinal i = 0; i < N; i++)
    newElements[i] = oldElements[i];

  std::vector<size_t> strData = map.getStridingData();
  return rcp(new StridedMap(map.lib(), map.getGlobalNumElements(), newElements, map.getIndexBase(), strData, map.getComm(), map.getStridedBlockId()));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
          GlobalOrdinal indexBase,
          std::vector<size_t>& stridingInfo,
          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          LocalOrdinal stridedBlockId,  // FIXME (mfh 03 Sep 2014) This breaks if LocalOrdinal is unsigned
          GlobalOrdinal /* offset */) {
  return rcp(new StridedMap(lib, numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId));
}

}  // namespace Xpetra

#endif  // XPETRA_STRIDEDMAPFACTORY_DEF_HPP

// TODO: removed unused methods
