// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_DEF_HPP_
#define PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_DEF_HPP_

#include "Xpetra_BlockedMap_decl.hpp"

#include "Xpetra_Exceptions.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_TpetraImport.hpp"

namespace Xpetra {

// Unwrap a vector of Xpetra maps into Tpetra maps.  A non-Tpetra (e.g. Epetra)
// map triggers Xpetra::Exceptions::BadCast inside toTpetra -- this is the
// intended Epetra guard in the Tpetra-only world.
namespace BlockedMapDetails {
template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>
toTpetraMaps(const std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>& maps) {
  std::vector<Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>> tmaps(maps.size());
  for (size_t i = 0; i < maps.size(); ++i)
    tmaps[i] = (maps[i] == Teuchos::null) ? Teuchos::null : toTpetra(maps[i]);
  return tmaps;
}
}  // namespace BlockedMapDetails

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap() {
  map_ = Teuchos::rcp(new tpetra_blockedmap_type());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const RCP<const Map>& fullmap, const std::vector<RCP<const Map>>& maps, bool bThyraMode) {
  map_ = Teuchos::rcp(new tpetra_blockedmap_type(toTpetra(fullmap),
                                                 BlockedMapDetails::toTpetraMaps<LocalOrdinal, GlobalOrdinal, Node>(maps),
                                                 bThyraMode));
  // Retain the original Xpetra sub-maps so getMap() preserves nested
  // BlockedMap identity (the Tpetra core keeps only flattened plain maps).
  // In Xpetra-style numbering these are the returned sub-maps; in Thyra-style
  // numbering the incoming maps are the Thyra-numbered sub-maps.
  if (bThyraMode)
    thyraMapsXpetra_ = maps;
  else
    mapsXpetra_ = maps;
  // Retain the original full map too, so getFullMap() preserves its dynamic
  // type (e.g. StridedMap).  The Tpetra core flattens it to a plain Map, so
  // re-wrapping map_->getFullMap() would lose that richer type -- which MueLu
  // factories (SubBlockAFactory, BlockedPFactory) rely on via rcp_dynamic_cast.
  fullmapXpetra_ = fullmap;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const std::vector<RCP<const Map>>& maps, const std::vector<RCP<const Map>>& thyramaps) {
  map_ = Teuchos::rcp(new tpetra_blockedmap_type(BlockedMapDetails::toTpetraMaps<LocalOrdinal, GlobalOrdinal, Node>(maps),
                                                 BlockedMapDetails::toTpetraMaps<LocalOrdinal, GlobalOrdinal, Node>(thyramaps)));
  // Retain the original Xpetra sub-maps so getMap() preserves nested
  // BlockedMap identity for both numbering styles.
  mapsXpetra_      = maps;
  thyraMapsXpetra_ = thyramaps;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const Teuchos::RCP<const tpetra_blockedmap_type>& map)
  : map_(map) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const BlockedMap& input) {
  map_ = Teuchos::rcp(new tpetra_blockedmap_type(*input.getTpetra_BlockedMap()));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedMap() {
  fullmapXpetra_ = Teuchos::null;
  mapsXpetra_.clear();
  thyraMapsXpetra_.clear();
  importersXpetra_.clear();
  map_ = Teuchos::null;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalNumElements() const {
  return map_->getGlobalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::getLocalNumElements() const {
  return map_->getLocalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getIndexBase() const {
  return map_->getIndexBase();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinLocalIndex() const {
  return map_->getMinLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxLocalIndex() const {
  return map_->getMaxLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinGlobalIndex() const {
  return map_->getMinGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxGlobalIndex() const {
  return map_->getMaxGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinAllGlobalIndex() const {
  return map_->getMinAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxAllGlobalIndex() const {
  return map_->getMaxAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalElement(GlobalOrdinal globalIndex) const {
  return map_->getLocalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalElement(LocalOrdinal localIndex) const {
  return map_->getGlobalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                       const Teuchos::ArrayView<int>& nodeIDList,
                       const Teuchos::ArrayView<LocalOrdinal>& LIDList) const {
  return toXpetra(map_->getRemoteIndexList(GIDList, nodeIDList, LIDList));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                       const Teuchos::ArrayView<int>& nodeIDList) const {
  return toXpetra(map_->getRemoteIndexList(GIDList, nodeIDList));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayView<const GlobalOrdinal>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalElementList() const {
  return map_->getLocalElementList();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Map<LocalOrdinal, GlobalOrdinal, Node>::global_indices_array_device_type
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMyGlobalIndicesDevice() const {
  return map_->getMyGlobalIndicesDevice();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isNodeLocalElement(LocalOrdinal localIndex) const {
  return map_->isNodeLocalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isNodeGlobalElement(GlobalOrdinal globalIndex) const {
  return map_->isNodeGlobalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isContiguous() const {
  return map_->isContiguous();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isDistributed() const {
  return map_->isDistributed();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isCompatible(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const {
  RCP<const Map> rcpMap         = Teuchos::rcpFromRef(map);
  RCP<const BlockedMap> rcpBMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rcpMap);
  if (!rcpBMap.is_null())
    return map_->isCompatible(*rcpBMap->getTpetra_BlockedMap());
  // A plain (non-blocked) map: dispatch to the Tpetra plain-map overload.
  return map_->isCompatible(*toTpetra(rcpMap));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isSameAs(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const {
  RCP<const Map> rcpMap         = Teuchos::rcpFromRef(map);
  RCP<const BlockedMap> rcpBMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rcpMap);
  if (!rcpBMap.is_null())
    return map_->isSameAs(*rcpBMap->getTpetra_BlockedMap());
  // A plain (non-blocked) map: dispatch to the Tpetra plain-map overload
  // (which handles the single-block special case).
  return map_->isSameAs(*toTpetra(rcpMap));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Teuchos::Comm<int>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getComm() const {
  return map_->getComm();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>&
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
operator=(const BlockedMap& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::getThyraMode() const {
  return map_->getThyraMode();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    removeEmptyProcesses() const {
  throw Xpetra::Exceptions::RuntimeError("BlockedMap::removeEmptyProcesses: routine not implemented.");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    replaceCommWithSubset(const Teuchos::RCP<const Teuchos::Comm<int>>& /* newComm */) const {
  throw Xpetra::Exceptions::RuntimeError("BlockedMap::replaceCommWithSubset: routine not implemented.");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
UnderlyingLib
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::lib() const {
  // The wrapped object is always Tpetra-backed.
  return Xpetra::UseTpetra;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return getFullMap();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getNumMaps() const {
  return map_->getNumMaps();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMap(size_t i,
           bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(i >= getNumMaps(),
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedMap::getMap: tried to access block " << i << ", but BlockedMap has only " << getNumMaps()
                                                                         << " blocks! Block indices must be between 0 and " << getNumMaps() - 1
                                                                         << ".");
  // Size the identity caches without discarding any entries the constructor
  // pre-populated: when this BlockedMap was built from Xpetra sub-maps, those
  // originals are cached here so nested BlockedMap identity survives
  // getMap()/rcp_dynamic_cast (the Tpetra core stores only flattened plain
  // maps and cannot represent a nested BlockedMap).  resize() preserves
  // existing slots; empty slots are filled lazily from the Tpetra core below.
  const size_t n = getNumMaps();
  if (mapsXpetra_.size() != n)
    mapsXpetra_.resize(n);
  if (thyraMapsXpetra_.size() != n)
    thyraMapsXpetra_.resize(n);

  if (map_->getThyraMode() == true && bThyraMode == true) {
    if (thyraMapsXpetra_[i].is_null())
      thyraMapsXpetra_[i] = toXpetra(map_->getMap(i, true));
    return thyraMapsXpetra_[i];
  }

  XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedMap::getMap: cannot return sub map in Thyra-style numbering if BlockedMap object is not created using "
                            "Thyra-style numbered submaps.");
  if (mapsXpetra_[i].is_null())
    mapsXpetra_[i] = toXpetra(map_->getMap(i, false));
  return mapsXpetra_[i];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getImporter(size_t i) const {
  XPETRA_TEST_FOR_EXCEPTION(i >= getNumMaps(),
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedMap::getImporter: tried to access block " << i << ", but BlockedMap has only " << getNumMaps()
                                                                              << " blocks! Block indices must be between 0 and " << getNumMaps() - 1
                                                                              << ".");
  if (importersXpetra_.size() != getNumMaps()) {
    const size_t n = getNumMaps();
    importersXpetra_.assign(n, Teuchos::null);
    for (size_t k = 0; k < n; ++k) {
      Teuchos::RCP<Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> imp = map_->getImporter(k);
      if (!imp.is_null())
        importersXpetra_[k] = Teuchos::rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(imp));
    }
  }
  return importersXpetra_[i];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getFullMap() const {
  if (fullmapXpetra_.is_null() && !map_->getFullMap().is_null())
    fullmapXpetra_ = toXpetra(map_->getFullMap());
  return fullmapXpetra_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMapIndexForGID(GlobalOrdinal gid) const {
  return map_->getMapIndexForGID(gid);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  return std::string("BlockedMap");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  map_->describe(out, verbLevel);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    assign(const BlockedMap& input) {
  map_           = Teuchos::rcp(new tpetra_blockedmap_type(*input.getTpetra_BlockedMap()));
  fullmapXpetra_ = Teuchos::null;
  mapsXpetra_.clear();
  thyraMapsXpetra_.clear();
  importersXpetra_.clear();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    concatenateMaps(const std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>& subMaps) {
  return toXpetra(tpetra_blockedmap_type::concatenateMaps(
      BlockedMapDetails::toTpetraMaps<LocalOrdinal, GlobalOrdinal, Node>(subMaps)));
}

}  // namespace Xpetra
#endif /* PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_DEF_HPP_ */
